/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "VMSystem.H"

#include <sstream>

#include "tbox/MathUtilities.H"
#include "LoadBalancer.H"
#include "Loki_Defines.H"
#include "TimerManager.H"
#include "RestartManager.H"
#include "RK4Integrator.H"
#include "RK6Integrator.H"

namespace Loki {

VMSystem::VMSystem(
   LokiInputParser& a_pp,
   int a_num_probes,
   int a_spatial_solution_order,
   int a_temporal_solution_order,
   bool a_coll_diag_on)
   : System(a_pp, a_spatial_solution_order, a_temporal_solution_order)
{
   // Create all the timers specific to a Vlasov-Maxwell system.
   TimerManager* timers(TimerManager::getManager());
   timers->createTimer("current driver");
   timers->createTimer("current density");
   timers->createTimer("BC (Maxwell)");
   timers->createTimer("Maxwell RHS");

   // Read all user input.
   parseParameters(a_pp);

   // Create the state for this system, essentially all the KineticSpecies and
   // a Maxwell object.
   createVMState(a_pp, a_coll_diag_on);

   // Now that all the loads are defined, load balance the system.
   loadBalance();

   // This is something of a hack.  We only want all the particles, especially
   // the noise source particles, to exist on the Maxwell object(s) on the
   // Maxwell processor(s).  But we only know which Maxwell object(s) are on
   // Maxwell processor(s) after load balancing.  So we read the particle data
   // here.  As the species need to know if there are particles in the problem
   // we tell them that here as well.
   KineticSpeciesPtrVect& state_kinetic_species = m_state.kineticSpecies();
   Maxwell& state_maxwell = *m_state.maxwell();
   int num_kinetic_species = static_cast<int>(state_kinetic_species.size());

   double vlo[2], vhi[2];
   state_kinetic_species.FindSpeciesVelocityBoundingBox(vlo, vhi);

   LokiInputParser ppp("maxwell");
   state_maxwell.readParticles(ppp, vlo, vhi);
   bool problem_has_particles =
      state_maxwell.numProblemTrackingParticles() > 0 ||
      state_maxwell.numProblemNoiseSourceParticles() > 0;
   if (problem_has_particles) {
      for (int s(0); s < num_kinetic_species; ++s) {
         state_kinetic_species[s]->problemHasParticles();
      }
   }

   // The system writes restart data so register it with the RestartManager
   // which will use the putToRestart/getFromRestart callbacks to get the
   // restart data written/read.
   RestartManager* restart_manager(RestartManager::getManager());
   restart_manager->registerRestart(this);

   // Allocate global arrays defined on configuration space for the net current
   // density and the charge and current densities for each species.  If there
   // are particles we need an array for the net external E field and the
   // external E field from each species' drivers.  Additionally, for certain
   // current drivers we need an array for the effective frequency at which its
   // generated EM wave propagates through the plasma.
   bool needs_disp_rel = state_maxwell.needsDispersionRelation();
   if (needs_disp_rel) {
      m_charge_density.resize(num_kinetic_species);
   }
   state_maxwell.newAuxVariable(m_net_x_current_density, 1);
   m_x_current_density.resize(num_kinetic_species);
   state_maxwell.newAuxVariable(m_net_y_current_density, 1);
   m_y_current_density.resize(num_kinetic_species);
   state_maxwell.newAuxVariable(m_net_z_current_density, 1);
   m_z_current_density.resize(num_kinetic_species);
   if (needs_disp_rel) {
     state_maxwell.newAuxVariable(m_omega_eff2, 1);
   }
   if (problem_has_particles) {
      state_maxwell.newAuxVariable(m_net_ext_efield, 2);
      m_ext_efield.resize(num_kinetic_species);
   }
   for (int s(0); s < num_kinetic_species; ++s) {
      if (needs_disp_rel) {
         state_maxwell.newAuxVariable(m_charge_density[s], 1);
      }
      state_maxwell.newAuxVariable(m_x_current_density[s], 1);
      state_maxwell.newAuxVariable(m_y_current_density[s], 1);
      state_maxwell.newAuxVariable(m_z_current_density[s], 1);
      if (problem_has_particles) {
         state_maxwell.newAuxVariable(m_ext_efield[s], 2);
      }
   }

   // Compute the number of individual time histories and allocate initial space
   // for the time histories and time stamps.
   m_num_seq =
      Maxwell::GLOBAL_TIME_HISTS +
      Maxwell::TIME_HISTS_PER_PROBE * a_num_probes +
      state_maxwell.numProblemTrackingParticles() * 4 +
      num_kinetic_species * 16;
   m_sequences.resize(m_num_seq);
   for (int i1 = 0; i1 < m_num_seq; ++i1) {
      vector<double>& this_seq = m_sequences[i1];
      this_seq.resize(m_length_seq);
      for (int i2 = 0; i2 < m_length_seq; ++i2) {
         this_seq[i2] = 0.0;
      }
   }
   m_time_seq.resize(m_length_seq);
   for (int i1(0); i1 < m_length_seq; ++i1) {
      m_time_seq[i1] = 0.0;
   }

   // Collision version
   m_num_seq_coll = Maxwell::NUM_COLL_TIME_HISTS;
   m_sequences_coll.resize(m_num_seq_coll);
   for (int i1(0); i1 < m_num_seq_coll; ++i1) {
      vector<double>& this_seq = m_sequences_coll[i1];
      this_seq.resize(m_length_seq_coll);
      for (int i2(0); i2 < m_length_seq_coll; ++i2) {
         this_seq[i2] = 0.0;
      }
   }
   m_time_seq_coll.resize(m_length_seq_coll);
   for (int i1(0); i1 < m_length_seq_coll; ++i1) {
      m_time_seq_coll[i1] = 0.0;
   }
}


VMSystem::~VMSystem()
{
}


void
VMSystem::initialize(
   bool a_is_from_restart,
   double a_time,
   double a_dt,
   bool a_coll_diag_on)
{
   // If this is an initial run (not from a restart file) then the
   // KineticSpecies' distribution functions must be initialized and the initial
   // condition stored.  Otherwise we just need to store the initial condition.
   Maxwell& state_maxwell = *m_state.maxwell();
   state_maxwell.initialize(static_cast<int>(m_species_names.size()),
                            m_plot_ke_vel_bdy_flux,
                            a_is_from_restart,
                            a_coll_diag_on);
   // Initialize each species distributed to this processor.
   KSPV::Iterator it_end(m_state.kineticSpecies().end_locals());
   for (KSPV::Iterator it(m_state.kineticSpecies().begin_locals());
        it != it_end; ++it) {
      (*it)->initialize(a_is_from_restart);
   }

   // Clone all the species and Maxwell from m_state into m_state_old.
   KineticSpeciesPtrVect& state_kinetic_species = m_state.kineticSpecies();
   PopulateInterspeciesData(state_kinetic_species);   
   KineticSpeciesPtrVect& old_state_kinetic_species =
      m_state_old.kineticSpecies();
   int num_kinetic_species = static_cast<int>(state_kinetic_species.size());
   old_state_kinetic_species.resize(num_kinetic_species);
   for (int s(0); s < num_kinetic_species; ++s) {
      old_state_kinetic_species[s] = state_kinetic_species[s]->clone();
   }
   m_state_old.maxwell() = state_maxwell.clone();

   // Build the time integrator of the requested order.
   if (m_temporal_solution_order == 4) {
      m_integrator = new RK4Integrator<VMState>(*this, m_state);
   }
   else {
      m_integrator = new RK6Integrator<VMState>(*this, m_state);
   }

   // Every processor participates in the computation of the current density of
   // each species.  If this loop was to be replaced by an iterator over
   // species then the Maxwell processor would do nothing and the final
   // communication to m_*_current_density which are owned by the Maxwell
   // processor would not be completed and the simulation would be wrong
   // however the code would not hang.  This is something of a mystery as one
   // would expect the senders to wait for their communication to be fulfilled
   // which it appears does not happen.
   // Additionally, embedded in this calculation is a communication of the 2D
   // transverse drift velocity from the Maxwell processor to each species so
   // that the current densities may be computed.  So the Maxwell processor
   // must be involved.
   bool needs_disp_rel = state_maxwell.needsDispersionRelation();
   for (int s(0); s < num_kinetic_species; ++s) {
      state_kinetic_species[s]->currentDensity(state_maxwell,
         m_x_current_density[s],
         m_y_current_density[s],
         m_z_current_density[s]);
      if (needs_disp_rel) {
         state_kinetic_species[s]->chargeDensity(m_charge_density[s]);
      }
   }

   // Only the Maxwell processor(s) are needed to apply any external currents.
   if (state_maxwell.isEMSolverProcessor()) {
      m_net_x_current_density = 0;
      m_net_y_current_density = 0;
      m_net_z_current_density = 0;
      if (needs_disp_rel) {
         m_omega_eff2 = 0;
      }
      for (int s(0); s < num_kinetic_species; ++s) {
         m_net_x_current_density += m_x_current_density[s];
         m_net_y_current_density += m_y_current_density[s];
         m_net_z_current_density += m_z_current_density[s];
         if (needs_disp_rel) {
            double factor = state_kinetic_species[s]->charge()/
               state_kinetic_species[s]->mass();
            m_charge_density[s] *= factor;
            m_omega_eff2 += m_charge_density[s];
         }
      }
      state_maxwell.computeAntennaSource(a_time, m_omega_eff2);
      state_maxwell.applyNoiseParticleCurrent(m_net_x_current_density,
         m_net_y_current_density,
         a_time);
   }

   // If this is a restarted run we need to compute the acceleration given the
   // restarted state so that the time step computation will be correct.  If
   // this is at t=0 the accelerations have been properly initialized.  This
   // computation requires that the 2D electricmagnetic field be communicated
   // from the Maxwell processor(s) to Vlasov processors so every processor must
   // participate.
   if (a_is_from_restart) {
      for (int s(0); s < num_kinetic_species; ++s) {
         state_kinetic_species[s]->computeAcceleration(state_maxwell,
            m_ext_efield[s],
            a_time,
            a_dt,
            1);
      }
   }

   // Make an initial call to evalRHS just to set reasonable values for
   // KineticSpecies::m_lambda_max, which will be used by
   // KineticSpecies::computeDt in determining the first timestep size.
   evalRHS(m_state_old, m_state, 0.0, 1.0e-16, false, false);
}


void
VMSystem::defineRHSData(
   VMState&       a_rhs,
   const VMState& a_prototype)
{
   // We need to define all species on all processors of a_rhs.
   KineticSpeciesPtrVect& rhs_kspv = a_rhs.kineticSpecies();
   const KineticSpeciesPtrVect& prototype_kspv = a_prototype.kineticSpecies();
   int num_kinetic_species = static_cast<int>(prototype_kspv.size());
   rhs_kspv.resize(num_kinetic_species);
   for (int s(0); s < num_kinetic_species; ++s) {
      rhs_kspv[s] = prototype_kspv[s]->clone();
   }
   a_rhs.maxwell() = a_prototype.maxwell()->clone();
}


void
VMSystem::copySolnData(
   VMState&       a_dst,
   const VMState& a_src)
{
   // Copy the solution data from a_src to a_dst only for the species
   // distributed to this processor.
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("CopySolnData");
   KineticSpeciesPtrVect& dst_kspv = a_dst.kineticSpecies();
   const KineticSpeciesPtrVect& src_kspv = a_src.kineticSpecies();
   KSPV::ConstIterator src_it_end(src_kspv.end_locals());
   KSPV::Iterator dst_it(dst_kspv.begin_locals());
   for (KSPV::ConstIterator src_it(src_kspv.begin_locals());
        src_it != src_it_end; ++src_it, ++dst_it) {
      (*dst_it)->copySolnData(*(*src_it));
   }
   tbox::Pointer<Maxwell>& a_dst_maxwell = a_dst.maxwell();
   if (a_dst_maxwell->isEMSolverProcessor()) {
      a_dst_maxwell->copySolnData(*a_src.maxwell());
   }
   timers->stopTimer("CopySolnData");
}


void
VMSystem::zeroSolnData(
   VMState& a_soln)
{
   // Zero the solution data of a_soln for the species distributed to this
   // processor.
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("ZeroSolnData");
   KineticSpeciesPtrVect& soln_kspv = a_soln.kineticSpecies();
   KSPV::Iterator it_end(soln_kspv.end_locals());
   for (KSPV::Iterator it(soln_kspv.begin_locals()); it != it_end; ++it) {
      (*it)->zeroData();
   }
   tbox::Pointer<Maxwell>& a_soln_maxwell = a_soln.maxwell();
   if (a_soln_maxwell->isEMSolverProcessor()) {
      a_soln_maxwell->zeroData();
   }
   timers->stopTimer("ZeroSolnData");
}


void
VMSystem::addSolnData(
   VMState&       a_soln,
   const VMState& a_increment,
   double         a_scale,
   bool           a_sum_reduce_inc)
{
   // Add the solution data to a_soln for the species distributed to this
   // processor.
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("AddSolnData");
   KineticSpeciesPtrVect& soln_kspv = a_soln.kineticSpecies();
   const KineticSpeciesPtrVect& increment_kspv = a_increment.kineticSpecies();
   KSPV::Iterator soln_it_end(soln_kspv.end_locals());
   KSPV::ConstIterator increment_it(increment_kspv.begin_locals());
   for (KSPV::Iterator soln_it(soln_kspv.begin_locals());
        soln_it != soln_it_end; ++soln_it, ++increment_it) {
      (*soln_it)->addData(*(*increment_it), a_scale);
   }
   tbox::Pointer<Maxwell>& a_soln_maxwell = a_soln.maxwell();
   if (a_soln_maxwell->isEMSolverProcessor()) {
      a_soln_maxwell->addData(*a_increment.maxwell(),
         a_scale,
         a_sum_reduce_inc);
   }
   timers->stopTimer("AddSolnData");
}


bool
VMSystem::validRHSData(
   const VMState& a_rhs,
   const VMState& a_protoSoln)
{
   const KineticSpeciesPtrVect& rhs_kspv = a_rhs.kineticSpecies();
   const KineticSpeciesPtrVect& protoSoln_kspv = a_protoSoln.kineticSpecies();
   bool valid = (rhs_kspv.size() == protoSoln_kspv.size());
   if (valid) {
      // Check the validity of only the species distributed to this processor
      // for a_rhs.
      KSPV::ConstIterator rhs_it_end(rhs_kspv.end_locals());
      KSPV::ConstIterator protoSoln_it(protoSoln_kspv.begin_locals());
      for (KSPV::ConstIterator rhs_it(rhs_kspv.begin_locals());
           rhs_it != rhs_it_end; ++rhs_it, ++protoSoln_it) {
         valid &= (*rhs_it)->conformsTo(*(*protoSoln_it), false);
      }
   }
   valid &= a_rhs.maxwell()->conformsTo(*a_protoSoln.maxwell(), false);
   return Loki_Utilities::reduceBoolean(valid);
}


bool
VMSystem::validSolnData(
   const VMState& a_soln,
   const VMState& a_protoSoln)
{
   const KineticSpeciesPtrVect& soln_kspv = a_soln.kineticSpecies();
   const KineticSpeciesPtrVect& protoSoln_kspv = a_protoSoln.kineticSpecies();
   bool valid = (soln_kspv.size() == protoSoln_kspv.size());
   if (valid) {
      // Check the validity of only the species distributed to this processor
      // for a_soln.
      KSPV::ConstIterator soln_it_end(soln_kspv.end_locals());
      KSPV::ConstIterator protoSoln_it(protoSoln_kspv.begin_locals());
      for (KSPV::ConstIterator soln_it(soln_kspv.begin_locals());
           soln_it != soln_it_end; ++soln_it, ++protoSoln_it) {
         valid &= (*soln_it)->conformsTo(*(*protoSoln_it));
      }
   }
   valid &= a_soln.maxwell()->conformsTo(*a_protoSoln.maxwell());
   return Loki_Utilities::reduceBoolean(valid);
}


void
VMSystem::evalRHS(
   VMState& a_rhs,
   VMState& a_state,
   double   a_time,
   double   a_dt,
   bool     a_first_rk_stage,
   bool     a_last_rk_stage)
{
   KineticSpeciesPtrVect& rhs_kspv = a_rhs.kineticSpecies();
   Maxwell& rhs_maxwell = *a_rhs.maxwell();
   KineticSpeciesPtrVect& state_kspv = a_state.kineticSpecies();
   Maxwell& state_maxwell = *a_state.maxwell();
   int num_kinetic_species = static_cast<int>(state_kspv.size());

   // Every processor participates in the computation of the current density of
   // each species.  If this loop was to be replaced by an iterator over
   // species then the Maxwell processor would do nothing and the final
   // communication to m_*_current_density which are owned by the Maxwell
   // processor would not be completed and the simulation would be wrong
   // however the code would not hang.  This is something of a mystery as one
   // would expect the senders to wait for their communication to be fulfilled
   // which it appears does not happen.
   // Additionally, embedded in this calculation is a communication of the 2D
   // transverse drift velocity from the Maxwell processor to each species so
   // that the current densities may be computed.  So the Maxwell processor
   // must be involved.
   bool needs_disp_rel = state_maxwell.needsDispersionRelation();
   for (int s(0); s < num_kinetic_species; ++s) {
      state_kspv[s]->currentDensity(state_maxwell,
         m_x_current_density[s],
         m_y_current_density[s],
         m_z_current_density[s]);
      if (needs_disp_rel) {
         state_kspv[s]->chargeDensity(m_charge_density[s]);
      }
   }

   // The external current application by the Maxwell processor(s) and the
   // advection computation by the Vlasov processors are done in parallel as
   // these computations are independent and the processor sets do not overlap.
   // The computation of the acceleration that follows requires parallel
   // communication between all of the processors.  If load imbalance exists
   // between the Maxwell and Vlasov processors the acceleration computation
   // will appear to be a bottleneck as its parallel communication can not
   // proceed until all processors have reached that point.
   if (state_maxwell.isEMSolverProcessor()) {
      state_maxwell.fillGhostCells(false);
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("current density");
      m_net_x_current_density = 0.0;
      m_net_y_current_density = 0.0;
      m_net_z_current_density = 0.0;
      timers->stopTimer("current density");
      if (needs_disp_rel) {
         m_omega_eff2 = 0;
      }
      for (int s(0); s < num_kinetic_species; ++s) {
         timers->startTimer("current density");
         m_net_x_current_density += m_x_current_density[s];
         m_net_y_current_density += m_y_current_density[s];
         m_net_z_current_density += m_z_current_density[s];
         timers->stopTimer("current density");
         if (needs_disp_rel) {
            double factor = state_kspv[s]->charge()/
               state_kspv[s]->mass();
            m_charge_density[s] *= factor;
            m_omega_eff2 += m_charge_density[s];
         }
      }
      state_maxwell.computeAntennaSource(a_time, m_omega_eff2);
      state_maxwell.applyNoiseParticleCurrent(m_net_x_current_density,
         m_net_y_current_density,
         a_time);
   }

   // Perform the advection computation for each species distributed to each
   // processor.
   KSPV::Iterator a_state_it_end(state_kspv.end_locals());
   KSPV::Iterator a_rhs_it(rhs_kspv.begin_locals());
   for (KSPV::Iterator a_state_it(state_kspv.begin_locals());
        a_state_it != a_state_it_end; ++a_state_it, ++a_rhs_it) {
      (*a_state_it)->fillAdvectionGhostCells();
      if (m_do_new_algorithm) {
         (*a_state_it)->evalAdvectionDerivatives(*(*a_rhs_it));
      }
      else {
         (*a_state_it)->evalAdvectionFluxes();
      }
   }

   // Compute the acceleration which is needed in the KineticSpecies RHS
   // evaluation below.  This computation requires that the 2D electricmagnetic
   // field be communicated from the Maxwell processor(s) to Vlasov processors
   // so every processor must participate.
   for (int s(0); s < num_kinetic_species; ++s) {
      state_kspv[s]->computeAcceleration(state_maxwell,
         m_ext_efield[s],
         a_time,
         a_dt,
         a_first_rk_stage);
   }

   // Now the Vlasov processors will fill the acceleration boundary conditions
   // (boundaries of velocity domain) of each species distributed on them and
   // compute the Vlasov equation RHS.  The Maxwell processor(s) will compute
   // the Maxwell equations' RHS in parallel.
   a_rhs_it.reset();
   for (KSPV::Iterator a_state_it(state_kspv.begin_locals());
        a_state_it != a_state_it_end; ++a_state_it, ++a_rhs_it) {
      (*a_state_it)->setAccelerationBCs();
      if (m_do_new_algorithm) {
         (*a_state_it)->evalAccelerationDerivatives(*(*a_rhs_it));
      }
      else {
         (*a_state_it)->evalAccelerationFluxes();
      }
      (*a_state_it)->completeRHS(*(*a_rhs_it),
         a_time,
         a_dt,
         m_do_new_algorithm,
         a_first_rk_stage);
   }
   if (state_maxwell.isEMSolverProcessor()) {
      // If there are particles then compute the net externally applied E field
      // which is needed to update them.
      if (state_maxwell.numProblemTrackingParticles() > 0 ||
          state_maxwell.numProblemNoiseSourceParticles() > 0) {
         m_net_ext_efield = 0.0;
         for (int s(0); s < num_kinetic_species; ++s) {
            m_net_ext_efield += m_ext_efield[s];
         }
         m_net_ext_efield.communicateGhostData();
      }
      state_maxwell.evalRHS(rhs_maxwell,
         rhs_kspv,
         m_net_x_current_density,
         m_net_y_current_density,
         m_net_z_current_density,
         m_net_ext_efield,
         a_time);
   }
}


void
VMSystem::postStageAdvance(
   VMState& a_soln,
   int a_stage)
{
   NULL_USE(a_soln);
   NULL_USE(a_stage);
}


double
VMSystem::stableDt()
{
   double dt_stable(tbox::MathUtilities<double>::getMax());
   KineticSpeciesPtrVect& state_kinetic_species = m_state.kineticSpecies();
   // Compute the time step contraint for each species distributed to this
   // processor and take the min of them.
   KSPV::Iterator it_end(state_kinetic_species.end_locals());
   for (KSPV::Iterator it(state_kinetic_species.begin_locals());
        it != it_end; ++it) {
      double dt_species = (*it)->computeDt();
      dt_stable = min(dt_stable, dt_species);
   }
   double dt_maxwell = m_state.maxwell()->computeDt();
   dt_stable = min(dt_stable, dt_maxwell);
   dt_stable = Loki_Utilities::getMinValue(dt_stable);
   
   return dt_stable;
}


double
VMSystem::advance(
   double a_cur_time,
   double a_dt)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("other");
   // Copy current to old and advance solution via RK integration.
   copySolnData(m_state_old, m_state);
   double new_time = m_integrator->advance(m_state,
      m_state_old,
      a_cur_time,
      a_dt);

   timers->stopTimer("other");

   return new_time;
}


void
VMSystem::plot(
   double a_time,
   double a_dt,
   const vector<vector<double> >& a_probes,
   int a_num_probes,
   int a_saved_seq,
   int& a_saved_save,
   string& a_time_hist_file_name)
{
   tbox::Pointer<Maxwell>& state_maxwell = m_state.maxwell();
   if (m_plot_ke_vel_bdy_flux) {
      // Compute the KE velocity boundary flux which will be a plotted quantity.
      // As there is parallel communication in this computation all processors
      // must participate for all species.
      KineticSpeciesPtrVect& state_kinetic_species = m_state.kineticSpecies();
      int num_kinetic_species = static_cast<int>(state_kinetic_species.size());
      for (int s(0); s < num_kinetic_species; ++s) {
         state_kinetic_species[s]->computeKEVelBdyFlux(*state_maxwell);
      }
   }

   // The Maxwell object actually writes the plot data.
   if (state_maxwell->isEMSolverProcessor()) {
      state_maxwell->plot(a_time,
         a_dt,
         m_sequences,
         m_species_names,
         m_time_seq,
         a_probes,
         a_num_probes,
         a_saved_seq,
         a_saved_save,
         m_plot_ke_vel_bdy_flux,
         a_time_hist_file_name);
   }
}


void
VMSystem::plotColl(
   double a_time,
   double a_dt,
   int a_saved_seq,
   int& a_saved_save,
   string& a_time_hist_file_name)
{
   // Save the plot data
   tbox::Pointer<Maxwell>& state_maxwell = m_state.maxwell();
   KineticSpeciesPtrVect& state_kinetic_species = m_state.kineticSpecies();

   state_maxwell->plotCollision(a_time,
      a_dt,
      state_kinetic_species,
      m_sequences_coll,
      m_species_names);
}


void
VMSystem::accumulateSequences(
   double a_time,
   const vector<vector<double> >& a_probes,
   int a_num_probes,
   int& a_saved_seq)
{
   // resize m_time_seq and m_sequences if necessary
   if (a_saved_seq >= m_length_seq) {
      m_length_seq *= 2;
      m_time_seq.resize(m_length_seq);
      for (int i = a_saved_seq; i < m_length_seq; ++i) {
         m_time_seq[i] = 0.0;
      }
      for (int i = 0; i < m_num_seq; ++i) {
         vector<double>& this_seq = m_sequences[i];
         this_seq.resize(m_length_seq);
         for (int j = a_saved_seq; j < m_length_seq; ++j) {
            this_seq[j] = 0.0;
         }
      }
   }
   m_time_seq[a_saved_seq] = a_time;

   // Get the Maxwell object's time histories.
   tbox::Pointer<Maxwell>& state_maxwell = m_state.maxwell();
   int seq_idx = 0;
   state_maxwell->accumulateSequences(m_sequences,
      a_probes,
      a_num_probes,
      a_saved_seq,
      seq_idx);

   // Get each kinetic species time history.  There is parallel communication
   // involved here so each processor must participate for each species.
   KineticSpeciesPtrVect& state_kinetic_species = m_state.kineticSpecies();
   int num_kinetic_species = static_cast<int>(state_kinetic_species.size());
   for (int s(0); s < num_kinetic_species; ++s) {
      state_kinetic_species[s]->accumulateSequences(*state_maxwell,
         m_sequences,
         a_saved_seq,
         seq_idx,
         a_time);
   }

   // We've saved another time history.
   ++a_saved_seq;
}

void
VMSystem::accumulateCollisionSequences(
   double a_time,
   double a_dt,
   int a_coll_op_idx,
   int& a_saved_seq)
{
   LOKI_ABORT("Unimplemented");
}

void
VMSystem::getFromRestart(
   RestartReader& a_reader)
{
   // Read all the species names from the restart file and ensure that they
   // match the names in the input deck.  If they don't then the restart file
   // and input deck are incompatible.
   int species_list_size;
   a_reader.readIntegerValue("species_list_size", species_list_size);
   a_reader.readDoubleValue("bz_const", m_bz_const); //IEO
   int num_species_names = static_cast<int>(m_species_names.size());
   if (species_list_size != num_species_names) {
      LOKI_ABORT("Kinetic species names list length mismatch!");
   }

   a_reader.pushSubDir("species_list");
   for (int s(0); s < species_list_size; ++s) {
      string tmp;
      stringstream tag;
      tag << "species." << s + 1;
      a_reader.readString(tag.str(), tmp);

      int i;
      for (i = 0; i < num_species_names; ++i) {
         if (m_species_names[i] == tmp) {
            break;
         }
      }
      if (i == num_species_names) {
         LOKI_ABORT("Kinetic species names from restart don't match names from input file!");
      }
   }

   a_reader.popSubDir();
}

  
void
VMSystem::putToRestart(
   RestartWriter& a_writer,
   double a_time)
{
   NULL_USE(a_time);

   bool write_data = Loki_Utilities::s_my_id == 0;

   int species_list_size(static_cast<int>(m_species_names.size()));
   a_writer.writeIntegerValue("species_list_size",
      species_list_size,
      write_data);
   a_writer.writeDoubleValue("bz_const", m_bz_const, write_data); //IEO
   a_writer.writeIntegerValue("plot_ke_vel_bdy_flux",
      static_cast<int>(m_plot_ke_vel_bdy_flux),
      write_data);

   // Write the species names to a directory in the restart file.
   a_writer.pushSubDir("species_list");
   for (int s(0); s < species_list_size; ++s) {
      stringstream tag;
      tag << "species." << s + 1;
      a_writer.writeString(tag.str(), m_species_names[s], write_data);
   }
   a_writer.popSubDir();
}

long int
VMSystem::problemSize() const
{
   // The species are 4D and hence determine the size of the problem.  Ask each
   // species for its size and take the max.  Species not distributed to this
   // processor have a size of 0.
   long int size(0);
   const KineticSpeciesPtrVect& state_kinetic_species =
      m_state.kineticSpecies();
   int num_kinetic_species = static_cast<int>(state_kinetic_species.size());
   for (int s(0); s < num_kinetic_species; ++s) {
      size = max(size, state_kinetic_species[s]->numberOfCells());
   }
   return size;
}

void
VMSystem::updateGhosts(
   bool a_particles_only)
{
   // Update the Maxwell ghost cells.
   Maxwell& state_maxwell = *m_state.maxwell();
   if (state_maxwell.isEMSolverProcessor()) {
      state_maxwell.fillGhostCells(a_particles_only);
   }

   // Update the advection ghost cells for each species distributed to this
   // processor.
   if (!a_particles_only) {
      KineticSpeciesPtrVect& state_kinetic_species = m_state.kineticSpecies();
      KSPV::Iterator species_it_end(state_kinetic_species.end_locals());
      for (KSPV::Iterator species_it(state_kinetic_species.begin_locals());
           species_it != species_it_end; ++species_it) {
         (*species_it)->fillAdvectionGhostCells();
      }
   }
}

void
VMSystem::printParameters() const
{
   m_state.maxwell()->printParameters();
   int num_kinetic_species = static_cast<int>(m_state.kineticSpecies().size());
   for (int s = 0; s < num_kinetic_species; ++s) {
      m_state.kineticSpecies()[s]->printParameters();
   }
   printDecomposition();
}

////// PRIVATE FUNCTIONS ////////////////////////////////////
void
VMSystem::createVMState(
   LokiInputParser& a_pp,
   bool a_coll_diag_on)
{
   int number_of_species(1);
   a_pp.query("number_of_species", number_of_species);

   string test_str = "false";
   a_pp.query("use_new_bcs", test_str);
   bool use_new_bcs = test_str.compare("true") == 0 ? true : false;

   // Make the Maxwell object's sub-database and construct it.
   int plot_times_per_file = 1;
   a_pp.query("plot_times_per_file", plot_times_per_file);
   LokiInputParser ppp("maxwell");
   m_state.maxwell() = new Maxwell(ppp,
                                   m_cfg_domain,
                                   number_of_species,
                                   m_spatial_solution_order,
                                   plot_times_per_file,
                                   m_plot_ke_vel_bdy_flux,
                                   a_coll_diag_on,
                                   m_bz_const); //IEO

   // For each species, make its sub-database and construct it.
   KineticSpeciesPtrVect& state_kinetic_species = m_state.kineticSpecies();
   state_kinetic_species.resize(number_of_species);

   for (int s(0); s < number_of_species; ++s) {

      ostringstream input_string;
      input_string << "kinetic_species." << s+1;
      LokiInputParser ppspecies(input_string.str().c_str());
      input_string.str("");

      state_kinetic_species[s] = new KineticSpecies(m_cfg_domain,
         ppspecies,
         s+1,
         number_of_species,
         m_spatial_solution_order,
         m_temporal_solution_order,
         m_plot_ke_vel_bdy_flux,
         use_new_bcs,
         true,
         m_bz_const);

      // Create list of species names
      m_species_names.push_back(state_kinetic_species[s]->name().c_str());
   }
}


void
VMSystem::loadBalance()
{
   // Every species, even those not distributed to this processor, is a load to
   // be balanced as is the Maxwell object.
   vector<Load*> load_vector;
   KineticSpeciesPtrVect& state_kinetic_species = m_state.kineticSpecies();
   int num_kinetic_species = static_cast<int>(state_kinetic_species.size());
   for (int s(0); s < num_kinetic_species; ++s) {
      load_vector.push_back(state_kinetic_species[s].getPointer());
   }
   load_vector.push_back(m_state.maxwell().getPointer());

   // Balance everything.
   LoadBalancer::balance(load_vector);
}


void
VMSystem::printDecomposition() const
{
   // Print the decompositions of the Maxwell and of each species.
   Loki_Utilities::printF("\n#*#*# VMSystem: Parallel Decomposition #*#*#\n");
   m_state.maxwell()->printDecomposition();

   const KineticSpeciesPtrVect& state_kinetic_species =
      m_state.kineticSpecies();
   int num_kinetic_species = static_cast<int>(state_kinetic_species.size());
   for (int s(0); s < num_kinetic_species; ++s) {
      state_kinetic_species[s]->printDecomposition();
   }
}

} // end namespace Loki
