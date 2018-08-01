/*************************************************************************
 *
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * Written by Jeffrey Banks banksj3@rpi.edu (Rensselaer Polytechnic Institute,
 * Amos Eaton 301, 110 8th St., Troy, NY 12180); Jeffrey Hittinger
 * hittinger1@llnl.gov, William Arrighi arrighi2@llnl.gov, Richard Berger
 * berger5@llnl.gov, Thomas Chapman chapman29@llnl.gov (LLNL, P.O Box 808,
 * Livermore, CA 94551); Stephan Brunner stephan.brunner@epfl.ch (Ecole
 * Polytechnique Federale de Lausanne, EPFL SB SPC-TH, PPB 312, Station 13,
 * CH-1015 Lausanne, Switzerland).
 * CODE-744849
 *
 * All rights reserved.
 *
 * This file is part of Loki.  For details, see.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 ************************************************************************/
#include "VMSystem.H"

#include <sstream>

#include "tbox/MathUtilities.H"
#include "LoadBalancer.H"
#include "TimerManager.H"
#include "ParallelUtility.h"
#include "RestartManager.H"
#include "RK4Integrator.H"
#include "RK6Integrator.H"

namespace Loki {

static const bool SKIP_DEEP_COPY = false;

VMSystem::VMSystem(
   ParmParse& a_pp,
   int a_num_probes,
   int a_spatial_solution_order,
   int a_temporal_solution_order)
   : m_cdim(CDIM),
     m_pdim(PDIM),
     m_cfg_domain(new ProblemDomain(m_cdim, a_spatial_solution_order, a_pp)),
     m_length_seq(8),
     m_spatial_solution_order(a_spatial_solution_order),
     m_temporal_solution_order(a_temporal_solution_order)
{
   // Create all the timers that will be active in a Vlasov-Maxwell system.
   TimerManager* timers(TimerManager::getManager());
   timers->createTimer("Vlasov");
   timers->createTimer("phys to phase");
   timers->createTimer("blowout");
   timers->createTimer("driver");
   timers->createTimer("contraction");
   timers->createTimer("reduction");
   timers->createTimer("reduction to phase");
   timers->createTimer("summation");
   timers->createTimer("Poisson");
   timers->createTimer("BC (Vlasov)");
   timers->createTimer("BC (Maxwell)");
   timers->createTimer("parallel ghost");
   timers->createTimer("collisions");
   timers->createTimer("krook");
   timers->createTimer("other");
   timers->createTimer("AddSolnData");
   timers->createTimer("CopySolnData");
   timers->createTimer("ZeroSolnData");
   timers->createTimer("Maxwell RHS");
   timers->createTimer("Tracking Particles");

   // Read all user input.
   parseParameters(a_pp);

   // Create the state for this system, essentially all the KineticSpecies and
   // a Maxwell object.
   createVMState(a_pp);

   // Now that all the loads are defined, load balance the system.
   loadBalance();

   // This is something of a hack.  We only want all the particles to exist on
   // the Maxwell object(s) on the Maxwell processor(s).  But we only know which
   // Maxwell object(s) are on Maxwell processor(s) after load balancing.  So we
   // read the particle data here.  As the species need to know if there are
   // tracking particles in the problem we tell them that here as well.
   KineticSpeciesPtrVect& state_kinetic_species = m_state.kineticSpecies();
   Maxwell& state_maxwell = *m_state.maxwell();
   int num_kinetic_species = static_cast<int>(state_kinetic_species.size());
   ParmParse ppp("maxwell");
   state_maxwell.readParticles(ppp);
   bool problem_has_tracking_particles =
      state_maxwell.numProblemTrackingParticles() > 0;
   if (problem_has_tracking_particles) {
      for (int s(0); s < num_kinetic_species; ++s) {
         state_kinetic_species[s]->problemHasParticles();
      }
   }

   // The system writes restart data so register it with the RestartManager
   // which will use the putToRestart/getFromRestart callbacks to get the
   // restart data written/read.
   RestartManager* restart_manager(RestartManager::getManager());
   restart_manager->registerRestart(this);

   // Allocate global arrays defined on configuration space for the net charge
   // and current densities and the charge and current densities for each
   // species.  In addition, if there are tracking particles we need an array
   // for the net external E field and the external E field from each species'
   // drivers.
   state_maxwell.newAuxVariable(m_net_charge_density, 1);
   m_charge_density.resize(num_kinetic_species);
   state_maxwell.newAuxVariable(m_net_x_current_density, 1);
   m_x_current_density.resize(num_kinetic_species);
   state_maxwell.newAuxVariable(m_net_y_current_density, 1);
   m_y_current_density.resize(num_kinetic_species);
   state_maxwell.newAuxVariable(m_net_z_current_density, 1);
   m_z_current_density.resize(num_kinetic_species);
   if (problem_has_tracking_particles) {
      state_maxwell.newAuxVariable(m_net_ext_efield, 2);
      m_ext_efield.resize(num_kinetic_species);
   }
   for (int s(0); s < num_kinetic_species; ++s) {
      state_maxwell.newAuxVariable(m_charge_density[s], 1);
      m_charge_density[s] = 0.0;
      state_maxwell.newAuxVariable(m_x_current_density[s], 1);
      m_x_current_density[s] = 0.0;
      state_maxwell.newAuxVariable(m_y_current_density[s], 1);
      m_y_current_density[s] = 0.0;
      state_maxwell.newAuxVariable(m_z_current_density[s], 1);
      m_z_current_density[s] = 0.0;
      if (problem_has_tracking_particles) {
         state_maxwell.newAuxVariable(m_ext_efield[s], 2);
         m_ext_efield[s] = 0.0;
      }
   }

   // Also need global arrays defined on configuration space for the species KE,
   // species KE flux on physical boundaries, and species KE flux on velocity
   // boundaries.
   state_maxwell.newAuxVariable(m_species_ke, 1);
   state_maxwell.newAuxVariable(m_species_phys_bdry_flux, 1);
   state_maxwell.newAuxVariable(m_species_vel_bdry_flux, 1);

   // NOTE:  Overture's plotStuff has a limit of 25 components per sequence.
   // The source of this limit is not clear but it looks like it may come from
   // the creation of the pull down menu of component names.  If there are
   // more than 25 components in a sequence plotStuff will core dump.  If more
   // than 25 components are needed it will be necessary to create multiple
   // sequences.
   // Compute the number of individual time histories and allocate initial space
   // for the time histories and time stamps.  This is not done efficiently.
   // The time history and time stamp arrays keep growing as the simulation
   // runs.  So we hold onto and repeatedly write old time history data.  These
   // arrays should be of size large enough to hold the time histories between
   // 2 plot cycles.  Then each plot file has that time intervals' worth of
   // time history data.  The post processor will need to change to accomodate
   // this change in handling time histories.
   m_num_seq =
      Maxwell::GLOBAL_TIME_HISTS +
      Maxwell::TIME_HISTS_PER_PROBE * a_num_probes +
      state_maxwell.numProblemTrackingParticles() * 4 +
      num_kinetic_species * 5;
   m_sequences.resize(m_length_seq, m_num_seq);
   m_time_seq.resize(m_length_seq);
   for (int i1(0); i1 < m_length_seq; ++i1) {
      m_time_seq(i1) = 0.0;
      for (int i2(0); i2 < m_num_seq; ++i2) {
         m_sequences(i1, i2) = 0.0;
      }
   }
}


VMSystem::~VMSystem()
{
}


void
VMSystem::initialize(
   bool a_is_from_restart,
   real a_time,
   real a_dt)
{
   // If this is an initial run (not from a restart file) then the
   // KineticSpecies' distribution functions must be initialized.
   Maxwell& state_maxwell = *m_state.maxwell();
   state_maxwell.initialize(m_species_names);
   if (!a_is_from_restart) {
      // Initialize each species distributed to this processor.
      KSPV::Iterator it_end(m_state.kineticSpecies().end_locals());
      for (KSPV::Iterator it(m_state.kineticSpecies().begin_locals());
           it != it_end; ++it) {
         (*it)->initialize(0.0);
      }
   }

   // Compute the static factorization for Poisson's equation.  This can
   // probably be removed as we don't solve Poisson's equation.
   state_maxwell.computeFactorization();

   // Clone all the species and Maxwell from m_state into m_state_old.
   KineticSpeciesPtrVect& state_kinetic_species = m_state.kineticSpecies();
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
   if (!a_is_from_restart) {
      state_maxwell.initializeVZ();
      state_maxwell.initializeEM();
   }
   for (int s(0); s < num_kinetic_species; ++s) {
      state_kinetic_species[s]->currentDensity(state_maxwell,
         m_x_current_density[s],
         m_y_current_density[s],
         m_z_current_density[s]);
   }

   // Only the Maxwell processor(s) are needed to apply any external currents.
   if (state_maxwell.isMaxwellProcessor()) {
      m_net_charge_density = 0;
      m_net_x_current_density = 0;
      m_net_y_current_density = 0;
      m_net_z_current_density = 0;
      for (int s(0); s < num_kinetic_species; ++s) {
         m_net_charge_density += m_charge_density[s];
         m_net_x_current_density += m_x_current_density[s];
         m_net_y_current_density += m_y_current_density[s];
         m_net_z_current_density += m_z_current_density[s];
      }
      state_maxwell.applyExternalCurrent(m_net_x_current_density,
         m_net_y_current_density,
         m_net_z_current_density,
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
      rhs_kspv[s] = prototype_kspv[s]->clone(tbox::IntVector::Zero(m_pdim),
         SKIP_DEEP_COPY);
   }
   a_rhs.maxwell() =
      a_prototype.maxwell()->clone(tbox::IntVector::Zero(m_cdim),
         SKIP_DEEP_COPY);
}


void
VMSystem::defineSolnData(
   VMState&       a_soln,
   const VMState& a_prototype)
{
   // We need to define all species on all processors of a_soln.
   KineticSpeciesPtrVect& soln_kspv = a_soln.kineticSpecies();
   const KineticSpeciesPtrVect& prototype_kspv = a_prototype.kineticSpecies();
   int num_kinetic_species = static_cast<int>(prototype_kspv.size());
   soln_kspv.resize(num_kinetic_species);
   for (int s(0); s < num_kinetic_species; ++s) {
      soln_kspv[s] = prototype_kspv[s]->clone(tbox::IntVector::Zero(m_pdim),
         SKIP_DEEP_COPY);
   }
   a_soln.maxwell() =
      a_prototype.maxwell()->clone(tbox::IntVector::Zero(m_cdim),
         SKIP_DEEP_COPY);
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
   if (a_dst_maxwell->isMaxwellProcessor()) {
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
   if (a_soln_maxwell->isMaxwellProcessor()) {
      a_soln_maxwell->zeroData();
   }
   timers->stopTimer("ZeroSolnData");
}


void
VMSystem::addSolnData(
   VMState&       a_soln,
   const VMState& a_increment,
   real           a_scale,
   bool           a_final_rk)
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
   if (a_soln_maxwell->isMaxwellProcessor()) {
      a_soln_maxwell->addData(*a_increment.maxwell(), a_scale, a_final_rk);
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
   real     a_time,
   real     a_dt,
   int      a_stage)
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
   for (int s(0); s < num_kinetic_species; ++s) {
      state_kspv[s]->currentDensity(state_maxwell,
         m_x_current_density[s],
         m_y_current_density[s],
         m_z_current_density[s]);
   }

   // The external current application by the Maxwell processor(s) and the
   // advection computation by the Vlasov processors are done in parallel as
   // these computations are independent and the processor sets do not overlap.
   // The computation of the acceleration that follows requires parallel
   // communication between all of the processors.  If load imbalance exists
   // between the Maxwell and Vlasov processors the acceleration computation
   // will appear to be a bottleneck as its parallel communication can not
   // proceed until all processors have reached that point.
   if (state_maxwell.isMaxwellProcessor()) {
      m_net_x_current_density = 0.0;
      m_net_y_current_density = 0.0;
      m_net_z_current_density = 0.0;
      for (int s(0); s < num_kinetic_species; ++s) {
         m_net_x_current_density += m_x_current_density[s];
         m_net_y_current_density += m_y_current_density[s];
         m_net_z_current_density += m_z_current_density[s];
      }
      state_maxwell.applyExternalCurrent(m_net_x_current_density,
         m_net_y_current_density,
         m_net_z_current_density,
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
         a_stage);
   }

   // Now the Vlasov processors will fill the acceleration ghost cells
   // (boundaries of velocity domain) of each species distributed on them and
   // compute the Vlasov equation RHS.  The Maxwell processor(s) will compute
   // the Maxwell equations' RHS in parallel.
   a_rhs_it.reset();
   for (KSPV::Iterator a_state_it(state_kspv.begin_locals());
        a_state_it != a_state_it_end; ++a_state_it, ++a_rhs_it) {
      (*a_state_it)->fillAccelerationGhostCells();
      if (m_do_new_algorithm) {
         (*a_state_it)->evalAccelerationDerivatives(*(*a_rhs_it));
      }
      else {
         (*a_state_it)->evalAccelerationFluxes();
      }
      (*a_state_it)->completeRHS(*(*a_rhs_it), a_time, m_do_new_algorithm);
   }
   if (state_maxwell.isMaxwellProcessor()) {
      // If there are tracking particles then compute the net externally applied
      // E field which is needed to update them.
     if (state_maxwell.numProblemTrackingParticles() > 0) {
         m_net_ext_efield = 0.0;
         for (int s(0); s < num_kinetic_species; ++s) {
            m_net_ext_efield += m_ext_efield[s];
         }
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
   // The hack found in VPSystem is not needed as the issue with the fields is
   // not present.
   NULL_USE(a_soln);
}


real
VMSystem::stableDt()
{
   real dt_stable(tbox::MathUtilities<real>::getMax());
   KineticSpeciesPtrVect& state_kinetic_species = m_state.kineticSpecies();
   // Compute the time step contraint for each species distributed to this
   // processor and take the min of them.
   KSPV::Iterator it_end(state_kinetic_species.end_locals());
   for (KSPV::Iterator it(state_kinetic_species.begin_locals());
        it != it_end; ++it) {
      real dt_species = (*it)->computeDt();
      dt_stable = std::min(dt_stable, dt_species);
   }
   real dt_maxwell = m_state.maxwell()->computeDt();
   dt_stable = std::min(dt_stable, dt_maxwell);
   dt_stable = ParallelUtility::getMinValue(dt_stable);
   
   return dt_stable;
}


real
VMSystem::advance(
   real a_cur_time,
   real a_dt)
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


int
VMSystem::getNumFrameSeries()
{
   return Maxwell::NUM_FRAME_SERIES;
}


void
VMSystem::plot(
   real a_time,
   real a_dt,
   const RealArray& a_probes,
   int a_num_probes,
   int a_saved_seq,
   int& a_saved_save,
   Ogshow& a_show)
{
   tbox::Pointer<Maxwell>& state_maxwell = m_state.maxwell();
   KineticSpeciesPtrVect& state_kinetic_species = m_state.kineticSpecies();
   // Compute the KE velocity boundary flux which will be a plotted quantity.
   // As there is parallel communication in this computation all processors
   // must participate for all species.
   int num_kinetic_species = static_cast<int>(state_kinetic_species.size());
   for (int s(0); s < num_kinetic_species; ++s) {
      state_kinetic_species[s]->computeKEVelBdyFlux(*state_maxwell,
         m_species_vel_bdry_flux);
   }

   // The Maxwell object actually writes the plot data.
   state_maxwell->plot(a_time,
      a_dt,
      m_sequences,
      m_species_names,
      m_time_seq,
      a_probes,
      a_num_probes,
      m_num_seq,
      a_saved_seq,
      a_saved_save,
      a_show);
}


void
VMSystem::accumulateSequences(
   real a_time,
   const RealArray& a_probes,
   int a_num_probes,
   int& a_saved_seq)
{
   // resize saved sequence if necessary
   if (a_saved_seq >= m_length_seq) {
      m_length_seq *= 2;
      m_time_seq.resize(m_length_seq);
      m_sequences.resize(m_length_seq, m_num_seq);

      for (int i1(a_saved_seq); i1 < m_length_seq; ++i1) {
         m_time_seq(i1) = 0.0;
         for (int i2(0); i2 < m_num_seq; ++i2) {
            m_sequences(i1, i2) = 0.0;
         }
      }
   }
   m_time_seq(a_saved_seq) = a_time;

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
   bool is_maxwell_proc = state_maxwell->isMaxwellProcessor();
   KineticSpeciesPtrVect& state_kinetic_species = m_state.kineticSpecies();
   int num_kinetic_species = static_cast<int>(state_kinetic_species.size());
   for (int s(0); s < num_kinetic_species; ++s) {
      state_kinetic_species[s]->accumulateSequences(is_maxwell_proc,
         m_sequences,
         a_saved_seq,
         seq_idx,
         m_species_ke,
         m_species_phys_bdry_flux);
   }

   // We've saved another time history.
   ++a_saved_seq;
}


void
VMSystem::getFromRestart(
   const HDF_DataBase& a_db)
{
   HDF_DataBase sub_db;
   a_db.locate(sub_db, "species_list");

   // Read all the species names from the restart file and ensure that they
   // match the names in the input deck.  If they don't then the restart file
   // and input deck are incompatible.
   int species_list_size;
   a_db.get(species_list_size, "species_list_size");
   int num_species_names = static_cast<int>(m_species_names.size());
   if (species_list_size != num_species_names) {
      OV_ABORT("Kinetic species names list length mismatch!");
   }

   for (int s(0); s < species_list_size; ++s) {
      aString tmp;
      std::stringstream tag;
      tag << "species." << s + 1;
      sub_db.get(tmp, (aString)tag.str());

      int i;
      for (i = 0; i < num_species_names; ++i) {
         if (m_species_names[i] == tmp) {  // TODO: Explicitly cast m_species_names to aString, or allow auto casting...?
            break;
         }
      }
      if (i == num_species_names) {
         OV_ABORT("Kinetic species names from restart don't match names from input file!");
      }
   }
}

  
void
VMSystem::putToRestart(
   HDF_DataBase& a_db,
   real a_time)
{
   NULL_USE(a_time);

   // Write the species names to the restart file.
   HDF_DataBase sub_db;
   a_db.create(sub_db, "species_list", "directory");

   int species_list_size(static_cast<int>(m_species_names.size()));
   a_db.put(species_list_size, "species_list_size");
   for (int s(0); s < species_list_size; ++s) {
      std::stringstream tag;
      tag << "species." << s + 1;
      sub_db.put(m_species_names[s], tag.str().c_str());
   }
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
      size = std::max(size, state_kinetic_species[s]->numberOfCells());
   }
   return size;
}

void
VMSystem::updateGhosts()
{
   KineticSpeciesPtrVect& state_kinetic_species = m_state.kineticSpecies();
   // Update the advection ghost cells for each species distributed to this
   // processor.
   KSPV::Iterator species_it_end(state_kinetic_species.end_locals());
   for (KSPV::Iterator species_it(state_kinetic_species.begin_locals());
        species_it != species_it_end; ++species_it) {
      (*species_it)->fillAdvectionGhostCells();
   }
}

////// PRIVATE FUNCTIONS ////////////////////////////////////
void
VMSystem::createVMState(
   ParmParse& a_pp)
{
   int number_of_species(1);
   a_pp.query("number_of_species", number_of_species);

   aString test_str = "false";
   a_pp.query("use_new_bcs", test_str);
   bool use_new_bcs = test_str.matches("true") ? true : false;

   // Make the Maxwell object's sub-database and construct it.
   ParmParse ppp("maxwell");
   m_state.maxwell() =
      new Maxwell(ppp, m_cfg_domain, number_of_species, m_spatial_solution_order);

   // For each species, make its sub-database and construct it.
   KineticSpeciesPtrVect& state_kinetic_species = m_state.kineticSpecies();
   state_kinetic_species.resize(number_of_species);

   for (int s(0); s < number_of_species; ++s) {

      char buffer[100];
      sprintf(buffer, "kinetic_species.%i", s+1);
      ParmParse ppspecies(buffer);

      state_kinetic_species[s] = new KineticSpecies(m_cfg_domain,
         ppspecies,
         s+1,
         number_of_species,
         m_spatial_solution_order,
         m_temporal_solution_order,
         use_new_bcs,
         true);

      // Create list of species names
      m_species_names.push_back(state_kinetic_species[s]->name().c_str());
   }
}


void
VMSystem::loadBalance()
{
   LoadBalancer load_balancer;

   // Every species, even those not distributed to this processor, is a load to
   // be balanced as is the Maxwell object.
   std::vector<Load*> load_vector;
   KineticSpeciesPtrVect& state_kinetic_species = m_state.kineticSpecies();
   int num_kinetic_species = static_cast<int>(state_kinetic_species.size());
   for (int s(0); s < num_kinetic_species; ++s) {
      load_vector.push_back(state_kinetic_species[s].getPointer());
   }
   load_vector.push_back(m_state.maxwell().getPointer());

   // Balance everything.
   load_balancer.balance(load_vector);

   printDecomposition();
}


void
VMSystem::printDecomposition() const
{
   // Print the decompositions of the Maxwell and of each species.
   printF("\n#*#*# VMSystem: Parallel Decomposition #*#*#\n");
   m_state.maxwell()->printDecomposition();

   const KineticSpeciesPtrVect& state_kinetic_species =
      m_state.kineticSpecies();
   int num_kinetic_species = static_cast<int>(state_kinetic_species.size());
   for (int s(0); s < num_kinetic_species; ++s) {
      state_kinetic_species[s]->printDecomposition();
   }
}

} // end namespace Loki
