/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "Maxwell.H"

#include "CurrentDriverFactory.H"
#include "EMICFactory.H"
#include "VELICFactory.H"
#include "MaxwellF.H"
#include "RestartManager.H"
#include "Loki_Defines.H"
#include "Loki_Utilities.H"
#include "Interpolator4.H"
#include "Interpolator6.H"

#include <sstream>

namespace Loki {

const int Maxwell::s_TAG_BATON = 4379;
int Maxwell::TIME_HISTS_PER_PROBE;
const int Maxwell::GLOBAL_TIME_HISTS = 12;
int Maxwell::NUM_COLL_TIME_HISTS;

Maxwell::Maxwell(LokiInputParser& a_pp,
   tbox::Pointer<ProblemDomain> a_domain,
   int a_num_kinetic_species,
   int a_solution_order,
   int a_plot_times_per_file,
   bool a_plot_ke_vel_bdy_flux,
   bool a_coll_diag_on,
   double a_bz_const)
   : EMSolverBase(a_domain,
                  a_num_kinetic_species,
                  a_solution_order,
                  NUM_EM_VARS,
                  a_plot_times_per_file,
                  a_plot_ke_vel_bdy_flux,
                  a_coll_diag_on,
                  a_bz_const),
     m_vz(a_num_kinetic_species),
     m_avWeak(0.0),
     m_avStrong(1.0),
     m_current_drivers(0),
     m_num_current_drivers(0),
     m_em_ics(0),
     m_vel_ics(0),
     m_ex_interp(0),
     m_ey_interp(0),
     m_bz_interp(0),
     m_jx_interp(0),
     m_jy_interp(0),
     m_supergrid_lo(m_domain->lower()),
     m_supergrid_hi(m_domain->upper())
{
   // 6 components of EM + vz for each species.
   TIME_HISTS_PER_PROBE = 6 + a_num_kinetic_species;

   // 4 moments (ke, ent, px, and py) plus change in these 4 moments due to
   //  collision operator.
   NUM_COLL_TIME_HISTS = 8*a_num_kinetic_species;

   // Read all input related to this object.
   parseParameters(a_pp);

   // The following code could all probably be placed into parseParameters but
   // isn't because it involves parsing separate sub-databases of this Maxwell.

   // Get the sub-database for any current drivers and create them.
   ostringstream input_string;
   for (int i = 0; i < m_num_current_drivers; ++i) {
      input_string << "maxwell.current_driver." << i+1;
      LokiInputParser cdf_pp(input_string.str().c_str());
      input_string.str("");
      m_current_drivers[i] =
         CurrentDriverFactory::create(cdf_pp, *m_domain, i+1);
   }

   // Get the sub-database for each EM initial condition and create them.  The
   // EM initial conditions are optional but if supplied there must be 2 and
   // only 2--one for E and one for B.
   if (a_pp.contains("em_ic.1.name")) {
      if (a_pp.contains("em_ic.3.name")) {
         LOKI_ABORT("There may be only 2 electromagnetic initial conditions.");
      }
      m_em_ics.resize(2, tbox::Pointer<EMICInterface>(0));
      int num_E_initializers = 0;
      int num_B_initializers = 0;
      for (int i = 0; i < 2; ++i) {
         input_string << "maxwell.em_ic." << i+1;
         LokiInputParser emic_pp(input_string.str().c_str());
         input_string.str("");
         m_em_ics[i] = EMICFactory::create(emic_pp);
         if (m_em_ics[i]->initializesE()) {
            ++num_E_initializers;
         }
         else {
            ++num_B_initializers;
         }
      }
      if (num_E_initializers != 1 || num_B_initializers != 1) {
         LOKI_ABORT("There must be 1 initial condition for each of E and B");
      }
      if (m_em_ics[0]->numWaves() != m_em_ics[1]->numWaves()) {
         LOKI_ABORT("Unequal number of E and B wave initializers.");
      }
   }

   // Get the sub-database for each species transverse drift velocity initial
   // condition and create them.  The transverse drift velocity initial
   // condition is optional but if supplied there must be one for each species.
   if (a_pp.contains("vel_ic.1.name")) {
      input_string << "vel_ic." << m_num_kinetic_species+1 << ".name";
      if (a_pp.contains(input_string.str().c_str())) {
         ostringstream error_msg;
         error_msg << "There must be only " << m_num_kinetic_species
                   << " transverse drift velocity initializers.";
         LOKI_ABORT(error_msg.str());
      }
      input_string.str("");
      m_vel_ics.resize(m_num_kinetic_species, tbox::Pointer<VELICInterface>(0));
      for (int i = 0; i < m_num_kinetic_species; ++i) {
         input_string << "maxwell.vel_ic." << i+1;
         LokiInputParser vel_ic_pp(input_string.str().c_str());
         input_string.str("");
         m_vel_ics[i] = VELICFactory::create(vel_ic_pp);
         if (i > 0 && (m_vel_ics[0]->numWaves() != m_vel_ics[i]->numWaves())) {
            LOKI_ABORT("Transverse drift velocity initializers have different number of waves.");
         }
      }
   }
   if (m_em_ics.size() > 0 && m_vel_ics.size() > 0 &&
       (m_em_ics[0]->numWaves() != m_vel_ics[0]->numWaves())) {
      LOKI_ABORT("Unequal number of EM and transverse drift velocity wave initializers.");
   }

   // The Maxwell object writes restart data so register this object with the
   // RestartManager which will use the putToRestart/getFromRestart callbacks to
   // get the restart data written/read.
   RestartManager* restart_manager(RestartManager::getManager());
   restart_manager->registerRestart(this);
}


Maxwell::Maxwell(
   const Maxwell& a_other)
   : EMSolverBase(a_other),
     m_vz(a_other.m_num_kinetic_species),
     m_ex_interp(0),
     m_ey_interp(0),
     m_bz_interp(0),
     m_jx_interp(0),
     m_jy_interp(0)
{
   if (m_dim != a_other.m_dim) {
      LOKI_ABORT("Attempt to copy incongruent Maxwells");
   }

   // If these 2 Maxwells are different, then copy their internals.
   if (&a_other != this) {
      for (int i = 0; i < m_num_kinetic_species; ++i) {
         m_vz[i] = a_other.m_vz[i];
      }

      m_avWeak = a_other.m_avWeak;
      m_avStrong = a_other.m_avStrong;
      m_current_drivers = a_other.m_current_drivers;
      m_num_current_drivers = a_other.m_num_current_drivers;
      m_antenna_source = a_other.m_antenna_source;
      m_em_ics = a_other.m_em_ics;
      m_vel_ics = a_other.m_vel_ics;
      m_supergrid_lo = a_other.m_supergrid_lo;
      m_supergrid_hi = a_other.m_supergrid_hi;
      buildInterpolators();
   }
}


Maxwell::~Maxwell()
{
   if (m_ex_interp) {
      for (int i = 0; i < m_solution_order; ++i) {
         delete [] m_ex_interp[i];
      }
      delete [] m_ex_interp;
   }
   if (m_ey_interp) {
      for (int i = 0; i < m_solution_order; ++i) {
         delete [] m_ey_interp[i];
      }
      delete [] m_ey_interp;
   }
   if (m_bz_interp) {
      for (int i = 0; i < m_solution_order; ++i) {
         delete [] m_bz_interp[i];
      }
      delete [] m_bz_interp;
   }
   if (m_jx_interp) {
      for (int i = 0; i < m_solution_order; ++i) {
         delete [] m_jx_interp[i];
      }
      delete [] m_jx_interp;
   }
   if (m_jy_interp) {
      for (int i = 0; i < m_solution_order; ++i) {
         delete [] m_jy_interp[i];
      }
      delete [] m_jy_interp;
   }
}


void
Maxwell::createPartition(
   int a_proc_lo,
   int a_proc_hi,
   const MPI_Comm& a_comm)
{
   // Part common to all EMSolvers.
   createPartitionCommon(false, a_proc_lo, a_proc_hi, a_comm);

   // Repeat the above process for the global array of each species' transverse
   // drift velocity.
   deque<bool> is_periodic(m_dim);
   vector<int> num_cells(m_dim);
   for (int dim = 0; dim < m_dim; ++dim) {
      is_periodic[dim] = m_domain->isPeriodic(dim);
      num_cells[dim] = m_domain->numberOfCells(dim);
   }
   for (int i = 0; i < m_num_kinetic_species; ++i) {
      m_vz[i].partition(m_dim,
         m_dim,
         m_proc_lo,
         m_proc_hi,
         m_n_ghosts,
         is_periodic,
         num_cells);
   }

   // Figure out the local and interior boxes.  Recall that all Maxwells exist
   // on all processors that only 1 (or possibly a few) processors have a
   // Maxwell that does work.
   if (isInRange(Loki_Utilities::s_my_id)) {
      for (int dir(0); dir < m_dim; ++dir) {
         if (num_cells[dir] < m_solution_order+1) {
            LOKI_ABORT("Too few interior points in decomposition!");
         }
      }
   }

   // If there are any current drivers (antennae), then allocate the antenna
   // source array.
   if (m_num_current_drivers > 0) {
      num_cells.push_back(m_num_em_vars);
      m_antenna_source.partition(m_dim+1,
         m_dim,
         m_proc_lo,
         m_proc_hi,
         m_n_ghosts,
         is_periodic,
         num_cells);
   }
}


bool
Maxwell::conformsTo(
   const Maxwell& a_other,
   bool a_include_ghost_cells) const
{
   // If the 2 Maxwells are different, check that they are defined on the same
   // piece of space and are partitioned the same.
   if (m_domain != a_other.m_domain ||
       m_proc_lo != a_other.m_proc_lo ||
       m_proc_hi != a_other.m_proc_hi ||
       m_partition_defined != a_other.m_partition_defined ||
       numTrackingParticles() != a_other.numTrackingParticles() ||
       numNoiseSourceParticles() != a_other.numNoiseSourceParticles()) {
      return false;
   }
   else if (!a_include_ghost_cells) {
      return (m_domain->box() == a_other.m_domain->box());
   }
   else {
      return (m_domain->box() == a_other.m_domain->box()) &&
         (m_n_ghosts == a_other.m_n_ghosts);
   }
}


void
Maxwell::addData(
   const Maxwell& a_increment_maxwell,
   double a_factor,
   bool a_sum_reduce_inc)
{
   // Check that we're copying between similar Maxwells.  Add a_factor times
   // a_increment_maxwell's EM fields, species transverse drift velocities to
   // this Maxwell's EM fields and species transverse drift velocities.
   if (conformsTo(a_increment_maxwell, false)) {
      // Handle the EM fields and the transverse drift velocities of each
      // species.
      FORT_XPBY_2D(*m_em_vars.getData(),
         *a_increment_maxwell.m_em_vars.getData(),
         a_factor,
         BOX2D_TO_FORT(m_em_vars.dataBox()),
         BOX2D_TO_FORT(m_em_vars.interiorBox()),
         NUM_EM_VARS);
      for (int i = 0; i < m_num_kinetic_species; ++i) {
         FORT_XPBY_2D(*m_vz[i].getData(),
            *a_increment_maxwell.m_vz[i].getData(),
            a_factor,
            BOX2D_TO_FORT(m_vz[i].dataBox()),
            BOX2D_TO_FORT(m_vz[i].interiorBox()),
            1);
      }

      // Add the tracking and noise source particle data.
      addParticleData(a_increment_maxwell.m_tracking_particles,
         a_increment_maxwell.m_noise_source_particles,
         a_factor,
         a_sum_reduce_inc);
   }
}


void
Maxwell::copySolnData(
   const Maxwell& a_rhs)
{
   if (m_dim != a_rhs.m_dim) {
      LOKI_ABORT("Attemtpt to copy incongruent Maxwells!");
   }

   // If the 2 Maxwells are different copy the EM fields, each species
   // transverse drift velocity, and the particles' x, y, vx, and vy.
   if (&a_rhs != this) {
      m_em_vars = a_rhs.m_em_vars;
      for (int i = 0; i < m_num_kinetic_species; ++i) {
         m_vz[i] = a_rhs.m_vz[i];
      }

      copyParticleData(a_rhs.m_tracking_particles,
                       a_rhs.m_noise_source_particles);
   }
}


void
Maxwell::initialize(
   int a_num_species,
   bool a_plot_ke_vel_bdy_flux,
   bool a_from_restart,
   bool a_coll_diag_on)
{
   // Do the initialization common to all EMSolvers.
  initializeCommon(a_num_species, a_plot_ke_vel_bdy_flux, a_coll_diag_on);

   // Set the initial conditions for the EM fields and the transverse drift
   // velocities at the beginning of the problem.
   if (!a_from_restart && isEMSolverProcessor()) {
      initializeEM();
      initializeVZ();
   }

   // Now we know enough to build the appropriate field interpolators if needed.
   buildInterpolators();
}


void
Maxwell::computeAntennaSource(
   double a_time,
   const ParallelArray a_omega_eff2)
{
   // Compute and accumulate all the antenna source terms.
   if (m_num_current_drivers > 0) {
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("current driver");
      m_antenna_source = 0.0;
      for (int i = 0; i < m_num_current_drivers; ++i) {
         m_current_drivers[i]->evaluate(m_antenna_source,
            a_omega_eff2,
            *m_domain,
            a_time);
      }
      timers->stopTimer("current driver");
   }
}


void
Maxwell::applyNoiseParticleCurrent(
   ParallelArray& a_Jx,
   ParallelArray& a_Jy,
   double a_time)
{
   int num_noise_source_particles = numNoiseSourceParticles();
   TimerManager* timers(TimerManager::getManager());
   if (num_noise_source_particles > 0) {
      timers->startTimer("Noisy Particles");
   }
   else {
      return;
   }

   // If there are noise particles, add its contribution of each to the x and y
   // net current densities.
   double x_lo = m_domain->lower(0);
   double x_hi = m_domain->upper(0);
   double dx = m_domain->dx(0);
   double nx = m_domain->numberOfCells(0);
   double y_lo = m_domain->lower(1);
   double y_hi = m_domain->upper(1);
   double dy = m_domain->dx(1);
   double ny = m_domain->numberOfCells(1);
   const ParallelArray::Box& interior_box = m_em_vars.interiorBox();
   int ib_x_lo = interior_box.lower(0);
   int ib_x_hi = interior_box.upper(0);
   int ib_y_lo = interior_box.lower(1);
   int ib_y_hi = interior_box.upper(1);
   const ParallelArray::Box& data_box = m_em_vars.dataBox();
   int db_x_lo = data_box.lower(0);
   int db_x_hi = data_box.upper(0);
   int db_y_lo = data_box.lower(1);
   int db_y_hi = data_box.upper(1);
   for (int i = 0; i < num_noise_source_particles; ++i) {
      // Check if this particle makes a contribution to the current density on
      // this Maxwell processor.
      const Particle& this_particle = m_noise_source_particles[i];
      double x_loc = this_particle.x();
      double y_loc = this_particle.y();
      int x_idx = (x_loc-x_lo)/dx;
      int y_idx = (y_loc-y_lo)/dy;
      double t_ramp = this_particle.startingTime();
      double envel;
      if (a_time < t_ramp) {
         envel = 0.5+0.5*tanh(8.0*a_time/t_ramp-4.0);
      }
      else {
         envel = 1.0;
      }

      // See if this particle is in the interior or ghost layer of this Maxwell
      // processor.  If so, then it contributes to the current density on the
      // processor.
      bool particle_contributes = x_idx >= db_x_lo && x_idx <= db_x_hi &&
         y_idx >= db_y_lo && y_idx <= db_y_hi;

      // If the particle is in the interior or the ghost layer of this domain
      // get its contribution.
      if (particle_contributes) {
         // Set the evaluation point to be the location of this particle.
         m_j_interpolator->setEvaluationPoint(x_loc, y_loc);

         // Now compute the current densities at the interpolating points that
         // is required to give the current densities of this noise source
         // particle.
         double noise_charge =
            envel*this_particle.charge()*this_particle.noiseSourceWeight();
         double Jx_noise = noise_charge*this_particle.vx();
         double Jy_noise = noise_charge*this_particle.vy();
         m_j_interpolator->deterpolate(Jx_noise, m_jx_interp);
         m_j_interpolator->deterpolate(Jy_noise, m_jy_interp);

         // Add the particle's contribution to the net x and y current density.
         for (int xidx = 0; xidx < m_solution_order; ++xidx) {
            int xGridIdx = m_j_interpolator->xInterpolatingGridIndex(xidx);
            // If this part of the stencil is in the ghost layer then it does
            // not contribute to this processor's current density.
            if (xGridIdx < ib_x_lo || xGridIdx > ib_x_hi) {
               continue;
            }
            for (int yidx = 0; yidx < m_solution_order; ++yidx) {
               int yGridIdx = m_j_interpolator->yInterpolatingGridIndex(yidx);
               // If thid part of the stencil is in the ghost layer then it does
               // not contribute to this processor's current density.
               if (yGridIdx < ib_y_lo || yGridIdx > ib_y_hi) {
                  continue;
               }
               a_Jx(xGridIdx, yGridIdx) += m_jx_interp[xidx][yidx];
               a_Jy(xGridIdx, yGridIdx) += m_jy_interp[xidx][yidx];
            }
         }
      }

      // If part of the particle's stencil is outside the domain then the
      // periodic image of that part of the stencil may contribute to the
      // current density on this processor.
      bool stencil_outside_domain =
         x_idx >= nx-(m_n_ghosts-1) || x_idx <= m_n_ghosts-1 ||
         y_idx >= ny-(m_n_ghosts-1) || y_idx <= m_n_ghosts-1;
      if (stencil_outside_domain) {
         // Set the evaluation point to be the location of this particle.
         m_j_interpolator->setEvaluationPoint(x_loc, y_loc);

         // Now compute the current densities at the interpolating points that
         // is required to give the current densities of this noise source
         // particle.
         double noise_charge =
            envel*this_particle.charge()*this_particle.noiseSourceWeight();
         double Jx_noise = noise_charge*this_particle.vx();
         double Jy_noise = noise_charge*this_particle.vy();
         m_j_interpolator->deterpolate(Jx_noise, m_jx_interp);
         m_j_interpolator->deterpolate(Jy_noise, m_jy_interp);

         // Add the particle's contribution to the net x and y current density.
         for (int xidx = 0; xidx < m_solution_order; ++xidx) {
            int xGridIdx = m_j_interpolator->xInterpolatingGridIndex(xidx);
            // If the x grid index of this stencil point is outside the domain
            // then shift it to it's periodic image in x.
            int x_shift = 0;
            if (xGridIdx < 0) {
               x_shift = nx;
            }
            else if (xGridIdx >= nx) {
               x_shift = -nx;
            }
            xGridIdx += x_shift;
            for (int yidx = 0; yidx < m_solution_order; ++yidx) {
               int yGridIdx = m_j_interpolator->yInterpolatingGridIndex(yidx);
               // If the y grid index of this stencil point is outside the
               // domain then shift it to it's periodic image in y.
               int y_shift = 0;
               if (yGridIdx < 0) {
                  y_shift = ny;
               }
               else if (yGridIdx >= ny) {
                  y_shift = -ny;
               }
               yGridIdx += y_shift;
               // If this stencil point is not outside the domain then ignore
               // it as we've already dealt with it.
               if (x_shift == 0 && y_shift == 0) {
                  continue;
               }
               // If this stencil point's periodic image is not in the interior
               // of this domain then it does not contribute to this processor's
               // current density.
               if (xGridIdx < ib_x_lo || xGridIdx > ib_x_hi ||
                   yGridIdx < ib_y_lo || yGridIdx > ib_y_hi) {
                  continue;
               }
               a_Jx(xGridIdx, yGridIdx) += m_jx_interp[xidx][yidx];
               a_Jy(xGridIdx, yGridIdx) += m_jy_interp[xidx][yidx];
            }
         }
      }
   }
   timers->stopTimer("Noisy Particles");
}


void
Maxwell::evalRHS(
   Maxwell& a_rhs,
   const KineticSpeciesPtrVect& a_kspv,
   const ParallelArray& a_Jx,
   const ParallelArray& a_Jy,
   const ParallelArray& a_Jz,
   const ParallelArray& a_net_ext_efield,
   double a_time)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("Maxwell RHS");

   // Evaluate the RHS of Maxwell's equations exclusive of any antenna source.
   FORT_MAXWELL_EVAL_RHS(BOX2D_TO_FORT(m_em_vars.dataBox()),
      BOX2D_TO_FORT(m_em_vars.interiorBox()),
      PROBLEMDOMAIN_TO_FORT((*m_domain)),
      Simulation::s_LIGHT_SPEED,
      m_avWeak,
      m_avStrong,
      m_solution_order,
      m_supergrid_lo[0],
      m_supergrid_hi[0],
      *(m_em_vars.getData()),
      *a_Jx.getData(),
      *a_Jy.getData(),
      *a_Jz.getData(),
      *(a_rhs.m_em_vars.getData()));

   // If there are antennae, then fold the antenna sources into the RHS of
   // Maxwell's equations.
   if (m_num_current_drivers > 0) {
      FORT_MAXWELL_ADD_ANTENNA_SOURCE(BOX2D_TO_FORT(m_em_vars.dataBox()),
         BOX2D_TO_FORT(m_em_vars.interiorBox()),
         PROBLEMDOMAIN_TO_FORT((*m_domain)),
         *(m_antenna_source.getData()),
         *(a_rhs.m_em_vars.getData()));
   }

   // For each species, evaluate the RHS of the transverse drift velocity
   // equations, dvz/dt.
   for (int i = 0; i < m_num_kinetic_species; ++i) {
      const KineticSpeciesPtr& ksp = a_kspv[i];
      FORT_MAXWELL_EVAL_VZ_RHS(BOX2D_TO_FORT(m_em_vars.dataBox()),
         BOX2D_TO_FORT(m_em_vars.interiorBox()),
         ksp->charge()/ksp->mass(),
         *(m_em_vars.getData()),
         *(a_rhs.m_vz[i].getData()));
   }

   if (numTrackingParticles() > 0) {
      timers->startTimer("Tracking Particles");
      evalRHSParticles(false, a_rhs, a_net_ext_efield, a_time);
      timers->stopTimer("Tracking Particles");
   }

   if (numNoiseSourceParticles() > 0) {
      timers->startTimer("Noisy Particles");
      evalRHSParticles(true, a_rhs, a_net_ext_efield, a_time);
      timers->stopTimer("Noisy Particles");
   }
   timers->stopTimer("Maxwell RHS");
}

void
Maxwell::plotCollision(
   double                         a_time,
   double                         a_dt,
   KineticSpeciesPtrVect const&   a_kinetic_species,
   const vector<vector<double> >& a_sequences,
   const vector<string>&          a_species_names)
{
   if (static_cast<int>(a_kinetic_species.size()) != m_num_kinetic_species) {
      LOKI_ABORT("Bug in Maxwell::plotCollision");
   }

   // TODO: implement this like Poisson
}


void
Maxwell::plot(
   double                         a_time,
   double                         a_dt,
   const vector<vector<double> >& a_sequences,
   const vector<string>&          a_species_names,
   const vector<double>&          a_time_seq,
   const vector<vector<double> >& a_probes, 
   int                            a_num_probes,
   int                            a_saved_seq,
   int&                           a_saved_save,
   bool                           a_plot_ke_vel_bdy_flux,
   string&                        a_time_hist_file_name)
{
   vector<string> frame_names;
   buildPlotNames(a_plot_ke_vel_bdy_flux, a_species_names, frame_names);
   vector<string> time_hist_names;
   buildTimeHistoryNames(a_num_probes,
                         numProblemTrackingParticles(),
                         a_species_names,
                         time_hist_names);
   bool is_lo_em_proc = Loki_Utilities::s_my_id == m_proc_lo;
   bool multiple_em_procs = m_proc_lo != m_proc_hi;

   if (multiple_em_procs && !is_lo_em_proc) {
      MPI_Status status;
      int info;
      int prev_writer = Loki_Utilities::s_my_id-1;
      int err = MPI_Recv(&info,
         1,
         MPI_INT,
         prev_writer,
         s_TAG_BATON,
         MPI_COMM_WORLD,
         &status);
      if (err != MPI_SUCCESS) {
         LOKI_ABORT("Could not get OK from previous writer.");
      }
   }

   m_field_writer->startTimeSlice(a_time,
      a_dt,
      numTrackingParticles(),
      a_num_probes,
      a_probes,
      m_domain->numberOfCells(),
      frame_names,
      is_lo_em_proc);
   plotCommon(a_sequences,
              a_time_seq,
              a_num_probes,
              a_saved_seq,
              a_saved_save,
              a_time_hist_file_name,
              frame_names,
              time_hist_names);

   // Write the plot data for each species following the same pattern as above.
   const ParallelArray::Box& data_box = m_vz[0].dataBox();
   const ParallelArray::Box& local_box = m_vz[0].localBox();
   int name_idx = NUM_EM_VARS;
   for (int is(0); is < m_num_kinetic_species; ++is) {
      // Write this species' vz to the next frame.
      m_field_writer->writeField(frame_names[name_idx++],
         m_vz[is].getData(),
         data_box,
         local_box,
         m_n_ghosts);
      if (a_plot_ke_vel_bdy_flux) {
         // Write this species' vx lo ke flux to the next frame.
         m_field_writer->writeField(frame_names[name_idx++],
            m_ke_flux_vx_lo[is].getData(),
            data_box,
            local_box,
            m_n_ghosts);

         // Write this species' vx hi ke flux to the next frame.
         m_field_writer->writeField(frame_names[name_idx++],
            m_ke_flux_vx_hi[is].getData(),
            data_box,
            local_box,
            m_n_ghosts);

         // Write this species' vy lo ke flux to the next frame.
         m_field_writer->writeField(frame_names[name_idx++],
            m_ke_flux_vy_lo[is].getData(),
            data_box,
            local_box,
            m_n_ghosts);

         // Write this species' vy hi ke flux to the next frame.
         m_field_writer->writeField(frame_names[name_idx++],
            m_ke_flux_vy_hi[is].getData(),
            data_box,
            local_box,
            m_n_ghosts);
      }
   }
   m_field_writer->endTimeSlice(is_lo_em_proc);

   if (multiple_em_procs && Loki_Utilities::s_my_id != m_proc_hi) {
      int info = 0;
      int next_writer = Loki_Utilities::s_my_id+1;
      MPI_Send(&info, 1, MPI_INT, next_writer, s_TAG_BATON, MPI_COMM_WORLD);
   }

   // We finished saving another plot cycle.
   ++a_saved_save;
}


void
Maxwell::accumulateSequences(
   vector<vector<double> >&       a_sequences,
   const vector<vector<double> >& a_probes,
   int                            a_num_probes,
   int                            a_saved_seq,
   int&                           a_seq_idx)
{
   // These are the scalar time history quantities.  Initialize them to 0.
   double e_sum_tot(0.0);
   double e_max(0.0);
   double e_tot(0.0);
   double ex_max(0.0);
   double ey_max(0.0);
   double ez_max(0.0);

   double b_sum_tot(0.0);
   double b_max(0.0);
   double b_tot(0.0);
   double bx_max(0.0);
   double by_max(0.0);
   double bz_max(0.0);

   // Do the same to the probes' time history data.
   vector<vector<double> > emProbe(TIME_HISTS_PER_PROBE);
   vector<vector<double> > emProbeTotal(TIME_HISTS_PER_PROBE);
   for (int ix(0); ix < TIME_HISTS_PER_PROBE; ++ix) {
      emProbe[ix].resize(a_num_probes, 0.0);
      emProbeTotal[ix].resize(a_num_probes, 0.0);
   }

   // Only processors onto which the Maxwell calculation is distributed need to
   // do this.
   if (isEMSolverProcessor()) {
      // Get the interior box bounds.
      const int x_lo = m_em_vars.interiorBox().lower(0);
      const int x_hi = m_em_vars.interiorBox().upper(0);
      const int y_lo = m_em_vars.interiorBox().lower(1);
      const int y_hi = m_em_vars.interiorBox().upper(1);

      // Get some other misc quantities we'll need along the way.
      const tbox::IntVector N(m_domain->numberOfCells());
      const vector<double>& dx(m_domain->dx());
      const double area(dx[X1] * dx[X2]);

      // Compute the local parts of the EM field time histories.
      for (int i1 = x_lo; i1 <= x_hi; ++i1) {
         for (int i2 = y_lo; i2 <= y_hi; ++i2) {
            const double ex = m_em_vars(i1, i2, EX);
            const double ey = m_em_vars(i1, i2, EY);
            const double ez = m_em_vars(i1, i2, EZ);
            const double bx = m_em_vars(i1, i2, BX);
            const double by = m_em_vars(i1, i2, BY);
            const double bz = m_em_vars(i1, i2, BZ);
            double tmp = ex*ex + ey*ey + ez*ez;
            double e_loc = sqrt(tmp);
            e_sum_tot += 0.5*tmp;
            e_max = max(e_max, e_loc);
            e_tot += e_loc;
            ex_max = max(ex_max, fabs(ex));
            ey_max = max(ey_max, fabs(ey));
            ez_max = max(ez_max, fabs(ez));

            tmp = bx*bx + by*by + bz*bz;
            double b_loc = sqrt(tmp);
            b_sum_tot += 0.5*tmp;
            b_max = max(b_max, b_loc);
            b_tot += b_loc;
            bx_max = max(bx_max, fabs(bx));
            by_max = max(by_max, fabs(by));
            bz_max = max(bz_max, fabs(bz));
         }
      }
      // These are integrated over configuration space.
      e_sum_tot *= area;
      e_tot     *= area;
      b_sum_tot *= area;
      b_tot     *= area;

      // Compute the local part of the probe time histories.
      for (int np(0); np < a_num_probes; ++np) {
         const int ip(int(floor(a_probes[0][np] * N[X1])));
         const int jp(int(floor(a_probes[1][np] * N[X2])));
         if (ip >= x_lo && ip <= x_hi && jp >= y_lo && jp <= y_hi) {
            emProbe[EX][np] = m_em_vars(ip, jp, EX);
            emProbe[EY][np] = m_em_vars(ip, jp, EY);
            emProbe[EZ][np] = m_em_vars(ip, jp, EZ);
            emProbe[BX][np] = m_em_vars(ip, jp, BX);
            emProbe[BY][np] = m_em_vars(ip, jp, BY);
            emProbe[BZ][np] = m_em_vars(ip, jp, BZ);
            for (int ns(0); ns < m_num_kinetic_species; ++ns) {
               emProbe[NUM_EM_VARS+ns][np] = m_vz[ns](ip, jp);
            }                
         }
      }
   }

   // Now get the appropriate max/sum of the EM field time histories from all
   // the processors and store to the time history.
   a_sequences[a_seq_idx++][a_saved_seq] =
      Loki_Utilities::getMaxValue(e_max, -1);
   a_sequences[a_seq_idx++][a_saved_seq] =
      Loki_Utilities::getSum(e_tot, -1);
   a_sequences[a_seq_idx++][a_saved_seq] =
      Loki_Utilities::getMaxValue(ex_max, -1);
   a_sequences[a_seq_idx++][a_saved_seq] =
      Loki_Utilities::getMaxValue(ey_max, -1);
   a_sequences[a_seq_idx++][a_saved_seq] =
      Loki_Utilities::getMaxValue(ez_max, -1);
   a_sequences[a_seq_idx++][a_saved_seq] =
      Loki_Utilities::getSum(e_sum_tot, -1);

   a_sequences[a_seq_idx++][a_saved_seq] =
      Loki_Utilities::getMaxValue(b_max, -1);
   a_sequences[a_seq_idx++][a_saved_seq] =
      Loki_Utilities::getSum(b_tot, -1);
   a_sequences[a_seq_idx++][a_saved_seq] =
      Loki_Utilities::getMaxValue(bx_max, -1);
   a_sequences[a_seq_idx++][a_saved_seq] =
      Loki_Utilities::getMaxValue(by_max, -1);
   a_sequences[a_seq_idx++][a_saved_seq] =
      Loki_Utilities::getMaxValue(bz_max, -1);
   a_sequences[a_seq_idx++][a_saved_seq] =
      Loki_Utilities::getSum(b_sum_tot, -1);

   // Get the sum of the probe time histories from all the processors and store
   // to the time history.
   for (int i = 0; i < TIME_HISTS_PER_PROBE; ++i) {
      Loki_Utilities::getSums(&emProbe[i][0],
         &emProbeTotal[i][0],
         a_num_probes,
         -1);
   }
   for (int ip(0); ip < a_num_probes; ++ip) {
      a_sequences[a_seq_idx++][a_saved_seq] = emProbeTotal[EX][ip];
      a_sequences[a_seq_idx++][a_saved_seq] = emProbeTotal[EY][ip];
      a_sequences[a_seq_idx++][a_saved_seq] = emProbeTotal[EZ][ip];
      a_sequences[a_seq_idx++][a_saved_seq] = emProbeTotal[BX][ip];
      a_sequences[a_seq_idx++][a_saved_seq] = emProbeTotal[BY][ip];
      a_sequences[a_seq_idx++][a_saved_seq] = emProbeTotal[BZ][ip];
      for (int is(0); is < m_num_kinetic_species;  ++is) {
         a_sequences[a_seq_idx++][a_saved_seq] =
            emProbeTotal[NUM_EM_VARS+is][ip];
      }
   }

   // Store the tracking particle positions and velocities to the time history.
   for (int ip = 0; ip < numTrackingParticles(); ++ip) {
      const Particle& tracking_particle = m_tracking_particles[ip];
      a_sequences[a_seq_idx++][a_saved_seq] = tracking_particle.x();
      a_sequences[a_seq_idx++][a_saved_seq] = tracking_particle.y();
      a_sequences[a_seq_idx++][a_saved_seq] = tracking_particle.vx();
      a_sequences[a_seq_idx++][a_saved_seq] = tracking_particle.vy();
   }
}


void
Maxwell::getFromRestart(
   RestartReader& a_reader)
{
   // Get the subdatabase containing this object's data.
   a_reader.pushSubDir("Maxwell");

   // Get the electromagnetic variables from the subdatabase.
   a_reader.readParallelArray("EMVars", m_em_vars);
   for (int i = 0; i < m_num_kinetic_species; ++i) {
      ostringstream name;
      name << "vz" << i;
      a_reader.readParallelArray(name.str(), m_vz[i]);
      name.str("");
   }

   // We don't read the physical boundary so we need to zero out that data.
   FORT_ZERO_GHOST_2D(*m_em_vars.getData(),
      BOX2D_TO_FORT(m_em_vars.interiorBox()),
      BOX2D_TO_FORT(m_em_vars.dataBox()),
      NUM_EM_VARS);
   for (int i = 0; i < m_num_kinetic_species; ++i) {
      FORT_ZERO_GHOST_2D(*m_vz[i].getData(),
         BOX2D_TO_FORT(m_vz[i].interiorBox()),
         BOX2D_TO_FORT(m_vz[i].dataBox()),
         1);
   }

   a_reader.popSubDir();
}


void
Maxwell::putToRestart(
   RestartWriter& a_writer,
   double a_time)
{
   NULL_USE(a_time);
   bool write_data = Loki_Utilities::s_my_id == m_proc_lo;
   putToRestartCommon(a_writer);
   a_writer.writeIntegerValue("isMaxwell", 1, write_data);

   // Make a subdirectory.
   a_writer.pushSubDir("Maxwell");

   // Put the electromagnetic variables into the subdirectory.
   a_writer.writeParallelArray("EMVars", m_em_vars, write_data);
   for (int i = 0; i < m_num_kinetic_species; ++i) {
      ostringstream name;
      name << "vz" << i;
      a_writer.writeParallelArray(name.str(), m_vz[i], write_data);
      name.str("");
   }

   a_writer.popSubDir();
}


void
Maxwell::buildTimeHistoryNames(
   int a_num_probes,
   int a_num_tracking_particles,
   const vector<string>& a_species_names,
   vector<string>& a_time_hist_names)
{
   ostringstream th_name;
   int num_species = static_cast<int>(a_species_names.size());
   a_time_hist_names.push_back("E_max");
   a_time_hist_names.push_back("norm E");
   a_time_hist_names.push_back("Ex_max");
   a_time_hist_names.push_back("Ey_max");
   a_time_hist_names.push_back("Ez_max");
   a_time_hist_names.push_back("E_tot");
   a_time_hist_names.push_back("B_max");
   a_time_hist_names.push_back("norm B");
   a_time_hist_names.push_back("Bx_max");
   a_time_hist_names.push_back("By_max");
   a_time_hist_names.push_back("Bz_max");
   a_time_hist_names.push_back("B_tot");
   for (int ip = 0; ip < a_num_probes; ++ip) {
      th_name << "Ex_probe" << ip;
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");

      th_name << "Ey_probe" << ip;
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");

      th_name << "Ez_probe" << ip;
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");

      th_name << "Bx_probe" << ip;
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");

      th_name << "By_probe" << ip;
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");

      th_name << "Bz_probe" << ip;
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");

      for (int is = 0; is < num_species; ++is) {
         th_name << "vz_probe" << ip << "_" << a_species_names[is];
         a_time_hist_names.push_back(th_name.str());
         th_name.str("");
      }
   }
   for (int ipart = 0; ipart < a_num_tracking_particles; ++ipart) {
      th_name << "particle" << ipart << "_x";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << "particle" << ipart << "_y";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << "particle" << ipart << "_vx";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << "particle" << ipart << "_vy";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
   }
   for (int is = 0; is < num_species; ++is) {
      th_name << a_species_names[is] << "_ke";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_ke_x";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_ke_y";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_px";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_py";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_xlo_flux";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_xhi_flux";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_ylo_flux";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_yhi_flux";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_vxlo_flux";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_vxhi_flux";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_vylo_flux";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_vyhi_flux";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_ke_e_dot";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_integrated_ke_e_dot";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_driver_time_envel";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
   }
}


void
Maxwell::buildCollTimeHistoryNames(
   const vector<string>& a_species_names,
   vector<string>& a_time_hist_names)
{
   ostringstream th_name;
   int num_species = static_cast<int>(a_species_names.size());
   for (int is = 0; is < num_species; ++is) {
      th_name << a_species_names[is] << "_dMomX";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_dMomY";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_dKE";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_dEnt";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_iMomX";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_iMomY";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_iKE";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
      th_name << a_species_names[is] << "_iEnt";
      a_time_hist_names.push_back(th_name.str());
      th_name.str("");
   }
}


void
Maxwell::buildPlotNames(
   bool a_plot_ke_vel_bdy_flux,
   const vector<string>& a_species_names,
   vector<string>& a_plot_names)
{
   int num_species = static_cast<int>(a_species_names.size());
   a_plot_names.push_back("EX");
   a_plot_names.push_back("EY");
   a_plot_names.push_back("EZ");
   a_plot_names.push_back("BX");
   a_plot_names.push_back("BY");
   a_plot_names.push_back("BZ");
   for (int i = 0; i < num_species; ++i) {
      ostringstream plot_name;
      plot_name << a_species_names[i] << " VZ";
      a_plot_names.push_back(plot_name.str());
      plot_name.str("");
      if (a_plot_ke_vel_bdy_flux) {
         plot_name << a_species_names[i] << " ke flux vx lo";
         a_plot_names.push_back(plot_name.str());
         plot_name.str("");
         plot_name << a_species_names[i] << " ke flux vx hi";
         a_plot_names.push_back(plot_name.str());
         plot_name.str("");
         plot_name << a_species_names[i] << " ke flux vy lo";
         a_plot_names.push_back(plot_name.str());
         plot_name.str("");
         plot_name << a_species_names[i] << " ke flux vy hi";
         a_plot_names.push_back(plot_name.str());
         plot_name.str("");
      }
   }
}


void
Maxwell::buildCollPlotNames(
   const vector<string>& a_species_names,
   vector<string>& a_plot_names)
{
   int num_species = static_cast<int>(a_species_names.size());
   for (int i = 0; i < num_species; ++i) {
      ostringstream plot_name;
      plot_name << a_species_names[i] << " diag MomX";
      a_plot_names.push_back(plot_name.str());
      plot_name.str("");
      plot_name << a_species_names[i] << " diag MomY";
      a_plot_names.push_back(plot_name.str());
      plot_name.str("");
      plot_name << a_species_names[i] << " diag KE";
      a_plot_names.push_back(plot_name.str());
      plot_name.str("");
      plot_name << a_species_names[i] << " diag Ent";
      a_plot_names.push_back(plot_name.str());
      plot_name.str("");
      plot_name << a_species_names[i] << " MomX";
      a_plot_names.push_back(plot_name.str());
      plot_name.str("");
      plot_name << a_species_names[i] << " MomY";
      a_plot_names.push_back(plot_name.str());
      plot_name.str("");
      plot_name << a_species_names[i] << " KE";
      a_plot_names.push_back(plot_name.str());
      plot_name.str("");
      plot_name << a_species_names[i] << " Ent";
      a_plot_names.push_back(plot_name.str());
      plot_name.str("");
   }
}


// Protected member functions.


void
Maxwell::addNoiseToChargeDist(
   TimerManager* a_timers,
   ParallelArray& a_charge_density,
   double a_time)
{
   // Noise particles contribute to the current density, not the charge density
   // for Maxwell.
   NULL_USE(a_timers);
   NULL_USE(a_charge_density);
   NULL_USE(a_time);
}


void
Maxwell::addExternalPotential(
   ParallelArray& a_phi_local,
   double a_time)
{
   // There are no external potential drivers for Maxwell.
   NULL_USE(a_phi_local);
   NULL_USE(a_time);
}


void
Maxwell::startBCTimer(
   TimerManager* a_timers)
{
   a_timers->startTimer("BC (Maxwell)");
}


void
Maxwell::stopBCTimer(
   TimerManager* a_timers)
{
   a_timers->stopTimer("BC (Maxwell)");
}


void
Maxwell::UGB()
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("parallel ghost");
   m_em_vars.communicateGhostData();
   for (int i = 0; i < m_num_kinetic_species; ++i) {
      m_vz[i].communicateGhostData();
   }
   timers->stopTimer("parallel ghost");
}


void
Maxwell::setPhysicalBCs(
   bool a_particles_only)
{
   if (!a_particles_only) {
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("BC (Maxwell)");
      // Cover both non-periodic and periodic possibilities.  The Fortran
      // handles non-periodic and the C is for periodic.
      FORT_SET_MAXWELL_EM_BCS(BOX2D_TO_FORT(m_em_vars.dataBox()),
         BOX2D_TO_FORT(m_em_vars.interiorBox()),
         *(m_em_vars.getData()),
         m_domain->numberOfCells(X1),
         m_domain->numberOfCells(X2),
         m_domain->isPeriodicAsInt(X1),
         m_domain->isPeriodicAsInt(X2),
         m_solution_order,
         Simulation::s_LIGHT_SPEED);
      m_em_vars.communicatePeriodicBoundaries();
      for (int i = 0; i < m_num_kinetic_species; ++i) {
         // Cover both non-periodic and periodic possibilities.  The Fortran
         // handles non-periodic and the C is for periodic.
         FORT_SET_MAXWELL_VZ_BCS(BOX2D_TO_FORT(m_vz[i].dataBox()),
            BOX2D_TO_FORT(m_vz[i].interiorBox()),
            *(m_vz[i].getData()),
            m_domain->numberOfCells(X1),
            m_domain->numberOfCells(X2),
            m_domain->isPeriodicAsInt(X1),
            m_domain->isPeriodicAsInt(X2),
            m_solution_order);
         m_vz[i].communicatePeriodicBoundaries();
      }
      timers->stopTimer("BC (Maxwell)");
   }
   setParticleBCs();
}


// Private member functions.


void
Maxwell::parseParameters(
   LokiInputParser& a_pp)
{
   // Get the Maxwell specific inputs.
   // See if these optional inputs have been specified.
   a_pp.query("avWeak", m_avWeak);
   a_pp.query("avStrong", m_avStrong);

   // See if there are any current drivers.
   a_pp.query("num_current_drivers", m_num_current_drivers);
   m_current_drivers.resize(m_num_current_drivers,
      tbox::Pointer<CurrentDriver>(0));

   // Get the super grid if specified.
   for (int d = 0; d < m_dim; ++d) {
      ostringstream name;
      name << "supergrid" << d+1 << "a";
      a_pp.query(name.str().c_str(), m_supergrid_lo[d]);
      name.str("");
      if (m_domain->lower(d) > m_supergrid_lo[d]) {
         ostringstream msg;
         msg << "Maxwell super grid must be inside physical domain (supergrid"
             << d+1 << "a)";
         LOKI_ABORT(msg.str());
      }

      name << "supergrid" << d+1 << "b";
      a_pp.query(name.str().c_str(), m_supergrid_hi[d]);
      name.str("");
      if (m_domain->upper(d) < m_supergrid_hi[d]) {
         ostringstream msg;
         msg << "Maxwell super grid must be inside physical domain (supergrid"
             << d+1 << "b)";
         LOKI_ABORT(msg.str());
      }
   }

   // Now get the parameters common to all EMSolvers.
   parseParametersCommon(a_pp);
}


void
Maxwell::printParameters() const
{
   // Print this Maxwell's parameters and the parameters of the entities that it
   // holds and are not accessible elsewhere like the current drivers and EM
   // initial conditions.
   Loki_Utilities::printF("\n#*#*# Maxwell #*#*#\n");
   Loki_Utilities::printF("  Number of Processsors = %i\n", m_number_of_procs);
   Loki_Utilities::printF("  avWeak                = %f\n", m_avWeak);
   Loki_Utilities::printF("  avStrong              = %f\n", m_avStrong);
   Loki_Utilities::printF("  super grid x bounds   = (%f, %f)\n",
      m_supergrid_lo[X1],
      m_supergrid_hi[X1]);
   Loki_Utilities::printF("  super grid y bounds   = (%f, %f)\n",
      m_supergrid_lo[X2],
      m_supergrid_hi[X2]);
   Loki_Utilities::printF("\n  num current drivers = %i\n",
      m_num_current_drivers);
   for (int i = 0; i < m_num_current_drivers; ++i) {
      m_current_drivers[i]->printParameters();
   }
   for (int i = 0; i < static_cast<int>(m_em_ics.size()); ++i) {
      m_em_ics[i]->printParameters();
   }
   for (int i = 0; i < static_cast<int>(m_vel_ics.size()); ++i) {
      m_vel_ics[i]->printParameters();
   }
}


void
Maxwell::evalRHSParticles(
   bool a_noise_source_particles,
   Maxwell& a_rhs,
   const ParallelArray& a_net_ext_efield,
   double a_time)
{
   int num_particles = a_noise_source_particles ?
      numNoiseSourceParticles() : numTrackingParticles();
   const vector<Particle>& particles = a_noise_source_particles ?
      m_noise_source_particles : m_tracking_particles;
   vector<Particle>& rhs_particles = a_noise_source_particles ?
      a_rhs.m_noise_source_particles : a_rhs.m_tracking_particles;

   // For each particle use an interpolating polynomial to compute the E and
   // B fields that it sees and update x, y, vx, and vy.
   double x_lo = m_domain->lower(0);
   double dx = m_domain->dx(0);
   double y_lo = m_domain->lower(1);
   double dy = m_domain->dx(1);
   const ParallelArray::Box& interior_box = m_em_vars.interiorBox();
   int ib_x_lo = interior_box.lower(0);
   int ib_x_hi = interior_box.upper(0);
   int ib_y_lo = interior_box.lower(1);
   int ib_y_hi = interior_box.upper(1);
   for (int i = 0; i < num_particles; ++i) {
      // Check if this particle is on this Maxwell processor and that it is
      // time to begin to move it.
      const Particle& this_particle = particles[i];
      double x_loc = this_particle.x();
      double y_loc = this_particle.y();
      int x_idx = (x_loc-x_lo)/dx;
      int y_idx = (y_loc-y_lo)/dy;
      double starting_time =
         a_noise_source_particles ? 0.0 : this_particle.startingTime();
      if (a_time < starting_time ||
          x_idx < ib_x_lo || x_idx > ib_x_hi ||
          y_idx < ib_y_lo || y_idx > ib_y_hi) {
         continue;
      }
      Particle& rhs_particle = rhs_particles[i];
      // Set the evaluation point to be the location of this particle.
      m_em_interpolator->setEvaluationPoint(x_loc, y_loc);

      // Figure out Ex, Ey, and Bz at each interpolating point.
      for (int xidx = 0; xidx < m_solution_order; ++xidx) {
         int xGridIdx = m_em_interpolator->xInterpolatingGridIndex(xidx);
         for (int yidx = 0; yidx < m_solution_order; ++yidx) {
            int yGridIdx = m_em_interpolator->yInterpolatingGridIndex(yidx);
            m_ex_interp[xidx][yidx] = m_em_vars(xGridIdx, yGridIdx, EX) +
               a_net_ext_efield(xGridIdx, yGridIdx, EX);
            m_ey_interp[xidx][yidx] = m_em_vars(xGridIdx, yGridIdx, EY) +
               a_net_ext_efield(xGridIdx, yGridIdx, EY);
            m_bz_interp[xidx][yidx] = m_em_vars(xGridIdx, yGridIdx, BZ);
         }
      }

      // Now compute Ex, Ey, Bz, for this particle.
      double Ex_particle = m_em_interpolator->interpolate(m_ex_interp);
      double Ey_particle = m_em_interpolator->interpolate(m_ey_interp);
      double Bz_particle = m_em_interpolator->interpolate(m_bz_interp);

      // The new x/y is the old vx/vy.  The new vx/vy is determined by the
      // Lorentz force.
      double t_ramp =
         a_noise_source_particles ? this_particle.startingTime() : 0.0;
      double envel;
      if (a_time < t_ramp) {
         envel = 0.5+0.5*tanh(8.0*a_time/t_ramp-4.0);
      }
      else {
         envel = 1.0;
      }
      double qoverm = envel*this_particle.charge()/this_particle.mass();
      double vx = this_particle.vx();
      double vy = this_particle.vy();
      rhs_particle.x() = vx;
      rhs_particle.y() = vy;
      rhs_particle.vx() =
         (Ex_particle+vy*(Bz_particle+m_bz_const))*qoverm; //IEO
      rhs_particle.vy() =
         (Ey_particle-vx*(Bz_particle+m_bz_const))*qoverm; //IEO
   }
}


void
Maxwell::buildInterpolators()
{
   // If this object is on the Maxwell processor and there are either tracking
   // or noise source particles then we need to construct the interpolators of
   // the appropriate order.
   if (isEMSolverProcessor() &&
       (m_problem_num_tracking_particles > 0 ||
        m_problem_num_noise_source_particles > 0)) {
      if (m_solution_order == 4) {
         m_em_interpolator = new Interpolator4(*m_domain);
      }
      else {
         m_em_interpolator = new Interpolator6(*m_domain);
      }
      m_ex_interp = new double* [m_solution_order];
      m_ey_interp = new double* [m_solution_order];
      m_bz_interp = new double* [m_solution_order];
      for (int i = 0; i < m_solution_order; ++i) {
         m_ex_interp[i] = new double [m_solution_order];
         m_ey_interp[i] = new double [m_solution_order];
         m_bz_interp[i] = new double [m_solution_order];
      }
      if (m_problem_num_noise_source_particles > 0) {
         if (m_solution_order == 4) {
            m_j_interpolator = new Interpolator4(*m_domain);
         }
         else {
            m_j_interpolator = new Interpolator6(*m_domain);
         }
         m_jx_interp = new double* [m_solution_order];
         for (int i = 0; i < m_solution_order; ++i) {
            m_jx_interp[i] = new double [m_solution_order];
         }
         m_jy_interp = new double* [m_solution_order];
         for (int i = 0; i < m_solution_order; ++i) {
            m_jy_interp[i] = new double [m_solution_order];
         }
      }
   }
}

} // end namespace Loki
