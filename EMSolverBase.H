/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _EMSOLVERBASE_H_
#define _EMSOLVERBASE_H_

#include <string>
#include "LokiInputParser.H"
#include "tbox/DescribedClass.H"
#include "tbox/Dimension.H"
#include "tbox/Pointer.H"
#include "tbox/Box.H"
#include "Directions.H"
#include "FieldWriter.H"
#include "Load.H"
#include "Loki_Defines.H"
#include "Loki_Utilities.H"
#include "PoissonSolver.H"
#include "Serializable.H"
#include "TimerManager.H"
#include "ParallelArray.H"
#include "ProblemDomain.H"
#include "Particle.H"

namespace Loki {

/**
 * Abstract base class for implementation of the common parts of the solver of
 * the electromagnetics of a Vlasov code.  Also implements much of the
 * abstract Load interface.
 */
class EMSolverBase : public tbox::DescribedClass, public Load,
                     public Serializable
{
public:
   /**
    * @brief Constructor
    *
    * @param[in] a_domain The global problem domain.
    * @param[in] a_num_kinetic_species Number of KineticSpecies in the problem.
    * @param[in] a_solution_order Order of solution.
    * @param[in] a_num_em_vars Number of EM variables this solver solves for.
    * @param[in] a_plot_times_per_file The number of time slices in a field
    *                                  file.
    * @param[in] a_plot_ke_vel_bdy_flux If true plot species' KE vel bdy fluxes.
    * @param[in] a_coll_diag_on If true, collision diagnostics are on.
    * @param[in] a_bz_const Constant Bz for particle tracking.
    */
   EMSolverBase(
      const tbox::Pointer<ProblemDomain>& a_domain,
      int a_num_kinetic_species,
      int a_solution_order,
      int a_num_em_vars,
      int a_plot_times_per_file,
      bool a_plot_ke_vel_bdy_flux,
      bool a_coll_diag_on,
      double a_bz_const);

   /**
    * @brief Copy constructor.
    *
    * @param[in] a_other The EMSolverBase to copy from.
    */
   EMSolverBase(
      const EMSolverBase& a_other);

   /**
    * @brief Destructor
    */
   virtual
   ~EMSolverBase();

   //@{

   //! @name Implementation of Load interface.

   /**
    * @brief Returns the computational cost of this load.
    */
   virtual
   float
   netCost() const;

   /**
    * @brief Returns the number of processors over which this load is
    *        distributed.
    */
   virtual
   int
   numberOfProcessors() const;

   /**
    * @brief Returns true if this load is distributed over a user specified
    *        fixed number of processors.
    */
   virtual
   bool
   fixedNumberOfProcessors() const;

   /**
    * @brief Returns true if this load is distributed to processor with
    *        supplied id.
    *
    * @param[in] a_proc_id ID of processor of interest.
    */
   virtual
   bool
   isInRange(
      int a_proc_id) const;

   /**
    * @brief Print information about processors this load is distributed over.
    */
   virtual
   void
   printDecomposition() const;

   //@}

   /**
    * @brief Allocate and distribute a variable defined on the 2D EM domain.
    *
    * @param[out] a_var The allocated variable.
    * @param[in] a_depth The number of components in var.
    * @param[in] a_n_ghosts The number of ghost cells var has in each
    *                       direction.
    */
   void
   newAuxVariable(
      ParallelArray& a_var,
      int a_depth,
      int a_n_ghosts) const;

   /**
    * @brief Allocate and distribute a variable defined on the 2D EM domain
    *        using the number of ghosts of this class.
    *
    * @param[out] a_var The allocated variable.
    * @param[in] a_depth The number of components in var.
    */
   void
   newAuxVariable(
      ParallelArray& a_var,
      int a_depth) const
   {
      newAuxVariable(a_var, a_depth, m_n_ghosts);
   }

   /**
    * @brief Accessor for the global kinetic energy flux for a given
    *        KineticSpecies, direction and side.
    *
    * @param[in] a_species_idx Index of the species of interest.
    * @param[in] a_dir x or y direction.
    * @param[in] a_side High or low side.
    */
   ParallelArray&
   getKEFluxVar(
      int a_species_idx,
      int a_dir,
      int a_side)
   {
      if (a_dir != V1 && a_dir != V2) {
         LOKI_ABORT("Illegal direction requested.");
      }
      if (a_side != LO && a_side != HI) {
         LOKI_ABORT("Illegal side requested.");
      }
      if (a_dir == V1) {
         if (a_side == LO) {
            return m_ke_flux_vx_lo[a_species_idx];
         }
         else {
            return m_ke_flux_vx_hi[a_species_idx];
         }
      }
      else {
         if (a_side == LO) {
            return m_ke_flux_vy_lo[a_species_idx];
         }
         else {
            return m_ke_flux_vy_hi[a_species_idx];
         }
      }
   }

   /**
    * @brief Accessor for the global electromagnetic fields.
    */
   const ParallelArray&
   getEMVars() const
   {
      return m_em_vars;
   }

   /**
    * @brief Returns the number of tracking particles in the problem.
    */
   int
   numProblemTrackingParticles() const
   {
      return m_problem_num_tracking_particles;
   }

   /**
    * @brief Returns the number of noise source particles in the problem.
    */
   int
   numProblemNoiseSourceParticles() const
   {
      return m_problem_num_noise_source_particles;
   }

   /**
    * @brief Read in tracking and noise source particles if any.
    *
    * @param[in] a_pp Where the input parameters will be found.
    * @param[in] a_vmin Minimal end of species maximal bounding box.
    * @param[in] a_vmax Maximal end of species maximal bounding box.
    */
   void
   readParticles(
      LokiInputParser& a_pp,
      const double* a_vmin,
      const double* a_vmax);

   /**
    * @brief Returns true if this load has been distributed to the processor
    *        on which the query is made.
    */
   bool
   isEMSolverProcessor() const
   {
      return m_is_emsolver_processor;
   }

   /**
    * @brief Compute the electric field resulting from the supplied charge
    *        density and any externally applied electric potential.
    *
    * @param[in] a_charge_density Charge density resulting from all
    *                             KineticSpecies.
    * @param[in] a_time Current simulation time to evaluate externally applied
    *                   potential.
    */
   void
   electricField(
      ParallelArray& a_charge_density,
      double a_time);

   /**
    * @brief Returns the problem domain.
    */
   const ProblemDomain&
   domain() const
   {
      return *m_domain;
   }

   /**
    * @brief Returns the solution order.
    */
   int
   order() const
   {
      return m_solution_order;
   }

   /**
    * @brief The number of EM processors.
    */
   int
   num_procs() const
   {
      return m_number_of_procs;
   }

   /**
    * @brief Returns the EM communicator.
    */
   MPI_Comm
   comm() const
   {
      return m_comm;
   }

protected:
   // Part of createPartition common to all EMSolvers.
   // Currently we can only enforce periodic boundary conditions on the Poisson
   // solve so we must force periodicity for this case.
   void
   createPartitionCommon(
      bool a_enforce_periodicity,
      int a_proc_lo,
      int a_proc_hi,
      const MPI_Comm& a_comm);

   // Part of putToRestart common to all EMSolvers.
   void
   putToRestartCommon(
      RestartWriter& a_writer)
   {
      bool write_data = Loki_Utilities::s_my_id == m_proc_lo;

      // We write the number of ghosts only so that the post processor knows it.
      a_writer.writeIntegerValue("nGhost", m_n_ghosts, write_data);
 
      // Now write the tracking and noise source particle data.
      if (write_data) {
         if (numTrackingParticles() > 0) {
            Particle::putParticlesToRestart(m_tracking_particle_file,
                                            m_tracking_particles,
                                            false);
         }
         if (numNoiseSourceParticles() > 0) {
            Particle::putParticlesToRestart(m_noise_source_particle_file,
                                            m_noise_source_particles,
                                            true);
         }
      }
   }

   // Zeros out all particle data.
   void
   zeroParticleData()
   {
      for (int i = 0; i < numTrackingParticles(); ++i) {
         m_tracking_particles[i].x() = 0.0;
         m_tracking_particles[i].y() = 0.0;
         m_tracking_particles[i].vx() = 0.0;
         m_tracking_particles[i].vy() = 0.0;
      }
      for (int i = 0; i < numNoiseSourceParticles(); ++i) {
         m_noise_source_particles[i].x() = 0.0;
         m_noise_source_particles[i].y() = 0.0;
         m_noise_source_particles[i].vx() = 0.0;
         m_noise_source_particles[i].vy() = 0.0;
      }
   }

   // Adds all particle data.
   void
   addParticleData(
      const vector<Particle>& a_tracking_particles_increment,
      const vector<Particle>& a_noise_source_particles_increment,
      double a_factor,
      bool a_sum_reduce_inc)
   {
      TimerManager* timers(TimerManager::getManager());

      // Handle the tracking particles.
      if (m_problem_num_tracking_particles > 0) {
         timers->startTimer("Tracking Particles");
         Particle::addParticleData(m_tracking_particles,
            a_tracking_particles_increment,
            a_factor,
            a_sum_reduce_inc && m_proc_lo != m_proc_hi,
            m_comm);
         timers->stopTimer("Tracking Particles");
      }

      // Handle the noise source particles.
      if (m_problem_num_noise_source_particles > 0) {
         timers->startTimer("Noisy Particles");
         Particle::addParticleData(m_noise_source_particles,
            a_noise_source_particles_increment,
            a_factor,
            a_sum_reduce_inc && m_proc_lo != m_proc_hi,
            m_comm);
         timers->stopTimer("Noisy Particles");
      }
   }

   // Copies all particle data from the rhs.
   void
   copyParticleData(
      const vector<Particle>& a_tracking_particles_rhs,
      const vector<Particle>& a_noise_source_particles_rhs)
   {
      TimerManager* timers(TimerManager::getManager());

      // Copy the tracking particles' data.
      if (m_problem_num_tracking_particles > 0) {
         timers->startTimer("Tracking Particles");
         Particle::copyParticleData(m_tracking_particles,
                                    a_tracking_particles_rhs);
         timers->stopTimer("Tracking Particles");
      }

      // Copy the noise source particles' data.
      if (m_problem_num_noise_source_particles > 0) {
         timers->startTimer("Noisy Particles");
         Particle::copyParticleData(m_noise_source_particles,
                                    a_noise_source_particles_rhs);
         timers->stopTimer("Noisy Particles");
      }
   }

   // Set boundary conditions on particles.
   void
   setParticleBCs()
   {
      TimerManager* timers(TimerManager::getManager());

      // Set tracking particle boundary conditions.
      if (m_problem_num_tracking_particles > 0) {
         timers->startTimer("Tracking Particles");
         Particle::setBC(m_tracking_particles,
            m_domain,
            m_vmin,
            m_vmax,
            false);
         timers->stopTimer("Tracking Particles");
      }

      // Set noise source particle boundary conditions.
      if (m_problem_num_noise_source_particles > 0) {
         timers->startTimer("Noisy Particles");
         Particle::setBC(m_noise_source_particles,
            m_domain,
            m_vmin,
            m_vmax,
            m_reflect_noise_source_particles);
         timers->stopTimer("Noisy Particles");
      }
   }

   // Part of initialize common to all EMSolvers.
   void
   initializeCommon(
      int a_num_species,
      bool a_plot_ke_vel_bdy_flux,
      bool a_coll_diag_on);

   // Part of plot common to all EMSolvers.
   void
   plotCommon(
      const vector<vector<double> >& a_sequences,
      const vector<double>& a_time_seq,
      int a_num_probes,
      int a_saved_seq,
      int& a_saved_save,
      string& a_time_hist_file_name,
      vector<string>& a_frame_names,
      vector<string>& a_time_hist_names);

   // Specialized version to write collision time histories to the time history
   // file.
   void
   writeTimeHistories(
      const vector<vector<double> >& a_sequences,
      const vector<string>& a_name,
      const vector<double>& a_time_seq,
      int a_saved_seq,
      int& a_saved_save,
      string& a_time_hist_file_name);

   // Write the time histories to the time history file.
   void
   writeTimeHistories(
      const vector<vector<double> >& a_sequences,
      const vector<string>& a_name,
      const vector<double>& a_time_seq,
      int a_num_probes,
      int a_saved_seq,
      int& a_saved_save,
      string& a_time_hist_file_name);

   // Part of parseParameters common to all EMSolvers.
   void
   parseParametersCommon(
      LokiInputParser& a_pp);

   // Adds noise to the charge distribution
   virtual
   void
   addNoiseToChargeDist(
      TimerManager* a_timers,
      ParallelArray& a_charge_density,
      double a_time) = 0;

   // Adds potential from any external potential drivers.
   virtual
   void
   addExternalPotential(
      ParallelArray& a_phi_local,
      double a_time) = 0;

   // Starts the right boundary condition timer.
   virtual
   void
   startBCTimer(
      TimerManager* a_timers) = 0;

   // Stops the right boundary condition timer.
   virtual
   void
   stopBCTimer(
      TimerManager* a_timers) = 0;

   // Returns the number of tracking particles in object.
   int
   numTrackingParticles() const
   {
      return static_cast<int>(m_tracking_particles.size());
   }

   // Returns the number of noise source particles in object.
   int
   numNoiseSourceParticles() const
   {
      return static_cast<int>(m_noise_source_particles.size());
   }

   // Configuration space dimension this class operates in.
   tbox::Dimension m_dim;

   // Global problem domain.
   tbox::Pointer<ProblemDomain> m_domain;

   // The number of processors this Poisson is distributed over.
   int m_number_of_procs;

   // The lowest rank processor this Poisson is distributed over.
   int m_proc_lo;

   // The highest rank processor this Poisson is distributed over.
   int m_proc_hi;

   // If true, the how this is partitioned across processors has been defined.
   bool m_partition_defined;

   // The number of ghost cells in each direction.
   int m_n_ghosts;

   // The communicator for the processors this is distributed over.
   MPI_Comm m_comm;

   // The electromagnetic fields.
   ParallelArray m_em_vars;

   // Number of EM variables being solved for.
   int m_num_em_vars;

   // Loki Poisson solver.
   tbox::Pointer<PoissonSolver> m_loki_solver;

   // Array containing potential.
   ParallelArray m_phi_solver;

   // Array containing the serialized potential.  This is only needed for
   // parallel EM.  For serial EM m_phi_solver and m_phi_global are the same.
   ParallelArray m_phi_global;

   // Array containing the serialized charge density.  This is only needed for
   // parallel EM.  For serial EM the charge density passed into electricField
   // is already seialized.
   ParallelArray m_rho_global;

   // True if this has been distributed to the proceesor owning it.
   bool m_is_emsolver_processor;

   // For each KineticSpecies, the flux of kinetic energy through each velocity
   // space boundary.
   vector<ParallelArray> m_ke_flux_vx_hi;
   vector<ParallelArray> m_ke_flux_vx_lo;
   vector<ParallelArray> m_ke_flux_vy_hi;
   vector<ParallelArray> m_ke_flux_vy_lo;

   // Diagnostic fields for each species.
   vector<ParallelArray> m_diag_plot;

   // Order of solution.
   int m_solution_order;

   // Tracking particle file.
   string m_tracking_particle_file;

   // Tracking particles
   vector<Particle> m_tracking_particles;

   // Number of tracking particles in the problem.
   int m_problem_num_tracking_particles;

   // Noise source particle file.
   string m_noise_source_particle_file;

   // Noise source particles
   vector<Particle> m_noise_source_particles;

   // Number of noise source particles in the problem.
   int m_problem_num_noise_source_particles;

   // If true noise source particles that cross the physical boundary will be
   // reflected.
   bool m_reflect_noise_source_particles;

   // Constant Bz component for particle tracking
   double m_bz_const; //IEO

   // Maximal velocity space bounding box for all species.
   double m_vmin[2];
   double m_vmax[2];

   // The object that writes plot data.
   tbox::Pointer<FieldWriter> m_field_writer;
   tbox::Pointer<FieldWriter> m_field_writer_coll;

   // The number of KineticSpecies in the problem.
   int m_num_kinetic_species;

private:
   // Unimplemented default constructor.
   EMSolverBase();

   // Unimplemented assignment operator.
   EMSolverBase&
   operator = (
      const EMSolverBase& rhs);
};

} // end namespace Loki

#endif
