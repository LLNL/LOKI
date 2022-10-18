/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _MAXWELL_H_
#define _MAXWELL_H_

#include "CurrentDriver.H"
#include "EMICInterface.H"
#include "EMSolverBase.H"
#include "VELICInterface.H"
#include "KineticSpeciesPtrVect.H"
#include "Interpolator.H"
#include "Simulation.H"

namespace Loki {

/**
 * Implementation of abstract Load interface for the solution of Maxwell's
 * equations.
 */
class Maxwell : public EMSolverBase
{
public:
   enum EM_Var_Dims{EX, EY, EZ, BX, BY, BZ, NUM_EM_VARS};

   /**
    * @brief The number of time histories to create for each probe.
    */
   static int TIME_HISTS_PER_PROBE;

   /**
    * @brief The number of time histories not related to probes.
    */
   static const int GLOBAL_TIME_HISTS;

   /**
    * @brief The number of collision time histories.
    */
   static int NUM_COLL_TIME_HISTS;

   /**
    * @brief Constructor.
    *
    * @param[in] a_pp Where the input parameters will be found.
    * @param[in] a_domain The global problem domain.
    * @param[in] a_num_kinetic_species Number of KineticSpecies in the problem.
    * @param[in] a_solution_order Order of solution.
    * @param[in] a_plot_times_per_file The number of time slices in a field
    *                                  file.
    * @param[in] a_plot_ke_vel_bdy_flux If true plot species' KE vel bdy fluxes.
    * @param[in] a_coll_diag_on If true, collision diagnostics are on.
    * @param[in] a_bz_const Constant Bz for particle tracking.
    */
   Maxwell(
      LokiInputParser& a_pp,
      tbox::Pointer<ProblemDomain> a_domain,
      int a_num_kinetic_species,
      int a_solution_order,
      int a_plot_times_per_file,
      bool a_plot_ke_vel_bdy_flux,
      bool a_coll_diag_on,
      double a_bz_const); //IEO

   /**
    * @brief Copy constructor.
    *
    * @param[in] a_other What is being copies.
    */
   Maxwell(
      const Maxwell& a_other);

   /**
    * @brief Destructor.
    */
   virtual
   ~Maxwell();

   //@{

   //! @name Implementation of Load interface.

   /**
    * @brief Partitions this load over the necessary number of processors.
    *
    * @param[in] a_proc_lo Lowest rank processor distributed over.
    * @param[in] a_proc_hi Highest rank processor distributed over.
    * @param[in] a_comm Corresponding communicator.
    */
   virtual
   void
   createPartition(
      int a_proc_lo,
      int a_proc_hi,
      const MPI_Comm& a_comm);

   //@}

   //@{

   //! @name Implementation of Serializable interface.

   /**
    * @brief Read the electromagnetic fields and transverse drift velocities
    *        from the restart database.
    *
    * @param[in] a_reader The restart file reader.
    */
   virtual
   void
   getFromRestart(
      RestartReader& a_reader);

   /**
    * @brief Write the electromagnetic fields and transverse drift velocities
    *        to the restart database.
    *
    * @param[out] a_writer The restart file writer.
    * @param[in] a_time Current simulation time.
    */
   virtual
   void
   putToRestart(
      RestartWriter& a_writer,
      double a_time);

   //@}

   /**
    * @brief Dynamically creates a new Maxwell with conforming layout and
    *        copies internal parameters.
    */
   tbox::Pointer<Maxwell>
   clone() const
   {
      tbox::Pointer<Maxwell> clone(new Maxwell(*this));
      return clone;
   }

   /**
    * @brief Checks that the rhs and this have a similar storage
    *        representations and returns true if so, false if not.
    *
    * @param[in] a_other Maxwell with which to compare.
    * @param[in] a_include_ghost_cells Include ghost cells in comparison.
    */
   bool
   conformsTo(
      const Maxwell& a_other,
      bool a_include_ghost_cells = true) const;

   /**
    * @brief Zeros out electromagnetic fields, transverse drift velocities, and
    *        tracking particles.
    */
   void
   zeroData()
   {
      m_em_vars = 0.0;
      for (int i = 0; i < m_num_kinetic_species; ++i) {
         m_vz[i] = 0.0;
      }
      zeroParticleData();
   }

   /**
    * @brief Add a multiple of another Maxwell's data to this.
    *
    * @param[in] a_increment_maxwell Maxwell to add.
    * @param[in] a_factor Multiple of input Maxwell to add.
    * @param[in] a_sum_reduce_inc If true, sum reduce the particle increment.
    *                             This is necessary at different points in the
    *                             RK schemes.
    */
   void
   addData(
      const Maxwell& a_increment_maxwell,
      double a_factor,
      bool a_sum_reduce_inc);

   /**
    * @brief Copies electromagnetic fields, transverse drift velocities, and
    *        particles of rhs into this.
    *
    * @param[in] a_rhs Maxwell from which to copy.
    */
   void
   copySolnData(
      const Maxwell& a_rhs);

   /**
    * @brief Compute and return stable time step limit for Maxwell's equations.
    */
   double
   computeDt()
   {
      return 1.0/(Simulation::s_LIGHT_SPEED*
                  (1.0/m_domain->dx(0) + 1.0/m_domain->dx(1)));
   }

   /**
    * @brief Initialize the electromagnetic fields and transverse drift
    *        velocities.  Set up all private data related to plotting and the
    *        Maxwell solver to compute EM fields.
    *
    * @param[in] a_num_species The number of kinetic species.
    * @param[in] a_plot_ke_vel_bdy_flux If true plot species' KE vel bdy fluxes.
    * @param[in] a_from_restart If true this is a restarted run so do not zero
    *                           out the EM vars or transverse drift velocities.
    * @param[in] a_coll_diag_on If true, collision diagnostics are on.
    */
   void
   initialize(
      int a_num_species,
      bool a_plot_ke_vel_bdy_flux,
      bool a_from_restart,
      bool a_coll_diag_on);

   /**
    * @brief Print the class' parameters.
    */
   void
   printParameters() const;

   /**
    * @brief Returns true if the dispersion relation is needed.  Currently that
    *        occurs only if there is a UnidirectionalCurrentDriver.
    */
   bool
   needsDispersionRelation()
   {
      bool needs_disp_rel = false;
      for (int i = 0; i < m_num_current_drivers; ++i) {
         if (m_current_drivers[i]->needsDispersionRelation()) {
            needs_disp_rel = true;
            break;
         }
      }
      return needs_disp_rel;
   }

   /**
    * @brief Compute antenna sources.
    *
    * @param[in] a_time Current simulation time.
    * @param[in] a_omega_eff2 Square of the effective plasma frequency.
    */
   void
   computeAntennaSource(
      double a_time,
      const ParallelArray a_omega_eff2);

   /**
    * @brief Apply current due to noise source particles.
    *
    * @param[out] a_Jx x component of current
    * @param[out] a_Jy y component of current
    * @param[in] a_time Current simulation time.
    */
   void
   applyNoiseParticleCurrent(
      ParallelArray& a_Jx,
      ParallelArray& a_Jy,
      double a_time);

   /**
    * @brief Evaluate the right hand side of Maxwell's equations and the
    *        particles' motion.  We are not enforcing the constraint on div(B)
    *        in this step.
    *
    * @param[out] a_rhs The right hand side.
    * @param[in] a_kspv KineticSpecies needed to get charge and mass.
    * @param[in] a_Jx x component of current
    * @param[in] a_Jy y component of current
    * @param[in] a_Jz z component of current
    * @param[in] a_net_ext_efield The total external efield from all drivers
    *                             from all species.
    * @param[in] a_time Current simulation time.
    */
   void
   evalRHS(
      Maxwell& a_rhs,
      const KineticSpeciesPtrVect& a_kspv,
      const ParallelArray& a_Jx,
      const ParallelArray& a_Jy,
      const ParallelArray& a_Jz,
      const ParallelArray& a_net_ext_efield,
      double a_time);

   /**
    * @brief Writes electromagnetic fields and accumulated sequences to the plot
    *        file.
    *
    * @param[in] a_time Simulation time.
    * @param[in] a_dt Simulation time step.
    * @param[in] a_sequences Sequences to save.
    * @param[in] a_species_names Names of KineticSpecies in the problem.
    * @param[in] a_time_seq Simulation times of sequences.
    * @param[in] a_probes Location of each probe.
    * @param[in] a_num_probes Number of probes.
    * @param[in] a_saved_seq Number of times sequences that have been saved.
    * @param[in] a_saved_save Number of plots that have been saved.
    * @param[in] a_plot_ke_vel_bdy_flux If true plot species KE vel bdy fluxes.
    * @param[in] a_time_hist_file_name Base name for time history files.
    */
   void
   plot(
      double a_time,
      double a_dt,
      const vector<vector<double> >& a_sequences,
      const vector<string>& a_species_names,
      const vector<double>& a_time_seq,
      const vector<vector<double> >& a_probes,
      int a_num_probes,
      int a_saved_seq,
      int& a_saved_save,
      bool a_plot_ke_vel_bdy_flux,
      string& a_time_hist_file_name);

   void
   plotCollision(
      double a_time,
      double a_dt,
      KineticSpeciesPtrVect const& a_kinetic_species,
      const vector<vector<double> >& a_sequences,
      const vector<string>& a_species_names);
  
   /**
    * @brief Accumulates sequences in a_sequences.
    *
    * @param[out] a_sequences Sequences to save.
    * @param[in] a_probes Location of each probe.
    * @param[in] a_num_probes Number of probes.
    * @param[in] a_saved_seq Number of times sequences that have been saved.
    * @param[in,out] a_seq_idx Starting index of sequences written by Maxwell.
    */
   void
   accumulateSequences(
      vector<vector<double> >& a_sequences,
      const vector<vector<double> >& a_probes,
      int a_num_probes,
      int a_saved_seq,
      int& a_seq_idx);

   /**
    * @brief Accessor for the global transverse drift velocity field for a
    *        given KineticSpecies.
    *
    * @param[in] a_species_idx Index of the species of interest.
    */
   const ParallelArray&
   getVZVar(
      int a_species_idx) const
   {
      return m_vz[a_species_idx];
   }

   /**
    * @brief Set ghost boundaries, parallel ghosts and physical boundary
    *        conditions.
    *
    * @param[in] a_particles_only If true only update particle boundary
    *                             conditions.
    */
   void
   fillGhostCells(
      bool a_particles_only)
   {
      // Perform an exchange of parallel ghost data.
      if (!a_particles_only) {
         UGB();
      }

      // Set the appropriate physical boundary conditions.
      setPhysicalBCs(a_particles_only);
   }

   /**
    * @brief Construct the list of time history names for a VM problem.
    *        This could be a member function and use Maxwell's knowledge of
    *        the number of tracking particles but it's useful to use this
    *        elsewhere (post-processor) so it's a static method which takes
    *        this info as input.
    *
    * @param[in] a_num_probes Number of probes in the problem.
    * @param[in] a_num_tracking_particles Number of tracking particles in the
    *                                     problem
    * @param[in] a_species_names Names of the kinetic species in the problem.
    * @param[out] a_time_hist_names The list of time history names.
    */
   static
   void
   buildTimeHistoryNames(
      int a_num_probes,
      int a_num_tracking_particles,
      const vector<string>& a_species_names,
      vector<string>& a_time_hist_names);

   /**
    * @brief Construct the list of collision time history names for a VM
    *        problem.  This could be a member function but it's useful to use
    *        this elsewhere (post-processor) so it's a static method.
    *
    * @param[in] a_species_names Names of the kinetic species in the problem.
    * @param[out] a_time_hist_names The list of time history names.
    */
   static
   void
   buildCollTimeHistoryNames(
      const vector<string>& a_species_names,
      vector<string>& a_time_hist_names);

   /**
    * @brief Construct the list of plot history names for a VM problem.
    *        This could be a member function and use Maxwell's knowledge of
    *        the kinetic species names' but it's useful to use this elsewhere
    *        (post-processor) so it's a static method which takes this info as
    *        input.
    *
    * @param[in] a_plot_ke_vel_bdy_flux If true KE flux on each velocity space
    *                                   boundary is plotted for each species.
    * @param[in] a_species_names Names of the kinetic species in the problem.
    * @param[out] a_plot_names The list of plot names.
    */
   static
   void
   buildPlotNames(
      bool a_plot_ke_vel_bdy_flux,
      const vector<string>& a_species_names,
      vector<string>& a_plot_names);

   /**
    * @brief Construct the list of collision plot names for a VM problem.
    *        This could be a member function and use Maxwell's knowledge of
    *        the kinetic species names' but it's useful to use this elsewhere
    *        (post-processor) so it's a static method which takes this info as
    *        input.
    *
    * @param[in] a_species_names Names of the kinetic species in the problem.
    * @param[out] a_plot_names The list of plot names.
    */
   static
   void
   buildCollPlotNames(
      const vector<string>& a_species_names,
      vector<string>& a_plot_names);

protected:
   // Adds noise to the charge distribution
   virtual
   void
   addNoiseToChargeDist(
      TimerManager* a_timers,
      ParallelArray& a_charge_density,
      double a_time);

   // Adds potential from any external potential drivers.
   virtual
   void
   addExternalPotential(
      ParallelArray& a_phi_local,
      double a_time);

   // Starts the right boundary condition timer.
   virtual
   void
   startBCTimer(
      TimerManager* a_timers);

   // Stops the right boundary condition timer.
   virtual
   void
   stopBCTimer(
      TimerManager* a_timers);

   // Update ghost boundaries.
   void
   UGB();

   // Sets physical boundary conditions.  Boundary condition depends on whether
   // periodicity is or is not on in a physical direction.
   void
   setPhysicalBCs(
      bool a_particles_only);

private:
   // Unimplemented default constructor.
   Maxwell();

   // Unimplemented assignment operator.
   Maxwell&
   operator = (
      const Maxwell& rhs);

   // Obtain control parameters from "maxwell" section of database.
   void
   parseParameters(
      LokiInputParser& a_pp);

   // Does the work of evalRHS for the supplied particles.
   void
   evalRHSParticles(
      bool a_noise_source_particles,
      Maxwell& a_rhs,
      const ParallelArray& a_net_ext_efield,
      double a_time);

   // Apply the initial conditions specified in the input deck to the
   // electromagnetic fields.
   void
   initializeEM()
   {
      for (int i = 0; i < static_cast<int>(m_em_ics.size()); ++i) {
         m_em_ics[i]->set(m_em_vars, *m_domain);
      }
   }

   // Apply the initial conditions specified in the input deck to vz for all
   // species.
   void
   initializeVZ()
   {
      for (int i = 0; i < static_cast<int>(m_vel_ics.size()); ++i) {
         m_vel_ics[i]->set(m_vz[i], *m_domain);
      }
   }

   // If this object is on the Maxwell process and there are either tracking or
   // noise source particles then this method builds the field interpolator
   // appropriate to the spatial order of the problem.
   void
   buildInterpolators();

   // Transverse drift velocity for each KineticSpecies.
   vector<ParallelArray> m_vz;

   // Undocumented artificial viscosity for RHS evaluation.
   double m_avWeak;
   double m_avStrong;

   // Antenna source term.
   ParallelArray m_antenna_source;

   // The current drivers.
   vector<tbox::Pointer<CurrentDriver> > m_current_drivers;

   // The number of current drivers.
   int m_num_current_drivers;

   // The electromagnetic field initial conditions.
   vector<tbox::Pointer<EMICInterface> > m_em_ics;

   // The z velocity initial conditions.
   vector<tbox::Pointer<VELICInterface> > m_vel_ics;

   // EM field data interpolator needed only for particle calculations.
   tbox::Pointer<Interpolator> m_em_interpolator;

   // Current density interpolator needed only for noise particle calculations.
   tbox::Pointer<Interpolator> m_j_interpolator;

   // Interpolation coefficient storage.
   double** m_ex_interp;
   double** m_ey_interp;
   double** m_bz_interp;
   double** m_jx_interp;
   double** m_jy_interp;

   // Coordinates of the lower end of the super grid.
   vector<double> m_supergrid_lo;

   // Coordinates of the upper end of the super grid.
   vector<double> m_supergrid_hi;

   // Tag for baton passing communication.
   static const int s_TAG_BATON;
};

} // end namespace Loki

#endif
