/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _POISSON_H_
#define _POISSON_H_

#include "EMSolverBase.H"
#include "KineticSpeciesPtrVect.H"
#include "Interpolator.H"

namespace Loki {

class ElectricPotentialDriver;

/**
 * Implementation of abstract Load interface describing computation of the
 * electric field from a charge distribution via Poisson's equation.
 */
class Poisson : public EMSolverBase
{
public:
   enum EM_Var_Dims{EX, EY, NUM_EM_VARS};

   /**
    * @brief The number of time histories to create for each probe.
    */
   static const int TIME_HISTS_PER_PROBE;

   /**
    * @brief The number of time histories not related to probes.
    */
   static const int GLOBAL_TIME_HISTS;

   /**
    * @brief The number of collision time histories.
    */
   static int NUM_COLL_TIME_HISTS;

   /**
    * @brief Constructor for a Poisson build from an input file.
    *
    * @param[in] a_pp Where the input parameters will be found.
    * @param[in] a_domain The global problem domain.
    * @param[in] a_num_kinetic_species Number of kinetic species in the problem.
    * @param[in] a_solution_order Order of solution.
    * @param[in] a_plot_times_per_file The number of time slices in a field
    *                                  file.
    * @param[in] a_plot_ke_vel_bdy_flux If true plot species' KE vel bdy fluxes.
    * @param[in] a_coll_diag_on If true, collision diagnostics are on.
    * @param[in] a_bz_const Constant Bz for particle tracking.
    */
   Poisson(
      LokiInputParser& a_pp,
      const tbox::Pointer<ProblemDomain>& a_domain,
      int a_num_kinetic_species,
      int a_solution_order,
      int a_plot_times_per_file,
      bool a_plot_ke_vel_bdy_flux,
      bool a_coll_diag_on,
      double a_bz_const); //IEO

   /**
    * @brief Copy constructor.
    *        Public function although it is only called from clone.
    *
    * @param[in] a_other The Poisson to copy from.
    */
   Poisson(
      const Poisson& a_other);

   /**
    * @brief Destructor.
    */
   virtual
   ~Poisson();

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
    * @brief Dynamically creates a new Poisson with conforming layout and
    *        copies internal parameters.
    */
   tbox::Pointer<Poisson>
   clone() const
   {
      tbox::Pointer<Poisson> clone(new Poisson(*this));
      return clone;
   }

   /**
    * @brief Checks that the rhs and this have a similar storage
    *        representations and returns true if so, false if not.
    *
    * @param[in] a_other Poisson with which to compare.
    */
   bool
   conformsTo(
      const Poisson& a_other) const;

   /**
    * @brief Zeros out tracking and noise source particles.
    */
   void
   zeroData()
   {
      zeroParticleData();
   }

   /**
    * @brief Add a multiple of another Poisson's data to this.
    *
    * @param[in] a_increment_poisson Poisson to add.
    * @param[in] a_factor Multiple of input Poisson to add.
    * @param[in] a_sum_reduce_inc If true, sum reduce the particle increment.
    *                             This is necessary at different points in the
    *                             RK schemes.
    */
   void
   addData(
      const Poisson& a_increment_poisson,
      double a_factor,
      bool a_sum_reduce_inc);

   /**
    * @brief Copies particles of rhs into this.
    *
    * @param[in] a_rhs Poisson from which to copy.
    */
   void
   copySolnData(
      const Poisson& a_rhs);

   /**
    * @brief Set up all private data related to plotting and the Poisson solver
    *        to compute the intrinsic electric field resulting from the charge
    *        distribution.
    *
    * @param[in] a_pp Where the input parameters will be found.
    * @param[in] a_num_species The number of kinetic species.
    * @param[in] a_plot_ke_vel_bdy_flux If true plot species' KE vel bdy fluxes.
    * @param[in] a_coll_diag_on If true, collision diagnostics are on.
    */
   void
   initialize(
      LokiInputParser& a_pp,
      int a_num_species,
      bool a_plot_ke_vel_bdy_flux,
      bool a_coll_diag_on);

   /**
    * @brief Print the class' parameters.
    */
   void
   printParameters() const;

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
      setParticleBCs();
   }

   /**
    * @brief Evaluate the right hand side of the particles' motion.
    *
    * @param[out] a_rhs The right hand side.
    * @param[in] a_net_ext_efield The total external efield from all drivers
    *                             from all species.
    * @param[in] a_time Current simulation time.
    */
   void
   evalRHS(
      Poisson& a_rhs,
      const ParallelArray& a_net_ext_efield,
      double a_time);

   /**
    * @brief Writes electrostatic fields and accumulated sequences to the plot
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
      const vector<string>& a_species_names,
      const vector<double>& a_time_seq,
      int a_saved_seq,
      int& a_saved_save,
      string& a_time_hist_file_name);
  
   /**
    * @brief Accumulates sequences in a_sequences.
    *
    * @param[out] a_sequences Sequences to save.
    * @param[in] a_probes Location of each probe.
    * @param[in] a_num_probes Number of probes.
    * @param[in] a_saved_seq Number of times sequences that have been saved.
    * @param[in,out] a_seq_idx Starting index of sequences written by Poisson.
    */
   void
   accumulateSequences(
      vector<vector<double> >& a_sequences,
      const vector<vector<double> >& a_probes, 
      int a_num_probes,
      int a_saved_seq,
      int& a_seq_idx);

   void
   accumulateCollisionSequences(
      vector<vector<double> >& a_sequences,
      int a_saved_seq,
      int& a_seq_idx);

   /**
    * @brief Construct the list of time history names for a VP problem.
    *        This could be a member function and use Poisson's knowledge of
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
    * @brief Construct the list of collision time history names for a VP
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
    * @brief Construct the list of plot names for a VP problem.
    *        This could be a member function and use Poisson's knowledge of
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
    * @brief Construct the list of collision plot names for a VP problem.
    *        This could be a member function and use Poisson's knowledge of
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

private:
   // Unimplemented default constructor.
   Poisson();

   // Unimplemented assignment operator.
   Poisson&
   operator = (
      const Poisson& rhs);

   // Obtain control parameters from "poisson" section of database.
   void
   parseParameters(
      LokiInputParser& a_pp);

   // Does work of evalRHS for the supplied particles.
   void
   evalRHSParticles(
      bool a_noise_source_particles,
      Poisson& a_rhs,
      const ParallelArray& a_net_ext_efield,
      double a_time);

   // If this object is on the Poisson process and there are either tracking or
   // noise source particles then this method builds the field interpolators
   // appropriate to the spatial order of the problem.
   void
   buildInterpolators();

   // The external potential driver.
   tbox::Pointer<ElectricPotentialDriver> m_ep_driver;

   // True if an external potential is to be applied.
   bool m_apply_external_potential;

   // EM field interpolator needed only for particle calculations.
   tbox::Pointer<Interpolator> m_em_interpolator;

   // Charge density interpolator needed only for noise particle calculations.
   tbox::Pointer<Interpolator> m_rho_interpolator;

   // Interpolation coefficient storage.
   double** m_ex_interp;
   double** m_ey_interp;
   double** m_rho_interp;

   // Tag for baton passing communication.
   static const int s_TAG_BATON;
};

} // end namespace Loki

#endif
