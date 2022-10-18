/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _KINETIC_SPECIES_H_
#define _KINETIC_SPECIES_H_

#include <string>
#include "tbox/Dimension.H"
#include "tbox/Box.H"
#include "tbox/Pointer.H"
#include "Directions.H"
#include "ProblemDomain.H"
#include "KrookLayer.H"
#include "ExternalDistKrook.H"
#include "ICInterface.H"
#include "Simulation.H"
#include "TZSourceInterface.H"
#include "TimerManager.H"
#include "ParallelArray.H"

#include "ContractionSchedule.H"
#include "ReductionSchedule.H"
#include "ExpansionSchedule.H"
#include "SummationSchedule.H"

#include "tbox/DescribedClass.H"
#include "Load.H"
#include "Serializable.H"

#include "ElectricFieldDriver.H"
#include "CollisionOperator.H"

#include "KineticSpeciesF.H"
#include "Loki_Utilities.H"

namespace Loki {

class EMSolverBase;
class Maxwell;
class Poisson;

/**
 * Implementation of abstract Load interface describing a distribution of
 * charged particles in 4 space.
 */
class KineticSpecies 
  : public tbox::DescribedClass, public Load, public Serializable
{
public:
   /**
    * @brief Constructor for a KineticSpecies built from a restart file.
    *        This is used by the post-processor.
    *
    * @param[in] a_reader Restart database reader.
    * @param[in] a_species_num Which species is being constructed.
    * @param[in] a_number_of_species Total number of KineticSpecies.
    * @param[in] a_spatial_solution_order Spatial order of solution.
    * @param[in] a_temporal_solution_order Temporal order of solution.
    * @param[in] a_plot_ke_vel_bdy_flux If true the machinery for computing the
    *                                   kinetic energy flux through the velocity
    *                                   boundaries must be computed.
    * @param[in] a_use_new_bcs Use Jeff Banks' new boundary condition
    *                          formulation.
    * @param[in] a_name Name of this species.
    */
   KineticSpecies(
      RestartReader& a_reader,
      int a_species_num,
      int a_number_of_species,
      int a_spatial_solution_order,
      int a_temporal_solution_order,
      bool a_plot_ke_vel_bdy_flux,
      bool a_use_new_bcs,
      const string& a_name);

   /**
    * @brief Constructor for a KineticSpecies built from an input file.
    *        This is used by the main code.
    *
    * @param[in] a_domain The domain of the problem.
    * @param[in] a_pp Where input parameters will be found.
    * @param[in] a_species_num Which species is being constructed.
    * @param[in] a_number_of_species Total number of KineticSpecies.
    * @param[in] a_spatial_solution_order Spatial order of solution.
    * @param[in] a_temporal_solution_order Temporal order of solution.
    * @param[in] a_plot_ke_vel_bdy_flux If true the machinery for computing the
    *                                   kinetic energy flux through the velocity
    *                                   boundaries must be computed.
    * @param[in] a_use_new_bcs Use Jeff Banks' new boundary condition
    *                          formulation.
    * @param[in] a_do_maxwell If true problem being run is electrodynamic.
    * @param[in] a_bz_const Constant magnetic field in z.
    */
   KineticSpecies(
      const tbox::Pointer<ProblemDomain>& a_domain,
      LokiInputParser& a_pp,
      int a_species_num,
      int a_number_of_species,
      int a_spatial_solution_order,
      int a_temporal_solution_order,
      bool a_plot_ke_vel_bdy_flux,
      bool a_use_new_bcs,
      bool a_do_maxwell,
      double a_bz_const);

   /**
    * @brief Copy constructor.
    *        Public function although it is only called by clone.
    *
    * @param[in] a_other KineticSpecies from which to initialize.
    */
   KineticSpecies(
      const KineticSpecies& a_other);

   /**
    * @brief Destructor.
    */
   virtual ~KineticSpecies();

   //@{

   //! @name Implementation of Load interface.

   /**
    * @brief Computes the computational cost of this object.
    *
    * @return The computational cost of this object.
    */
   virtual
   float
   netCost() const;

   /**
    * @brief Returns the number of processors over which this species is
    *        distributed.
    */
   virtual
   int
   numberOfProcessors() const;

   /**
    * @brief Returns true if this species is distributed over a user specified
    *        fixed number of processors.
    */
   virtual
   bool
   fixedNumberOfProcessors() const;

   /**
    * @brief Partitions this species over the necessary number of processors.
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

   /**
    * @brief Returns true if this KineticSpecies is distributed to processor
    *        with supplied id.
    *
    * @param[in] a_proc_id ID of processor of interest.
    */
   virtual
   bool
   isInRange(
      int a_proc_id) const;

   /**
    * @brief Print information about processors this KineticSpecies is
    *        distributed over.
    */
   virtual
   void
   printDecomposition() const;

   //@}

   //@{

   //! @name Implementation of Serializable interface.

   /**
    * @brief Read distribution function from restart.
    *
    * @param[in] a_reader Restart file reader.
    */
   virtual
   void
   getFromRestart(
      RestartReader& a_reader);

   /**
    * @brief Write distribution function and Krook layer information to
    *        restart.
    *
    * @param[out] a_writer Restart file writer.
    * @param[in] a_time Current simulation time.
    */
   virtual
   void
   putToRestart(
      RestartWriter& a_writer,
      double a_time);

   //@}

   /**
    * @brief Write only distribution function to restart.
    *
    * @param[out] a_db Restart database
    * @param[in] a_time Current simulation time.
    * @param[in] a_return_to_root If true, on completion return to root
    *                             directory in restart file.
    */
   void
   putToRestart_SkipKrook(
      RestartWriter& a_db,
      double a_time,
      bool a_return_to_root = true);

   /**
    * @brief Returns species name.
    */
   const string&
   name() const
   {
      return m_name;
   }

   /**
    * @brief Returns species mass.
    */
   double
   mass() const
   {
      return m_mass;
   }

   /**
    * @brief Returns species charge.
    */
   double
   charge() const
   {
      return m_charge;
   }

   /**
    * @brief Returns species charge density.
    *
    * @param[out] a_rho Computed charge density.
    */
   void
   chargeDensity(
      ParallelArray& a_rho)
   {
      m_moment_schedule->execute(a_rho);
   }

   /**
    * @brief Computes species current density.
    *
    * @param[in] a_maxwell The entity owning the transverse drift velocity.
    * @param[out] a_Jx The x component of the current density.
    * @param[out] a_Jy The y component of the current density.
    * @param[out] a_Jz The z component of the current density.
    */
   void
   currentDensity(
      const Maxwell& a_maxwell,
      ParallelArray& a_Jx,
      ParallelArray& a_Jy,
      ParallelArray& a_Jz);

   /**
    * @brief Returns dimension.
    */
   tbox::Dimension
   dim() const
   {
      return m_pdim;
   }

   /**
    * @brief Zeros out distribution function and, if there are any external
    *        field drivers, the integrated ke_e_dot.
    */
   void
   zeroData()
   {
      m_dist_func = 0.0;
      if (m_num_external_drivers > 0) {
         m_integrated_ke_e_dot = 0.0;
      }
   }

   /**
    * @brief Add a multiple of another KineticSpecies' data to this
    *        KineticSpecies.
    *
    * @param[in] a_rhs KineticSpecies to add.
    * @param[in] a_factor Multiple of input KineticSpecies to add.
    */
   void
   addData(
      const KineticSpecies& a_rhs,
      double a_factor);

   /**
    * @brief Checks that the rhs and this KineticSpecies have a similar storage
    *        representations and returns true if so, false if not.
    *
    * @param[in] a_rhs KineticSpecies with which to compare.
    * @param[in] a_include_ghost_cells Include ghost cells in comparison.
    */
   bool
   conformsTo(
      const KineticSpecies& a_rhs,
      bool a_include_ghost_cells = true) const;

   /**
    * @brief Dynamically creates a new KineticSpecies with conforming layout
    *        and copies internal parameters.
    *
    * @param[in] a_ghost_cells Number of ghost cells in cloned version
    */
   tbox::Pointer<KineticSpecies>
   clone(
      int a_ghost_cells) const
   {
      tbox::Pointer<KineticSpecies> clone(new KineticSpecies(*this));
      clone->m_n_ghosts = a_ghost_cells;
      return clone;
   }

   /**
    * @brief Dynamically creates a new KineticSpecies with conforming layout
    *        and copies internal parameters.
    */
   tbox::Pointer<KineticSpecies>
   clone() const
   {
      tbox::Pointer<KineticSpecies> clone(new KineticSpecies(*this));
      return clone;
   }

   /**
    * @brief Copies only distribution function of rhs into this.
    *
    * @param[in] a_rhs KineticSpecies from which to copy.
    */
   void
   copySolnData(
      const KineticSpecies& a_rhs);

   /**
    * @brief Initialize the distribution function and velocities.
    *
    * @param[in] a_from_restart True if running from a restart file.
    */
   void
   initialize(
      bool a_from_restart)
   {
      if (m_krook_layer->hasKrookLayer()) {
         m_krook_layer->initialize(m_dist_func, m_spatial_solution_order);
      }
      if (m_external_dist_krook->hasKrookLayer()) {
         m_external_dist_krook->initialize(m_dist_func,
            m_spatial_solution_order);
      }
      m_initial_condition->cache(m_dist_func);
      if (!a_from_restart) {
         m_initial_condition->set(m_dist_func);
         initializeVelocity();
      }
   }

   /**
    * @brief Print this class' parameters.
    */
   void
   printParameters() const;

   /**
    * @brief Set ghost boundaries for advection boundaries on the physical
    *        domain.
    *
    * Part of evaluation of rhs to allow interleving of rhs evaluation with
    * electric field computation.
    */
   void
   fillAdvectionGhostCells()
   {
      // Perform an exchange of parallel ghost data.
      UGB();

      // Set the appropriate physical boundary conditions.
      setPhysicalBCs();
   }

   /**
    * @brief Set boundary condition for acceleration boundaries on the physical
    *        domain.
    *
    * Part of evaluation of rhs to allow interleving of rhs evaluation with
    * electric field computation.
    */
   void
   setAccelerationBCs()
   {
      // set boundary conditions for acceleration boundaries
      // (acceleration boundaries of physical domain)
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("BC (Vlasov)");

      if (m_use_new_bcs) {
         FORT_SET_ACCELERATION_BCS_4D_JB(*m_dist_func.getData(),
            BOX4D_TO_FORT(m_global_box),
            BOX4D_TO_FORT(dataBox()),
            BOX4D_TO_FORT(interiorBox()),
            m_spatial_solution_order,
            m_n_ghosts,
            *(m_vel_face[V1].getData()),
            *(m_vel_face[V2].getData()),
            PROBLEMDOMAIN_TO_FORT((*m_domain)),
            int64_t(m_initial_condition.getPointer()));
      }
      else {
         FORT_SET_ACCELERATION_BCS_4D(*m_dist_func.getData(),
            BOX4D_TO_FORT(m_global_box),
            BOX4D_TO_FORT(dataBox()),
            BOX4D_TO_FORT(interiorBox()),
            m_spatial_solution_order,
            *(m_vel_face[V1].getData()),
            *(m_vel_face[V2].getData()),
            int64_t(m_initial_condition.getPointer()));
      }

      timers->stopTimer("BC (Vlasov)");
   }

   /**
    * @brief Compute advection fluxes.
    *
    * Part of evaluation of rhs to allow interleving of rhs evaluation with
    * electric field computation.
    */
   void
   evalAdvectionFluxes()
   {
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("Vlasov");

      FORT_COMPUTE_ADVECTION_FLUXES_4D(*(m_flux[X1].getData()),
         *(m_flux[X2].getData()),
         BOX4D_TO_FORT(dataBox()),
         *(m_vel_face[X1].getData()),
         *(m_vel_face[X2].getData()),
         *(m_u_face[X1].getData()),
         *(m_u_face[X2].getData()),
         *(m_dist_func.getData()),
         m_domain->dx()[0],
         m_spatial_solution_order);

      timers->stopTimer("Vlasov");
   }

   /**
    * @brief Compute advection derivatives.
    *
    * @param[out] a_rhs The right hand side advection derivatives computed for.
    *
    * Part of evaluation of rhs to allow interleving of rhs evaluation with
    * electric field computation.
    */
   void
   evalAdvectionDerivatives(
      KineticSpecies& a_rhs ) const
   {
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("Vlasov");

      FORT_COMPUTE_ADVECTION_DERIVATIVES_4D(
         *(a_rhs.m_dist_func.getData()),
         *(m_dist_func.getData()),
         BOX4D_TO_FORT(dataBox()),
         BOX4D_TO_FORT(interiorBox()),
         *m_vel_face[X1].getData(),
         *m_vel_face[X2].getData(),
         m_domain->dx()[0],
         m_spatial_solution_order);

      timers->stopTimer("Vlasov");
   }

   /**
    * @brief Compute acceleration fluxes.
    *
    * Part of evaluation of rhs to allow interleving of rhs evaluation with
    * electric field computation.
    */
   void
   evalAccelerationFluxes()
   {
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("Vlasov");

      FORT_COMPUTE_ACCELERATION_FLUXES_4D(*(m_flux[V1].getData()),
         *(m_flux[V2].getData()),
         BOX4D_TO_FORT(dataBox()),
         *(m_vel_face[V1].getData()),
         *(m_vel_face[V2].getData()),
         *(m_u_face[V1].getData()),
         *(m_u_face[V2].getData()),
         *(m_dist_func.getData()),
         m_domain->dx()[0],
         m_spatial_solution_order);

      timers->stopTimer("Vlasov");
   }

   /**
    * @brief Compute advection derivatives.
    *
    * @param[out] a_rhs The right hand side acceleration derivatives computed
    * for.
    *
    * Part of evaluation of rhs to allow interleving of rhs evaluation with
    * electric field computation.
    */
   void
   evalAccelerationDerivatives(
      KineticSpecies& a_rhs ) const
   {
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("Vlasov");

      FORT_COMPUTE_ACCELERATION_DERIVATIVES_4D(
         *(a_rhs.m_dist_func.getData()),
         *(m_dist_func.getData()),
         BOX4D_TO_FORT(dataBox()),
         BOX4D_TO_FORT(interiorBox()),
         *(m_vel_face[V1].getData()),
         *(m_vel_face[V2].getData()),
         m_domain->dx()[0],
         m_spatial_solution_order);

      timers->stopTimer("Vlasov");
   }

   /**
    * @brief Adds flux divergence to rhs, adds in collision and krook layer
    *        operators.  Also adds in any ExternalDistKrookLayer.
    *
    * Part of evaluation of rhs to allow interleving of rhs evaluation with
    * electric field computation.
    *
    * @param[out] a_rhs The right hand side.
    * @param[in] a_time Current simulation time.
    * @param[in] a_dt Current simulation time step.
    * @param[in] a_use_new_alg If true use new time saving formulation.
    * @param[in] a_last_rk_stage True if called from the last RK stage.
    */
   void
   completeRHS(
      KineticSpecies& a_rhs,
      double a_time,
      double a_dt,
      bool a_use_new_alg,
      bool a_last_rk_stage) const;

   void
   copyCollisionDiagnostics(
      const int a_collOperIndex,
      vector<ParallelArray>& a_diags);

   void
   copyMomentDiagnostics(
      vector<ParallelArray>& a_diags);

   /**
    * @brief Compute and return stable time step limit for this KineticSpecies.
    */
   double
   computeDt();

   /**
    * @brief Compute a and y accelerations on the faces perpendicular to vx and
    *        vy respectively for electrostatic case.
    *
    * Includes evaluation of externally applied electric fields.
    *
    * @param[in] a_poisson Object holding electrodynamic fields.
    * @param[in] a_ext_efield The electric field from all external drivers from
    *                         this species.
    * @param[in] a_time The time at which velocity is computed.
    * @param[in] a_dt The time step at which velocity is computed.
    * @param[in] a_first_rk_stage True if called from the first RK stage.
    */
   void
   computeAcceleration(
      Poisson& a_poisson,
      ParallelArray& a_ext_efield,
      double a_time,
      double a_dt,
      bool a_first_rk_stage);

   /**
    * @brief Compute a and y accelerations on the faces perpendicular to vx and
    *        vy respectively for electrodynamic case.
    *
    * Includes evaluation of externally applied electric fields.
    *
    * @param[in] a_maxwell Object holding electrodynamic fields.
    * @param[in] a_ext_efield The electric field from all external drivers from
    *                         this species.
    * @param[in] a_time The time at which velocity is computed.
    * @param[in] a_dt The time step at which velocity is computed.
    * @param[in] a_first_rk_stage True if called from the first RK stage.
    */
   void
   computeAcceleration(
      Maxwell& a_maxwell,
      ParallelArray& a_ext_efield,
      double a_time,
      double a_dt,
      bool a_first_rk_stage);

   /**
    * @brief Return the global number of cells.
    */
   long int
   numberOfCells() const
   {
      return m_global_box.size();
   }

   /**
    * @brief Compute the kinetic energy flux through the specified velocity
    *        space boundary.
    *
    * @param[in] a_em_solver EMSolver object that contains velocity boundary
    *                        flux arrays.
    */
   void
   computeKEVelBdyFlux(
      EMSolverBase& a_em_solver);

   /**
    * @brief Accumulates sequences in a_sequences for a VPSystem.
    *
    * @param[in] a_sequences Sequences to save.
    * @param[in] a_saved_seq Number of times sequences that have been saved.
    * @param[in,out] a_seq_idx Index of next sequence member being accumulated.
    * @param[in] a_time  Current simulation time.
    */
   void
   accumulateSequences(
      vector<vector<double> >& a_sequences,
      int a_saved_seq,
      int& a_seq_idx,
      double a_time);

   /**
    * @brief Accumulates sequences in a_sequences for a VMSystem.
    *
    * @param[in] a_maxwell The entity owning the transverse drift velocity.
    * @param[in] a_sequences Sequences to save.
    * @param[in] a_saved_seq Number of times sequences that have been saved.
    * @param[in,out] a_seq_idx Index of next sequence member being accumulated.
    * @param[in] a_time  Current simulation time.
    */
   void
   accumulateSequences(
      const Maxwell& a_maxwell,
      vector<vector<double> >& a_sequences,
      int a_saved_seq,
      int& a_seq_idx,
      double a_time);

   void
   accumulateMomentSequences(
      bool a_is_2d_proc,
      vector<vector<double> >& a_sequences,
      int a_saved_seq,
      int& a_seq_idx,
      ParallelArray& a_momx_global,
      ParallelArray& a_momy_global,
      ParallelArray& a_ke_global,
      ParallelArray& a_ent_global);

   void
   accumulateCollisionSequences(
      double a_dt,
      int a_coll_op_idx,
      vector<vector<double> >& a_sequences,
      int a_saved_seq,
      int& a_seq_idx);

   /**
    * @brief Indicates that the problem has particles.
    */
   void
   problemHasParticles()
   {
      m_problem_has_particles = true;
      defineExtEfieldContractionSchedule();
   }

   /**
    * @brief Returns the lowest rank processor this species is distributed over.
    */
   int
   procLo() const
   {
      return m_proc_lo;
   }

   /**
    * @brief Returns the highest rank processor this species is distributed
    *        over.
    */
   int
   procHi() const
   {
      return m_proc_hi;
   }

   /**
    * @brief Returns the MPI communicator for this species.
    */
   const MPI_Comm&
   communicator() const
   {
      return m_comm;
   }

   /**
    * @brief Returns the box containing the interior and ghost/boundary layer
    *        for this species on this processor.
    */
   const ParallelArray::Box&
   dataBox() const
   {
      return m_dist_func.dataBox();
   }

   /**
    * @brief Returns the interior box for this species on this processor.
    */
   const ParallelArray::Box&
   interiorBox() const
   {
      return m_dist_func.interiorBox();
   }

   /**
    * @brief Returns the global box for the problem.
    */
   const tbox::Box&
   globalBox() const
   {
     return m_global_box;
   }

   /**
    * @brief Returns the number of ghosts cells in each direction.
    */
   int
   numGhosts() const
   {
      return m_n_ghosts;
   }

   /**
    * @brief Returns the array of velocities for the velocity space of this
    *        species on this processor.
    */
   const ParallelArray&
   velocities() const
   {
      return m_velocities;
   }

   /**
    * @brief Returns the masses of all the species in the problem.
    */
   const vector<double>&
   speciesMasses() const
   {
      return m_speciesMass;
   }

   /**
    * @brief Returns the head processors of all the species in the problem.
    */
   const vector<int>&
   speciesHeads() const
   {
      return m_speciesHeads;
   }

   /**
    * @brief Returns the index of this species.
    */
   int
   speciesIndex() const
   {
      return m_species_index;
   }

   /**
    * @brief Returns the number of species in the problem.
    */
   int
   numSpecies() const
   {
      return m_number_of_species;
   }

   /**
    * @brief Returns the box defining the problem domain.
    */
   const ProblemDomain&
   domain() const
   {
     return *m_domain;
   }

   /**
    * @brief Returns the box defining the problem domain.
    */
   const tbox::Box&
   domainBox() const
   {
     return m_domain->box();
   }

   /**
    * @brief Returns the spatial solution order.
    */
   int
   spatialSolutionOrder() const
   {
      return m_spatial_solution_order;
   }

   /**
    * @brief Returns the initial x flow velocity.
    */
   double
   vflowinitx() const
   {
      return m_initial_condition->vflowinitx();
   }

   /**
    * @brief Returns the initial y flow velocity.
    */
   double
   vflowinity() const
   {
      return m_initial_condition->vflowinity();
   }

   /**
    * @brief Returns the distribution function.
    */
   ParallelArray&
   distribution()
   {
      return m_dist_func;
   }

   /**
    * @brief Returns the distribution function.
    */
   const ParallelArray&
   distribution() const
   {
      return m_dist_func;
   }

   int
   headRank() const
   {
      return m_proc_lo;
   }

   void
   SetSpeciesHeadsAndMass(
      const int numSpecies,
      const int *heads,
      const double *mass);

   const vector<double>&
   GetGlobalDomainLower()
   {
      return m_domain->lower();
   }

   const vector<double>&
   GetGlobalDomainUpper()
   {
      return m_domain->upper();
   }

private:
   // Unimplemented default constructor.
   KineticSpecies();

   // Unimplemented assignment operator.
   KineticSpecies&
   operator = (
      const KineticSpecies&);

   // Parse the input parameters for this KineticSpecies.
   void
   parseParameters(
      LokiInputParser& a_pp,
      const tbox::Pointer<ProblemDomain> a_cfg_domain);

   // Copies internals of rhs into this.
   void
   copy(
      const KineticSpecies& a_rhs);

   // Allocate private array data members and associated info.
   void
   allocateAuxArrays();

   // Initialize face velocities.
   void
   initializeVelocity();

   // Enforce periodicity of distribution function in period directions.
   void
   setPeriodicBCs();

   // Adds flux divergence to RHS.
   void
   addFluxDivergence(
      KineticSpecies& a_rhs) const;

   // Defines external efield contraction schedule.
   void
   defineExtEfieldContractionSchedule();

   // Defines charge density reduction schedule.
   void
   defineChargeDensityReductionSchedule();

   // Defines kinetic energy summation schedule.
   void
   defineKineticEnergySummationSchedule();

   // Defines reduction schedule for collision operator momentum, kinetic
   // energy, and entropy diagnostics.
   void
   defineMomentReductionSchedule();

   // Defines collision operator diagnostic field summation schedule.
   void
   defineDiagnosticSummationSchedule();

   // Defines current density reduction schedules when electrodynamic problem
   // is run.
   void
   defineCurrentDensityReductionSchedules();

   // Set boundary conditions for advection boundaries on the physical domain.
   void
   UGB()
   {
      // set boundary conditions for advection boundaries (boundaries of
      // physical domain)
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("parallel ghost");
      m_dist_func.communicateGhostData();
      timers->stopTimer("parallel ghost");
   }

   // Sets physical boundary conditions.  Boundary condition depends on whether
   // periodicity is or is not on in a physical direction.
   void
   setPhysicalBCs()
   {
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("BC (Vlasov)");
      if (m_use_new_bcs) {
         FORT_SET_ADVECTION_BCS_4D_JB(*m_dist_func.getData(),
            BOX4D_TO_FORT(dataBox()),
            BOX4D_TO_FORT(interiorBox()),
            m_spatial_solution_order,
            m_n_ghosts,
            *m_vel_face[X1].getData(),
            *m_vel_face[X2].getData(),
            PROBLEMDOMAIN_TO_FORT((*m_domain)),
            m_domain->isPeriodicAsInt(X1),
            m_domain->isPeriodicAsInt(X2),
            int64_t(m_initial_condition.getPointer()));
      }
      else {
         FORT_SET_ADVECTION_BCS_4D(*m_dist_func.getData(),
            BOX4D_TO_FORT(m_global_box),
            BOX4D_TO_FORT(dataBox()),
            BOX4D_TO_FORT(interiorBox()),
            m_spatial_solution_order,
            *m_vel_face[X1].getData(),
            *m_vel_face[X2].getData(),
            m_domain->isPeriodicAsInt(X1),
            m_domain->isPeriodicAsInt(X2),
            int64_t(m_initial_condition.getPointer()));
      }

      setPeriodicBCs();
      timers->stopTimer("BC (Vlasov)");
   }

   // Build arrays on the "velocity" space containing velocities.  Most useful
   // for relativistic problems.
   void
   buildVelocityArrays();

   // Part of accumulateSequences common to both VP and VM systems.
   void
   accumulateSequencesCommon(
      vector<vector<double> >& a_sequences,
      int a_saved_seq,
      int& a_seq_idx,
      double a_time);

   // Phase space dimension of this KineticSpecies.
   const tbox::Dimension m_pdim;

   // Configuration space dimension of this KineticSpecies.
   const tbox::Dimension m_cdim;

   // Name of this KineticSpecies.
   string m_name;

   // This KineticSpecies' ID.
   int m_species_index;

   // Total number of KineticSpecies.
   int m_number_of_species;

   // Mass of this KineticSpecies.
   double m_mass;

   // Charge of this KineticSpecies.
   double m_charge;

   // The global problem domain.
   tbox::Pointer<ProblemDomain> m_domain;

   // The distribution function for this species.
   ParallelArray m_dist_func;

   // The x, y, and z components of the current.
   ParallelArray m_Jx;
   ParallelArray m_Jy;
   ParallelArray m_Jz;

   // The kinetic energy, momentum components, entropy
   ParallelArray m_ke;
   ParallelArray m_momx;
   ParallelArray m_momy;
   ParallelArray m_ent;

   // The kinetic energy velocity boundary flux.
   ParallelArray m_ke_vel_bdry_flux;

   // Spatial order of solution
   int m_spatial_solution_order;

   // Temporal order of solution.
   int m_temporal_solution_order;

   // The number of ghost cells in each direction.
   int m_n_ghosts;

   // Algorithm stencil width.
   int m_stencil_width;

   // Box describing the global extent including ghosts of this KineticSpecies.
   tbox::Box m_global_box;

   // How to initialize the distribution function.
   tbox::Pointer<ICInterface> m_initial_condition;

   // Any Krook layers in the problem.
   tbox::Pointer<KrookLayer> m_krook_layer;

   // Any ExternalDistKrookLayer in the problem.
   tbox::Pointer<ExternalDistKrookLayer> m_external_dist_krook;

   // The number of processors this KineticSpecies is distributed over.
   int m_number_of_procs;

   // If true the user has specified the number of processors.
   bool m_fixed_number_of_procs;

   // The lowest rank processor this KineticSpecies is distributed over.
   int m_proc_lo;

   // The highest rank processor this KineticSpecies is distributed over.
   int m_proc_hi;

   // If true, the distribution function has been partitioned.
   bool m_partition_defined;

   // The communicator for the processors this KineticSpecies is distributed
   // over.
   MPI_Comm m_comm;

   // The external electric field drivers.
   vector<tbox::Pointer<ElectricFieldDriver> > m_ef_drivers;

   // Number of external electric field drivers.
   int m_num_external_drivers;

   // If true, external electric field drivers use old syntax.
   bool m_old_driver_syntax;

   // Contraction schedule for the external efield.
   tbox::Pointer<ContractionSchedule> m_ext_efield_schedule;

   // Reduction schedule for the charge density.
   tbox::Pointer<ReductionSchedule> m_moment_schedule;

   // Expansion schedule for configuration space electric field.
   tbox::Pointer<ExpansionSchedule> m_efield_expansion_schedule;

   // Expansion schedule for configuration space electromagnetic fields.
   tbox::Pointer<ExpansionSchedule> m_em_expansion_schedule;

   // Expansion schedule for configuration space transverse drift velocity
   // field.
   tbox::Pointer<ExpansionSchedule> m_vz_expansion_schedule;

   // Reduction schedules for the x, y, and z components of the current.
   tbox::Pointer<ReductionSchedule> m_jx_schedule;
   tbox::Pointer<ReductionSchedule> m_jy_schedule;
   tbox::Pointer<ReductionSchedule> m_jz_schedule;

   // Reduction schedules for the momentum components, kinetic energy and
   // entropy collision diagnostics.
   tbox::Pointer<ReductionSchedule> m_ke_schedule;
   tbox::Pointer<ReductionSchedule> m_momy_schedule;
   tbox::Pointer<ReductionSchedule> m_momx_schedule;
   tbox::Pointer<ReductionSchedule> m_ent_schedule;

   // Reduction schedule for the kinetic energy velocity boundary flux.
   tbox::Pointer<SummationSchedule> m_ke_vel_bdry_flux_schedule;

   // Reduction schedule for diagnostic fields.
   tbox::Pointer<SummationSchedule> m_diagnostic_schedule;

   // Acceleration.
   ParallelArray m_accel;

   // Local external electric field.
   ParallelArray m_ext_efield;

   // Electromagnetic variables.
   ParallelArray m_em_vars;

   // Transverse drift velocity.
   ParallelArray m_vz;

   // Maximum x and y velocities/accelerations.
   vector<double> m_lambda_max;

   // For each phase space dimension, the distribution function on those faces.
   vector<ParallelArray> m_u_face;

   // Something of a misnomer as its meaning is different for configuration
   // space dimensions and velocity space dimensions.  For the configuration
   // space dimensions, X1 and X2, this is a zone centered velocity in that
   // dimension for all the faces perpendicular to that dimension.  Since these
   // are x and y faces, the velocities are only dependent on vx and vy.  Hence
   // m_vel_face[X1] is the x velocity for the velocity cell that each face
   // perpendicular to x appears in.  This means that for dimensions X1 and X2
   // m_vel_face for a given i3 and i4 index, m_vel_face[X1/X2] are the same as
   // corresponding velocities in m_velocities for the same velocity space
   // indices.  (Recall though that m_velocities contain both x and y
   // velocities while m_vel_face contains only 1.)  This make for a bit of
   // ambiguity and the decision of whether to use m_vel_face or
   // m_velocities is somewhat arbitrary for X1/X2.  For the velocity
   // space dimensions, V1 and V2, this is the acceleration in that dimension
   // at the faces perpendicular to that dimension.  Hence m_vel_face[V1] is the
   // x acceleration at faces perpendicular to vx.  For these 2 dimensions
   // m_vel_face does depend on all 4 indices as the EM forces are spatially
   // dependent and v is clearly velocity space dependent.  There's really no
   // ambiguity about these 2 dimensions.  The reason m_vel_face is 4D for all 4
   // dimensions is that both advection and acceleration use the same undelying
   // WENO scheme which must assume dependence on all 4 dimensions for the
   // acceleration calculations.
   vector<ParallelArray> m_vel_face;

   // For each phase space dimension, the fluxes on those faces.
   vector<ParallelArray> m_flux;

   // Array of collision operator diagnostics.
   vector<ParallelArray> m_diagnostics;

   // If true then an electrodynamic problem is being run.
   bool m_do_maxwell;

   // The collision operators.
   vector<tbox::Pointer<CollisionOperator> > m_collision_operators;

   // Number of collision drivers.
   int m_num_collision_operators;

   // Twilight zone source.
   tbox::Pointer<TZSourceInterface> m_tz_source;

   // True if the problem contains particles.
   bool m_problem_has_particles;

   // Data for interspecies communication.
   vector<int> m_speciesHeads;
   vector<double> m_speciesMass;

   // If true use Jeff Banks' new acceleration and advection BC treatment.
   bool m_use_new_bcs;

   // If true the machinery for computing the kinetic energy flux through the
   // velocity boundaries must be computed.
   bool m_plot_ke_vel_bdy_flux;

   // Constant magnetic field in z
   double m_bz_const; //IEO

   // For each "velocity cell", the velocities at the "velocity cell"
   // center.  It's easy to compute these when needed in a non-relativistic
   // setting.  However, in a relativistic setting that computation is
   // expensive.  Instead, in both contexts we compute the "velocity cell"
   // centered velocities once and look them up wherever needed.
   ParallelArray m_velocities;

   // For each "velocity cell's" face perpendicular to a direction, the face
   // centered velocities (both vx and vy) at that face.  Again this is
   // easy to compute when needed in a non-relativistic setting.  However, in a
   // relativistic setting that computation is expensive.  Instead, in both
   // contexts we compute these face centered velocities once and look them up
   // wherever needed.
   ParallelArray m_vxface_velocities;
   ParallelArray m_vyface_velocities;

   // The ke_e_dot time history integrated over time.
   double m_integrated_ke_e_dot;
};

} // end namespace Loki

#endif
