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
#ifndef _KINETIC_SPECIES_H_
#define _KINETIC_SPECIES_H_

#include <string>
#include "Overture.h"
#include "tbox/Dimension.H"
#include "tbox/Box.H"
#include "tbox/DescribedClass.H"
#include "tbox/Pointer.H"
#include "Directions.H"
#include "ProblemDomain.H"
#include "KrookLayer.H"
#include "ICInterface.H"
#include "TZSourceInterface.H"
#include "Load.H"
#include "TimerManager.H"

#include "SumOp.H"
#include "Reduction.H"
#include "ContractionSchedule.H"
#include "ReductionSchedule.H"
#include "ExpansionSchedule.H"
#include "SummationSchedule.H"

#include "Serializable.H"

#include "ElectricFieldDriver.H"
#include "CollisionOperator.H"

#include "KineticSpeciesF.H"
#include "Loki_Utilities.H"

namespace Loki {

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
    * @param[in] a_db Restart database.
    * @param[in] a_species_num Which species is being constructed.
    * @param[in] a_number_of_species Total number of KineticSpecies.
    * @param[in] a_spatial_solution_order Spatial order of solution.
    * @param[in] a_temporal_solution_order Temporal order of solution.
    * @param[in] a_use_new_bcs Use Jeff Banks' new boundray condition
    *                          formulation.
    * @param[in] a_name Name of this species.
    */
   KineticSpecies(
      HDF_DataBase& a_db,
      int a_species_num,
      int a_number_of_species,
      int a_spatial_solution_order,
      int a_temporal_solution_order,
      bool a_use_new_bcs,
      const aString& a_name);

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
    * @param[in] a_use_new_bcs Use Jeff Banks' new boundray condition
    *                          formulation.
    * @param[in] a_do_maxwell If true problem being run is electrodynamic.
    */
   KineticSpecies(
      const tbox::Pointer<ProblemDomain>& a_domain,
      ParmParse& a_pp,
      int a_species_num,
      int a_number_of_species,
      int a_spatial_solution_order,
      int a_temporal_solution_order,
      bool a_use_new_bcs,
      bool a_do_maxwell);

   /**
    * @brief Copy constructor.
    *        Public function although it is only called by clone.
    *
    * @param[in] a_other KineticSpecies from which to initialize.
    * @param[in] a_deep_copy If true copy the distribution function.
    */
   KineticSpecies(
      const KineticSpecies& a_other,
      bool a_deep_copy = true);

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
    * distributed.
    */
   virtual
   int
   numberOfProcessors() const;

   /**
    * @brief Returns true if this species is distributed over a user specified
    * fixed number of processors.
    */
   virtual
   bool
   fixedNumberOfProcessors() const;

   /**
    * @brief Partitions this species over the necessary number of processors.
    *
    * @param[in] a_range Processor range distributed over.
    * @param[in] a_comm Corresponding communicator.
    */
   virtual
   void
   createPartition(
      const Range& a_range,
      const MPI_Comm& a_comm);

   /**
    * @brief Returns true if this KineticSpecies is distributed to processor
    * with supplied id.
    *
    * @param[in] a_proc_id ID of processor of interest.
    */
   virtual
   bool
   isInRange(
      int a_proc_id) const;

   /**
    * @brief Print information about processors this KineticSpecies is
    * distributed over.
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
    * @param[in] a_db Restart database
    */
   virtual
   void
   getFromRestart(
      const HDF_DataBase& a_db);

   /**
    * @brief Write distribution function and Krook layer information to
    * restart.
    *
    * @param[out] a_db Restart database
    * @param[in] a_time Current simulation time.
    */
   virtual
   void
   putToRestart(
      HDF_DataBase& a_db,
      real a_time);

   //@}

   /**
    * @brief Write only distribution function to restart.
    *
    * @param[out] a_db Restart database
    * @param[in] a_time Current simulation time.
    */
   void
   putToRestart_SkipKrook(
      HDF_DataBase& a_db,
      real a_time);

   /**
    * @brief Returns species name.
    */
   const std::string&
   name() const
   {
      return m_name;
   }

   /**
    * @brief Returns species mass.
    */
   real
   mass() const
   {
      return m_mass;
   }

   /**
    * @brief Returns species charge.
    */
   real
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
      realArray& a_rho) const
   {
      m_moment_schedule->execute(a_rho, DensityKernel(m_charge));
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
      realArray& a_Jx,
      realArray& a_Jy,
      realArray& a_Jz);

   /**
    * @brief Returns dimension.
    */
   tbox::Dimension
   dim() const
   {
      return m_pdim;
   }

   /**
    * @brief Zeros out distribution function.
    */
   void
   zeroData()
   {
      m_dist_func_local = 0.0;
   }

   /**
    * @brief Add a multiple of another KineticSpecies' data to this
    * KineticSpecies.
    *
    * @param[in] a_rhs KineticSpecies to add.
    * @param[in] a_factor Multiple of input KineticSpecies to add.
    */
   void
   addData(
      const KineticSpecies& a_rhs,
      real                  a_factor = 1.0);

   /**
    * @brief Checks that the rhs and this KineticSpecies have a similar storage
    * representations and returns true if so, false if not.
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
    * and copies internal parameters.
    *
    * @param[in] a_ghost_cells IntVector of ghost cells in cloned version
    * @param[in] a_deep_copy If true copy the distribution function.
    */
   tbox::Pointer<KineticSpecies>
   clone(
      const tbox::IntVector& a_ghost_cells,
      bool a_deep_copy = true) const
   {
      tbox::Pointer<KineticSpecies> clone(
         new KineticSpecies(*this, a_deep_copy));
      clone->m_n_ghosts = a_ghost_cells;
      return clone;
   }

   /**
    * @brief Dynamically creates a new KineticSpecies with conforming layout
    * and copies internal parameters.
    *
    * @param[in] a_deep_copy If true copy the distribution function.
    */
   tbox::Pointer<KineticSpecies>
   clone(
      bool a_deep_copy = true) const
   {
      tbox::Pointer<KineticSpecies> clone(
         new KineticSpecies(*this, a_deep_copy));
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
    * @brief Initialize the distribution function and velocities at the
    * supplied time.
    *
    * @param[in] a_time Simulation time to initialize.
    */
   void
   initialize(
      real a_time)
   {
      tbox::Box grown_global_box(m_global_box);
      grown_global_box.grow(m_n_ghosts);
      m_initial_condition->set(m_dist_func_local,
                               *m_domain,
                               grown_global_box,
                               a_time);
      initializeVelocity();
   }

   /**
    * @brief Set ghost boundaries for advection boundaries on the physical
    * domain.
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
    * @brief Set ghost boundaries for advection acceleration boundaries on the
    * physical domain.
    *
    * Part of evaluation of rhs to allow interleving of rhs evaluation with
    * electric field computation.
    */
   void
   fillAccelerationGhostCells()
   {
      // set boundary conditions for advection acceleration boundaries
      // (boundaries of physical domain)
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("BC (Vlasov)");

      if (m_use_new_bcs) {
         FORT_SET_ACCELERATION_BCS_4D_JB(*m_dist_func_local.getDataPointer(),
            BOX4D_TO_FORT(m_global_box),
            BOX4D_TO_FORT(m_local_box),
            BOX4D_TO_FORT(m_interior_box),
            m_spatial_solution_order,
            m_n_ghosts[2],
            *(m_vel_face[V1].getDataPointer()),
            *(m_vel_face[V2].getDataPointer()),
            PROBLEMDOMAIN_TO_FORT((*m_domain)),
            *(m_initial_condition->parameters()).getDataPointer());
      }
      else {
         FORT_SET_ACCELERATION_BCS_4D(*m_dist_func_local.getDataPointer(),
            BOX4D_TO_FORT(m_global_box),
            BOX4D_TO_FORT(m_local_box),
            BOX4D_TO_FORT(m_interior_box),
            m_spatial_solution_order,
            *(m_vel_face[V1].getDataPointer()),
            *(m_vel_face[V2].getDataPointer()),
            PROBLEMDOMAIN_TO_FORT((*m_domain)),
            *(m_initial_condition->parameters()).getDataPointer());
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

      FORT_COMPUTE_ADVECTION_FLUXES_4D(*(m_flux[X1].getDataPointer()),
         *(m_flux[X2].getDataPointer()),
         BOX4D_TO_FORT(m_local_box),
         *(m_vel_face[X1].getDataPointer()),
         *(m_vel_face[X2].getDataPointer()),
         *(m_u_face[X1].getDataPointer()),
         *(m_u_face[X2].getDataPointer()),
         *(m_dist_func_local.getDataPointer()),
         *(m_domain->dx()).dataPtr(),
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
         *(a_rhs.m_dist_func_local.getDataPointer()),
         *(m_dist_func_local.getDataPointer()),
         BOX4D_TO_FORT(m_local_box),
         BOX4D_TO_FORT(m_interior_box),
         *(m_vel_face[X1].getDataPointer()),
         *(m_vel_face[X2].getDataPointer()),
         *(m_domain->dx()).dataPtr(),
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

      FORT_COMPUTE_ACCELERATION_FLUXES_4D(*(m_flux[V1].getDataPointer()),
         *(m_flux[V2].getDataPointer()),
         BOX4D_TO_FORT(m_local_box),
         *(m_vel_face[V1].getDataPointer()),
         *(m_vel_face[V2].getDataPointer()),
         *(m_u_face[V1].getDataPointer()),
         *(m_u_face[V2].getDataPointer()),
         *(m_dist_func_local.getDataPointer()),
         *(m_domain->dx()).dataPtr(),
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
         *(a_rhs.m_dist_func_local.getDataPointer()),
         *(m_dist_func_local.getDataPointer()),
         BOX4D_TO_FORT(m_local_box),
         BOX4D_TO_FORT(m_interior_box),
         *(m_vel_face[V1].getDataPointer()),
         *(m_vel_face[V2].getDataPointer()),
         *(m_domain->dx()).dataPtr(),
         m_spatial_solution_order);

      timers->stopTimer("Vlasov");
   }

   /**
    * @brief Adds flux divergence to rhs and adds in collision and krook layer
    * operators.
    *
    * Part of evaluation of rhs to allow interleving of rhs evaluation with
    * electric field computation.
    *
    * @param[out] a_rhs The right hand side.
    * @param[in] a_time Current simulation time.
    * @param[in] a_use_new_alg If true use new time saving formulation.
    */
   void
   completeRHS(
      KineticSpecies& a_rhs,
      real a_time,
      bool a_use_new_alg) const
   {
      if (!a_use_new_alg) {
         addFluxDivergence(a_rhs);
      }

      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("collisions");
      for (int i = 0; i < m_num_collision_operators; ++i) {
         m_collision_operators[i]->evaluate(a_rhs.m_dist_func_local,
            m_dist_func_local,
            m_local_box,
            m_interior_box,
            *m_domain,
            m_vflowx,
            m_vflowy);
      }
      timers->stopTimer("collisions");
      timers->startTimer("krook");
      if (m_krook_layer->hasKrookLayer()) {
         real krook_power = m_krook_layer->power();
         real krook_coef = m_krook_layer->coefficient();
         FORT_APPEND_KROOK(
            BOX4D_TO_FORT(m_local_box),
            BOX4D_TO_FORT(m_interior_box),
            PROBLEMDOMAIN_TO_FORT((*m_domain)),
            *(m_krook_layer->lower()).dataPtr(),
            *(m_krook_layer->upper()).dataPtr(),
            krook_power,
            krook_coef,
            *m_dist_func_local.getDataPointer(),
            *(a_rhs.m_dist_func_local.getDataPointer()));
      }
      if (m_tz_source) {
         m_tz_source->set(a_rhs.m_dist_func_local, *m_domain, a_time);
      }
      timers->stopTimer("krook");
   }

   /**
    * @brief Compute and return stable time step limit for this KineticSpecies.
    */
   real
   computeDt();

   /**
    * @brief Compute a and y accelerations on the faces perpendicular to vx and
    *        vy respectively for electrostatic case.
    *
    * Includes evaluation of externally applied electric fields.
    *
    * @param[in] a_efield Electric field resulting from all KineticSpecies.
    * @param[in] a_ext_efield_global The global electric field from all external
    *                                drivers from this species.
    * @param[in] a_time The time at which velocity is computed.
    * @param[in] a_dt The time step at which velocity is computed.
    * @param[in] a_stage RK stage at which velocity is computed.
    */
   void
   computeAcceleration(
      const realArray& a_efield,
      realArray&       a_ext_efield_global,
      real             a_time,
      real             a_dt,
      int              a_stage);

   /**
    * @brief Compute a and y accelerations on the faces perpendicular to vx and
    *        vy respectively for electrodynamic case.
    *
    * Includes evaluation of externally applied electric fields.
    *
    * @param[in] a_maxwell Object holding electrodynamic fields.
    * @param[in] a_ext_efield_global The electric field from all external
    *                                drivers from this species.
    * @param[in] a_time The time at which velocity is computed.
    * @param[in] a_dt The time step at which velocity is computed.
    * @param[in] a_stage RK stage at which velocity is computed.
    */
   void
   computeAcceleration(
      Maxwell&   a_maxwell,
      realArray& a_ext_efield_global,
      real       a_time,
      real       a_dt,
      int        a_stage);

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
    * space boundary.
    *
    * @param[in] a_poisson Poisson object that contains velocity boundary flux
    * grid function.
    * @param[in] a_ke_vel_bdry_flux Total velocity boundary flux.
    */
   void
   computeKEVelBdyFlux(
      Poisson& a_poisson,
      realArray& a_ke_vel_bdry_flux);

   /**
    * @brief Compute the kinetic energy flux through the specified velocity
    * space boundary.
    *
    * @param[in] a_maxwell Maxwell object that contains velocity boundary flux
    * grid function.
    * @param[in] a_ke_vel_bdry_flux Total velocity boundary flux.
    */
   void
   computeKEVelBdyFlux(
      Maxwell& a_maxwell,
      realArray& a_ke_vel_bdry_flux);

   /**
    * @brief Accumulates sequences in a_sequences.
    *
    * @param[in] a_is_2d_proc True if this is a KineticSpecies on a Poisson
    *                         or Maxwell processor.
    * @param[in] a_sequences Sequences to save.
    * @param[in] a_saved_seq Number of times sequences that have been saved.
    * @param[in,out] a_seq_idx Index of next sequence member being accumulated.
    * @param[in] a_ke_global The global 2D kinetic energy.
    * @param[in] a_ke_flux_global The global 2D kinetic energy flux.
    */
   void
   accumulateSequences(
      bool a_is_2d_proc,
      RealArray& a_sequences,
      int a_saved_seq,
      int& a_seq_idx,
      realArray& a_ke_global,
      realArray& a_ke_flux_global);

   /**
    * @brief Indicates that the problem has particles.
    */
   void
   problemHasParticles()
   {
      m_problem_has_particles = true;
      m_ext_efield_local.redim(Range(m_accel_box.lower(X1),
                                     m_accel_box.upper(X1)),
                               Range(m_accel_box.lower(X2),
                                     m_accel_box.upper(X2)),
                               Range(0, 1));
      defineExtEfieldContractionSchedule();
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
      ParmParse& a_pp,
      const tbox::Pointer<ProblemDomain> a_cfg_domain);

   // Print this class' parameters.
   void
   printParameters() const;

   // Copies internals of rhs into this.
   void
   copy(
      const KineticSpecies& a_rhs,
      bool a_deep_copy = true);

   // Allocate private array data members and associated info.
   void
   allocateLocalAuxArrays();

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

   //  Defines charge density reduction schedule.
   void
   defineChargeDensityReductionSchedule();

   //  Defines kinetic energy reduction schedule.
   void
   defineKineticEnergyReductionSchedule();

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
      m_dist_func_global.updateGhostBoundaries();
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
         FORT_SET_ADVECTION_BCS_4D_JB(*m_dist_func_local.getDataPointer(),
            BOX4D_TO_FORT(m_global_box),
            BOX4D_TO_FORT(m_local_box),
            BOX4D_TO_FORT(m_interior_box),
            m_spatial_solution_order,
            m_n_ghosts[0],
            *(m_vel_face[X1].getDataPointer()),
            *(m_vel_face[X2].getDataPointer()),
            PROBLEMDOMAIN_TO_FORT((*m_domain)),
            *(m_initial_condition->parameters()).getDataPointer(),
            m_domain->isPeriodicAsInt(X1),
            m_domain->isPeriodicAsInt(X2));
      }
      else {
         FORT_SET_ADVECTION_BCS_4D(*m_dist_func_local.getDataPointer(),
            BOX4D_TO_FORT(m_global_box),
            BOX4D_TO_FORT(m_local_box),
            BOX4D_TO_FORT(m_interior_box),
            m_spatial_solution_order,
            *(m_vel_face[X1].getDataPointer()),
            *(m_vel_face[X2].getDataPointer()),
            PROBLEMDOMAIN_TO_FORT((*m_domain)),
            *(m_initial_condition->parameters()).getDataPointer(),
            m_domain->isPeriodicAsInt(X1),
            m_domain->isPeriodicAsInt(X2));
      }

      setPeriodicBCs();
      timers->stopTimer("BC (Vlasov)");
   }

   // Phase space dimension of this KineticSpecies.
   const tbox::Dimension m_pdim;

   // Configuration space dimension of this KineticSpecies.
   const tbox::Dimension m_cdim;

   // Name of this KineticSpecies.
   std::string m_name;

   // This KineticSpecies' ID.
   int m_species_index;

   // Total number of KineticSpecies.
   int m_number_of_species;

   // Mass of this KineticSpecies.
   real m_mass;

   // Charge of this KineticSpecies.
   real m_charge;

   // The global problem domain.
   tbox::Pointer<ProblemDomain> m_domain;

   // The global and local distribution functions.
   realArray m_dist_func_global;
#ifdef USE_PPP
   RealArray m_dist_func_local;
#else
   RealArray& m_dist_func_local;
#endif

   // The local x, y, and z components of the current.
   RealArray m_Jx_local;
   RealArray m_Jy_local;
   RealArray m_Jz_local;

   // The local kinetic energy.
   RealArray m_ke_local;

   // The local kinetic energy physical boundary flux.
   RealArray m_ke_phys_bdry_flux_local;

   // The local kinetic energy velocity boundary flux.
   RealArray m_ke_vel_bdry_flux_local;

   // Spatial order of solution
   int m_spatial_solution_order;

   // Temporal order of solution.
   int m_temporal_solution_order;

   // The number of ghost cells in each direction.
   tbox::IntVector m_n_ghosts;

   // Algorithm stencil width.
   int m_stencil_width;

   // Box describing the local extent including ghosts of this KineticSpecies.
   tbox::Box m_local_box;

   // Box describing the local extent not including ghosts of this
   // KineticSpecies.
   tbox::Box m_interior_box;

   // Box describing the global extent including ghosts of this KineticSpecies.
   tbox::Box m_global_box;

   // Box describing the extents of m_accel.
   tbox::Box m_accel_box;

   // Box describing the extents of m_em_vars.
   tbox::Box m_em_vars_box;

   // Box describing the extents of m_vz.
   tbox::Box m_vz_box;

   // How to initialize the distribution function.
   tbox::Pointer<ICInterface> m_initial_condition;

   // Any Krook layers in the problem.
   tbox::Pointer<KrookLayer> m_krook_layer;

   // The number of processors this KineticSpecies is distributed over.
   int m_number_of_procs;

   // If true the user has specified the number of processors.
   bool m_fixed_number_of_procs;

   // The range of processors this KineticSpecies is distributed over.
   Range m_processor_range;

   // If true, the distribution function has been partitioned.
   bool m_partition_defined;

   // The communicator for the processors this KineticSpecies is distributed
   // over.
   MPI_Comm m_comm;

   // The external electric field drivers.
   std::vector<tbox::Pointer<ElectricFieldDriver> > m_ef_drivers;

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

   // Reduction schedule for the kinetic energy.
   tbox::Pointer<ReductionSchedule> m_ke_schedule;

   // Reduction schedule for the kinetic energy physical boundray flux.
   tbox::Pointer<ReductionSchedule> m_ke_phys_bdry_flux_schedule;

   // Reduction schedule for the kinetic energy velocity boundray flux.
   tbox::Pointer<SummationSchedule> m_ke_vel_bdry_flux_schedule;

   // Acceleration.
   RealArray m_accel;

   // Local external electric field.
   RealArray m_ext_efield_local;

   // Electromagnetic variables.
   RealArray m_em_vars;

   // Transverse drift velocity.
   RealArray m_vz;

   // Maximum x and y accelerations.
   std::vector<real> m_lambda_max;

   // For each phase space dimension, the distribution function on those faces.
   std::vector<RealArray> m_u_face;

   // Something of a misnomer as its meaning is different for configuration
   // space dimensions and velocity space dimensions.  For the configuration
   // space dimensions, X1 and X2, this is the velocity in that dimension at
   // the faces perpendicular to that dimension.  Hence m_vel_face[X1] is the
   // x velocity at faces perpendicular to x.  For the velocity space
   // dimensions, V1 and V2, this is the acceleration in that dimension at the
   // faces perpendicular to that dimension.  Hence m_vel_face[V1] is the
   // x acceleration at faces perpendicular to vx.
   std::vector<RealArray> m_vel_face;

   // For each phase space dimension, the fluxes on those faces.
   std::vector<RealArray> m_flux;

   // If true then an electrodynamic problem is being run.
   bool m_do_maxwell;

   // The collision operators.
   std::vector<tbox::Pointer<CollisionOperator> > m_collision_operators;

   // Number of collision drivers.
   int m_num_collision_operators;

   // Flow velocity in x.
   real m_vflowx;

   // Flow velocity in y.
   real m_vflowy;

   // Twilight zone source.
   tbox::Pointer<TZSourceInterface> m_tz_source;

   // True if the problem contains particles.
   bool m_problem_has_particles;

   // If true use Jeff Banks' new acceleration and advection BC treatment.
   bool m_use_new_bcs;
};

} // end namespace Loki

#endif
