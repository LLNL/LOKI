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
#ifndef _MAXWELL_H_
#define _MAXWELL_H_

#include <string>
#include "Overture.h"
#include "CompositeGridOperators.h"
#include "ParallelUtility.h"
#include "ParmParse.H"
#include "Oges.h"
#include "Ogshow.h"
#include "SquareMapping.h"

#include "tbox/DescribedClass.H"
#include "tbox/Dimension.H"
#include "tbox/Pointer.H"
#include "tbox/Box.H"
#include "BoxOps.H"
#include "CurrentDriver.H"
#include "Directions.H"
#include "EMICInterface.H"
#include "VELICInterface.H"
#include "KineticSpeciesPtrVect.H"
#include "Load.H"
#include "ProblemDomain.H"
#include "Particle.H"
#include "MaxwellF.H"
#include "Serializable.H"
#include "TimerManager.H"
#include "Interpolator.H"

namespace Loki {

class ElectricPotentialDriver;

/**
 * Implementation of abstract Load interface for the solution of Maxwell's
 * equations.
 */
class Maxwell : public tbox::DescribedClass, public Load, public Serializable
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
    * @brief The number of grid functions.
    */
   static int NUM_FRAME_SERIES;

   /**
    * @brief Constructor.
    *
    * @param[in] a_pp Where the input parameters will be found.
    * @param[in] a_domain The global problem domain.
    * @param[in] a_num_kinetic_species Number of KineticSpecies in the problem.
    * @param[in] a_solution_order Order of solution.
    */
   Maxwell(
      ParmParse& a_pp,
      tbox::Pointer<ProblemDomain> a_domain,
      int a_num_kinetic_species,
      int a_solution_order);

   /**
    * @brief Copy constructor.
    *
    * @param[in] a_other What is being copies.
    * @param[in] a_deep_copy If true also copy electromagnetic fields and
    *                        transverse drift velocities.
    */
   Maxwell(
      const Maxwell& a_other,
      bool a_deep_copy = true);

   /**
    * @brief Destructor.
    */
   virtual
   ~Maxwell();

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
    * distributed.
    */
   virtual
   int
   numberOfProcessors() const;

   /**
    * @brief Returns true if this load is distributed over a user specified
    * fixed number of processors.
    */
   virtual
   bool
   fixedNumberOfProcessors() const;

   /**
    * @brief Partitions this load over the necessary number of processors.
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
    * @brief Returns true if this load is distributed to processor with
    * supplied id.
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

   //@{

   //! @name Implementation of Serializable interface.

   /**
    * @brief Read the electromagnetic fields and transverse drift velocities
    * from the restart database.
    *
    * @param[in] a_db The restart database.
    */
   virtual
   void
   getFromRestart(
      const HDF_DataBase& a_db);

   /**
    * @brief Write the electromagnetic fields and transverse drift velocities
    * to the restart database.
    *
    * @param[out] a_db The restart database.
    * @param[in] a_time Current simulation time.
    */
   virtual
   void
   putToRestart(
      HDF_DataBase& a_db,
      real a_time);

   //@}

   /**
    * @brief Returns true if this load has been distributed to the processor
    * on which the query is made.
    */
   bool
   isMaxwellProcessor() const
   {
      return m_is_maxwell_processor;
   }

   /**
    * @brief Dynamically creates a new Maxwell with conforming layout and
    * copies internal parameters.
    *
    * @param[in] a_ghost_cells IntVector of ghost cells in cloned version
    * @param[in] a_deep_copy If true copy the distribution function.
    */
   tbox::Pointer<Maxwell>
   clone(
      const tbox::IntVector& a_ghost_cells,
      bool a_deep_copy = true) const
   {
      tbox::Pointer<Maxwell> clone(new Maxwell(*this, a_deep_copy));
      clone->m_n_ghosts = a_ghost_cells;
      return clone;
   }

   /**
    * @brief Dynamically creates a new Maxwell with conforming layout and
    * copies internal parameters.
    *
    * @param[in] a_deep_copy If true copy the distribution function.
    */
   tbox::Pointer<Maxwell>
   clone(
      bool a_deep_copy = true) const
   {
      tbox::Pointer<Maxwell> clone(new Maxwell(*this, a_deep_copy));
      return clone;
   }

   /**
    * @brief Checks that the rhs and this have a similar storage
    * representations and returns true if so, false if not.
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
    * tracking particles.
    */
   void
   zeroData()
   {
      m_em_vars_local = 0.0;
      for (int i = 0; i < m_num_kinetic_species; ++i) {
         m_vz_local[i] = 0.0;
      }
      for (int i = 0; i < numTrackingParticles(); ++i) {
         m_tracking_particles[i].x() = 0.0;
         m_tracking_particles[i].y() = 0.0;
         m_tracking_particles[i].vx() = 0.0;
         m_tracking_particles[i].vy() = 0.0;
      }
   }

   /**
    * @brief Add a multiple of another Maxwell's data to this.
    *
    * @param[in] a_increment_maxwell Maxwell to add.
    * @param[in] a_factor Multiple of input Maxwell to add.
    * @param[in] a_final_rk If true, this has been called as part of the final
    *                       RK update and tracking particle positions need to
    *                       be limited.
    */
   void
   addData(
      const Maxwell& a_increment_maxwell,
      real a_factor = 1.0,
      bool a_final_rk = false);

   /**
    * @brief Copies electromagnetic fields and transverse drift velocities of
    * rhs into this.
    *
    * @param[in] a_rhs Maxwell from which to copy.
    */
   void
   copySolnData(
      const Maxwell& a_rhs);

   /**
    * @brief Compute and return stable time step limit for Maxwell's equations.
    */
   real
   computeDt();

   /**
    * @brief Initialize the electromagnetic fields and transverse drift
    * velocities.  Set up all private data related to plotting and the Poisson
    * solver to compute the intrinsic electric field resulting from the charge
    * distribution.
    *
    * @param[in] a_species_names The name of each kinetic species.
    */
   void
   initialize(
      const std::vector<string> a_species_names);

   /**
    * @brief Apply the initial conditions specified in the input deck to the
    * electromagnetic fields.
    */
   void
   initializeEM();

   /**
    * @brief Apply the initial conditions specified in the input deck to vz
    *        for all species.
    */
   void
   initializeVZ();

   /**
    * @brief Computes static factorization of Poisson solve matrix.
    */
   void
   computeFactorization();

   /**
    * @brief Apply external current drivers to compute applied currents.
    *
    * @param[out] a_Jx x component of external current
    * @param[out] a_Jy y component of external current
    * @param[out] a_Jz z component of external current
    * @param[in] a_time Simulation time to compute currents.
    */
   void
   applyExternalCurrent(
      realArray& a_Jx,
      realArray& a_Jy,
      realArray& a_Jz,
      real a_time);

   /**
    * @brief Compute the electric field resulting from the supplied charge
    * density.
    *
    * @param[in] a_net_charge_density Charge density resulting from all
    *                                 KineticSpecies.
    */
   void
   electricField(
      const realArray& a_net_charge_density);

   /**
    * @brief Allocate and distribute a variable defined on the 2D
    * electromagnetic domain.
    *
    * @param[out] a_var The allocated variable.
    * @param[in] a_depth The number of components in var.
    * @param[in] a_n_ghosts The number of ghost cells var has in each
    *                       direction.
    */
   void
   newAuxVariable(
      realArray& a_var,
      int a_depth,
      const tbox::IntVector& a_n_ghosts) const;

   /**
    * @brief Allocate and distribute a variable defined on the 2D
    * electromagnetic domain using the number of ghosts of this class.
    *
    * @param[out] a_var The allocated variable.
    * @param[in] a_depth The number of components in var.
    */
   void
   newAuxVariable(
      realArray& a_var,
      int a_depth) const
   {
      newAuxVariable(a_var, a_depth, m_n_ghosts);
   }

   /**
    * @brief Evaluate the right hand side of Maxwell's equations.  We are not
    * enforcing the constraint on div(B) in this step.
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
      const realArray& a_Jx,
      const realArray& a_Jy,
      const realArray& a_Jz,
      const realArray& a_net_ext_efield,
      real a_time);

   /**
    * @brief Writes electromagnetic fields and accumulated sequences to an
    * Overture plot file.
    *
    * @param[in] a_time Simulation time.
    * @param[in] a_dt Simulation time step.
    * @param[in] a_sequences Sequences to save.
    * @param[in] a_species_names Names of KineticSpecies in the problem.
    * @param[in] a_time_seq Simulation times of sequences.
    * @param[in] a_probes Location of each probe.
    * @param[in] a_num_probes Number of probes.
    * @param[in] a_num_seq Number of sequences at each saved time.
    * @param[in] a_saved_seq Number of times sequences that have been saved.
    * @param[in] a_saved_save Number of plots that have been saved.
    * @param[out] a_show Overture plot file.
    */
   void
   plot(
      real a_time,
      real a_dt,
      const RealArray& a_sequences,
      const std::vector<string>& a_species_names,
      const RealArray& a_time_seq,
      const RealArray& a_probes,
      int a_num_probes,
      int a_num_seq,
      int a_saved_seq,
      int& a_saved_save,
      Ogshow& a_show);

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
      RealArray& a_sequences,
      const RealArray& a_probes,
      int a_num_probes,
      int a_saved_seq,
      int& a_seq_idx);

   /**
    * @brief Accessor for the global electromagnetic fields.
    */
   const realArray&
   getGlobalEMVars() const
   {
      return m_em_vars_global;
   }

   /**
    * @brief Accessor for the global transverse drift velocity field for a
    *        given KineticSpecies.
    *
    * @param[in] a_species_idx Index of the species of interest.
    */
   const realArray&
   getGlobalVZVar(
      int a_species_idx) const
   {
      return m_vz_global[a_species_idx];
   }

   /**
    * @brief Accessor for the global kinetic energy flux for a given
    *        KineticSpecies, direction and side.
    *
    * @param[in] a_species_idx Index of the species of interest.
    * @param[in] a_dir x or y direction.
    * @param[in] a_side High or low side.
    */
   realArray&
   getKEFluxVar(
      int a_species_idx,
      int a_dir,
      int a_side)
   {
      if (a_dir != V1 && a_dir != V2) {
         OV_ABORT("Illegal direction requested.");
      }
      if (a_side != LO && a_side != HI) {
         OV_ABORT("Illegal side requested.");
      }
      if (a_dir == V1) {
         if (a_side == LO) {
            return m_ke_flux_vx_lo[a_species_idx][0];
         }
         else {
            return m_ke_flux_vx_hi[a_species_idx][0];
         }
      }
      else {
         if (a_side == LO) {
            return m_ke_flux_vy_lo[a_species_idx][0];
         }
         else {
            return m_ke_flux_vy_hi[a_species_idx][0];
         }
      }
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
   * @brief Read in tracking particles if any.
   *
   * @param[in] a_pp Where the input parameters will be found.
   */
  void
  readParticles(
      ParmParse& a_pp);

private:
   // Unimplemented default constructor.
   Maxwell();

   // Unimplemented assignment operator.
   Maxwell&
   operator = (
      const Maxwell& rhs);

   // Different possible solver methods.
   enum SolverMethods {
      ORIGINAL_ITERATIVE,
      OVERTURE_BEST_ITERATIVE,
      SUPERLU_DIRECT
   };

   // Obtain control parameters from "maxwell" section of database.
   void
   parseParameters(
      ParmParse& a_pp);

   // Print the class' parameters.
   void
   printParameters() const;

   // Copies internals of a_rhs into this.  If a_deep_copy is true also copies
   // electromagnetic fields and transverse drift velocities.
   void
   copy(
      const Maxwell& a_rhs,
      bool a_deep_copy);

   // Creates mapped grid for Poisson solver, sets solver parameters and sets
   // boundary conditions.
   void
   defineSolver();

   // Returns the number of tracking particles in object.
   int
   numTrackingParticles() const
   {
      return static_cast<int>(m_tracking_particles.size());
   }

   // Creates a SquareMapping which essentially defines x and y coordinates.
   Mapping*
   createMapping(
      const int& a_sign)
   {
      const Array<double>& x_hi(m_domain->upper());
      const Array<double>& x_lo(m_domain->lower());
      const Array<double>& dx(m_domain->dx());
      const tbox::IntVector& ncells(m_domain->numberOfCells());

      Mapping *mapping;
      mapping = new SquareMapping(x_lo[X1] + 0.5 * dx[X1],
         x_hi[X1] + a_sign * 0.5 * dx[X1],
         x_lo[X2] + 0.5 * dx[X2],
         x_hi[X2] + a_sign * 0.5 * dx[X2]);
      mapping->setGridDimensions(X1, ncells[X1] + (1 + a_sign)/2);
      mapping->setGridDimensions(X2, ncells[X2] + (1 + a_sign)/2);
      mapping->setIsPeriodic(X1, Mapping::functionPeriodic);
      mapping->setIsPeriodic(X2, Mapping::functionPeriodic);
      mapping->setNumberOfGhostPoints(Start, X1, m_n_ghosts[X1]);
      mapping->setNumberOfGhostPoints(End,   X1, m_n_ghosts[X1]);
      mapping->setNumberOfGhostPoints(Start, X2, m_n_ghosts[X2]);
      mapping->setNumberOfGhostPoints(End,   X2, m_n_ghosts[X2]);
      return mapping;
   }

   // Creates mapped grid for the 2D configuration space with the specified
   // stencil width.
   void
   setupMappedGrid(
      MappedGrid& a_mapped_grid,
      int stencil_width)
   {
      a_mapped_grid.setDiscretizationWidth(X1, stencil_width);
      a_mapped_grid.setDiscretizationWidth(X2, stencil_width);
      a_mapped_grid.setNumberOfGhostPoints(Start, X1, m_n_ghosts[X1]);
      a_mapped_grid.setNumberOfGhostPoints(End,   X1, m_n_ghosts[X1]);
      a_mapped_grid.setNumberOfGhostPoints(Start, X2, m_n_ghosts[X2]);
      a_mapped_grid.setNumberOfGhostPoints(End,   X2, m_n_ghosts[X2]);
      a_mapped_grid.specifyProcesses(m_processor_range);
      a_mapped_grid.updateReferences();
   }

   // Sets solver parameters given requested solver type and other user input.
   void
   setupSolver(
      MappedGrid& a_mapped_grid)
   {
      m_solver = new Oges(a_mapped_grid);

      m_solver->setCommunicator(m_comm);

      if (m_solver_method == ORIGINAL_ITERATIVE) {
         m_solver->set(OgesParameters::THEsolverType,
            OgesParameters::PETScNew);
 
         m_solver->set(OgesParameters::THEsolverMethod,
            OgesParameters::biConjugateGradient);

         m_solver->set(OgesParameters::THEpreconditioner,
            OgesParameters::luPreconditioner);

         real tol(std::max(m_tol, 10.0*std::numeric_limits<real>::epsilon()));
         m_solver->set(OgesParameters::THErelativeTolerance, tol);

         m_solver->set(OgesParameters::THEmaximumNumberOfIterations,
            m_max_iterations);
      }
      else if (m_solver_method == OVERTURE_BEST_ITERATIVE) {
         m_solver->set(OgesParameters::THEbestIterativeSolver);

         real tol(std::max(m_tol, 10.0*std::numeric_limits<real>::epsilon()));
         m_solver->set(OgesParameters::THErelativeTolerance, tol);

         m_solver->set(OgesParameters::THEmaximumNumberOfIterations,
            m_max_iterations);
      }
      else {
         m_solver->set(OgesParameters::THEsolverType,
            OgesParameters::PETScNew);

         m_solver->set(OgesParameters::THEsolverMethod,
            OgesParameters::preonly);

         m_solver->set(OgesParameters::THEparallelSolverMethod,
            OgesParameters::preonly);

         m_solver->set(OgesParameters::THEexternalSolver,
            OgesParameters::superlu_dist);

         m_solver->set(OgesParameters::THEparallelExternalSolver,
            OgesParameters::superlu_dist);

         m_solver->set(OgesParameters::THEpreconditioner,
            OgesParameters::luPreconditioner);

         m_solver->set(OgesParameters::THEparallelPreconditioner,
            OgesParameters::luPreconditioner);
      }

      m_solver_composite_grid.add(a_mapped_grid);
      m_solver_composite_grid.updateReferences();
      m_solver_composite_grid.update(MappedGrid::THEvertex | MappedGrid::THEcenter);

      setupBoundaryConditions(m_solver_composite_grid);

      m_rho_solver = realCompositeGridFunction(m_solver_composite_grid);
      m_rho_solver = 0.0;
      real singular_offset(0.0);
      m_solver->setExtraEquationValues(m_rho_solver, &singular_offset);

      printF("\n\n ==== Solver:\n %s\n =====\n\n",
         static_cast<const char*>(m_solver->parameters.getSolverName()));
   }

   // Sets up Poisson solver boundary conditions.
   void
   setupBoundaryConditions(
      CompositeGrid& a_composite_grid)
   {
      CompositeGridOperators poisson_op(a_composite_grid);
      poisson_op.setStencilSize(m_stencil_size);
      poisson_op.setOrderOfAccuracy(m_solution_order);

      IntegerArray boundary_conditions(2, 3, 1);   // for Oges
      boundary_conditions = 0;

      RealArray boundary_condition_data(2, 2, 3, 1); // for Oges
      boundary_condition_data = 0.0;

      m_solver->setEquationAndBoundaryConditions(OgesParameters::laplaceEquation,
         poisson_op,
         boundary_conditions,
         boundary_condition_data);
   }

   // If this object is on the Maxwell processor and there are tracking
   // particles then this method builds the field interpolator appropriate to
   // the spatial order of the problem.
   void
   buildInterpolator();

   // Configuration space dimension this class operates in.
   tbox::Dimension m_dim;

   // Global problem domain.
   tbox::Pointer<ProblemDomain> m_domain;

   // The number of processors this Maxwell is distributed over.
   int m_number_of_procs;

   // The range of processors this Maxwell is distributed over.
   Range m_processor_range;

   // If true, the electromagnetic fields and transverse drift velocities have
   // been partitioned.
   bool m_partition_defined;

   // The number of ghost cells in each direction.
   tbox::IntVector m_n_ghosts;

   // Box describing the local extent including ghosts of this Maxwell.
   tbox::Box m_local_box;

   // Box describing the local extent not including ghosts of this Maxwell.
   tbox::Box m_interior_box;

   // Box describing the global extent including ghosts of this Maxwell.
   tbox::Box m_global_box;

   // True if this Maxwell has been distributed to the proceesor owning it.
   bool m_is_maxwell_processor;

   // The communicator for the processors this Maxwell is distributed
   // over.
   MPI_Comm m_comm;

   // The partitioning of this Maxwell.
   Partitioning_Type m_partition;

   // The number of KineticSpecies in the problem.
   int m_num_kinetic_species;

   // The electromagnetic fields.
   realArray m_em_vars_global;
#ifdef USE_PPP
   RealArray m_em_vars_local;
#else
   RealArray& m_em_vars_local;
#endif

   // Transverse drift velocity for each KineticSpecies.
   std::vector<realArray> m_vz_global;
#ifdef USE_PPP
   std::vector<RealArray> m_vz_local;
#else
   std::vector<RealArray>& m_vz_local;
#endif

   // Light speed.
   real m_c;

   // Undocumented artificial viscosity for RHS evaluation.
   real m_avWeak;
   real m_avStrong;

   // Poisson solver.
   tbox::Pointer<Oges> m_solver;

   // Grid Poisson solve takes place on.
   CompositeGrid m_solver_composite_grid;

   // Grid plottable grid functions are defined on.
   CompositeGrid m_composite_grid;

   // Grid function containing charge density.
   realCompositeGridFunction m_rho_solver;

   // Grid function containing potential.
   realCompositeGridFunction m_phi_solver;

   // Grid function containing potential.
   realCompositeGridFunction m_phi;

   // X, Y, and Z components of electric fields as plottable grid functions.
   realCompositeGridFunction m_ex_plot;
   realCompositeGridFunction m_ey_plot;
   realCompositeGridFunction m_ez_plot;

   // X, Y, and Z components of magnetic fields as plottable grid functions.
   realCompositeGridFunction m_bx_plot;
   realCompositeGridFunction m_by_plot;
   realCompositeGridFunction m_bz_plot;

   // VZ of each species as plottable grid functions.
   std::vector<realCompositeGridFunction> m_vz_plot;

   // For each KineticSpecies, the flux of kinetic energy through each velocity
   // space boundary as plottable grid functions.
   std::vector<realCompositeGridFunction> m_ke_flux_vx_hi;
   std::vector<realCompositeGridFunction> m_ke_flux_vx_lo;
   std::vector<realCompositeGridFunction> m_ke_flux_vy_hi;
   std::vector<realCompositeGridFunction> m_ke_flux_vy_lo;

   // Solver method.
   SolverMethods m_solver_method;

   // Solver tolerance.
   real m_tol;

   // Boundary condition stencil size.
   int m_stencil_size;

   // Mapped grid stencil width.
   int m_stencil_width;

   // Maximum number of solver iterations.
   int m_max_iterations;

   // The current drivers.
   std::vector<tbox::Pointer<CurrentDriver> > m_current_drivers;

   // The number of current drivers.
   int m_num_current_drivers;

   // Order of solution.
   int m_solution_order;

   // The electromagnetic field initial conditions.
   std::vector<tbox::Pointer<EMICInterface> > m_em_ics;

   // The z velocity initial conditions.
   std::vector<tbox::Pointer<VELICInterface> > m_vel_ics;

   // Number of tracking particles in the problem.
   int m_problem_num_tracking_particles;

   // Tracking particles
   std::vector<Particle> m_tracking_particles;

   // Field data interpolator.  Currently needed only for particle calculations.
   tbox::Pointer<Interpolator> m_interpolator;

   // Interpolation coefficient storage.
   double** m_ex_interp;
   double** m_ey_interp;
   double** m_bz_interp;
};

} // end namespace Loki

#endif
