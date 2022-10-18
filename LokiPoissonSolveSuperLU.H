/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _LOKIPOISSONSOLVESUPERLU_H_
#define _LOKIPOISSONSOLVESUPERLU_H_

#include "PoissonSolver.H"
#include "petscksp.h"

namespace Loki {

/**
 * Class to solve Poisson's equation using SuperLU via PETSc.
 */
class LokiPoissonSolveSuperLU : public PoissonSolver
{
public:
   /**
    * @brief Constructor
    *
    * @param[in] a_domain Configuration space domain.
    * @param[in] a_order Order of accuracy.
    * @param[in] a_comm MPI communicator context for this solver.
    */
   LokiPoissonSolveSuperLU(
      const ProblemDomain& a_domain,
      int a_order,
      MPI_Comm a_comm);

   /**
    * @brief Destructor
    */
   ~LokiPoissonSolveSuperLU();

   //@{

   //! @name Implementation of PoissonSolver interface.

   /**
    * @brief Solve the linear system.
    *
    * @param[out] a_solution Solution as a ParallelArray.
    * @param[out] a_rhs RHS as a ParallelArray.
    */
   virtual
   void
   solve(
      ParallelArray& a_solution,
      const ParallelArray& a_rhs);

   /**
    * @brief Print the solver's parameters.
    */
   virtual
   void
   printParameters() const;

   //@}

   /**
    * @brief Returns true if name describes this class' type.
    *
    * @param[in] a_name Name of a type of Poisson solver.
    *
    * @return true if name describes this class' type.
    */
   static
   bool
   isType(
      const string& a_name);

private:
   // Unimplemented default constructor.
   LokiPoissonSolveSuperLU();

   // Unimplemented copy constructor.
   LokiPoissonSolveSuperLU(
      const LokiPoissonSolveSuperLU& a_other);

   // Unimplemented assignment operator.
   LokiPoissonSolveSuperLU&
   operator = (
      const LokiPoissonSolveSuperLU& a_rhs);

   // Build the linear system matrix, the RHS and solution vectors.  Fill in
   // terms of matrix given problem solution order and dimensions.
   void
   buildAndFillSystemMat(
      const ProblemDomain& a_domain,
      int a_order);

   // The name of this type of Poisson solver.
   static const string s_CLASS_NAME;

   // Create superlu linear system context.
   void
   createLinearSystemContext();

   // PETSc solution vector.
   Vec m_x;

   // PETSc RHS vector.
   Vec m_b;

   // PETSc linear system matrix
   Mat m_A;

   // PETSc linear solver context
   KSP m_ksp;

   // Interior size of problem in x.
   int m_nx;

   // Interior size of problem in y.
   int m_ny;

   // Size of problem in x including boundary.
   int m_nx_w_bdy;

   // Size of problem in y including boundary.
   int m_ny_w_bdy;

   // Total size of problem in x including boundary and extra equation.
   int m_nx_tot;

   // Total size of problem in y including boundary and extra equation.
   int m_ny_tot;

   // Number of ghosts.
   int m_nghosts;

   // Scale factor, the central value in the Laplacian stencil.
   double m_scale;
};

} // end namespace Loki

#endif
