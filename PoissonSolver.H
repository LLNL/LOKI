/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _POISSONSOLVERS_H_
#define _POISSONSOLVERS_H_

#include "ParallelArray.H"
#include "ProblemDomain.H"

namespace Loki {

/**
 * Abstract base class for solvers of Poisson's equation.
 */
class PoissonSolver
{
public:
   /**
    * @brief Default constructor
    */
   PoissonSolver();

   /**
    * @brief Destructor
    */
   ~PoissonSolver();

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
      const ParallelArray& a_rhs) = 0;

   /**
    * @brief Print the solver's parameters.
    */
   virtual
   void
   printParameters() const = 0;

private:
   // Unimplemented copy constructor.
   PoissonSolver(
      const PoissonSolver& a_other);

   // Unimplemented assignment operator.
   PoissonSolver&
   operator = (
      const PoissonSolver& a_rhs);
};

} // end namespace Loki

#endif
