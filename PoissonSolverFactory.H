/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _POISSONSOLVER_FACTORY_H_
#define _POISSONSOLVER_FACTORY_H_

#include "LokiInputParser.H"

namespace Loki {

class Poisson;
class PoissonSolver;

/**
 * Factory for creating Poisson solver objects.  Only contains static create
 * method.
 */
class PoissonSolverFactory
{
public:
   /**
    * @brief Creates the requested PoissonSolver.
    *
    * @param[in] a_pp Where the solver's input parameters will be found.
    * @param[in] a_poisson The Poisson object the solver works for.
    *
    * @return The requested PoissonSolver.
    */
   static
   PoissonSolver*
   create(
      LokiInputParser& a_pp,
      const Poisson& a_poisson);

private:
   // Unimplemented default constructor.
   PoissonSolverFactory();

   // Unimplemented destructor.
   ~PoissonSolverFactory();

   // Unimplemented copy constructor.
   PoissonSolverFactory(
      const PoissonSolverFactory& other);

   // Unimplemented assignment operator.
   PoissonSolverFactory&
   operator = (
      const PoissonSolverFactory& rhs);
};

} // end namespace Loki

#endif
