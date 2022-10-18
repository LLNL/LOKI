/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "PoissonSolverFactory.H"
#include "Poisson.H"

// Add new derived PoissonSolver headers here
#include "LokiPoissonSolveSuperLU.H"
#include "LokiPoissonSolveFFT.H"

namespace Loki {

PoissonSolver*
PoissonSolverFactory::create(
   LokiInputParser& a_pp,
   const Poisson& a_poisson)
{
   // Get the name of the Poisson solver.  LokiPoissonSolveSuperLU is the
   // default for backward compatibility.
   string solver_name("superlu");
   a_pp.query("name", solver_name);
   PoissonSolver* solver = 0;

   // Build the requested Poisson solver.
   if (LokiPoissonSolveSuperLU::isType(solver_name)) {
      MPI_Comm comm(MPI_COMM_NULL);
      if (a_poisson.num_procs() > 1) {
         // Currently, "parallel" solves of Poisson's equation are done by each
         // Poisson processor getting the serialized charge density and doing
         // the serial solve.  To make PETSc happy, the communicator in the
         // serial solver must contain only one processor.  So for "parallel"
         // Poisson the serial solver can't just be the communicator generated
         // by the load balancer.  In this case, each Poisson processor must
         // create a communicator containing only itself to pass to the serial
         // solver.
         int comm_id;
         MPI_Comm_rank(a_poisson.comm(), &comm_id);
         const int status = MPI_Comm_split(a_poisson.comm(),
            Loki_Utilities::s_my_id,
            comm_id,
            &comm);
         if (status != MPI_SUCCESS) {
            LOKI_ABORT("Communicator split failed.");
         }
      }
      else {
         comm = a_poisson.comm();
      }
      solver = new LokiPoissonSolveSuperLU(a_poisson.domain(),
         a_poisson.order(),
         comm);
   }
   else if (LokiPoissonSolveFFT::isType(solver_name)) {
      solver = new LokiPoissonSolveFFT(a_pp,
         a_poisson.domain(),
         a_poisson.order());
   }
   // Add new cases here in this form:
   // else if (NewSolver::isType(solver_name)) {
   //    solver = new NewSolver(a_pp);
   // }
   else {
      ostringstream error_msg;
      error_msg << "Unknown poisson solver \"" << solver_name
                << "\" ... quitting";
      LOKI_ABORT(error_msg);
   }

   return solver;
}

} // end namespace Loki
