/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _LOKIPOISSONSOLVEFFT_H_
#define _LOKIPOISSONSOLVEFFT_H_

#include "PoissonSolver.H"
#include <fftw3.h>

namespace Loki {

/**
 * Class to solve Poisson's equation using FFT.
 */
class LokiPoissonSolveFFT : public PoissonSolver
{
public:
   /**
    * @brief Constructor
    *
    * @param[in] a_pp The input parameter database.
    * @param[in] a_domain Configuration space domain.
    * @param[in] a_order Order of accuracy.
    */
   LokiPoissonSolveFFT(
      LokiInputParser& a_pp,
      const ProblemDomain& a_domain,
      int a_order);

   /**
    * @brief Destructor
    */
   ~LokiPoissonSolveFFT();

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
   LokiPoissonSolveFFT();

   // Unimplemented copy constructor.
   LokiPoissonSolveFFT(
      const LokiPoissonSolveFFT& a_other);

   // Unimplemented assignment operator.
   LokiPoissonSolveFFT&
   operator = (
      const LokiPoissonSolveFFT& a_rhs);

   // Parse the input parameters.
   void
   parseParameters(
      LokiInputParser& a_pp);

   // The name of this type of Poisson solver.
   static const string s_CLASS_NAME;

   // Problem size in x.
   int m_nx;

   // Problem size in y.
   int m_ny;

   // The number of ghosts.
   int m_nghosts;

   // The number problem size in x including the boundary layer.
   int m_nx_w_bdy;

   // The number problem size in y including the boundary layer.
   int m_ny_w_bdy;

   // Block of memory to do transform in-place.  Hold real valued input,
   // complex valued forward transform, and real valued inverse transform which
   // is the output.
   fftw_complex* m_mem;

   // Pointer to m_mem for real valued input.
   double* m_in;

   // Pointer to m_mem for complex valued transform.
   fftw_complex* m_transform;

   // Pointer to m_mem for real valued inverse transform, the output.
   double* m_out;

   // Forward transform.
   fftw_plan m_pf;

   // Inverse transform.
   fftw_plan m_pb;

   // X symbols to divide transform by.
   vector<double> m_sx;

   // Y symbols to divide transform by.
   vector<double> m_sy;

   // Solution order.
   int m_order;
};

} // end namespace Loki

#endif
