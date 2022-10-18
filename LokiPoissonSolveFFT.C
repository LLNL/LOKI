/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "LokiPoissonSolveFFT.H"
#include "Directions.H"
#include "Loki_Utilities.H"
#include <math.h>

namespace Loki {

const string LokiPoissonSolveFFT::s_CLASS_NAME("fft");


bool
LokiPoissonSolveFFT::isType(
   const string& a_name)
{
   if (a_name.compare(s_CLASS_NAME) == 0) {
      return true;
   }
   return false;
}


LokiPoissonSolveFFT::LokiPoissonSolveFFT(
   LokiInputParser& a_pp,
   const ProblemDomain& a_domain,
   int a_order)
   : m_nx(a_domain.numberOfCells(X1)),
     m_ny(a_domain.numberOfCells(X2)),
     m_nghosts(a_order == 4 ? 2 : 3),
     m_sx(a_domain.numberOfCells(X1)),
     m_sy(a_domain.numberOfCells(X2)/2+1),
     m_order(a_order)
{
   parseParameters(a_pp);
   m_nx_w_bdy = m_nx + 2*m_nghosts;
   m_ny_w_bdy = m_ny + 2*m_nghosts;
   const double pi = 4.0*atan(1.0);
   double Lx = a_domain.upper(X1)-a_domain.lower(X1);
   double dx = Lx/m_nx;
   double Ly = a_domain.upper(X2)-a_domain.lower(X2);
   double dy = Ly/m_ny;

   // Construct memory block and assign pointers.
   m_mem = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m_nx*(m_ny/2+1));
   m_in = m_mem[0];
   m_transform = m_mem;
   m_out = m_mem[0];

   // Construct the transforms.
   m_pf = fftw_plan_dft_r2c_2d(m_nx, m_ny, m_in, m_transform, FFTW_MEASURE);
   m_pb = fftw_plan_dft_c2r_2d(m_nx, m_ny, m_transform, m_out, FFTW_MEASURE);

   // Execute the transforms once so fftw can optimize them.
   fftw_execute(m_pf);
   fftw_execute(m_pb);

   // Compute the symbols m_sx and m_sy.
   // First compute kx and ky.
   vector<double> kx(m_nx);
   for (int i = 1; i < m_nx; ++i) {
      if (2*i < m_nx) {
         kx[i] = (2.0*pi/Lx)*i;
      }
      else {
         kx[i] = (2.0*pi/Lx)*(m_nx-i);
      }
   }
   vector<double> ky(m_ny/2+1);
   for (int i = 1; i < m_ny/2+1; ++i) {
      ky[i] = (2.0*pi/Ly)*i;
   }

   // Now store m_sx and m_sy.
   for (int i = 0; i < m_nx; ++i) {
      double dpdm = (2.0*cos(dx*kx[i])-2.0)/pow(dx, 2.0);
      if (m_order == 4) {
         m_sx[i] = dpdm - pow(dx, 2.0)/12.0*pow(dpdm, 2.0);
      }
      else if (m_order == 6) {
         m_sx[i] = dpdm - pow(dx, 2.0)/12.0*pow(dpdm, 2.0) +
            pow(dx, 4.0)/90.0*pow(dpdm, 3.0);
      }
      else if (m_order == -1) {
         m_sx[i] = -kx[i]*kx[i];
      }
      else {
         LOKI_ABORT("Order not supported");
      }
      m_sx[i] *= m_nx*m_ny;
   }
   for (int i = 0; i < m_ny/2+1; ++i) {
      double dpdm = (2.0*cos(dy*ky[i])-2.0)/pow(dy, 2.0);
      if (m_order == 4) {
         m_sy[i] = dpdm - pow(dy, 2.0)/12.0*pow(dpdm, 2.0);
      }
      else if (m_order == 6) {
         m_sy[i] = dpdm - pow(dy, 2.0)/12.0*pow(dpdm, 2.0) +
            pow(dy, 4.0)/90.0*pow(dpdm, 3.0);
      }
      else if (m_order == -1) {
         m_sy[i] = -ky[i]*ky[i];
      }
      else {
         LOKI_ABORT("Order not supported");
      }
      m_sy[i] *= m_nx*m_ny;
   }
}


LokiPoissonSolveFFT::~LokiPoissonSolveFFT()
{
   // Destroy the forward and inverse transforms and the complex data buffer.
   fftw_destroy_plan(m_pf);
   fftw_destroy_plan(m_pb);
   fftw_free(m_mem);
}


void
LokiPoissonSolveFFT::solve(
   ParallelArray& a_solution,
   const ParallelArray& a_rhs)
{
   // Dump the interior of a_rhs into m_in.
   int incr = 2*(m_ny/2+1);
   int idx;
   for (int j = 0; j < m_ny; ++j) {
      idx = j;
      for (int i = 0; i < m_nx; ++i) {
         m_in[idx] = a_rhs(i, j);
         idx += incr;
      }
   }

   // Perform forward transform.
   fftw_execute(m_pf);

   // Add symbols to transform.
   idx = 0;
   for (int i = 0; i < m_nx; ++i) {
      for (int j = 0; j < m_ny/2+1; ++j) {
         if (m_sx[i] != 0.0 || m_sy[j] != 0.0) {
            m_transform[idx][0] /= m_sx[i]+m_sy[j];
            m_transform[idx][1] /= m_sx[i]+m_sy[j];
         }
         ++idx;
      }
   }

   // Perform inverse transform.
   fftw_execute(m_pb);

   // Dump result into a_solution.
   for (int j = 0; j < m_ny; ++j) {
      idx = j;
      for (int i = 0; i < m_nx; ++i) {
         a_solution(i, j) = m_out[idx];
         idx += incr;
      }
   }
}


void
LokiPoissonSolveFFT::printParameters() const
{
   // Write the solver's parameters.
   Loki_Utilities::printF("  Using FFT Poisson solver:\n" );
   if (m_order == -1) {
      Loki_Utilities::printF("    discretization = spectral\n");
   }
   else {
      Loki_Utilities::printF("    discretization = solution order\n");
   }
}


void
LokiPoissonSolveFFT::parseParameters(
   LokiInputParser& a_pp)
{
   if (a_pp.contains("discretization")) {
      string tmp;
      a_pp.get("discretization", tmp);
      if (tmp.compare("spectral") == 0) {
         m_order = -1;
      }
      else if (tmp.compare("order") != 0) {
         LOKI_ABORT("Unknown discretization method.");
      }
   }
}

} // end namespace Loki
