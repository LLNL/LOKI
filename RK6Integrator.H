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
#ifndef _RK6_INTEGRATOR_H_
#define _RK6_INTEGRATOR_H_

#include "RKIntegrator.H"
#include "ODEInterface.H"

namespace Loki {

/**
 * This implements the interface of RKIntegrator for a sixth-order Runge-Kutta
 * time advance of an ODE.  ODE is d(TSoln)/dt = RHS
 *
 * Template types:
 *
 *   TSoln is the datatype for the solution
 */
template <class TSoln>
class RK6Integrator : public RKIntegrator<TSoln>
{
public:
   /**
    * @brief Constructor.
    *
    * Initializes the integrator, primarily by associating an operator.
    * Note that the operator must outlive the integrator.
    *
    * @param[in] a_ode The ODE operator.
    * @param[in] a_protoSoln Prototype soln.
    */
   RK6Integrator(
      ODEInterface<TSoln>& a_ode,
      const TSoln& a_protoSoln)
      : m_op(a_ode)
   {
      m_k = new TSoln [8];
      for (int i = 0; i < 8; ++i) {
         m_op.defineRHSData(m_k[i], a_protoSoln);
      }
      m_tmp = new TSoln(a_protoSoln);
   }

   /**
    * @brief Destructor.
    */
   ~RK6Integrator()
   {
      delete [] m_k;
      delete m_tmp;
   }

   /**
    * @brief Advance by one time step using 6th-order Runge-Kutta and return
    * the new time (time + dt).
    *
    * @param[out] a_newSoln The solution at the new time.
    * @param[in] a_oldSoln The solution at the old time.
    * @param[in,out] a_time The time centering of oldSoln.
    * @param[in] a_dt The time step.
    *
    * @return New simulation time, a_time + a_dt.
    */
   real
   advance(
      TSoln& a_newSoln,
      TSoln& a_oldSoln,
      real   a_time,
      real   a_dt)
   {
      // Useful definitions
      static const real A[8][8] =
         {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
          {1.0/9.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
          {1.0/24.0, 1.0/8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
          {1.0/6.0, -1.0/2.0, 2.0/3.0, 0.0, 0.0, 0.0, 0.0, 0.0},
          {935.0/2536.0, -2781.0/2536.0, 309.0/317.0, 321.0/1268.0, 0.0, 0.0, 0.0, 0.0},
          {-12710.0/951.0, 8287.0/317.0, -40.0/317.0, -6335.0/317.0, 8.0, 0.0, 0.0, 0.0},
          {5840285.0/3104064.0, -7019.0/2536.0, -52213.0/86224.0, 1278709.0/517344.0, -433.0/2448.0, 33.0/1088.0, 0.0, 0.0},
          {-5101675.0/1767592.0, 112077.0/25994.0, 334875.0/441898.0, -973617.0/883796.0, -1421.0/1394.0, 333.0/5576.0, 36.0/41.0, 0.0}};
      
      static const real b[8] = {41.0/840.0, 
                                0.0, 
                                9.0/35.0, 
                                9.0/280.0,
                                34.0/105.0, 
                                9.0/280.0, 
                                9.0/35.0, 
                                41/840.0};

      static const real c[8] = {0.0,
                                1.0/9.0,
                                1.0/6.0,
                                1.0/3.0,
                                1.0/2.0,
                                2.0/3.0,
                                5.0/6.0,
                                1.0};
      // Initialize the results for the 8 stages.
      for (int i = 0; i < 8; ++i) {
         m_op.zeroSolnData(m_k[i]);
      }

      // Initialize temp with the old solution.
      m_op.copySolnData(*m_tmp, a_oldSoln);

      // RK Stage 1
      m_op.evalRHS(m_k[0], *m_tmp, a_time+c[0]*a_dt, a_dt, 1);

      // RK Stages 2-8
      for (int i = 1; i < 8; ++i) {
         m_op.copySolnData(*m_tmp, a_oldSoln);
         for (int j = 0; j < i; ++j) {
            m_op.addSolnData(*m_tmp, m_k[j], a_dt*A[i][j]);
         }
         m_op.evalRHS(m_k[i], *m_tmp, a_time+c[i]*a_dt, a_dt, i+1);
         m_op.postStageAdvance(*m_tmp, i+1);
      }

      // Form solution at end of time step from each RK stage result.
      m_op.copySolnData(a_newSoln, a_oldSoln);
      for (int i = 0; i < 8; ++i) {
         m_op.addSolnData(a_newSoln, m_k[i], a_dt*b[i], i == 7);
      }

      return (a_time + a_dt);
   }

private:
   // Unimplemented default constructor.
   RK6Integrator();

   // Unimplemented copy constructor.
   RK6Integrator(
      const RK6Integrator& other);

   // Unimplemented assignment operator.
   RK6Integrator&
   operator = (
      const RK6Integrator& rhs);

   // The ODE operator.
   ODEInterface<TSoln>& m_op;

   // The right hand side of the ODE computed for each stage advance.
   TSoln* m_k;

   // Temp solution needed for each stage advance.
   TSoln* m_tmp;
};

} // end namespace Loki

#endif
