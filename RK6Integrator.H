/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
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
   }

   /**
    * @brief Destructor.
    */
   ~RK6Integrator()
   {
      delete [] m_k;
   }

   /**
    * @brief Advance by one time step using 6th-order Runge-Kutta and return
    *        the new time (time + dt).
    *
    * @param[out] a_newSoln The solution at the new time.
    * @param[in] a_oldSoln The solution at the old time.
    * @param[in,out] a_time The time centering of oldSoln.
    * @param[in] a_dt The time step.
    *
    * @return New simulation time, a_time + a_dt.
    */
   double
   advance(
      TSoln& a_newSoln,
      TSoln& a_oldSoln,
      double   a_time,
      double   a_dt)
   {
      // Useful definitions
      static const double A[8][8] =
         {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
          {1.0/9.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
          {1.0/24.0, 1.0/8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
          {1.0/6.0, -1.0/2.0, 2.0/3.0, 0.0, 0.0, 0.0, 0.0, 0.0},
          {935.0/2536.0, -2781.0/2536.0, 309.0/317.0, 321.0/1268.0, 0.0, 0.0, 0.0, 0.0},
          {-12710.0/951.0, 8287.0/317.0, -40.0/317.0, -6335.0/317.0, 8.0, 0.0, 0.0, 0.0},
          {5840285.0/3104064.0, -7019.0/2536.0, -52213.0/86224.0, 1278709.0/517344.0, -433.0/2448.0, 33.0/1088.0, 0.0, 0.0},
          {-5101675.0/1767592.0, 112077.0/25994.0, 334875.0/441898.0, -973617.0/883796.0, -1421.0/1394.0, 333.0/5576.0, 36.0/41.0, 0.0}};
      
      static const double b[8] = {41.0/840.0, 
                                  0.0, 
                                  9.0/35.0, 
                                  9.0/280.0,
                                  34.0/105.0, 
                                  9.0/280.0, 
                                  9.0/35.0, 
                                  41/840.0};

      static const double c[8] = {0.0,
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
      // note I've used a_newSoln for the temporary storage
      m_op.copySolnData(a_newSoln, a_oldSoln);

      // RK Stage 1
      m_op.evalRHS(m_k[0], a_newSoln, a_time+c[0]*a_dt, a_dt, true, false);

      // RK Stages 2-8
      for (int i = 1; i < 8; ++i) {
         m_op.copySolnData(a_newSoln, a_oldSoln);
         for (int j = 0; j < i; ++j) {
            m_op.addSolnData(a_newSoln, m_k[j], a_dt*A[i][j], true);
         }
         m_op.evalRHS(m_k[i], a_newSoln, a_time+c[i]*a_dt, a_dt, false, i == 7);
         m_op.postStageAdvance(a_newSoln, i+1);
      }

      // Form solution at end of time step from each RK stage result.
      m_op.copySolnData(a_newSoln, a_oldSoln);
      for (int i = 0; i < 8; ++i) {
         m_op.addSolnData(a_newSoln, m_k[i], a_dt*b[i], true);
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
};

} // end namespace Loki

#endif
