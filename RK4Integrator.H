/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _RK4_INTEGRATOR_H_
#define _RK4_INTEGRATOR_H_

#include "RKIntegrator.H"
#include "ODEInterface.H"

namespace Loki {

/**
 * This implements the interface of RKIntegrator for a fourth-order Runge-Kutta
 * time advance of an ODE.  ODE is d(TSoln)/dt = RHS
 *
 * Template types:
 *
 *   TSoln is the datatype for the solution
 */
template <class TSoln>
class RK4Integrator : public RKIntegrator<TSoln>
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
   RK4Integrator(
      ODEInterface<TSoln>& a_ode,
      const TSoln& a_protoSoln)
      : m_op(a_ode)
   {
      m_op.defineRHSData(m_rhs, a_protoSoln);
      m_op.defineRHSData(m_delta, a_protoSoln);
   }

   /**
    * @brief Destructor.
    */
   ~RK4Integrator()
   {
   }

   /**
    * @brief Advance by one time step using 4th-order Runge-Kutta and return
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
      static const double THIRD = 1.0 / 3.0;
      double dtOn2 = 0.5 * a_dt;
      double dtOn3 = THIRD * a_dt;
      double dtOn6 = 0.5 * dtOn3;

      // Initialize change accumulator
      m_op.zeroSolnData(m_delta);

      // RK Stage 1:
      stageAdvance(a_newSoln,
         a_oldSoln,
         a_oldSoln,
         a_time,
         a_dt,
         dtOn6,
         dtOn2,
         1);

      // RK Stage 2:
      stageAdvance(a_newSoln,
         a_oldSoln,
         a_newSoln,
         a_time + dtOn2,
         a_dt,
         dtOn3,
         dtOn2,
         2);

      // RK Stage 3:
      stageAdvance(a_newSoln,
         a_oldSoln,
         a_newSoln,
         a_time + dtOn2,
         a_dt,
         dtOn3,
         a_dt,
         3);

      // RK Stage 4:
      stageAdvance(a_newSoln,
         a_oldSoln,
         a_newSoln,
         a_time + a_dt,
         a_dt,
         dtOn6,
         1.0,
         4);

      return (a_time + a_dt);
   }

private:
   // Unimplemented default constructor.
   RK4Integrator();

   // Unimplemented copy constructor.
   RK4Integrator(
      const RK4Integrator& other);

   // Unimplemented assignment operator.
   RK4Integrator&
   operator = (
      const RK4Integrator& rhs);

   /**
    * Conduct a single stage advance.
    *
    * @param[in,out] a_predSoln The predicted solution.
    * @param[in] a_oldSoln The solution at the old time.
    * @param[in, out] a_evalSoln Same as predSoln, the predicted solution.
    * @param[in] a_time The time centering of predSoln.
    * @param[in] a_dt_eval The time step for the RHS evaluation.
    * @param[in] a_dt_update The time step for the update.
    * @param[in] a_stage RK state.
    */
   void
   stageAdvance(
      TSoln& a_predSoln,
      const TSoln& a_oldSoln,
      TSoln& a_evalSoln,
      double a_time,
      double a_dt,
      double a_dt_eval,
      double a_dt_update,
      int a_stage)
   {
      m_op.zeroSolnData(m_rhs);
      m_op.evalRHS(m_rhs, a_evalSoln, a_time, a_dt, a_stage == 1, a_stage == 4);
      m_op.addSolnData(m_delta, m_rhs, a_dt_eval, false);
      m_op.copySolnData(a_predSoln, a_oldSoln);
      if (a_stage < 4) {
         m_op.addSolnData(a_predSoln, m_rhs, a_dt_update, true);
      }
      else {
         m_op.addSolnData(a_predSoln, m_delta, a_dt_update, true);
      }
      m_op.postStageAdvance(a_predSoln, a_stage);
   }

   // The ODE operator.
   ODEInterface<TSoln>& m_op;

   // The right hand side of the ODE computed for a single stage advance.
   TSoln m_rhs;

   // The accumulation of the right hand side stage advances.
   TSoln m_delta;
};

} // end namespace Loki

#endif
