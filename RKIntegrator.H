/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _RK_INTEGRATOR_H_
#define _RK_INTEGRATOR_H_

#include "tbox/DescribedClass.H"

namespace Loki {

/**
 * This defines interface for Runge-Kutta time advance of an ODE.
 * ODE is d(TSoln)/dt = RHS
 *
 * Template types:
 *
 *   TSoln is the datatype for the solution
 */
template <class TSoln>
class RKIntegrator : public tbox::DescribedClass
{
public:
   /**
    * @brief Constructor.
    */
   RKIntegrator()
   {
   }

   /**
    * @brief Destructor.
    */
   ~RKIntegrator()
   {
   }

   /**
    * @brief Advance by one time step using Runge-Kutta and return the new
    *        time (time + dt).
    *
    * @param[out] a_newSoln The solution at the new time.
    * @param[in] a_oldSoln The solution at the old time.
    * @param[in,out] a_time The time centering of oldSoln.
    * @param[in] a_dt The time step.
    *
    * @return New simulation time, teim + dt.
    */
   virtual
   double
   advance(
      TSoln& a_newSoln,
      TSoln& a_oldSoln,
      double a_time,
      double a_dt) = 0;

private:
   // Unimplemented copy constructor.
   RKIntegrator(
      const RKIntegrator& other);

   // Unimplemented assignment operator.
   RKIntegrator&
   operator = (
      const RKIntegrator& rhs);
};

} // end namespace Loki

#endif
