/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _SHAPED_RAMPED_COSINE_CURRENT_DRIVER_H_
#define _SHAPED_RAMPED_COSINE_CURRENT_DRIVER_H_

#include "CurrentDriver.H"
#include "LokiInputParser.H"

namespace Loki {

/**
 * Implementation of the CurrentDriver interface for a shaped ramped cosine
 * driver with the following functional form:
 *
 * @verbatim
 *  J_ext = J0*f(t)*g(x)*h(y)*cos(omega*(t-t0))
 *
 *          / 0.5+0.5*tanh(4.0*(2.0*(t-t0)/t_rampup-1.0))     for t0 <= t < t0+t_rampup
 *         /  0.5-0.5*tanh(4.0)                               for t0+t_rampup <= t < t0+t_rampup+t_hold
 *  f(t) =|
 *         \  0.5-0.5*tanh(4.0*(2.0*(t-t0-t_rampup-t_hold)/t_rampdown-1.0)) for t0+t_rampup+t_hold <= t < t0+t_rampup+t_hold+t_rampdown
 *          \ 0                                               else
 *
 *         / 1                                          x == x0
 *  g(x) =|
 *         \ 0                                          else
 *
 *         / 1-shape*(sin(pi*y/width))**2               for |y| < width/2
 *  h(y) =|
 *         \ 1-shape                                    else
 * @endverbatim
 */
class ShapedRampedCosineCurrentDriver : public CurrentDriver
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_pp Where the input parameters are to be found.
    */
   ShapedRampedCosineCurrentDriver(
      LokiInputParser& a_pp);

   /**
    * @brief Destructor.
    */
   virtual
   ~ShapedRampedCosineCurrentDriver();

   //@{

   //! @name Implementation of CurrentDriver interface.

   /**
    * @brief Compute and accumulate the antenna source at the supplied time.
    *
    * @param[out] a_antenna_source The antenna source accumulator.
    * @param[in] a_omega_eff2 Square of the effective plasma frequency.
    * @param[in] a_domain Information about the problem domain.
    * @param[in] a_time The time at which the current is to be computed.
    */
   virtual
   void
   evaluate(
      ParallelArray& a_antenna_source,
      const ParallelArray& a_omega_eff2,
      const ProblemDomain& a_domain,
      double a_time) const;

   /**
    * @brief Print the driver's parameters.
    */
   virtual
   void
   printParameters() const;

   /**
    * @brief Returns true if the dispersion relation is needed.
    */
   virtual
   bool
   needsDispersionRelation();

   //@}

   /**
    * @brief Returns true if name describes this class' type.
    *
    * @param[in] a_name Name of a type of current driver.
    *
    * @return true if name describes this class' type.
    */
   static
   bool
   isType(
      const string& a_name);

private:
   // Unimplemented default constructor.
   ShapedRampedCosineCurrentDriver();

   // Unimplemented copy constructor.
   ShapedRampedCosineCurrentDriver(
      const ShapedRampedCosineCurrentDriver& other);

   // Unimplemented assignment operator.
   ShapedRampedCosineCurrentDriver&
   operator = (
      const ShapedRampedCosineCurrentDriver& rhs);

   // Parse the input parameters.
   void
   parseParameters(
      LokiInputParser& a_pp);

   // The name of this type of current driver.
   static const string s_CLASS_NAME;

   // Indices of specific parameters in m_parameters.
   enum {
      WIDTH,
      APPLY_DIR,
      SHAPE,
      OMEGA,
      J0,
      T0,
      TRAMPUP,
      THOLD,
      TRAMPDOWN,
      X0,
      PLANE,
      NUM_PARAMS
   };

   // The parameters necessary for the evaluation of this driver.
   vector<double> m_parameters;
};

} // end namespace Loki

#endif
