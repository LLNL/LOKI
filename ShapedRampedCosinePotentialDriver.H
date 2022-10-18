/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _SHAPED_RAMPED_COSINE_POTENTIAL_DRIVER_H_
#define _SHAPED_RAMPED_COSINE_POTENTIAL_DRIVER_H_

#include "ElectricPotentialDriver.H"
#include "LokiInputParser.H"
#include <vector>

namespace Loki {

/**
 * Implementation of ElectricPotentialDriver interface for a shaped rampled
 * cosine driver with the following functional form:
 *
 * @verbatim
 *  phi = phi_0*f(t)*g(x)*h(y)*cos(pi*x/xwidth-omega*(t-t0))
 *
 *          / 0.5+0.5*tanh(4.0*(2.0*(t-t0)/t_rampup-1.0))     for t0 <= t < t0+t_rampup
 *         /  0.5-0.5*tanh(4.0)                               for t0+t_rampup <= t < t0+t_rampup+t_hold
 *  f(t) =|
 *         \  0.5-0.5*tanh(4.0*(2.0*(t-t0-t_rampup-t_hold)/t_rampdown-1.0)) for t0+t_rampup+t_hold <= t < t0+t_rampup+t_hold+t_rampdown
 *          \ 0                                               else
 *
 *         / 1-x_shape*(sin(pi*(x-x0)/lwidth))**2             for |x-x0| < lwidth/2
 *  g(x) =|
 *         \ 1-x_shape                                        else
 *
 *         / 1-shape*(sin(pi*y/ywidth))**2                    for |y| < ywidth/2
 *  h(y) =|
 *         \ 1-shape                                          else
 * @endverbatim
 */
class ShapedRampedCosinePotentialDriver : public ElectricPotentialDriver
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_pp Where the input parameters are to be found.
    */
   ShapedRampedCosinePotentialDriver(
      LokiInputParser& a_pp);

   /**
    * @brief Destructor.
    */
   virtual
   ~ShapedRampedCosinePotentialDriver();

   //@{

   //! @name Implementation of ElectricPotentialDriver interface.

   /**
    * @brief Compute the value of the electric potential at the supplied time.
    *
    * @param[out] a_phi The electric potential.
    * @param[in] a_domain Information about the problem domain.
    * @param[in] a_time The time at which the potential is to be computed.
    */
   virtual
   void
   evaluate(
      ParallelArray& a_phi,
      const ProblemDomain& a_domain,
      double a_time) const;

   /**
    * @brief Print the driver's parameters.
    */
   virtual
   void
   printParameters() const;

   //@}

   /**
    * @brief Returns true if name describes this class' type.
    *
    * @param[in] a_name Name of a type of electric potential driver.
    *
    * @return true if name describes this class' type.
    */
   static
   bool
   isType(
      const string& a_name);

private:
   // Unimplemented default constructor.
   ShapedRampedCosinePotentialDriver();

   // Unimplemented copy constructor.
   ShapedRampedCosinePotentialDriver(
      const ShapedRampedCosinePotentialDriver& other);

   // Unimplemented assignment operator.
   ShapedRampedCosinePotentialDriver&
   operator = (
      const ShapedRampedCosinePotentialDriver& rhs);

   // Parse the input parameters.
   void
   parseParameters(
      LokiInputParser& a_pp);

   // The name of this type of electric potential driver.
   static const string s_CLASS_NAME;

   // Indices of specific parameters in m_parameters.
   enum {
      XWIDTH,
      YWIDTH,
      SHAPE,
      OMEGA,
      AMP,
      T0,
      TRAMPUP,
      THOLD,
      TRAMPDOWN,
      XSHAPE,
      LWIDTH,
      X0,
      NUM_PARAMS
   };

   // The parameters necessary for the evaluation of this driver.
   vector<double> m_parameters;
};

} // end namespace Loki

#endif
