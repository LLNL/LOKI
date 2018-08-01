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
#ifndef _SHAPED_RAMPED_COSINE_CURRENT_DRIVER_H_
#define _SHAPED_RAMPED_COSINE_CURRENT_DRIVER_H_

#include "CurrentDriver.H"
#include "ParmParse.H"

namespace Loki {

/**
 * Implementation of the CurrentDriver interface for a shaped ramped cosine
 * driver with the following functional form:
 *
 * @verbatim
 *  J_ext = J0*f(t)*g(x)*h(y)*cos(omega*(t-t0))
 *
 *          / 0.5+0.5*tanh(4.0*(2.0*(t-t0)/t_ramp-1.0))       for t0 <= t < t0+t_ramp
 *
 *         /
 *  f(t) =|   0.5-0.5*tanh(4.0*(2.0*(t-t0-t_ramp)/t_off-1.0)) for t0+t_ramp <= t < t0+t_ramp+t_off
 *         \
 *          \ 0
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
      ParmParse& a_pp);

   /**
    * @brief Destructor.
    */
   virtual
   ~ShapedRampedCosineCurrentDriver();

   //@{

   //! @name Implementation of CurrentDriver interface.

   /**
    * @brief Compute the values of the 3 components of the current at the
    * supplied time.
    *
    * @param[out] a_Jx The x component of the current.
    * @param[out] a_Jy The y component of the current.
    * @param[out] a_Jz The z component of the current.
    * @param[in] a_local_box Box describing the bounds of Jx, Jy, and Jz.
    * @param[in] a_domain Information about the problem domain.
    * @param[in] a_time The time at which the current is to be computed.
    */
   virtual
   void
   evaluate(
      RealArray& a_Jx,
      RealArray& a_Jy,
      RealArray& a_Jz,
      const tbox::Box& a_local_box,
      const ProblemDomain& a_domain,
      real a_time) const;

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
    * @param[in] a_name Name of a type of current driver.
    *
    * @return true if name describes this class' type.
    */
   static
   bool
   isType(
      const aString& a_name);

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
      ParmParse& a_pp);

   // The name of this type of current driver.
   static const aString s_CLASS_NAME;

   // Indices of specific parameters in m_parameters.
   enum {
      WIDTH,
      APPLY_DIR,
      SHAPE,
      OMEGA,
      J0,
      T0,
      TRAMP,
      TOFF,
      X0,
      NUM_PARAMS
   };

   // The parameters necessary for the evaluation of this driver.
   RealArray m_parameters;
};

} // end namespace Loki

#endif
