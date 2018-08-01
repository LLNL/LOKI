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
#ifndef _SHAPED_RAMPED_COSINE_DRIVER_H_
#define _SHAPED_RAMPED_COSINE_DRIVER_H_

#include "ElectricFieldDriver.H"
#include "ParmParse.H"
#include "Serializable.H"

#include <string>

namespace Loki {

/**
 * Implementation of the ElectricFieldDriver interface for a shaped rampled
 * cosine drive with the following functional form:
 *
 * @verbatim
 *  E_ext = E0*f(t)*g(x)*h(y)*cos(pi*x/xwidth-omega*(t-t0)+phase-0.5*alpha*(t-t0-t_res)^2)
 *
 *          / 0.5+0.5*tanh(4.0*(2.0*(t-t0)/t_ramp-1.0))       for t0 <= t < t0+t_ramp
 *         /
 *  f(t) =|   0.5-0.5*tanh(4.0*(2.0*(t-t0-t_ramp)/t_off-1.0)) for t0+t_ramp <= t < t0+t_ramp+t_off
 *         \
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
 *
 * This class needs to save state to and restore it from a restart file to
 * handle random, noisy phases.  We use a random number generator to compute
 * the phases.  So on restart we need to know how many times the random number
 * generator has already been called during this simulation.  We call the
 * random number generator that many times before computing the next random
 * phase so that they still have the right properties.
 */
class ShapedRampedCosineDriver : public ElectricFieldDriver, public Serializable
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_pp Where the input parameters are to be found.
    * @param[in] a_driver_num Which driver this one is.
    */
   ShapedRampedCosineDriver(
      ParmParse& a_pp,
      int a_driver_num);

   /**
    * @brief Destructor.
    */
   virtual
   ~ShapedRampedCosineDriver();

   //@{

   //! @name Implementation of ElectricFieldDriver interface.

   /**
    * @brief Compute the values of the 2 components of the electric field at
    * the supplied time.
    *
    * @param[out] a_Ex The x component of the electric field.
    * @param[out] a_Ey The y component of the electric field.
    * @param[in] a_fill_box Box describing the bounds of Ex and Ey.
    * @param[in] a_domain Information about the problem domain.
    * @param[in] a_time The time for this electric field evaluation.
    * @param[in] a_dt The time step for this electric field evaluation.
    * @param[in] a_stage The RK stage for this electric field evaluation.
    */
   virtual
   void
   evaluate(
      RealArray& a_Ex,
      RealArray& a_Ey,
      const tbox::Box& a_fill_box,
      const ProblemDomain& a_domain,
      real a_time,
      real a_dt,
      int a_stage) const;

   /**
    * @brief Print the driver's parameters.
    */
   virtual
   void
   printParameters() const;

   //@}

   //@{

   //! @name Implementation of Serializable interface.

   /**
    * @brief Read m_num_phase_evals from the restart database.
    *
    * @param[in] a_db The restart database.
    */
   virtual
   void
   getFromRestart(
      const HDF_DataBase& a_db);

   /**
    * @brief Write m_num_phase_evals to the restart database.
    *
    * @param[out] a_db The restart database.
    * @param[in] a_time Current simulation time.
    */
   virtual
   void
   putToRestart(
      HDF_DataBase& a_db,
      real a_time);

   //@}

   /**
    * @brief Returns true if name describes this class' type.
    *
    * @param[in] a_name Name of a type of electric field driver.
    *
    * @return true if name describes this class' type.
    */
   static
   bool
   isType(
      const aString& a_name);

private:
   // Unimplemented default constructor.
   ShapedRampedCosineDriver();

   // Unimplemented copy constructor.
   ShapedRampedCosineDriver(
      const ShapedRampedCosineDriver& other);

   // Unimplemented assignment operator.
   ShapedRampedCosineDriver&
   operator = (
      const ShapedRampedCosineDriver& rhs);

   // Parse the input parameters.
   void
   parseParameters(
      ParmParse& a_pp);

   // Computes the phase of the cos dependent part of the expression.
   void
   computePhase(
      real a_dt,
      int a_stage) const;

   // Compute the next random number.
   double
   rand() const
   {
      m_next = m_next * 1103515245 + 12345;
      double random_val = (unsigned)(m_next/65536) % 32768;
      random_val = random_val/32768.0;
      return random_val;
   }

   // The name of this type of electric field driver.
   static const aString s_CLASS_NAME;

   // Indices of specific parameters in m_parameters.
   // Note that PHASE_DECAY_TIME_STEPS and FWHM are not used in the fortran
   // to which m_parameters is passed.  We only need them in order to write
   // all the parameters in printParameters.
   enum {
      XWIDTH,
      YWIDTH,
      SHAPE,
      OMEGA,
      E0,
      T0,
      TRAMP,
      TOFF,
      XSHAPE,
      LWIDTH,
      X0,
      ALPHA,
      TRES,
      PHASE_DECAY_TIME_STEPS,
      FWHM,
      NUM_PARAMS
   };

   // The parameters necessary for the evaluation of this driver.
   RealArray m_parameters;

   // The name of this driver.  This is for restart capability.  There can
   // be multiple of these drivers for a given KineticSpecies so we construct
   // a unique name for each which is the name of the sub-database in the
   // restart file containing this object's state.
   std::string m_name;

   // If true code will compute a random, noisy phase.
   bool m_do_random_phase;

   // Phase in cos varying part of expression.
   mutable double m_phase;

   // Half time step phase in cos varying part of expression.
   mutable double m_phase_h;

   // The number of evaluations of the phase that have occurred.  We need to
   // track this and write it to restart so that we can continue the sequence
   // of random phases properly.
   mutable int m_num_phase_evals;

   // The underlying random value for the computation of the random phase.
   mutable unsigned long m_next;

   // Constant factors used in the computation of the random phase.
   double m_factor_1;
   double m_factor_2;
};

} // end namespace Loki

#endif
