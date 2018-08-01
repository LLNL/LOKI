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
#ifndef _TRIG_TZ_SOURCE_H_
#define _TRIG_TZ_SOURCE_H_

#include "TZSourceInterface.H"
#include "ParmParse.H"

namespace Loki {

/**
 * Implements the interface of TZSourceInterface for a trigonometric twilight
 * zone source having the following functional form:
 *
 * @verbatim
 *  f = f + amp*cos(x)*cos(y)*cos(t)
 * @endverbatim
 */
class TrigTZSource : public TZSourceInterface
{
public:
   /**
    * @brief Default Constructor.
    *
    * @param[in] a_pp Where the input parameters are to be found.
    */
   TrigTZSource(
      ParmParse& a_pp);

   /**
    * @brief Destructor.
    */
   virtual
   ~TrigTZSource();

   //@{

   //! @name Implementation of TZSourceInterface interface.

   /**
    * @brief Compute the twilight zone source for this class and add it to the
    *        distribution function.
    *
    * @param[out] a_dist_func Distribution function.
    * @param[in] a_domain Information about the problem domain.
    * @param[in] a_time Current simulation time.
    */
   virtual
   void
   set(
      RealArray& a_dist_func,
      const ProblemDomain& a_domain,
      real a_time) const;

   /**
    * @brief Compute the error in the twilight zone source solution.
    *
    * @param[out] a_tz_error Errir in the twilight zone source solution.
    * @param[in] a_dist_func Distribution function.
    * @param[in] a_domain Information about the problem domain.
    * @param[in] a_time Current simulation time.
    */
   virtual
   void
   computeError(
      RealArray& a_tz_error,
      const RealArray& a_dist_func,
      const ProblemDomain& a_domain,
      real a_time) const;

   /**
    * @brief Print the parameters defining the initial condition.
    */
   virtual
   void
   printParameters() const;

   //@}

   /**
    * @brief Returns true if name describes this class' type.
    *
    * @param[in] a_name Name of a type of twilight zone source.
    *
    * @return true if name describes this class' type.
    */
   static
   bool
   isType(
      const aString& a_name);

private:
   // Unimplemented copy constructor.
   TrigTZSource(
      const TrigTZSource& other);

   // Unimplemented assignment operator.
   TrigTZSource&
   operator = (
      const TrigTZSource& rhs);

   // The name of this type of twilight zone source.
   static const aString s_CLASS_NAME;

   // Indices of specific parameters in m_dparameters.
   enum {
      AMP,
      NUM_DPARAMS};

   // The double valued parameters.
   RealArray m_dparameters;
};

} // end namespace Loki
#endif
