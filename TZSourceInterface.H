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
#ifndef _TZ_SOURCE_INTERFACE_H_
#define _TZ_SOURCE_INTERFACE_H_

#include "OvertureTypes.h"
#include "ProblemDomain.H"

namespace Loki {

/**
 * TZSourceInterface is an abstract base class for twilight zone Vlasov RHS
 * source implementations.  As noted elsewhere twilight zone means Method of
 * Manufactured Solutions.
 */
class TZSourceInterface
{
public:
   /**
    * @brief Destructor.
    */
   virtual
   ~TZSourceInterface();

   /**
    * @brief Compute the twilight zone source and add it to a_dist_func.
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
      real a_time) const = 0;

   /**
    * @brief Compute the error in the twilight zone source solution.
    *
    * @param[out] a_tz_error Error in the twilight zone source solution.
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
      real a_time) const = 0;

   /**
    * @brief Print the parameters defining the twilight zone source.
    */
   virtual
   void
   printParameters() const = 0;

protected:
   /**
    * @brief Default Constructor.
    */
   TZSourceInterface();

private:
   // Unimplemented copy constructor.
   TZSourceInterface(
      const TZSourceInterface& other);

   // Unimplemented assignment operator.
   TZSourceInterface&
   operator = (
      const TZSourceInterface& rhs);
};

} // end namespace Loki

#endif
