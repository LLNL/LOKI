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
#ifndef _IC_INTERFACE_H_
#define _IC_INTERFACE_H_

#include "OvertureTypes.h"
#include "ProblemDomain.H"

namespace Loki {

/**
 * ICInterface is an abstract base class for all initial condition
 * implementations.
 */
class ICInterface
{
public:
   /**
    * @brief Default Constructor.
    */
   ICInterface();

   /**
    * @brief Destructor.
    */
   virtual
   ~ICInterface();

   /**
    * @brief Compute the initial conditions of u at the supplied time.
    *
    * @param[out] a_u The thing to compute.
    * @param[in] a_domain Information about the problem domain.
    * @param[in] a_grown_global_box The global box for the entire problem with
    *                               an additional set of ghosts because that's
    *                               how Overture does things.
    * @param[in] a_time The time at which u is to be computed.
    */
   virtual
   void
   set(
      RealArray& a_u,
      const ProblemDomain& a_domain,
      const tbox::Box& a_grown_global_box,
      real a_time) const = 0;

   /**
    * @brief Print the driver's parameters.
    */
   virtual
   void
   printParameters() const = 0;

   // FIXME: temporary accessor until we weed out other references
   /**
    * @brief Returns the user supplied initial condition input parameters.
    */
   virtual
   const RealArray&
   parameters() const = 0;

private:
   // Unimplemented copy constructor.
   ICInterface(
      const ICInterface& other);

   // Unimplemented assignment operator.
   ICInterface&
   operator = (
      const ICInterface& rhs);
};

} // end namespace Loki

#endif
