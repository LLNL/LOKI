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
#ifndef _ELECTRIC_FIELD_DRIVER_H_
#define _ELECTRIC_FIELD_DRIVER_H_

#include "Overture.h"
#include "OvertureTypes.h"

#include "tbox/Box.H"
#include "ProblemDomain.H"

namespace Loki {

/**
 * ElectricFieldDriver is an abstract base class for all ElectricFieldDriver
 * implementations.
 */
class ElectricFieldDriver
{
public:
   /**
    * @brief Default Constructor.
    */
   ElectricFieldDriver();

   /**
    * @brief Destructor.
    */
   virtual
   ~ElectricFieldDriver();

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
      int a_stage) const = 0;

   /**
    * @brief Print the driver's parameters.
    */
   virtual
   void
   printParameters() const = 0;

private:
   // Unimplemented copy constructor.
   ElectricFieldDriver(
      const ElectricFieldDriver& other);

   // Unimplemented assignment operator.
   ElectricFieldDriver&
   operator = (
      const ElectricFieldDriver& rhs);
};

} // end namespace Loki

#endif
