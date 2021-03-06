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
#include "CurrentDriverFactory.H"

#include "OvertureTypes.h"
#include "ParmParse.H"
#include "aString.H"

// Add new derived CurrentDriver headers here
#include "ShapedRampedCosineCurrentDriver.H"

namespace Loki {

CurrentDriver*
CurrentDriverFactory::create(
   ParmParse& a_pp,
   int a_driver_num)
{
   aString name("none");
   a_pp.query("name", name);
   CurrentDriver* driver = 0;

   if (name.matches("none")) {
      char buffer[100];
      sprintf(buffer,
         "Missing current driver number %i ... quitting",
         a_driver_num);
      OV_ABORT(buffer);
   }
   else {
      if (ShapedRampedCosineCurrentDriver::isType(name)) {
         driver = new ShapedRampedCosineCurrentDriver(a_pp);
      }
      // Add new cases here in this form:
      // else if (NewDriver::isType(name)) {
      //    driver = new NewDriver(a_pp);
      // }
      else {
         aString msg("Unknown external current \"" + name + "\" ... quitting");
         OV_ABORT(msg);
      }
   }
   return driver;
}

} // end namespace Loki
