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
#include "CollisionOperatorFactory.H"

#include "OvertureTypes.h"
#include "ParmParse.H"
#include "aString.H"

// Add new derived CollisionOperator headers here
#include "PitchAngleCollisionOperator.H"
#include "RosenbluthCollisionOperator.H"

namespace Loki {

CollisionOperator*
CollisionOperatorFactory::create(
   ParmParse& a_pp,
   int a_solution_order)
{
   // Get the name of the collision operator.  This is required.
   aString name("none");
   a_pp.query("name", name);
   CollisionOperator* coll_op = 0;

   // If the user forgot the name, error.  Otherwise try to build what they
   // asked for.
   if (name.matches("none")) {
      OV_ABORT("Missing collision operator ... quitting");
   }
   else {
      if (PitchAngleCollisionOperator::isType(name)) {
         coll_op = new PitchAngleCollisionOperator(a_pp, a_solution_order);
      }
      else if (RosenbluthCollisionOperator::isType(name)) {
         coll_op = new RosenbluthCollisionOperator(a_pp, a_solution_order);
      }
      // Add new cases here in this form:
      // else if (NewOperator::isType(name)) {
      //    coll_op = new NewOperator(a_pp);
      // }
      else {
         aString msg("Unknown collision operator \"" + name + "\" ... quitting");
         OV_ABORT(msg);
      }
   }
   return coll_op;
}

} // end namespace Loki
