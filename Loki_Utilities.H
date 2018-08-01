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
#ifndef _LOKI_UTILITIES_H
#define _LOKI_UTILITIES_H

#include "ProblemDomain.H"

/*!
 * A null use of a variable, use to avoid GNU compiler
 * warnings about unused variables.
 */
#define NULL_USE(variable)                               \
   do {                                                  \
      if (0) { char* temp = (char *)&variable; temp++; } \
   } while (0)

namespace Loki {

/**
 * Application wide utility functions.
 */
class Loki_Utilities
{
public:
   /**
    * @brief Enforce periodicity on a_f in periodic directions.
    *
    * As only the configuration space directions are periodic this only operates
    * on X1 and X2.  It's also written to deal with 2D fields with the specified
    * number of components.
    *
    * @param[out] a_f Array being made periodic.
    * @param[in] a_domain The global problem domain.
    * @param[in] a_n_ghosts The number of ghosts in each direction.
    * @param[in] a_depth The number of components in a_f.
    */
   static
   void
   fixPeriodicity(
      realArray& a_f,
      const ProblemDomain& a_domain,
      const tbox::IntVector& a_n_ghosts,
      int a_depth);

   /**
    * @brief Turn a boolean into an int, perform a logical and reduction on it,
    * and return the result as a boolean.
    *
    * @param[in] a_bool The boolean to be reduced.
    */
   static
   bool
   reduceBoolean(
      bool a_bool);

private:
   // Unimplemented default constructor.
   Loki_Utilities();

   // Unimplemented destructor.
   ~Loki_Utilities();

   // Unimplemented copy constructor.
   Loki_Utilities(
      const Loki_Utilities& other);

   // Unimplemented assignment operator.
   Loki_Utilities&
   operator = (
      const Loki_Utilities& rhs);

};

} // end namespace Loki

#endif
