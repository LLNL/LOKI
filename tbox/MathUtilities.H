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
#ifndef included_tbox_MathUtilities
#define included_tbox_MathUtilities

#include <cmath>
#include <math.h>

#include <float.h>
#include <limits.h>
#include <stdlib.h>

namespace Loki {
namespace tbox {

/**
 * Class MathUtilities is a utility that provides some basic math-related
 * functions and routines for initializing data to POSIX constants like
 * INT_MAX, FLT_MAX, DBL_MAX, etc.
 */
template<class TYPE>
class MathUtilities
{
public:
   /*!
    * @brief Return max value for the template type.
    *
    * For boolean type, will return "true".  For dcomplex type,
    * will return a dcomplex value with both real and imaginary
    * parts set to the POSIX max value for type double.  For
    * other types, will return the POSIX max value for the type.
    */
   static TYPE
   getMax()
   {
      return s_max;
   }

   /*!
    * @brief Return min value for the template type.
    *
    * For boolean type, will return "false".  For dcomplex type,
    * will return a dcomplex value with both real and imaginary
    * parts set to the POSIX min value for type double.  For
    * other types, will return the POSIX min value for the type.
    */
   static TYPE
   getMin()
   {
      return s_min;
   }

private:
   static TYPE s_max;
   static TYPE s_min;
};

}
}

#endif
