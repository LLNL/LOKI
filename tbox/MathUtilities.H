/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
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
