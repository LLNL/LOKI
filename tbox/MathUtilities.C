/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "MathUtilities.H"

namespace Loki {
namespace tbox {

template<>
 bool MathUtilities<bool>::s_max = true;
template<>
 bool MathUtilities<bool>::s_min = false;

template<>
 char MathUtilities<char>::s_max = CHAR_MAX;
template<>
 char MathUtilities<char>::s_min = CHAR_MIN;

template<>
 int MathUtilities<int>::s_max = INT_MAX;
template<>
 int MathUtilities<int>::s_min = INT_MIN;

template<>
 float MathUtilities<float>::s_max = FLT_MAX;
template<>
 float MathUtilities<float>::s_min = FLT_MIN;

template<>
 double MathUtilities<double>::s_max = DBL_MAX;
template<>
 double MathUtilities<double>::s_min = DBL_MIN;

}
}
