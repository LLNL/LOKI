/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "Dimension.H"

namespace Loki {
namespace tbox {

/* Setup Invalid Dimension values */
unsigned short Dimension::INVALID_DIMENSION_VALUE =
   std::numeric_limits<unsigned short>::max();
const Dimension Dimension::INVALID_DIMENSION = Dimension(
   Dimension::INVALID_DIMENSION_VALUE);

Dimension::Dimension()
   : d_dim(Dimension::INVALID_DIMENSION_VALUE)
{
}

Dimension::Dimension(
   const unsigned short& dim)
   : d_dim(dim)
{
   TBOX_DIM_ASSERT((!isValid()) ||
      (d_dim > 0 && d_dim <= Dimension::MAXIMUM_DIMENSION));
}

Dimension::Dimension(
   const Dimension& dimension)
   : d_dim(dimension.d_dim)
{
   TBOX_DIM_ASSERT((!isValid()) ||
      (d_dim > 0 && d_dim <= Dimension::MAXIMUM_DIMENSION));
}

std::ostream& operator << (
   std::ostream& s,
   const Dimension& dim)
{
   s << (unsigned short)dim << 'D';
   return s;
}

}
}
