/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "ICInterface.H"

namespace Loki {

ICInterface::ICInterface(
   const ProblemDomain& a_domain)
  : m_domain(a_domain),
    m_vflowinitx(0.0),
    m_vflowinity(0.0),
    m_factorable(false)
{
}


ICInterface::~ICInterface()
{
}

} // end namespace Loki


// Prototype for Fortran callback to get the species initial condition at a
// point.  This must have C linkage to avoid name mangling.
extern "C"
{
   double
   initialconditionatpoint_(Loki::ICInterface** a_ic,
      int* a_i1,
      int* a_i2,
      int* a_i3,
      int* a_i4);
}


// Fortran callback to get species initial condition at a point.  Just gets the
// pointer to the initial condition object that has been laundered through
// Fortran and invokes the virtual function to get the initial condition value.
double
initialconditionatpoint_(Loki::ICInterface** a_ic,
   int* a_i1,
   int* a_i2,
   int* a_i3,
   int* a_i4)
{
   Loki::ICInterface* ic = *a_ic;
   return ic->getIC_At_Pt(*a_i1, *a_i2, *a_i3, *a_i4);
}
