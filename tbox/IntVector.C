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
#include "IntVector.H"

namespace Loki {
namespace tbox {

bool IntVector::s_registered_shutdown = IntVector::initializeStatics();
IntVector * IntVector::s_zeros[Dimension::MAXIMUM_DIMENSION];
IntVector * IntVector::s_ones[Dimension::MAXIMUM_DIMENSION];

std::istream& operator >> (
   std::istream& s,
   IntVector& rhs)
{
   while (s.get() != '(') ;

   for (int i = 0; i < rhs.getDim(); i++) {
      s >> rhs(i);
      if (i < rhs.getDim() - 1)
         while (s.get() != ',') ;
   }

   while (s.get() != ')') ;

   return s;
}

std::ostream& operator << (
   std::ostream& s,
   const IntVector& rhs)
{
   s << '(';

   for (int i = 0; i < rhs.getDim(); i++) {
      s << rhs(i);
      if (i < rhs.getDim() - 1)
         s << ",";
   }
   s << ')';

   return s;
}

IntVector::IntVector():
   d_dim(Dimension::INVALID_DIMENSION)
{
//#ifdef DEBUG_INITIALIZE_UNDEFINED
   for (int i = 0; i < Dimension::MAXIMUM_DIMENSION; i++) {
      d_vector[i] = MathUtilities<int>::getMin();
   }
//#endif
}

IntVector::IntVector(
   const Dimension& dim):
   d_dim(dim)
{
   // an explicit setting Invalid is allowed.
   TBOX_DIM_ASSERT((!d_dim.isValid()) ||
      (d_dim > 0 && d_dim <= Dimension::MAXIMUM_DIMENSION));

//#ifdef DEBUG_INITIALIZE_UNDEFINED
   for (int i = 0; i < Dimension::MAXIMUM_DIMENSION; i++) {
      d_vector[i] = MathUtilities<int>::getMin();
   }
//#endif

}

IntVector::IntVector(
   const Dimension& dim,
   const int value):
   d_dim(dim)
{
   // an explicit setting Invalid is allowed.
   TBOX_DIM_ASSERT((!d_dim.isValid()) ||
      (d_dim > 0 && d_dim <= Dimension::MAXIMUM_DIMENSION));

   if (d_dim.isValid()) {
      for (int i = 0; i < d_dim; i++)
         d_vector[i] = value;

//#ifdef DEBUG_INITIALIZE_UNDEFINED
      for (int i = d_dim; i < Dimension::MAXIMUM_DIMENSION;
           i++) {
         d_vector[i] = MathUtilities<int>::getMin();
      }
//#endif
   } else {
      for (int i = 0; i < Dimension::MAXIMUM_DIMENSION; i++) {
         d_vector[i] = MathUtilities<int>::getMin();
      }
   }
}

IntVector::IntVector(
   const std::vector<int>& a):
   d_dim(static_cast<unsigned short>(a.size()))
{
//   TBOX_DIM_ASSERT(a.getSize() > 1 &&
//      a.getSize() <= Dimension::MAXIMUM_DIMENSION);
   TBOX_DIM_ASSERT(a.size() > 1 &&
      a.size() <= Dimension::MAXIMUM_DIMENSION);

   for (int i = 0; i < d_dim; i++)
      d_vector[i] = a[i];

//#ifdef DEBUG_INITIALIZE_UNDEFINED
   for (int i = d_dim; i < Dimension::MAXIMUM_DIMENSION; i++) {
      d_vector[i] = MathUtilities<int>::getMin();
   }
//#endif
}

IntVector::IntVector(
   const Dimension& dim,
   const int array[]):
   d_dim(dim)
{
   TBOX_DIM_ASSERT_CHECK_DIM(dim);

   for (int i = 0; i < d_dim; i++)
      d_vector[i] = array[i];
}

IntVector::IntVector(
   const IntVector& rhs):
   d_dim(rhs.getDim())
{
   /*
    * STL needs to be able to copy invalid values.
    */
   if (rhs.getDim().isValid()) {
      TBOX_DIM_ASSERT_CHECK_DIM(rhs.getDim());

      for (int i = 0; i < d_dim; i++)
         d_vector[i] = rhs.d_vector[i];
   }
}

IntVector::~IntVector()
{
}

void IntVector::putToDatabase(
   HDF_DataBase& database,
   const std::string& name) const
{
   database.put(d_vector, name, d_dim);
}

void IntVector::getFromDatabase(
   const HDF_DataBase& database,
   const std::string& name)
{
   database.get(d_vector, name, d_dim);
}
   
bool IntVector::initializeStatics()
{
   for (unsigned short d = 0; d < Dimension::MAXIMUM_DIMENSION; ++d) {
      s_zeros[d] = new IntVector(Dimension(static_cast<unsigned short>(d + 1)), 0);
   }

   for (unsigned short d = 0; d < Dimension::MAXIMUM_DIMENSION; ++d) {
      s_ones[d] = new IntVector(Dimension(static_cast<unsigned short>(d + 1)), 1);
   }

   return true;
}

void IntVector::freeStatics()
{
   for (int d = 0; d < Dimension::MAXIMUM_DIMENSION; ++d) {
      delete s_zeros[d];
   }

   for (int d = 0; d < Dimension::MAXIMUM_DIMENSION; ++d) {
      delete s_ones[d];
   }
}

}
}
