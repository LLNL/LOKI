/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "IntVector.H"
#include "../RestartReader.H"
#include "../RestartWriter.H"

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
   RestartWriter& writer,
   const std::string& name,
   bool a_write_data) const
{
   writer.writeIntegerArray(name, d_vector, d_dim, a_write_data);
}

void IntVector::getFromDatabase(
   RestartReader& reader,
   const std::string& name)
{
   reader.readIntegerArray(name, d_vector, d_dim);
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
