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
#include "Box.H"

namespace Loki {
namespace tbox {

bool Box::s_registered_shutdown = Box::initializeStatics();
Box * Box::s_emptys[Dimension::MAXIMUM_DIMENSION];
Box * Box::s_universes[Dimension::MAXIMUM_DIMENSION];

Box::Box():
   d_lo(Dimension::INVALID_DIMENSION, MathUtilities<int>::getMax()),
   d_hi(Dimension::INVALID_DIMENSION, MathUtilities<int>::getMin())
{
}

Box::Box(
   const Dimension& dim):
   d_lo(dim, MathUtilities<int>::getMax()),
   d_hi(dim, MathUtilities<int>::getMin())
{
   TBOX_DIM_ASSERT_CHECK_DIM(dim);
}

Box::Box(
   const IntVector& lower,
   const IntVector& upper):
   d_lo(lower),
   d_hi(upper)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(lower, upper);
}

Box::Box(
   const Box& box):
   d_lo(box.d_lo),
   d_hi(box.d_hi)
{
}

Box::~Box()
{
}

/*
 *************************************************************************
 *                                                                       *
 * Return the dimension of the box that is the longest.                  *
 *                                                                       *
 *************************************************************************
 */

int Box::longestDimension() const
{
   int max = upper(0) - lower(0);
   int dim = 0;

   for (int i = 1; i < getDim(); i++)
      if ((upper(i) - lower(i)) > max) {
         max = upper(i) - lower(i);
         dim = i;
      }
   return dim;
}

/*
 *************************************************************************
 *                                                                       *
 * Stream input/output operators: [(l0,...,ln),(u0,...,un)].             *
 *                                                                       *
 *************************************************************************
 */

#ifdef __INSURE__
#define NULL_STATEMENT if (0) int nullstatement = 0
#else
#define NULL_STATEMENT
#endif

std::istream& operator >> (
   std::istream& s,
   Box& box)
{
   TBOX_DIM_ASSERT_CHECK_DIM(box.getDim());

   while (s.get() != '[') ;
   s >> box.lower();
   while (s.get() != ',') NULL_STATEMENT;
   s >> box.upper();
   while (s.get() != ']') NULL_STATEMENT;
   return s;
}

std::ostream& operator << (
   std::ostream& s,
   const Box& box)
{
   TBOX_DIM_ASSERT_CHECK_DIM(box.getDim());

   if (box.empty()) {
      s << "[(),()]";
   } else {
      s << '[' << box.lower() << ',' << box.upper() << ']';
   }
   return s;
}

Box& Box::operator += (
   const Box& box)
{

   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   if (!box.empty()) {
      if (empty()) {
         *this = box;
      } else {
         d_lo.min(box.d_lo);
         d_hi.max(box.d_hi);
      }
   }
   return *this;
}

/*
 *************************************************************************
 *                                                                       *
 * Static member function called from coalesceWith().  It attempts to    *
 * recursively coalesce intervals individual dimensions in index space.  *
 * If it is possible to coalesce two intervals (defined by a proper      *
 * overlap or adjacency relationship), the value true is returned.       *
 * If this is impossible, false is returned.                             *
 *                                                                       *
 *************************************************************************
 */

bool Box::coalesceIntervals(
   const int* lo1,
   const int* hi1,
   const int* lo2,
   const int* hi2,
   const int dim)
{
   bool retval = false;
   if (dim == 1) {
      // interval 1 to the right of interval 2.
      if ((lo1[0] <= hi2[0] + 1) && (hi2[0] <= hi1[0])) {
         retval = true;
         return retval;
      }
      // interval 1 to the left of interval 2.
      if ((lo1[0] <= lo2[0]) && (lo2[0] <= hi1[0] + 1)) {
         retval = true;
         return retval;
      }
   } else {
      for (int id = 0; id < dim; id++) {
         if ((lo1[id] == lo2[id]) && (hi1[id] == hi2[id])) {
            int id2;
            int low1[Dimension::MAXIMUM_DIMENSION];
            int high1[Dimension::MAXIMUM_DIMENSION];
            int low2[Dimension::MAXIMUM_DIMENSION];
            int high2[Dimension::MAXIMUM_DIMENSION];
            for (id2 = 0; id2 < id; id2++) {
               low1[id2] = lo1[id2];
               high1[id2] = hi1[id2];
               low2[id2] = lo2[id2];
               high2[id2] = hi2[id2];
            }
            for (id2 = id + 1; id2 < dim; id2++) {
               int id1 = id2 - 1;
               low1[id1] = lo1[id2];
               high1[id1] = hi1[id2];
               low2[id1] = lo2[id2];
               high2[id1] = hi2[id2];
            }
            if (coalesceIntervals(low1, high1, low2, high2, dim - 1)) {
               retval = true;
               return retval;
            }
         }
      }
   }

   return retval;
}

/*
 *************************************************************************
 *                                                                       *
 * Return true if this box can be coalesced with the argument box,       *
 * and set this box to the union of the boxes.  Otherwise, return false  *
 * and leave this box as is.  Two boxes may be coalesced if their union  *
 * is a box.  This routine attempts to coalesce the boxes along          *
 * each coordinate direction using the coalesceIntervals() function.     *
 *                                                                       *
 *************************************************************************
 */

bool Box::coalesceWith(
   const Box& box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   bool retval = false;

   if (empty() || box.empty()) {
      retval = true;
      *this += box;
   } else {
      int id;
      const int* box_lo = &box.lower()[0];
      const int* box_hi = &box.upper()[0];
      int me_lo[Dimension::MAXIMUM_DIMENSION];
      int me_hi[Dimension::MAXIMUM_DIMENSION];
      for (id = 0; id < getDim(); id++) {
         me_lo[id] = d_lo(id);
         me_hi[id] = d_hi(id);
      }
      if (coalesceIntervals(box_lo, box_hi, me_lo, me_hi, getDim())) {
         retval = true;
      } else { // test for one box containing the other...
         // test whether me contains box.
         retval = true;
         id = 0;
         while (retval && (id < getDim())) {
            retval = ((me_lo[id] <= box_lo[id]) && (me_hi[id] >= box_hi[id]));
            id++;
         }
         if (!retval) { // me doesn't contain box; check other way around...
            retval = true;
            id = 0;
            while (retval && (id < getDim())) {
               retval = ((box_lo[id] <= me_lo[id])
                         && (box_hi[id] >= me_hi[id]));
               id++;
            }
         }
      }
   }

   if (retval) *this += box;

   return retval;
}

/*
 *************************************************************************
 *                                                                       *
 * Rotates a 3-Dimensional box 45*num_rotations degrees around the given *
 * and set this box to the union of the boxes.                           *
 *                                                                       *
 *************************************************************************
 */

void Box::rotateAboutAxis(
   const int axis,
   const int num_rotations)
{
   TBOX_DIM_ASSERT_CHECK_DIM(getDim());
   TBOX_ASSERT(axis < getDim());
   TBOX_ASSERT(getDim() == 3);

   const Dimension& dim(getDim());

   const int a = (axis + 1) % dim;
   const int b = (axis + 2) % dim;

   IntVector tmp_lo(dim);
   IntVector tmp_hi(dim);

   for (int j = 0; j < num_rotations; j++) {
      tmp_lo = d_lo;
      tmp_hi = d_hi;
      d_lo(a) = tmp_lo(b);
      d_lo(b) = -tmp_hi(a) - 1;
      d_hi(a) = tmp_hi(b);
      d_hi(b) = -tmp_lo(a) - 1;
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Rotate a box in the manner determined by the rotation number          *
 *                                                                       *
 *************************************************************************
 */

void Box::rotate(
   const int rotation_number)
{
   TBOX_DIM_ASSERT_CHECK_DIM(getDim());
   TBOX_ASSERT(getDim() == 2 || getDim() == 3);

   if (getDim() == 2) {
      for (int j = 0; j < rotation_number; j++) {
         IntVector tmp_lo(d_lo);
         IntVector tmp_hi(d_hi);

         d_lo(0) = tmp_lo(1);
         d_lo(1) = -tmp_hi(0) - 1;
         d_hi(0) = tmp_hi(1);
         d_hi(1) = -tmp_lo(0) - 1;
      }
   } else {

      if (getDim() == 3) {
         if (rotation_number == 0) {
            return;
         } else if (rotation_number == 1) {
            rotateAboutAxis(0, 3);
            rotateAboutAxis(2, 3);
         } else if (rotation_number == 2) {
            rotateAboutAxis(1, 1);
            rotateAboutAxis(2, 1);
         } else if (rotation_number == 3) {
            rotateAboutAxis(1, 2);
            rotateAboutAxis(0, 3);
         } else if (rotation_number == 4) {
            rotateAboutAxis(1, 3);
         } else if (rotation_number == 5) {
            rotateAboutAxis(2, 1);
         } else if (rotation_number == 6) {
            rotateAboutAxis(1, 1);
         } else if (rotation_number == 7) {
            rotateAboutAxis(0, 3);
         } else if (rotation_number == 8) {
            rotateAboutAxis(0, 2);
            rotateAboutAxis(2, 3);
         } else if (rotation_number == 9) {
            rotateAboutAxis(0, 3);
            rotateAboutAxis(2, 1);
         } else if (rotation_number == 10) {
            rotateAboutAxis(1, 2);
         } else if (rotation_number == 11) {
            rotateAboutAxis(0, 3);
            rotateAboutAxis(1, 3);
         } else if (rotation_number == 12) {
            rotateAboutAxis(2, 3);
         } else if (rotation_number == 13) {
            rotateAboutAxis(0, 1);
         } else if (rotation_number == 14) {
            rotateAboutAxis(0, 2);
            rotateAboutAxis(1, 1);
         } else if (rotation_number == 15) {
            rotateAboutAxis(0, 1);
            rotateAboutAxis(1, 3);
         } else if (rotation_number == 16) {
            rotateAboutAxis(0, 2);
            rotateAboutAxis(1, 2);
         } else if (rotation_number == 17) {
            rotateAboutAxis(0, 1);
            rotateAboutAxis(2, 1);
         } else if (rotation_number == 18) {
            rotateAboutAxis(0, 3);
            rotateAboutAxis(1, 1);
         } else if (rotation_number == 19) {
            rotateAboutAxis(0, 1);
            rotateAboutAxis(2, 3);
         } else if (rotation_number == 20) {
            rotateAboutAxis(0, 2);
         } else if (rotation_number == 21) {
            rotateAboutAxis(0, 2);
            rotateAboutAxis(2, 1);
         } else if (rotation_number == 22) {
            rotateAboutAxis(0, 2);
            rotateAboutAxis(1, 3);
         } else if (rotation_number == 23) {
            rotateAboutAxis(1, 2);
            rotateAboutAxis(0, 1);
         }
      }
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
bool Box::initializeStatics()
{
   TBOX_ASSERT(s_registered_shutdown == false);

   for (unsigned short d = 0; d < Dimension::MAXIMUM_DIMENSION; ++d) {
      Dimension dim(static_cast<unsigned short>(d + 1));
      s_emptys[d] = new Box(dim);

      /*
       * Note we can't use IntVector getMin, getMax here as that
       * would create a dependency between static initializers
       */
      s_universes[d] = new Box(
            IntVector(dim, MathUtilities<int>::getMin()),
            IntVector(dim, MathUtilities<int>::getMax()));
   }

   return true;
}

/*
 *************************************************************************
 *************************************************************************
 */
void Box::freeStatics()
{
   TBOX_ASSERT(s_registered_shutdown == true);

   for (int d = 0; d < Dimension::MAXIMUM_DIMENSION; ++d) {
      delete s_emptys[d];
      delete s_universes[d];
   }

   s_registered_shutdown = false;
}

}
}
