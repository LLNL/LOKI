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
#ifndef included_hier_Box
#define included_hier_Box

#include "IntVector.H"

#include <iostream>

namespace Loki {
namespace tbox {

/**
 * Class Box represents a n-dimensional box in the index space.  It is defined
 * by lower and upper bounds given by index objects.  The box semantics assume
 * that the box is cell-centered.  A cell-centered convention implies that the
 * index set covered by the box includes both the lower and upper bounds.
 */
class Box
{
public:
   /**
    * Creates an ``empty'' box.
    */
   explicit Box( const Dimension& dim);

   /**
    * Create a box describing the index space between lower and upper.  The
    * box is assumed to be cell centered and include all elements between lower
    * and upper, including the end points.
    */
   Box( const IntVector& lower,
        const IntVector& upper);

   /**
    * The copy constructor copies the index space of the argument box.
    */
   Box( const Box& box);

   /**
    * The destructor for Box.
    */
   ~Box();

   /**
    * The assignment operator copies the index space of the argument box.
    *
    * An assignment to an uninitialized box is allowed but assigning
    * from an uninitialized box will result in an assert.
    */
   Box&
   operator = (
      const Box& box)
   {
      /*
       * Allow assignment of to an uninitialized box but
       * not from an uninitialized box.
       */
      if (getDim().isValid()) {
         TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
      } else {
         TBOX_DIM_ASSERT_CHECK_DIM(box.getDim());
      }

      d_lo = box.d_lo;
      d_hi = box.d_hi;
      return *this;
   }

   /**
    * Return a non-const lower index of the box.
    */
   IntVector&
   lower()
   {
      return d_lo;
   }

   /**
    * Return a non-const upper index of the box.
    */
   IntVector&
   upper()
   {
      return d_hi;
   }

   /**
    * Return a const lower index of the box.
    */
   const IntVector&
   lower() const
   {
      return d_lo;
   }

   /**
    * Return a const upper index of the box.
    */
   const IntVector&
   upper() const
   {
      return d_hi;
   }

   /**
    * Return the i'th component (non-const) of the lower index of the box.
    */
   int&
   lower(
      const int i)
   {
      return d_lo(i);
   }

   /**
    * Return the i'th component (non-const) of the upper index of the box.
    */
   int&
   upper(
      const int i)
   {
      return d_hi(i);
   }

   /**
    * Return the i'th component (const) of the lower index of the box.
    */
   const int&
   lower(
      const int i) const
   {
      return d_lo(i);
   }

   /**
    * Return the i'th component (const) of the upper index of the box.
    */
   const int&
   upper(
      const int i) const
   {
      return d_hi(i);
   }

   /**
    * Set the index space represented by the box to empty.
    */
   void
   setEmpty()
   {
      const Dimension& dim(getDim());

      d_lo = IntVector(dim, MathUtilities<int>::getMax());
      d_hi = IntVector(dim, MathUtilities<int>::getMin());
   }

   /**
    * @brief Return whether the box is ``empty''.
    *
    * isEmpty() is preferred to match "is" standard syntax for
    * boolean methods.
    *
    * @see isEmpty()
    */
   bool
   empty() const
   {
      for (int i = 0; i < getDim(); i++) {
         if (d_hi(i) < d_lo(i)) return true;
      }

      return false;
   }

   /**
    * @brief Return whether the box is ``empty''.
    *
    * A box is empty if any of the lower bounds is greater than the
    * corresponding upper bound.  An empty box has a size of zero.
    */
   bool
   isEmpty() const
   {
      for (int i = 0; i < getDim(); i++) {
         if (d_hi(i) < d_lo(i)) return true;
      }

      return false;
   }

   /**
    * Return the number of cells (an integer) represented by the box in
    * the given coordinate direction.
    */
   int
   numberCells(
      const int i) const
   {
      if (empty()) {
         return 0;
      } else {
         return d_hi(i) - d_lo(i) + 1;
      }
   }

   /**
    * Return the number of cells (a vector of integers) represented by
    * the box in every coordinate direction.
    */
   IntVector
   numberCells() const
   {
      if (empty()) {
         return IntVector::Zero(getDim());
      } else {
         return d_hi - d_lo + 1;
      }
   }

   /**
    * Calculate the number of indices represented by the box.  If the box
    * is empty, then the number of index points within the box is zero.
    */
   long int
   size() const
   {
      long int mysize = 0;
      if (!empty()) {
         mysize = 1;
         for (int i = 0; i < getDim(); i++)
            mysize *= (d_hi(i) - d_lo(i) + 1);
      }

      return mysize;
   }

   /**
    *  Return the dimension of the box that is longest.
    */
   int
   longestDimension() const;

   /**
    * Given an index, calculate the offset of the index into the box.
    * This function assumes column-major (e.g., Fortran) ordering of
    * the indices within the box.  This operation is a convenience
    * function for later array indexing operations.
    */
   int
   offset(
      const IntVector& p) const
   {
      int myoffset = 0;

      for (int i = getDim() - 1; i > 0; i--) {
         myoffset =
            (d_hi(i - 1) - d_lo(i - 1) + 1) * (p(i) - d_lo(i) + myoffset);
      }
      myoffset += p(0) - d_lo(0);

      return myoffset;
   }

   /**
    * Given an offset, calculate the index of the offset into the box.
    * This function assumes column-major (e.g., Fortran) ordering of
    * the indices within the box.  This operation is a convenience
    * function for later array indexing operations.
    */
   IntVector
   index(
      const int offset) const
   {
      TBOX_ASSERT(offset >= 0);
      TBOX_ASSERT(offset <= size());

      IntVector n(getDim());
      IntVector index(getDim());

      n = numberCells();

      int remainder = offset;

      for (int d = getDim() - 1; d > -1; d--) {
         /* Compute the stride for indexing */
         int stride = 1;
         for (int stride_dim = 0; stride_dim < d; stride_dim++) {
            stride *= n[stride_dim];
         }

         /* Compute the local index */
         index[d] = remainder / stride;
         remainder -= index[d] * stride;

         /* Compute the global index */
         index[getDim()] += lower(d);
      }

      IntVector idx(index);

      return idx;
   }

   /**
    * Check whether an index lies within the bounds of the box.
    */
   bool
   contains(
      const IntVector& p) const
   {
      for (int i = 0; i < getDim(); i++) {
         if ((p(i) < d_lo(i)) || (p(i) > d_hi(i))) return false;
      }
      return true;
   }

   /**
    * Check whether a given box lies within the bounds of the box.
    *
    * If @c b is empty, always return true.
    */
   bool
   contains(
      const Box& b) const
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, b);

      if (b.empty()) return true;

      if (!contains(b.lower())) return false;

      if (!contains(b.upper())) return false;

      return true;
   }

   /**
    * Check whether two boxes represent the same portion of index space.
    */
   int
   operator == (
      const Box& box) const
   {
      return ((d_lo == box.d_lo) && (d_hi == box.d_hi)) ||
              (empty() && box.empty());
   }

   /**
    * Check whether two boxes cover different portions of index space.
    */
   int
   operator != (
      const Box& box) const
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

      return !(*this == box);
   }

   /**
    * Calculate the intersection of the index spaces of two boxes.  The
    * intersection with an empty box always yields an empty box.
    */
   Box
   operator * (
      const Box& box) const
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

      Box both(*this);
      both.d_lo.max(box.d_lo);
      both.d_hi.min(box.d_hi);
      return both;
   }

   /**
    * Return true if two boxes have a non-empty intersection.
    * Otherwise, return false.
    */
   bool
   intersects(
      const Box& box) const
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

      for (int i = 0; i < getDim(); i++) {
         if (std::max(d_lo(i), box.d_lo(i)) > std::min(d_hi(i), box.d_hi(i))) {
            return false;
         }
      }
      return true;
   }

   /**
    * Calculate the bounding box of two boxes.  Note that this is not
    * the union of the two boxes (since union is not closed over boxes),
    * but rather the smallest box that contains both boxes.
    */
   Box
   operator + (
      const Box& box) const
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

      Box bbox(*this);
      bbox += box;
      return bbox;
   }

   /**
    * Increase the bounding box to include the argument box.
    */
   Box&
   operator += (
      const Box& box);

   /**
    * Return true if this box can be coalesced with the argument box,
    * and set this box to the union of the boxes.  Otherwise, return false
    * and leave boxes as is.  Two boxes may be coalesced if their
    * union is a box (recall that index set union is not closed over boxes).
    * If either box is empty, then the return value is true and this box
    * becomes the union of the two.
    */
   bool
   coalesceWith(
      const Box& box);

   /**
    * Grow a box by the specified ghost cell width.  The lower bound is
    * decremented by the width, and the upper bound is incremented by the
    * width.  All dimensions are grown by the corresponding component in
    * the IntVector; ghost cell widths may be different in each dimension.
    * Negative ghost cell widths will shrink the box.
    */
   void
   grow(
      const IntVector& ghosts)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);

      if (!empty()) {
         d_lo -= ghosts;
         d_hi += ghosts;
      }
   }

   /**
    * Grow a box by the specified ghost cell width in the given coordinate
    * direction in index space.  The lower bound is decremented by the
    * width, and the upper bound is incremented by the width.  Note that
    * negative ghost cell widths will shrink the box.
    */
   void
   grow(
      const int direction,
      const int ghosts)
   {
      TBOX_ASSERT((direction >= 0) && (direction < getDim()));

      if (!empty()) {
         d_lo(direction) -= ghosts;
         d_hi(direction) += ghosts;
      }
   }

   /**
    * Similar to grow() functions. However, box is only grown in lower
    * directions (i.e., only lower index is changed).
    */
   void
   growLower(
      const IntVector& ghosts)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);

      if (!empty()) {
         d_lo -= ghosts;
      }
   }

   /**
    * Similar to grow() functions. However, box is only grown in lower
    * bound of given direction in index space.
    */
   void
   growLower(
      const int direction,
      const int ghosts)
   {
      TBOX_ASSERT((direction >= 0) && (direction < getDim()));

      if (!empty()) {
         d_lo(direction) -= ghosts;
      }
   }

   /**
    * Similar to grow() function. However, box is only grown in upper
    * directions (i.e., only upper index is changed).
    */
   void
   growUpper(
      const IntVector& ghosts)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);

      if (!empty()) {
         d_hi += ghosts;
      }
   }

   /**
    * Similar to grow() functions. However, box is only grown in upper
    * bound of given direction in index space.
    */
   void
   growUpper(
      const int direction,
      const int ghosts)
   {
      TBOX_ASSERT((direction >= 0) && (direction < getDim()));

      if (!empty()) {
         d_hi(direction) += ghosts;
      }
   }

   /**
    * Similar to growUpper() and growLower() functions. However, box is
    * lengthened (never shortened).  The sign of @c ghosts refer to whether
    * the box is lengthened in the upper or lower side.
    */
   void
   lengthen(
      const int direction,
      const int ghosts)
   {
      TBOX_ASSERT((direction >= 0) && (direction < getDim()));

      if (!empty()) {
         if (ghosts > 0) {
            d_hi(direction) += ghosts;
         } else {
            d_lo(direction) += ghosts;
         }
      }
   }

   /**
    * Similar to growUpper() and growLower() functions. However, box is
    * shortened (never lengthened).  The sign of @c ghosts refer to whether
    * the box is shortened in the upper or lower side.
    */
   void
   shorten(
      const int direction,
      const int ghosts)
   {
      TBOX_ASSERT((direction >= 0) && (direction < getDim()));

      if (!empty()) {
         if (ghosts > 0) {
            d_hi(direction) -= ghosts;
         } else {
            d_lo(direction) -= ghosts;
         }
      }
   }

   /**
    * Shift a box by the specified amount (a vector of integers).
    * The new box is located at (lower+offset, upper+offset).
    */
   void
   shift(
      const IntVector& offset)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, offset);

      d_lo += offset;
      d_hi += offset;
   }

   /**
    * Similar to shift() function above, but shift occurs only in specified
    * direction in index space.  The new box is located at (lower+offset,
    * upper+offset) in that direction.
    */
   void
   shift(
      const int direction,
      const int offset)
   {
      TBOX_ASSERT((direction >= 0) && (direction < getDim()));

      d_lo(direction) += offset;
      d_hi(direction) += offset;
   }

   /**
    * Rotate 90 degrees around origin.  Currently works for 2D only.
    */
   void
   rotate(
      const int rotation_number);

   /**
    * Refine the index space of a box by specified vector ratio.  Each
    * component of the box is multiplied by the refinement ratio,
    * then @c (ratio-1) is added to the upper corner.
    */
   void
   refine(
      const IntVector& ratio)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ratio);

      bool negative_ratio = false;
      for (int d = 0; d < getDim(); d++) {
         if (ratio(d) < 0) {
            negative_ratio = true;
            break;
         }
      }

      if (!negative_ratio) {
         d_lo *= ratio;
         d_hi = d_hi * ratio + (ratio - 1);
      } else {
         for (int i = 0; i < getDim(); i++) {
            if (ratio(i) > 0) {
               d_lo(i) *= ratio(i);
               d_hi(i) = d_hi(i) * ratio(i) + (ratio(i) - 1);
            } else {
               d_lo(i) = coarsen(d_lo(i), -ratio(i));
               d_hi(i) = coarsen(d_hi(i), -ratio(i));
            }
         }
      }
   }

   /**
    * Coarsen the index space of a box by specified vector ratio.  Each
    * component is divided by the specified coarsening ratio and rounded
    * (if necessary) such that the coarsened box contains the cells that
    * are the parents of the refined box.  In other words, refining a
    * coarsened box will always yield a box that is equal to or larger
    * than the original box.
    */
   void
   coarsen(
      const IntVector& ratio)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ratio);

      for (int i = 0; i < getDim(); i++) {
         d_lo(i) = coarsen(d_lo(i), ratio(i));
         d_hi(i) = coarsen(d_hi(i), ratio(i));
      }
   }

   /**
    * Utility function to grow a box by the specified vector ghost cell
    * width.  A new box is returned and the argument is not changed.
    */
   static Box
   grow(
      const Box& box,
      const IntVector& ghosts)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(box, ghosts);

      Box tmp = box;
      tmp.grow(ghosts);
      return tmp;
   }

   /**
    * Utility function to shift a box by the specified offset.  A new
    * box is returned and the argument is not changed.
    */
   static Box
   shift(
      const Box& box,
      const IntVector& offset)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(box, offset);

      return Box(box.lower() + offset, box.upper() + offset);
   }

   /**
    * Utility function to refine the index space of a box by the specified
    * refinement ratio.  A new box is returned and the argument is not changed.
    */
   static Box
   refine(
      const Box& box,
      const IntVector& ratio)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(box, ratio);

      Box tmp = box;
      tmp.refine(ratio);
      return tmp;
   }

   /**
    * Utility function to coarsen the index space of a box by the specified
    * coarsening ratio.  A new box is returned and the argument is not changed.
    */
   static Box
   coarsen(
      const Box& box,
      const IntVector& ratio)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(box, ratio);

      Box tmp = box;
      tmp.coarsen(ratio);
      return tmp;
   }

   /**
    * Return the dimension of this object.
    */
   const Dimension&
   getDim() const
   {
      return d_lo.getDim();
   }

   /**
    * Read the box description in the form [L,U], where L and U are the
    * lower and upper bounds of the box.
    */
   friend std::istream&
   operator >> (
      std::istream& s,
      Box& box);

   /**
    * Output the box description in the form [L,U], where L and U are the
    * lower and upper bounds of the box.
    */
   friend std::ostream&
   operator << (
      std::ostream& s,
      const Box& box);

   /*!
    * @brief Return an empty Box of the specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const Box&
   getEmptyBox(
      const Dimension& dim)
   {
      TBOX_DIM_ASSERT_CHECK_DIM(dim);
      return *(s_emptys[dim - 1]);
   }

   /**
    * Returns a Box that represents the maximum allowed index extents
    * for a given dimension.   The "universe" that can be represented.
    */
   static const Box&
   getUniverse(
      const Dimension& dim)
   {
      TBOX_DIM_ASSERT_CHECK_DIM(dim);
      return *(s_universes[dim - 1]);
   }

   void
   static
   shutdown()
   {
      TBOX_ASSERT(s_registered_shutdown == true);
      for (int i = 0; i <Dimension::MAXIMUM_DIMENSION; ++i) {
         delete s_emptys[i];
         delete s_universes[i];
      }
      s_registered_shutdown = false;
   }

   friend class std::vector<Box>;

private:
   /**
    * The default constructor creates an uninitialized box.
    *
    * This should never be invoked, it will cause assertions
    */
   Box();

   static int
   coarsen(
      const int index,
      const int ratio)
   {
      return index < 0 ? (index + 1) / ratio - 1 : index / ratio;
   }

   static bool
   coalesceIntervals(
      const int* lo1,
      const int* hi1,
      const int* lo2,
      const int* hi2,
      const int dim);

   void
   rotateAboutAxis(
      const int axis,
      const int num_rotations);

   /*!
    * @brief Initialize static objects and register shutdown routine.
    */
   static bool
   initializeStatics();

   /*!
    * @brief Method registered with ShutdownRegister to cleanup statics.
    */
   static void
   freeStatics();

   IntVector d_lo;
   IntVector d_hi;

   static bool s_registered_shutdown;

   /*
    * Array of empty boxes for each dimension.  Preallocated
    * as a performance enhancement to avoid constructing
    * them in multiple places.
    */

   static Box* s_emptys[Dimension::MAXIMUM_DIMENSION];

   /*
    * Box that represents the maximum allowed index extents,
    * the "universe" that can be represented.
    */
   static Box* s_universes[Dimension::MAXIMUM_DIMENSION];
};

}
}


#define BOX1D_TO_FORT(__box)                       \
                   __box.lower(0), __box.upper(0)

#define BOX2D_TO_FORT(__box)                       \
                   __box.lower(0), __box.upper(0), \
                   __box.lower(1), __box.upper(1)

#define BOX3D_TO_FORT(__box)                       \
                   __box.lower(0), __box.upper(0), \
                   __box.lower(1), __box.upper(1), \
                   __box.lower(2), __box.upper(2)

#define BOX4D_TO_FORT(__box) \
                   __box.lower(0), __box.upper(0), \
                   __box.lower(1), __box.upper(1), \
                   __box.lower(2), __box.upper(2), \
                   __box.lower(3), __box.upper(3)

#endif
