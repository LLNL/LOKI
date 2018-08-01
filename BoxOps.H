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
#ifndef __BOX_OPS__
#define __BOX_OPS__

#include "Overture.h"
#include "ParallelUtility.h"

#include "Directions.H"
#include "tbox/Dimension.H"
#include "tbox/IntVector.H"
#include "tbox/Box.H"

namespace Loki {

/**
 * BoxOps provides a collection of static functions for Box construction and
 * information extraction using Overture objects.
 */
class BoxOps
{
public:
   /**
    * @brief Returns the span of the Box, a_box, in the direction, a_dir, as an
    * Overture Range.
    *
    * @param[in] a_box The Box of interest.
    * @param[in] a_dir The direction of interest.
    *
    * @return The span of a_box in direction a_dir as an Overture Range.
    */
   static
   Range
   range(
      const tbox::Box& a_box,
      const int& a_dir)
   {
      return Range(a_box.lower(a_dir), a_box.upper(a_dir));
   }

   /**
    * @brief Returns the Box corresponding to the local part of the Array a_u.
    *
    * @param[in] a_u The Array of interest.
    *
    * @return The Box corresponding to the local part of a_u.
    */
   static
   tbox::Box
   getLocalBox(
      const RealArray& a_u)
   {
      const tbox::Dimension dim(
         static_cast<unsigned short>(a_u.numberOfDimensions()));
      tbox::IntVector lo(dim);
      tbox::IntVector hi(dim);
      for (int dir(0); dir<dim; dir++) {
         lo.operator[](dir) = a_u.getBase(dir);
         hi[dir] = a_u.getBound(dir);
      }
      return tbox::Box(lo, hi);
   }

   /**
    * @brief Returns the Box corresponding to the restriction of the local array
    *        to the computational domain with no boundary ghosts.
    *
    * @param[in] a_u The local array.
    * @param[in] a_uGlobal The global array.
    * @param[in] a_domain_box Box representing the entire domain.
    * @param[in] a_global_box Box representing the entire domain plus any
    *            ghosts.
    *
    * @return The Box corresponding to the restriction of the local array to the
    *         computational domain with no boundary ghosts.
    */
   static
   tbox::Box
   getOperationalBox(
      const RealArray& a_u,
      const realArray& a_uGlobal,
      const tbox::Box& a_domain_box,
      const tbox::Box& a_global_box)
   {
      const tbox::Dimension dim(a_domain_box.getDim());
      if (dim > 4) {
         OV_ABORT("Not defined for dim>4!");
      }
      bool includeGhost(false); // set to 1 to include parallel ghost
      Index J[4];
      for (int dir(0); dir < dim; ++dir) {
         J[dir] = range(a_global_box, dir);
      }
      for (int dir(dim); dir < 4; ++dir) {
         J[dir] = Range(0, 0);
      }
      const tbox::Dimension max_dim(4);
      tbox::IntVector lo(max_dim);
      tbox::IntVector hi(max_dim);

      ParallelUtility::getLocalArrayBounds(a_uGlobal,
         a_u,
         J[0], J[1], J[2], J[3],
         lo[0], hi[0],
         lo[1], hi[1],
         lo[2], hi[2],
         lo[3], hi[3],
         includeGhost);

      tbox::IntVector op_box_lo(dim);
      tbox::IntVector op_box_hi(dim);
      for (int dir(0); dir < dim; ++dir) {
         op_box_lo[dir] = max(a_domain_box.lower(dir), lo[dir]);
         op_box_hi[dir] = min(a_domain_box.upper(dir), hi[dir]);
      }
      return tbox::Box(op_box_lo, op_box_hi);
   }

   /**
    * @brief Return the IntVector representing only the directions in a_iv
    * that are true in a_remove_dir.
    *
    * @param[in] a_iv The original IntVector.
    * @param[in] a_remove_dir The directions to remove from a_iv.
    *
    * @return The IntVector corresponding to a_iv without the directions that
    * are true in a_remove_dir.
    */
   static
   tbox::IntVector
   restrictIntVector(
      const tbox::IntVector& a_iv,
      const Array<bool>& a_remove_dir)
   {
      unsigned short count(0);
      for (int dir(0); dir < a_iv.getDim(); ++dir) {
         if (a_remove_dir[dir]) {
            ++count;
         }
      }
      tbox::Dimension new_dim(count);

      count = 0;
      tbox::IntVector new_iv(new_dim);
      for (int dir(0); dir < a_iv.getDim(); ++dir) {
         if (a_remove_dir[dir]) {
            new_iv[count] = a_iv[dir];
            ++count;
         }
      }

      return new_iv;
   }

   /**
    * @brief Return the Box representing only the directions in a_box that are
    * true in a_remove_dir.
    *
    * @param[in] a_box The original Box.
    * @param[in] a_remove_dir The directions to remove from a_box.
    *
    * @return The Box corresponding to a_box without the directions that are
    * true in a_remove_dir.
    */
   static
   tbox::Box
   restrictBox(
      const tbox::Box& a_box,
      const Array<bool>& a_remove_dir)
   {
      tbox::IntVector new_lo(restrictIntVector(a_box.lower(), a_remove_dir));
      tbox::IntVector new_hi(restrictIntVector(a_box.upper(), a_remove_dir));
      return tbox::Box(new_lo, new_hi);
   }

   /**
    * @brief Convert an Array of Indexes into the Box represented by them.
    *
    * @param[in] a_index_array The Array of Indexes.
    *
    * @return The Box represented by the Indexes in a_index_array.
    */
   static
   tbox::Box
   indexArrayToBox(
      const Array<Index>& a_index_array)
   {
      const tbox::Dimension dim(
         static_cast<unsigned short>(a_index_array.length()));
      tbox::IntVector lo(dim);
      tbox::IntVector hi(dim);
      for (int dir(0); dir < dim; ++dir) {
         lo[dir] = a_index_array[dir].getBase();
         hi[dir] = a_index_array[dir].getBound();
      }
      return tbox::Box(lo, hi);
   }

   /**
    * @brief Convert a Box into an IndexBox.
    *
    * @param[in] a_box The Box to convert.
    *
    * @return The IndexBox corresponding to a_box.
    */
   static
   IndexBox
   boxToIndexBox(
      const tbox::Box& a_box)
   {
      switch (a_box.getDim()) {
         case 1:
            return IndexBox(a_box.lower(0), a_box.upper(0));
         case 2:
            return IndexBox(a_box.lower(0), a_box.upper(0),
               a_box.lower(1), a_box.upper(1));
         case 3:
            return IndexBox(a_box.lower(0), a_box.upper(0),
               a_box.lower(1), a_box.upper(1),
               a_box.lower(2), a_box.upper(2));
         case 4:
            return IndexBox(a_box.lower(0), a_box.upper(0),
               a_box.lower(1), a_box.upper(1),
               a_box.lower(2), a_box.upper(2),
               a_box.lower(3), a_box.upper(3));
         default:
            OV_ABORT("Not defined for dim>4!");
            return IndexBox();
      }
   }

private:
   // Unimplemented default constructor.
   BoxOps();

   // Unimplemented destructor.
   ~BoxOps();

   // Unimplemented copy constructor.
   BoxOps(
      const BoxOps& other);

   // Unimplemented assignment operator.
   BoxOps&
   operator = (
      const BoxOps& rhs);
};

} // end namespace Loki


#endif
