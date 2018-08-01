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
#include "SumOp.H"
#include "BoxOps.H"
#include "Loki_Utilities.H"

namespace Loki {

void
checkInput(
   const tbox::Dimension& a_dim,
   const RealArray& a_source,
   const Array<bool>& a_collapse_dir)
{
   // If we're about to do something that won't work, stop.
   if (a_dim > MAX_DIMENSION) {
      OV_ABORT("a_dim > MAX_DIMENSION");
   }
   if (a_dim != a_source.numberOfDimensions()) {
      OV_ABORT("a_dim != a_source dimension");
   }
   if (a_dim != a_collapse_dir.length()) {
      OV_ABORT("a_dim != a_collapse_dir length");
   }
}


tbox::Dimension
computeNewDimension(
   const tbox::Dimension& a_old_dim,
   const Array<bool>& a_collapse_dir)
{
   // Basically how many dimensions of a_old_dim are not being collapsed.
   unsigned short n_lost_dims(0);
   for (int dir(0); dir < a_old_dim; ++dir) {
      if (a_collapse_dir[dir]) {
         ++n_lost_dims;
      }
   }
   tbox::Dimension new_dim(
      static_cast<unsigned short>(a_old_dim - n_lost_dims));
   if (new_dim >= MAX_DIMENSION) {
      OV_ABORT("new dim >= MAX_DIMENSION");
   }
   return new_dim;
}


void
resultingDirections(
   tbox::IntVector& a_rdir,
   const Array<bool>& a_collapse_dir)
{
   // Which directions are not collapsed.
   const tbox::Dimension dim(
      static_cast<unsigned short>(a_collapse_dir.length()));
   int count(0);
   for (int dir(0); dir < dim; ++dir) {
      if (!a_collapse_dir[dir]) {
         a_rdir[count] = dir;
         ++count;
      }
   }
}


void
resizeResultArray(
   RealArray& a_result,
   const tbox::Box& a_box,
   const tbox::IntVector& a_rdir)
{
   // Forms ranges from a_box for each uncollapsed direction.
   const tbox::Dimension& dim(a_rdir.getDim());
   if (dim == 0) {
      a_result.resize(1);
   }
   else if (dim == 1) {
      Range r0(a_box.lower(a_rdir[0]), a_box.upper(a_rdir[0]));
      a_result.redim(r0);
   }
   else if (dim == 2) {
      Range r0(a_box.lower(a_rdir[0]), a_box.upper(a_rdir[0]));
      Range r1(a_box.lower(a_rdir[1]), a_box.upper(a_rdir[1]));
      a_result.redim(r0, r1);
   }
   else if (dim == 3) {
      Range r0(a_box.lower(a_rdir[0]), a_box.upper(a_rdir[0]));
      Range r1(a_box.lower(a_rdir[1]), a_box.upper(a_rdir[1]));
      Range r2(a_box.lower(a_rdir[2]), a_box.upper(a_rdir[2]));
      a_result.redim(r0, r1, r2);
   }
}


void
sum_reduce_1d_to_0d(
   RealArray& a_result,
   const RealArray& a_source,
   const tbox::Box& a_box)
{
   a_result = 0;
   tbox::IntVector i(a_box.lower());
   for (i[0] = a_box.lower(0); i[0] <= a_box.upper(0); ++i[0]) {
      a_result += a_source(i[0]);
   }
}


void
sum_reduce_2d_to_0d(
   RealArray& a_result,
   const RealArray& a_source,
   const tbox::Box& a_box)
{
   a_result = 0;
   tbox::IntVector i(a_box.lower());
   for (i[1] = a_box.lower(1); i[1] <= a_box.upper(1); ++i[1]) {
      for (i[0] = a_box.lower(0); i[0] <= a_box.upper(0); ++i[0]) {
         a_result += a_source(i[0], i[1]);
      }
   }
}


void
sum_reduce_2d_to_1d(
   RealArray& a_result,
   const RealArray& a_source,
   const tbox::Box& a_box,
   const tbox::IntVector& a_rdir)
{
   a_result = 0;
   tbox::IntVector i(a_box.lower());
   for (i[1] = a_box.lower(1); i[1] <= a_box.upper(1); ++i[1]) {
      for (i[0] = a_box.lower(0); i[0] <= a_box.upper(0); ++i[0]) {
         a_result(i[a_rdir[0]]) += a_source(i[0], i[1]);
      }
   }
}


void
sum_reduce_3d_to_0d(
   RealArray& a_result,
   const RealArray& a_source,
   const tbox::Box& a_box)
{
   a_result = 0;
   tbox::IntVector i(a_box.lower());
   for (i[2] = a_box.lower(2); i[2] <= a_box.upper(2); ++i[2]) {
      for (i[1] = a_box.lower(1); i[1] <= a_box.upper(1); ++i[1]) {
         for (i[0] = a_box.lower(0); i[0] <= a_box.upper(0); ++i[0]) {
            a_result += a_source(i[0], i[1], i[2]);
         }
      }
   }
}


void
sum_reduce_3d_to_1d(
   RealArray& a_result,
   const RealArray& a_source,
   const tbox::Box& a_box,
   const tbox::IntVector& a_rdir)
{
   a_result = 0;
   tbox::IntVector i(a_box.lower());
   for (i[2] = a_box.lower(2); i[2] <= a_box.upper(2); ++i[2]) {
      for (i[1] = a_box.lower(1); i[1] <= a_box.upper(1); ++i[1]) {
         for (i[0] = a_box.lower(0); i[0] <= a_box.upper(0); ++i[0]) {
            a_result(i[a_rdir[0]]) += a_source(i[0], i[1], i[2]);
         }
      }
   }
}


void
sum_reduce_3d_to_2d(
   RealArray& a_result,
   const RealArray& a_source,
   const tbox::Box& a_box,
   const tbox::IntVector& a_rdir)
{
   a_result = 0;
   tbox::IntVector i(a_box.lower());
   for (i[2] = a_box.lower(2); i[2] <= a_box.upper(2); ++i[2]) {
      for (i[1] = a_box.lower(1); i[1] <= a_box.upper(1); ++i[1]) {
         for (i[0] = a_box.lower(0); i[0] <= a_box.upper(0); ++i[0]) {
            a_result(i[a_rdir[0]], i[a_rdir[1]]) += a_source(i[0], i[1], i[2]);
         }
      }
   }
}


void
sum_reduce_4d_to_0d(
   RealArray& a_result,
   const RealArray& a_source,
   const tbox::Box& a_box)
{
   a_result = 0;
   tbox::IntVector i(a_box.lower());
   for (i[3] = a_box.lower(3); i[3] <= a_box.upper(3); ++i[3]) {
      for (i[2] = a_box.lower(2); i[2] <= a_box.upper(2); ++i[2]) {
         for (i[1] = a_box.lower(1); i[1] <= a_box.upper(1); ++i[1]) {
            for (i[0] = a_box.lower(0); i[0] <= a_box.upper(0); ++i[0]) {
               a_result += a_source(i[0], i[1], i[2], i[3]);
            }
         }
      }
   }
}


void
sum_reduce_4d_to_1d(
   RealArray& a_result,
   const RealArray& a_source,
   const tbox::Box& a_box,
   const tbox::IntVector& a_rdir)
{
   a_result = 0;
   tbox::IntVector i(a_box.lower());
   for (i[3] = a_box.lower(3); i[3] <= a_box.upper(3); ++i[3]) {
      for (i[2] = a_box.lower(2); i[2] <= a_box.upper(2); ++i[2]) {
         for (i[1] = a_box.lower(1); i[1] <= a_box.upper(1); ++i[1]) {
            for (i[0] = a_box.lower(0); i[0] <= a_box.upper(0); ++i[0]) {
               a_result(i[a_rdir[0]]) += a_source(i[0], i[1], i[2], i[3]);
            }
         }
      }
   }
}


void
sum_reduce_4d_to_2d(
   RealArray& a_result,
   const RealArray& a_source,
   const tbox::Box& a_box,
   const tbox::IntVector& a_rdir)
{
   a_result = 0;
   tbox::IntVector i(a_box.lower());
   for (i[3] = a_box.lower(3); i[3] <= a_box.upper(3); ++i[3]) {
      for (i[2] = a_box.lower(2); i[2] <= a_box.upper(2); ++i[2]) {
         for (i[1] = a_box.lower(1); i[1] <= a_box.upper(1); ++i[1]) {
            for (i[0] = a_box.lower(0); i[0] <= a_box.upper(0); ++i[0]) {
               a_result(i[a_rdir[0]], i[a_rdir[1]]) +=
                  a_source(i[0], i[1], i[2], i[3]);
            }
         }
      }
   }
}


void
sum_reduce_4d_to_3d(
   RealArray& a_result,
   const RealArray& a_source,
   const tbox::Box& a_box,
   const tbox::IntVector& a_rdir)
{
   a_result = 0;
   tbox::IntVector i(a_box.lower());
   for (i[3] = a_box.lower(3); i[3] <= a_box.upper(3); ++i[3]) {
      for (i[2] = a_box.lower(2); i[2] <= a_box.upper(2); ++i[2]) {
         for (i[1] = a_box.lower(1); i[1] <= a_box.upper(1); ++i[1]) {
            for (i[0] = a_box.lower(0); i[0] <= a_box.upper(0); ++i[0]) {
               a_result(i[a_rdir[0]], i[a_rdir[1]], i[a_rdir[2]]) +=
                  a_source(i[0], i[1], i[2], i[3]);
            }
         }
      }
   }
}


SumOp::SumOp()
{
}


SumOp::~SumOp()
{
}


void
SumOp::applyLocal(
   RealArray& a_result,
   const RealArray& a_source,
   const tbox::Box& a_box,
   const Array<bool>& a_collapse_dir) const
{
   // Figure out what we start with and what we collapse.  Then call the
   // specific method to sum over the collapsed directions.
   tbox::Dimension old_dim(a_box.getDim());
   checkInput(old_dim, a_source, a_collapse_dir);
   tbox::Dimension new_dim(computeNewDimension(old_dim, a_collapse_dir));
   tbox::IntVector rdir(new_dim);
   resultingDirections(rdir, a_collapse_dir);
   resizeResultArray(a_result, a_box, rdir);

   switch (old_dim) {
      case 1:
         switch (new_dim) {
            case 0:
               sum_reduce_1d_to_0d(a_result, a_source, a_box);
               break;
         }
         break;
      case 2:
         switch (new_dim) {
            case 0:
               sum_reduce_2d_to_0d(a_result, a_source, a_box);
               break;
            case 1:
               sum_reduce_2d_to_1d(a_result, a_source, a_box, rdir);
               break;
         }
         break;
      case 3:
         switch (new_dim) {
            case 0:
               sum_reduce_3d_to_0d(a_result, a_source, a_box);
               break;
            case 1:
               sum_reduce_3d_to_1d(a_result, a_source, a_box, rdir);
               break;
            case 2:
               sum_reduce_3d_to_2d(a_result, a_source, a_box, rdir);
               break;
         }
         break;
      case 4:
         switch (new_dim) {
            case 0:
               sum_reduce_4d_to_0d(a_result, a_source, a_box);
               break;
            case 1:
               sum_reduce_4d_to_1d(a_result, a_source, a_box, rdir);
               break;
            case 2:
               sum_reduce_4d_to_2d(a_result, a_source, a_box, rdir);
               break;
            case 3:
               sum_reduce_4d_to_3d(a_result, a_source, a_box, rdir);
               break;
         }
         break;
   }
}


void
SumOp::applyGlobal(
   RealArray& a_result,
   const RealArray& a_source,
   MPI_Comm& a_comm,
   ReductionTargets a_reduction_target) const
{
   // FIXME: this won't work if config space is distributed unless the
   // communicator is for processors that share exactly the same config
   // space domain.
   ParallelUtility::getSums(&(*a_source.getDataPointer()),
      &(*a_result.getDataPointer()),
      a_source.elementCount(),
      a_reduction_target,
      a_comm);
}

} // end namespace Loki
