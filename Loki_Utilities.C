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
#include "Loki_Utilities.H"

#include "BoxOps.H"

namespace Loki {

void
Loki_Utilities::fixPeriodicity(
   realArray& a_f,
   const ProblemDomain& a_domain,
   const tbox::IntVector& a_n_ghosts,
   int a_depth)
{
   const tbox::Box& domain_box(a_domain.box());
   tbox::Box outer_box(domain_box);
   outer_box.grow(a_n_ghosts);
   tbox::Box inner_box(domain_box);
   inner_box.grow(-a_n_ghosts);

   Index src[4];
   Index dest[4];
   if (a_domain.isPeriodic(X1)) {
      // Copy the data from the lower X1 boundary to the ghosts beyond the
      // upper X1 boundary.
      src[0]  = Range(domain_box.lower(X1), inner_box.lower(X1) - 1);
      src[1]  = Range(BoxOps::range(outer_box, X2));
      src[2]  = Range(0, a_depth-1);
      src[3]  = Range(0, 0);
      dest[0] = Range(domain_box.upper(X1) + 1, outer_box.upper(X1));
      dest[1] = Range(BoxOps::range(outer_box, X2));
      dest[2] = Range(0, a_depth-1);
      dest[3] = Range(0, 0);
      ParallelUtility::copy(a_f, dest, a_f, src, 4);

      // Copy the data from the upper X1 boundary to the ghosts beyond the
      // lower X1 boundary.
      src[0]  = Range(inner_box.upper(X1) + 1, domain_box.upper(X1));
      dest[0] = Range(outer_box.lower(X1), domain_box.lower(X1) - 1);
      ParallelUtility::copy(a_f, dest, a_f, src, 4);
   }
   if (a_domain.isPeriodic(X2)) {
      // Copy the data from the lower X2 boundary to the ghosts beyond the
      // upper X2 boundary.
      src[0] = Range(BoxOps::range(outer_box, X1));
      src[1] = Range(domain_box.lower(X2), inner_box.lower(X2) - 1);
      src[2] = Range(0, a_depth-1);
      src[3] = Range(0, 0);
      dest[0] = Range(BoxOps::range(outer_box, X1));
      dest[1] = Range(domain_box.upper(X2) + 1, outer_box.upper(X2));
      dest[2] = Range(0, a_depth-1);
      dest[3] = Range(0, 0);
      ParallelUtility::copy(a_f, dest, a_f, src, 4);

      // Copy the data from the upper X2 boundary to the ghosts beyond the
      // lower X2 boundary.
      src[1]  = Range(inner_box.upper(X2) + 1, domain_box.upper(X2));
      dest[1] = Range(outer_box.lower(X2), domain_box.lower(X2) - 1);
      ParallelUtility::copy(a_f, dest, a_f, src, 4);
   }
}

bool
Loki_Utilities::reduceBoolean(
   bool a_bool)
{
   // Booleans are some implementation dependent form of an int.  Assign the
   // boolean's value to a literal int and reduce it as there is no MPI type
   // for a bool.
   int local(a_bool ? 1 : 0);
   int global(1);
   MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
   return (global == 1);
}

}
