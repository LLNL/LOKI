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
#ifndef _COLLISION_OPERATOR_FORT_H_
#define _COLLISION_OPERATOR_FORT_H_

#define FORT_COMPUTE_DIST_FUNC_AVG computedistfuncavg_
#define FORT_COMPUTE_VTH_LOCAL     computevthlocal_
#define FORT_COMPUTE_VTH_GLOBAL    computevthglobal_

/**
 * Prototypes of fortran functions called by CollisionOperator.
 */
extern "C" {
   void
   FORT_COMPUTE_DIST_FUNC_AVG(
      const int& nd1a,
      const int& nd1b,
      const int& nd2a,
      const int& nd2b,
      const int& nd3a,
      const int& nd3b,
      const int& nd4a,
      const int& nd4b,
      const int& n1a,
      const int& n1b,
      const int& n2a,
      const int& n2b,
      const int& n3a,
      const int& n3b,
      const int& n4a,
      const int& n4b,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const int& vel_space_size,
      const int& config_space_size,
      const real& u,
      real& uxybar);

   void
   FORT_COMPUTE_VTH_LOCAL(
      const int& nd1a,
      const int& nd1b,
      const int& nd2a,
      const int& nd2b,
      const int& nd3a,
      const int& nd3b,
      const int& nd4a,
      const int& nd4b,
      const int& n1a,
      const int& n1b,
      const int& n2a,
      const int& n2b,
      const int& n3a,
      const int& n3b,
      const int& n4a,
      const int& n4b,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const real& vflowx,
      const real& vflowy,
      const int& comm,
      const real& u,
      real& vth);

   void
   FORT_COMPUTE_VTH_GLOBAL(
      const int& nd1a,
      const int& nd1b,
      const int& nd2a,
      const int& nd2b,
      const int& nd3a,
      const int& nd3b,
      const int& nd4a,
      const int& nd4b,
      const int& n1a,
      const int& n1b,
      const int& n2a,
      const int& n2b,
      const int& n3a,
      const int& n3b,
      const int& n4a,
      const int& n4b,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const int& vel_space_size,
      const real& vflowx,
      const real& vflowy,
      const real& uxybar,
      const int& is_config_space_head_node,
      const int& comm,
      real& vth);
};

#endif
