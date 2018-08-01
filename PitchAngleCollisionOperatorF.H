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
#ifndef _PITCH_ANGLE_COLLISION_OPERATOR_FORT_H_
#define _PITCH_ANGLE_COLLISION_OPERATOR_FORT_H_

#define FORT_APPEND_PITCH_ANGLE_COLLISION appendpitchanglecollision_

/**
 * Prototypes of fortran functions called by PitchAngleCollisionOperator.
 */
extern "C" {

   void
   FORT_APPEND_PITCH_ANGLE_COLLISION(
      real &rhs,
      const real &u,
      const int &nd1lo,
      const int &nd1hi,
      const int &nd2lo,
      const int &nd2hi,
      const int &nd3lo,
      const int &nd3hi,
      const int &nd4lo,
      const int &nd4hi,
      const int &n1lo,
      const int &n1hi,
      const int &n2lo,
      const int &n2hi,
      const int &n3lo,
      const int &n3hi,
      const int &n4lo,
      const int &n4hi,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const real& dparams,
      const int& iparams);

};

#endif
