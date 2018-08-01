c
c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c
c Written by Jeffrey Banks banksj3@rpi.edu (Rensselaer Polytechnic Institute,
c Amos Eaton 301, 110 8th St., Troy, NY 12180); Jeffrey Hittinger
c hittinger1@llnl.gov, William Arrighi arrighi2@llnl.gov, Richard Berger
c berger5@llnl.gov, Thomas Chapman chapman29@llnl.gov (LLNL, P.O Box 808,
c Livermore, CA 94551); Stephan Brunner stephan.brunner@epfl.ch (Ecole
c Polytechnique Federale de Lausanne, EPFL SB SPC-TH, PPB 312, Station 13,
c CH-1015 Lausanne, Switzerland).
c CODE-744849
c
c All rights reserved.
c
c This file is part of Loki.  For details, see.
c
c Permission is hereby granted, free of charge, to any person obtaining a
c copy of this software and associated documentation files (the "Software"),
c to deal in the Software without restriction, including without limitation
c the rights to use, copy, modify, merge, publish, distribute, sublicense,
c and/or sell copies of the Software, and to permit persons to whom the
c Software is furnished to do so, subject to the following conditions:
c
c The above copyright notice and this permission notice shall be included in
c all copies or substantial portions of the Software.
c
c THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
c OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
c FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
c THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
c LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
c FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
c DEALINGS IN THE SOFTWARE.
c
c Fortran functions called by SimpleEMIC.
c
      subroutine setsimpleemic(nd1a, nd1b, nd2a, nd2b,
     *                         u,
     *                         xlo, xhi, dx,
     *                         dparams,
     *                         iparams)
c
c.. function to set solution
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b
      real u(nd1a:nd1b, nd2a:nd2b, 1:6)
      real xlo(1:2), xhi(1:2), dx(1:2)
      real dparams(*)
      integer iparams(*)
c
c.. declarations of local variables
      integer i1, i2
      real x1, x2, env
      real xamp, yamp, zamp, kx, ky, phi
      integer field

      xamp     = dparams(1)
      yamp     = dparams(2)
      zamp     = dparams(3)
      kx       = dparams(4)
      ky       = dparams(5)
      phi      = dparams(6)
      field    = iparams(1)
      do i2 = nd2a,nd2b
        x2 = xlo(2)+(i2+0.5)*dx(2)
        do i1 = nd1a,nd1b
          x1 = xlo(1)+(i1+0.5)*dx(1)
          env = cos(kx*x1+ky*x2+phi)
          u(i1, i2, field)   = xamp*env
          u(i1, i2, field+1) = yamp*env
          u(i1, i2, field+2) = zamp*env
        end do
      end do
c
      return
      end
