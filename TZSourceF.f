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
c Fortran functions called by twilight zone sources.
c
      subroutine settrigtzsource(f,
     *                           nd1a, nd1b, nd2a, nd2b,
     *                           nd3a, nd3b, nd4a, nd4b,
     *                           xlo, xhi, dx,
     *                           time,
     *                           dparams)
c
c.. function to set solution
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      real f(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real xlo(1:4), xhi(1:4), dx(1:4)
      real time
      real dparams(*)
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      real x,y,vx,vy
      real kx,ky,kt,alpha,A,t
      real h
      real amp
      real pi

      amp = dparams(1)
c
      kx = 1.d0
      ky = 1.d0
      kt = 1.d0
      alpha = 1.d0
      A = amp
      t = time
      pi = 4.d0*atan(1.d0)
c
      do i4 = nd4a,nd4b
        vy = xlo(4)+(i4+0.5)*dx(4)
        do i3 = nd3a,nd3b
          vx = xlo(3)+(i3+0.5)*dx(3)
          do i2 = nd2a,nd2b
            y = xlo(2)+(i2+0.5)*dx(2)
            do i1 = nd1a,nd1b
              x = xlo(1)+(i1+0.5)*dx(1)

              h = 
     #-0.1D1 / (kx ** 2 + ky ** 2) * A * sin(kx * x) * kx * cos(ky 
     #* y) * sin(kt * t) * alpha ** 2 / pi * vx * exp(-alpha * (vx ** 2 
     #+ vy ** 2) / 0.2D1) * (0.1D1 + A * cos(kx * x) * cos(ky * y) * sin
     #(kt * t)) / 0.2D1 - 0.1D1 / (kx ** 2 + ky ** 2) * A * cos(kx * x) 
     #* sin(ky * y) * ky * sin(kt * t) * alpha ** 2 / pi * vy * exp(-alp
     #ha * (vx ** 2 + vy ** 2) / 0.2D1) * (0.1D1 + A * cos(kx * x) * cos
     #(ky * y) * sin(kt * t)) / 0.2D1 - alpha / pi * exp(-alpha * (vx **
     # 2 + vy ** 2) / 0.2D1) * A * sin(kx * x) * kx * cos(ky * y) * sin(
     #kt * t) * vx / 0.2D1 - alpha / pi * exp(-alpha * (vx ** 2 + vy ** 
     #2) / 0.2D1) * A * cos(kx * x) * sin(ky * y) * ky * sin(kt * t) * v
     #y / 0.2D1 + alpha / pi * exp(-alpha * (vx ** 2 + vy ** 2) / 0.2D1)
     # * A * cos(kx * x) * cos(ky * y) * cos(kt * t) * kt / 0.2D1

              f(i1, i2, i3, i4) = f(i1, i2, i3, i4) + h
            end do
          end do
        end do
      end do

      return
      end
c
c **************
c
      subroutine computetrigtzsourceerror(error, soln,
     *                           nd1a, nd1b, nd2a, nd2b,
     *                           nd3a, nd3b, nd4a, nd4b,
     *                           xlo, xhi, dx,
     *                           time,
     *                           dparams)
c
c.. function to set solution
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      real error(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real soln(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real xlo(1:4), xhi(1:4), dx(1:4)
      real time
      real dparams(*)
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      real x,y,vx,vy
      real kx,ky,kt,alpha,A,t
      real fexact
      real amp
      real pi

      amp = dparams(1)
c
      kx = 1.d0
      ky = 1.d0
      kt = 1.d0
      alpha = 1.d0
      A = amp
      t = time
      pi = 4.d0*atan(1.d0)
c
      do i4 = nd4a,nd4b
        vy = xlo(4)+(i4+0.5)*dx(4)
        do i3 = nd3a,nd3b
          vx = xlo(3)+(i3+0.5)*dx(3)
          do i2 = nd2a,nd2b
            y = xlo(2)+(i2+0.5)*dx(2)
            do i1 = nd1a,nd1b
              x = xlo(1)+(i1+0.5)*dx(1)

              fexact = 
     #alpha / pi * exp(-alpha * (vx ** 2 + vy ** 2) / 0.2D1) * (0.
     #1D1 + A * cos(kx * x) * cos(ky * y) * sin(kt * t)) / 0.2D1

              error(i1, i2, i3, i4) = soln(i1, i2, i3, i4) - fexact

            end do
          end do
        end do
      end do
c
      return
      end
