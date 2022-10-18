c
c Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
c See the top-level LICENSE file for details.
c Produced at the Lawrence Livermore National Laboratory
c
c SPDX-License-Identifier: MIT
c
c Fortran functions called by twilight zone sources.
c
      subroutine setelectrontrigtzsource(f,
     *                           nd1a, nd1b, nd2a, nd2b,
     *                           nd3a, nd3b, nd4a, nd4b,
     *                           xlo, xhi, dx,
     *                           time,
     *                           velocities,
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
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real dparams(*)
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      real x,y,vx,vy
      real kx,ky,kt,alpha,A,t
      real h
      real amp
      real pi
      real E1 
      real E2

      amp = dparams(1)
c
      kx = 0.4D1
      ky = 0.4D1
      kt = 0.1D1
      alpha = 1.d0
      A = amp
      t = time
      pi = 4.d0*atan(1.d0)
c
      do i4 = nd4a,nd4b
        do i3 = nd3a,nd3b
          vx = velocities(i3, i4, 0)
          vy = velocities(i3, i4, 1)
          do i2 = nd2a,nd2b
            y = xlo(2)+(i2+0.5)*dx(2)
            do i1 = nd1a,nd1b
              x = xlo(1)+(i1+0.5)*dx(1)

      h = alpha / pi * exp(-alpha * (vx ** 2 + vy ** 2) / 0.2D1) * A *
     # cos(kx * x) * cos(ky * y) * kt * cos(kt * t) / 0.2D1 - vx * alpha
     # / pi * exp(-alpha * (vx ** 2 + vy ** 2) / 0.2D1) * A * kx * sin(k
     #x * x) * cos(ky * y) * sin(kt * t) / 0.2D1 - vy * alpha / pi * exp
     #(-alpha * (vx ** 2 + vy ** 2) / 0.2D1) * A * cos(kx * x) * ky * si
     #n(ky * y) * sin(kt * t) / 0.2D1 - A * kx * sin(kx * x) * cos(ky * 
     #y) * sin(kt * t) / (kx ** 2 + ky ** 2) * alpha ** 2 / pi * vx * ex
     #p(-alpha * (vx ** 2 + vy ** 2) / 0.2D1) * (0.1D1 + A * cos(kx * x)
     # * cos(ky * y) * sin(kt * t)) / 0.2D1 - A * cos(kx * x) * ky * sin
     #(ky * y) * sin(kt * t) / (kx ** 2 + ky ** 2) * alpha ** 2 / pi * v
     #y * exp(-alpha * (vx ** 2 + vy ** 2) / 0.2D1) * (0.1D1 + A * cos(k
     #x * x) * cos(ky * y) * sin(kt * t)) / 0.2D1





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
      subroutine computeelectrontrigtzsourceerror(error, soln,
     *                           nd1a, nd1b, nd2a, nd2b,
     *                           nd3a, nd3b, nd4a, nd4b,
     *                           xlo, xhi, dx,
     *                           time,
     *                           velocities,
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
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
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
      kx = 0.4D1
      ky = 0.4D1
      kt = 0.1D1
      alpha = 1.d0
      A = amp
      t = time
      pi = 4.d0*atan(1.d0)
c
      do i4 = nd4a,nd4b
        do i3 = nd3a,nd3b
          vx = velocities(i3, i4, 0)
          vy = velocities(i3, i4, 1)
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
