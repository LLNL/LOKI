c
c Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
c See the top-level LICENSE file for details.
c Produced at the Lawrence Livermore National Laboratory
c
c SPDX-License-Identifier: MIT
c
c Fortran functions called by twilight zone sources.
c
      subroutine settwoiontrigtzsource(f,
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
      real kxI,kyI,ktI,kxE,kyE,ktE 
      real alphaE,alphaI,A,t
      real h1,h2,h3
      real cg1
      real amp
      real pi
      real me, mi 

      amp = dparams(1)
      me  = dparams(2)
      mi  = dparams(3)
c
      kxI = 0.2D1
      kyI = 0.4D1
      ktI = 0.1D1
      kxE = 0.4D1
      kyE = 0.2D1
      ktE = 0.1D1

c      me = 0.1D1
c      mi = 0.1D2      
       

      alphaE = sqrt(me) ;  
      alphaI = sqrt(mi) ; 

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

c     Test 1 - Decoupled

c      h1 = alphaI ** 2 * exp(-alphaI ** 2 * (vx ** 2 + vy ** 2) / 0.2D1
c     #) * A * (-vx * kxI * sin(kxI * x) * cos(kyI * y) * sin(ktI * t) - vy *
c     # cos(kxI * x) * kyI * sin(kyI * y) * sin(ktI * t) + cos(kxI * x) * cos(
c     #kyI * y) * ktI * cos(ktI * t)) / pi / 0.2D1

c     Test 2 - Coupled, same frequency

c      h2 = A * alphaI ** 2 * exp(-alphaI ** 2 * (vx ** 2 + vy ** 2) / 0
c     #.2D1) * (A * vy * alphaI ** 2 * cos(kyI * y) * sin(kyI * y) * kyI * (
c     #cos(ktI * t) - 0.1D1) * (cos(ktI * t) + 0.1D1) * cos(kxI * x) ** 2 + 
c     #(A * alphaI ** 2 * sin(kxI * x) * kxI * vx * (cos(ktI * t) - 0.1D1) *
c     # (cos(ktI * t) + 0.1D1) * cos(kyI * y) ** 2 + mi * cos(ktI * t) * ktI 
c     #* (kxI ** 2 + kyI ** 2) * cos(kyI * y) - (mi * kxI ** 2 + mi * kyI ** 2
c     # + alphaI ** 2 / 0.2D1) * kyI * sin(ktI * t) * sin(kyI * y) * vy) * c
c     #os(kxI * x) - (mi * kxI ** 2 + mi * kyI ** 2 + alphaI ** 2 / 0.2D1) *
c     # vx * sin(ktI * t) * kxI * cos(kyI * y) * sin(kxI * x)) / mi / (kxI ** 
c     #2 + kyI ** 2) / pi

c     Test 3 - Coupled, different frequency 

c     h3 = 
c     #0.2D1 * (A * kyI * alphaI ** 2 * cos(kyI * y) * sin(kyI * y)
c     # * vy * (cos(ktI * t) - 0.1D1) * (cos(ktI * t) + 0.1D1) * (kxE ** 
c     #2 + kyE ** 2) * cos(kxI * x) ** 2 + (A * kxI * alphaI ** 2 * sin(k
c     #xI * x) * vx * (cos(ktI * t) - 0.1D1) * (cos(ktI * t) + 0.1D1) * (
c     #kxE ** 2 + kyE ** 2) * cos(kyI * y) ** 2 + (A * alphaI ** 2 * sin(
c     #ktE * t) * (kxE * cos(kyE * y) * sin(kxE * x) * vx + kyE * cos(kxE
c     # * x) * sin(kyE * y) * vy) * sin(ktI * t) + ktI * cos(ktI * t) * m
c     #i * (kxE ** 2 + kyE ** 2)) * (kxI ** 2 + kyI ** 2) * cos(kyI * y) 
c     #/ 0.2D1 - (alphaI ** 2 + (kxI ** 2 + kyI ** 2) * mi) * sin(kyI * y
c     #) * kyI * (kxE ** 2 + kyE ** 2) * sin(ktI * t) * vy / 0.2D1) * cos
c     #(kxI * x) - vx * (alphaI ** 2 + (kxI ** 2 + kyI ** 2) * mi) * (kxE
c     # ** 2 + kyE ** 2) * sin(ktI * t) * sin(kxI * x) * kxI * cos(kyI * 
c     #y) / 0.2D1 + alphaI ** 2 * sin(ktE * t) * (kxI ** 2 + kyI ** 2) * 
c     #(kxE * cos(kyE * y) * sin(kxE * x) * vx + kyE * cos(kxE * x) * sin
c     #(kyE * y) * vy) / 0.4D1) * A * exp(-alphaI ** 2 * (vx ** 2 + vy **
c     # 2) / 0.2D1) * alphaI ** 2 / mi / (kxI ** 2 + kyI ** 2) / (kxE ** 
c     #2 + kyE ** 2) / pi

      cg1 = alphaI ** 2 / pi * exp(-alphaI ** 2 * (vx ** 2 + vy ** 2) / 
     #0.2D1) * A * cos(kxI * x) * cos(kyI * y) * ktI * cos(ktI * t) - vx
     # * alphaI ** 2 / pi * exp(-alphaI ** 2 * (vx ** 2 + vy ** 2) / 0.2
     #D1) * A * kxI * sin(kxI * x) * cos(kyI * y) * sin(ktI * t) - vy * 
     #alphaI ** 2 / pi * exp(-alphaI ** 2 * (vx ** 2 + vy ** 2) / 0.2D1)
     # * A * cos(kxI * x) * kyI * sin(kyI * y) * sin(ktI * t) + 0.1D1 / 
     #mi * (sin(ktE * t) * sin(kxE * x) * kxE * (kxI ** 2 + kyI ** 2) * 
     #cos(kyE * y) - 0.2D1 * cos(kyI * y) * sin(ktI * t) * sin(kxI * x) 
     #* kxI * (kxE ** 2 + kyE ** 2)) * A / (kxI ** 2 + kyI ** 2) / (kxE 
     #** 2 + kyE ** 2) * alphaI ** 4 / pi * vx * exp(-alphaI ** 2 * (vx 
     #** 2 + vy ** 2) / 0.2D1) * (0.1D1 + 0.2D1 * A * cos(kxI * x) * cos
     #(kyI * y) * sin(ktI * t)) / 0.2D1 + 0.1D1 / mi * A * (sin(ktE * t)
     # * sin(kyE * y) * kyE * (kxI ** 2 + kyI ** 2) * cos(kxE * x) - 0.2
     #D1 * kyI * cos(kxI * x) * sin(ktI * t) * sin(kyI * y) * (kxE ** 2 
     #+ kyE ** 2)) / (kxI ** 2 + kyI ** 2) / (kxE ** 2 + kyE ** 2) * alp
     #haI ** 4 / pi * vy * exp(-alphaI ** 2 * (vx ** 2 + vy ** 2) / 0.2D
     #1) * (0.1D1 + 0.2D1 * A * cos(kxI * x) * cos(kyI * y) * sin(ktI * 
     #t)) / 0.2D1
 
              f(i1, i2, i3, i4) = f(i1, i2, i3, i4) + cg1
            end do
          end do
        end do
      end do

      return
      end
c
c **************
c
      subroutine computetwoiontrigtzsourceerror(error, soln,
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
      real kxI,kyI,ktI,kxE,kyE,ktE 
      real alphaE,alphaI,A,t
      real fexact1,fexact2,fexact3
      real cg2
      real amp
      real pi
      real me, mi

      amp = dparams(1)
      me  = dparams(2)
      mi  = dparams(3) 
c
      kxI = 0.2D1
      kyI = 0.4D1
      ktI = 0.1D1
      kxE = 0.4D1
      kyE = 0.2D1
      ktE = 0.1D1

 
c      me = 0.1D1
c      mi = 0.1D2      
       

      alphaE = sqrt(me) ;  
      alphaI = sqrt(mi) ; 

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

c     Test 1 - Decoupled

c      fexact1 = 
c     # alphaI ** 2 / pi * exp(-alphaI ** 2 * (vx ** 2 + vy ** 2) / 
c     #0.2D1) * (0.1D1 + A * cos(kxI * x) * cos(kyI * y) * sin(ktI * t)) / 0
c     #.2D1

c     Test 2 - Coupled, same frequency

c      fexact2 = 
c     # alphaI ** 2 / pi * exp(-alphaI ** 2 * (vx ** 2 + vy ** 2) / 
c     #0.2D1) * (0.1D1 + 0.2D1 * A * cos(kxI * x) * cos(kyI * y) * sin(ktI *
c     # t)) / 0.2D1

c     Test 3 - Coupled, different frequency 

c      fexact3 = 
c     #alphaI ** 2 / pi * exp(-alphaI ** 2 * (vx ** 2 + vy ** 2) / 0.2D1)
c     #* (0.1D1 + 0.2D1 * A * cos(kxI * x) * cos(kyI * y) * sin(ktI * t))
c     #/ 0.2D1

      cg2 = alphaI ** 2 / pi * exp(-alphaI ** 2 * (vx ** 2 + vy ** 2) / 
     #0.2D1) * (0.1D1 + 0.2D1 * A * cos(kxI * x) * cos(kyI * y) * sin(kt
     #I * t)) / 0.2D1

              error(i1, i2, i3, i4) = soln(i1, i2, i3, i4) - cg2

            end do
          end do
        end do
      end do
c
      return
      end
