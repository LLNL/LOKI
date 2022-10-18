c
c Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
c See the top-level LICENSE file for details.
c Produced at the Lawrence Livermore National Laboratory
c
c SPDX-License-Identifier: MIT
c
c Fortran functions called by twilight zone sources.
c
      subroutine settwoelectrontrigtzsource(f,
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

      real cg
      real amp
      real pi
      real E1 
      real E2 
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

      A = amp
      t = time
      pi = 4.d0*atan(1.d0)

c      me = 0.1D1
c      mi = 0.1D2 
       

      alphaE = sqrt(me) ;  
      alphaI = sqrt(mi) ; 

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

c      h1 = alphaE ** 2 * exp(-alphaE ** 2 * (vx ** 2 + vy ** 2) / 0.2D1)
c     # * A * (-vx * kxI * sin(kxI * x) * cos(kyI * y) * sin(ktI * t) - vy * 
c     #cos(kxI * x) * kyI * sin(kyI * y) * sin(ktI * t) + cos(kxI * x) * cos(k
c     #y * y) * ktI * cos(ktI * t)) / pi / 0.2D1

c     Test 2 - Coupled, same frequency

c      h2 = -(A * vy * alphaE ** 2 * cos(kyI * y) * sin(kyI * y) * kyI * (co
c     #s(ktI * t) - 0.1D1) * (cos(ktI * t) + 0.1D1) * cos(kxI * x) ** 2 + (A
c     # * alphaE ** 2 * sin(kxI * x) * kxI * vx * (cos(ktI * t) - 0.1D1) * (
c     #cos(ktI * t) + 0.1D1) * cos(kyI * y) ** 2 - me * cos(ktI * t) * ktI * 
c     #(kxI ** 2 + kyI ** 2) * cos(kyI * y) + vy * sin(ktI * t) * sin(kyI * y)
c     # * kyI * (kxI ** 2 * me + kyI ** 2 * me - alphaE ** 2)) * cos(kxI * x)
c     # + cos(kyI * y) * sin(ktI * t) * sin(kxI * x) * kxI * vx * (kxI ** 2 * 
c     #me + kyI ** 2 * me - alphaE ** 2)) * A * alphaE ** 2 * exp(-alphaE 
c     #** 2 * (vx ** 2 + vy ** 2) / 0.2D1) / me / (kxI ** 2 + kyI ** 2) / p
c     #i / 0.2D1

c     Test 3 - Coupled, different frequency 

c      h3 = (A * sin(kyE * y) * cos(kyE * y) * alphaE ** 2 * kyE * vy * 
c     #( kxI ** 2 + kyI ** 2) * (cos(ktE * t) - 0.1D1) * (cos(ktE * t) +
c     #0.1D1) * cos(kxE * x) ** 2 / 0.2D1 + (A * sin(kxE * x) * alphaE 
c     #** 2 * kxE * vx * (kxI ** 2 + kyI ** 2) * (cos(ktE * t) - 0.1D1) 
c     #* (cos(ktE * t) + 0.1D1) * cos(kyE * y) ** 2 / 0.2D1 + (A * alphaE
c     # ** 2 * sin(ktI * t) * (vy * cos(kxI * x) * sin(kyI * y) * kyI + 
c     #cos(kyI * y) * sin(kxI * x) * kxI * vx) * sin(ktE * t) + me * ktE 
c     #* cos(ktE * t) * (kxI ** 2 + kyI ** 2) / 0.2D1) * (kxE ** 2 + 
c     #kyE ** 2) * cos(kyE * y ) - sin(ktE * t) * (alphaE ** 2 + me * 
c     #(kxE ** 2 + kyE ** 2)) * kyE * sin(kyE * y) * vy * (kxI ** 2 +
c     #kyI ** 2) / 0.2D1) * cos(kxE * x ) - sin(ktE * t) * vx * 
c     #(alphaE ** 2 + me * (kxE ** 2 + kyE ** 2)) * (kxI ** 2 + kyI ** 2)
c     #* kxE * sin(kxE * x) * cos(kyE * y) / 0.2D1 + alphaE ** 2 * 
c     #sin(ktI * t) * (kxE ** 2 + kyE ** 2) * (vy * cos(kxI * x) * 
c     #sin(kyI * y) * kyI + cos(kyI * y) * sin(kxI * x) * kxI * vx)) 
c     #* A * exp(-alphaE ** 2 * (vx ** 2 + vy ** 2) / 0.2D1) * alphaE ** 
c     #2 / me / (kxI ** 2 + kyI ** 2) / (kxE ** 2 + kyE ** 2) / pi


c      h3 = A * exp(-alphaE ** 2 * (vx ** 2 + vy ** 2) / 0.2D1) * (A * co
c     #s(kyE * y) * sin(kyE * y) * alphaE ** 2 * kyE * vy * (cos(ktE * t)
c     # - 0.1D1) * (cos(ktE * t) + 0.1D1) * (kxI ** 2 + kyI ** 2) * cos(k
c     #xE * x) ** 2 / 0.2D1 + (A * sin(kxE * x) * alphaE ** 2 * kxE * vx 
c     #* (cos(ktE * t) - 0.1D1) * (cos(ktE * t) + 0.1D1) * (kxI ** 2 + ky
c     #I ** 2) * cos(kyE * y) ** 2 / 0.2D1 + (kxE ** 2 + kyE ** 2) * (A *
c     # alphaE ** 2 * sin(ktI * t) * (kxI * cos(kyI * y) * sin(kxI * x) *
c     # vx + kyI * cos(kxI * x) * sin(kyI * y) * vy) * sin(ktE * t) + cos
c     #(ktE * t) * me * ktE * (kxI ** 2 + kyI ** 2) / 0.2D1) * cos(kyE * 
c     #y) - (alphaE ** 2 + me * (kxE ** 2 + kyE ** 2)) * kyE * sin(kyE * 
c     #y) * sin(ktE * t) * (kxI ** 2 + kyI ** 2) * vy / 0.2D1) * cos(kxE 
c     #* x) - sin(kxE * x) * vx * (alphaE ** 2 + me * (kxE ** 2 + kyE ** 
c     #2)) * kxE * sin(ktE * t) * (kxI ** 2 + kyI ** 2) * cos(kyE * y) / 
c     #0.2D1 + alphaE ** 2 * sin(ktI * t) * (kxE ** 2 + kyE ** 2) * (kxI 
c     #* cos(kyI * y) * sin(kxI * x) * vx + kyI * cos(kxI * x) * sin(kyI 
c     #* y) * vy)) * alphaE ** 2 / me / (kxI ** 2 + kyI ** 2) / (kxE ** 2
c     # + kyE ** 2) / pi
 
      cg = alphaE ** 2 / pi * exp(-alphaE ** 2 * (vx ** 2 + vy ** 2) / 0
     #.2D1) * A * cos(kxE * x) * cos(kyE * y) * ktE * cos(ktE * t) / 0.2
     #D1 - vx * alphaE ** 2 / pi * exp(-alphaE ** 2 * (vx ** 2 + vy ** 2
     #) / 0.2D1) * A * kxE * sin(kxE * x) * cos(kyE * y) * sin(ktE * t) 
     #/ 0.2D1 - vy * alphaE ** 2 / pi * exp(-alphaE ** 2 * (vx ** 2 + vy
     # ** 2) / 0.2D1) * A * cos(kxE * x) * kyE * sin(kyE * y) * sin(ktE 
     #* t) / 0.2D1 - 0.1D1 / me * (sin(ktE * t) * sin(kxE * x) * kxE * (
     #kxI ** 2 + kyI ** 2) * cos(kyE * y) - 0.2D1 * cos(kyI * y) * sin(k
     #tI * t) * sin(kxI * x) * kxI * (kxE ** 2 + kyE ** 2)) * A / (kxI *
     #* 2 + kyI ** 2) / (kxE ** 2 + kyE ** 2) * alphaE ** 4 / pi * vx * 
     #exp(-alphaE ** 2 * (vx ** 2 + vy ** 2) / 0.2D1) * (0.1D1 + A * cos
     #(kxE * x) * cos(kyE * y) * sin(ktE * t)) / 0.2D1 - 0.1D1 / me * A 
     #* (sin(ktE * t) * sin(kyE * y) * kyE * (kxI ** 2 + kyI ** 2) * cos
     #(kxE * x) - 0.2D1 * kyI * cos(kxI * x) * sin(ktI * t) * sin(kyI * 
     #y) * (kxE ** 2 + kyE ** 2)) / (kxI ** 2 + kyI ** 2) / (kxE ** 2 + 
     #kyE ** 2) * alphaE ** 4 / pi * vy * exp(-alphaE ** 2 * (vx ** 2 + 
     #vy ** 2) / 0.2D1) * (0.1D1 + A * cos(kxE * x) * cos(kyE * y) * sin
     #(ktE * t)) / 0.2D1
 

              f(i1, i2, i3, i4) = f(i1, i2, i3, i4) + cg
            end do
          end do
        end do
      end do

      return
      end
c
c **************
c
      subroutine computetwoelectrontrigtzsourceerror(error, soln,
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
      real fexact1, fexact2, fexact3 
      real cg0
      real amp
      real pi
      real me, mi

      amp = dparams(1)
      me  = dparams(2)
      mi  = dparams(3)
c
      kxE = 0.4D1
      kyE = 0.2D1
      ktE = 0.1D1
      kxI = 0.2D1
      kyI = 0.4D1
      ktI = 0.1D1

 
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
c     # alphaE ** 2 / pi * exp(-alphaE ** 2 * (vx ** 2 + vy ** 2) / 
c     #0.2D1) * (0.1D1 + A * cos(kxI * x) * cos(kyI * y) * sin(ktI * t)) / 0
c     #.2D1

c     Test 2 - Coupled, same frequency

c      fexact2 =
c     # alphaE ** 2 / pi * exp(-alphaE ** 2 * (vx ** 2 + vy ** 2) / 
c     #0.2D1) * (0.1D1 + A * cos(kxI * x) * cos(kyI * y) * sin(ktI * t)) / 0
c     #.2D1

c     Test 3 - Coupled, different frequency 

c      fexact3 = 
c     # alphaE ** 2 / pi * exp(-alphaE ** 2 * (vx ** 2 + vy ** 2) / 
c     #0.2D1) * (0.1D1 + A * cos(kxE * x) * cos(kyE * y) * sin(ktE * t)) 
c     #/ 0.2D1

c      fexact3 = 
c     # alphaE ** 2 / pi * exp(-alphaE ** 2 * (vx ** 2 + vy ** 2) / 
c     #0.2D1) * (0.1D1 + A * cos(kxE * x) * cos(kyE * y) * sin(ktE * t)) 
c     #/ 0.2D1

      cg0 = alphaE ** 2 / pi * exp(-alphaE ** 2 * (vx ** 2 + vy ** 2) / 
     #0.2D1) * (0.1D1 + A * cos(kxE * x) * cos(kyE * y) * sin(ktE * t)) 
     #/ 0.2D1

              error(i1, i2, i3, i4) = soln(i1, i2, i3, i4) - cg0

            end do
          end do
        end do
      end do
c
      return
      end
