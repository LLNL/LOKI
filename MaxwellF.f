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
c Fortran functions called by Maxwell.
c
      subroutine zeroghost2d(
     & u,
     & n1a,n1b,n2a,n2b,
     & nd1a,nd1b,nd2a,nd2b,
     & dim )
c
      implicit none
      integer nd1a,nd1b,nd2a,nd2b
      integer n1a,n1b,n2a,n2b
      integer dim
      real u(nd1a:nd1b,nd2a:nd2b,1:dim)
      integer i1,i2,i3
c
c .. i1 left
      do i3=1,dim
      do i2=nd2a,nd2b
      do i1=nd1a,n1a-1
        u(i1,i2,i3) = 0.0
      end do
      end do
      end do
c .. i1 right
      do i3=1,dim
      do i2=nd2a,nd2b
      do i1=n1b+1,nd1b
        u(i1,i2,i3) = 0.0
      end do
      end do
      end do
c
c .. i2 left
      do i3=1,dim
      do i2=nd2a,n2a-1
      do i1=nd1a,nd1b
        u(i1,i2,i3) = 0.0
      end do
      end do
      end do
c .. i2 right
      do i3=1,dim
      do i2=n2b+1,nd2b
      do i1=nd1a,nd1b
        u(i1,i2,i3) = 0.0
      end do
      end do
      end do

      return
      end
c
c ++++++++++++++
c
      subroutine xpby2d(
     &     x,
     &     nd1a, nd1b, nd2a, nd2b,
     &     y,
     &     n1a, n1b, n2a, n2b,
     &     b,
     &     dim)
c
c.. compute x = x + by
      implicit none
c
c.. declaration of incoming variables
      integer nd1a, nd1b, nd2a, nd2b
      integer n1a, n1b, n2a, n2b
      integer dim
      real x(nd1a:nd1b, nd2a:nd2b, 1:dim)
      real y(nd1a:nd1b, nd2a:nd2b, 1:dim)
      real b
c
c.. declaration of local variables
      integer i1, i2, i3
c
      do i3 = 1, dim
        do i2 = n2a, n2b
          do i1 = n1a, n1b
            x(i1, i2, i3) = x(i1, i2, i3) + b * y(i1, i2, i3)
          end do
        end do
      end do
c
      return
      end
c
c ++++++++++++++
c
      subroutine maxwellevalrhs(
     &     md1a, md1b, md2a, md2b,
     &     m1a, m1b, m2a,  m2b,
     &     xlo, xhi, dx,
     &     c, avWeak, avStrong, solution_order,
     &     EMvars,
     &     Jx, Jy, Jz,
     &     dEMvars)
c
c.. compute rhs of Maxwell's equations
      implicit none
c
c.. declaration of incoming variables
      integer md1a, md1b, md2a, md2b
      integer m1a, m1b, m2a, m2b
      real xlo(1:4), xhi(1:4), dx(1:4)
      real c, avWeak, avStrong
      integer solution_order
      real EMvars(md1a:md1b, md2a:md2b, 1:6)
      real Jx(md1a:md1b, md2a:md2b)
      real Jy(md1a:md1b, md2a:md2b)
      real Jz(md1a:md1b, md2a:md2b)
      real dEMvars(md1a:md1b, md2a:md2b, 1:6)
c
c.. declaration of local variables
      integer i1, i2, comp
      real csquared
      real Exdy, Eydx, Ezdx, Ezdy
      real Bxdy, Bydx, Bzdx, Bzdy
      real uxxxx, uyyyy, uxxxxxx, uyyyyyy
c
      csquared = c**2
c
      if( solution_order .eq. 4 ) then
        ! 4th order
        do i2 = m2a, m2b
        do i1 = m1a, m1b
          Exdy = (
     *         EMvars(i1, i2-2, 1)
     *         -8.0*EMvars(i1, i2-1, 1)
     *         +8.0*EMvars(i1, i2+1, 1)
     *         -EMvars(i1, i2+2, 1))/(12.0*dx(2))
          Eydx = (
     *         EMvars(i1-2, i2, 2)
     *         -8.0*EMvars(i1-1, i2, 2)
     *         +8.0*EMvars(i1+1, i2, 2)
     *         -EMvars(i1+2, i2, 2))/(12.0*dx(1))
          Ezdx = (
     *         EMvars(i1-2, i2, 3)
     *         -8.0*EMvars(i1-1, i2, 3)
     *         +8.0*EMvars(i1+1, i2, 3)
     *         -EMvars(i1+2, i2, 3))/(12.0*dx(1))
          Ezdy = (
     *         EMvars(i1, i2-2, 3)
     *         -8.0*EMvars(i1, i2-1, 3)
     *         +8.0*EMvars(i1, i2+1, 3)
     *         -EMvars(i1, i2+2, 3))/(12.0*dx(2))

          Bxdy = (
     *         EMvars(i1, i2-2, 4)
     *         -8.0*EMvars(i1, i2-1, 4)
     *         +8.0*EMvars(i1, i2+1, 4)
     *         -EMvars(i1, i2+2, 4))/(12.0*dx(2))
          Bydx = (
     *         EMvars(i1-2, i2, 5)
     *         -8.0*EMvars(i1-1, i2, 5)
     *         +8.0*EMvars(i1+1, i2, 5)
     *         -EMvars(i1+2, i2, 5))/(12.0*dx(1))
          Bzdx = (
     *         EMvars(i1-2, i2, 6)
     *         -8.0*EMvars(i1-1, i2, 6)
     *         +8.0*EMvars(i1+1, i2, 6)
     *         -EMvars(i1+2, i2, 6))/(12.0*dx(1))
          Bzdy = (
     *         EMvars(i1, i2-2, 6)
     *         -8.0*EMvars(i1, i2-1, 6)
     *         +8.0*EMvars(i1, i2+1, 6)
     *         -EMvars(i1, i2+2, 6))/(12.0*dx(2))
          
          dEMvars(i1, i2, 1) =  csquared*(Bzdy)-Jx(i1, i2)
          dEMvars(i1, i2, 2) = -csquared*(Bzdx)-Jy(i1, i2)
          dEMvars(i1, i2, 3) =  csquared*(Bydx-Bxdy)-Jz(i1, i2)
          
          dEMvars(i1, i2, 4) = -Ezdy
          dEMvars(i1, i2, 5) =  Ezdx
          dEMvars(i1, i2, 6) =  Exdy-Eydx
          
        end do
        end do
      else
        ! 6th order
        do i2 = m2a, m2b
        do i1 = m1a, m1b
          Exdy = (
     *         -1.0 *EMvars(i1, i2-3, 1)
     *         +9.0 *EMvars(i1, i2-2, 1)
     *         -45.0*EMvars(i1, i2-1, 1)
     *         +45.0*EMvars(i1, i2+1, 1)
     *         -9.0 *EMvars(i1, i2+2, 1)
     *         +1.0 *EMvars(i1, i2+3, 1))/(60.0*dx(2))
          Eydx = (
     *         -1.0 *EMvars(i1-3, i2, 2)
     *         +9.0 *EMvars(i1-2, i2, 2)
     *         -45.0*EMvars(i1-1, i2, 2)
     *         +45.0*EMvars(i1+1, i2, 2)
     *         -9.0 *EMvars(i1+2, i2, 2)
     *         +1.0 *EMvars(i1+3, i2, 2))/(60.0*dx(1))
          Ezdx = (
     *         -1.0 *EMvars(i1-3, i2, 3)
     *         +9.0 *EMvars(i1-2, i2, 3)
     *         -45.0*EMvars(i1-1, i2, 3)
     *         +45.0*EMvars(i1+1, i2, 3)
     *         -9.0 *EMvars(i1+2, i2, 3)
     *         +1.0 *EMvars(i1+3, i2, 3))/(60.0*dx(1))

          Ezdy = (
     *         -1.0 *EMvars(i1, i2-3, 3)
     *         +9.0 *EMvars(i1, i2-2, 3)
     *         -45.0*EMvars(i1, i2-1, 3)
     *         +45.0*EMvars(i1, i2+1, 3)
     *         -9.0 *EMvars(i1, i2+2, 3)
     *         +1.0 *EMvars(i1, i2+3, 3))/(60.0*dx(2))

          Bxdy = (
     *         -1.0 *EMvars(i1, i2-3, 4)
     *         +9.0 *EMvars(i1, i2-2, 4)
     *         -45.0*EMvars(i1, i2-1, 4)
     *         +45.0*EMvars(i1, i2+1, 4)
     *         -9.0 *EMvars(i1, i2+2, 4)
     *         +1.0 *EMvars(i1, i2+3, 4))/(60.0*dx(2))
          Bydx = (
     *         -1.0 *EMvars(i1-3, i2, 5)
     *         +9.0 *EMvars(i1-2, i2, 5)
     *         -45.0*EMvars(i1-1, i2, 5)
     *         +45.0*EMvars(i1+1, i2, 5)
     *         -9.0 *EMvars(i1+2, i2, 5)
     *         +1.0 *EMvars(i1+3, i2, 5))/(60.0*dx(1))
          Bzdx = (
     *         -1.0 *EMvars(i1-3, i2, 6)
     *         +9.0 *EMvars(i1-2, i2, 6)
     *         -45.0*EMvars(i1-1, i2, 6)
     *         +45.0*EMvars(i1+1, i2, 6)
     *         -9.0 *EMvars(i1+2, i2, 6)
     *         +1.0 *EMvars(i1+3, i2, 6))/(60.0*dx(1))
          Bzdy = (
     *         -1.0 *EMvars(i1, i2-3, 6)
     *         +9.0 *EMvars(i1, i2-2, 6)
     *         -45.0*EMvars(i1, i2-1, 6)
     *         +45.0*EMvars(i1, i2+1, 6)
     *         -9.0 *EMvars(i1, i2+2, 6)
     *         +1.0 *EMvars(i1, i2+3, 6))/(60.0*dx(2))
          
          dEMvars(i1, i2, 1) =  csquared*(Bzdy)-Jx(i1, i2)
          dEMvars(i1, i2, 2) = -csquared*(Bzdx)-Jy(i1, i2)
          dEMvars(i1, i2, 3) =  csquared*(Bydx-Bxdy)-Jz(i1, i2)
          
          dEMvars(i1, i2, 4) = -Ezdy
          dEMvars(i1, i2, 5) =  Ezdx
          dEMvars(i1, i2, 6) =  Exdy-Eydx
          
        end do
        end do
      end if
c
      if ((avWeak .gt. 0.0) .or. (avStrong .gt. 0.0)) then
        if( solution_order .eq. 4 ) then
          ! 4th order artificial dissipation
          do comp = 1, 6
          do i2 = m2a, m2b
          do i1 = m1a, m1b
            uxxxx =
     *           (1.0*EMvars(i1-2, i2, comp)
     *           -4.0*EMvars(i1-1, i2, comp)
     *           +6.0*EMvars(i1,   i2, comp)
     *           -4.0*EMvars(i1+1, i2, comp)
     *           +1.0*EMvars(i1+2, i2, comp))/(dx(1)**4)
            uyyyy =
     *           (1.0*EMvars(i1, i2-2, comp)
     *           -4.0*EMvars(i1, i2-1, comp)
     *           +6.0*EMvars(i1, i2,   comp)
     *           -4.0*EMvars(i1, i2+1, comp)
     *           +1.0*EMvars(i1, i2+2, comp))/(dx(2)**4)

            dEMvars(i1, i2, comp) = dEMvars(i1, i2, comp)
     *           -(avWeak*dx(1)**4+avStrong*dx(1)**3)*uxxxx
     *           -(avWeak*dx(2)**4+avStrong*dx(2)**3)*uyyyy
          end do
          end do
          end do
        else
          ! 6th order artificial dissipation
          do comp = 1, 6
          do i2 = m2a, m2b
          do i1 = m1a, m1b
            uxxxxxx =
     *           (1.0 *EMvars(i1-3, i2, comp)
     *           -6.0 *EMvars(i1-2, i2, comp)
     *           +15.0*EMvars(i1-1, i2, comp)
     *           -20.0*EMvars(i1,   i2, comp)
     *           +15.0*EMvars(i1+1, i2, comp)
     *           -6.0 *EMvars(i1+2, i2, comp)
     *           +1.0 *EMvars(i1+3, i2, comp))/(dx(1)**6)
            uyyyyyy =
     *           (1.0 *EMvars(i1, i2-3, comp)
     *           -6.0 *EMvars(i1, i2-2, comp)
     *           +15.0*EMvars(i1, i2-1, comp)
     *           -20.0*EMvars(i1, i2,   comp)
     *           +15.0*EMvars(i1, i2+1, comp)
     *           -6.0 *EMvars(i1, i2+2, comp)
     *           +1.0 *EMvars(i1, i2+3, comp))/(dx(2)**6)
            
            dEMvars(i1, i2, comp) = dEMvars(i1, i2, comp)
     *           +(avWeak*dx(1)**6+avStrong*dx(1)**5)*uxxxxxx
     *           +(avWeak*dx(2)**6+avStrong*dx(2)**5)*uyyyyyy
          end do
          end do
          end do
        end if
      end if
c
      return
      end
c
c ++++++++++++++
c
      subroutine maxwellevalvzrhs(
     &     md1a, md1b, md2a, md2b,
     &     m1a, m1b, m2a,  m2b,
     &     charge_per_mass,
     &     EMvars,
     &     dvz)
c
c.. compute rhs of Maxwell's equations
      implicit none
c
c.. declaration of incoming variables
      integer md1a, md1b, md2a, md2b
      integer m1a, m1b, m2a, m2b
      real charge_per_mass
      real EMvars(md1a:md1b, md2a:md2b, 1:6)
      real dvz(md1a:md1b, md2a:md2b)
c
c.. declaration of local variables
      integer i1, i2
c
      do i2 = m2a, m2b
        do i1 = m1a, m1b
          dvz(i1, i2) = charge_per_mass*EMvars(i1, i2, 3)
        end do
      end do
c
      return
      end
c
c+++++++++++
c
      subroutine computeEFieldFromPotentialMaxwell(
     *     nd1a, nd1b, nd2a, nd2b,
     *     n1a, n1b, n2a, n2b,
     *     solution_order,
     *     dx,
     *     E, phi)
c
c.. function to neutralize charge densities
      implicit none
c
      integer nd1a, nd1b, nd2a, nd2b
      integer n1a, n1b, n2a, n2b
      integer solution_order
      real dx(*)
      real E(nd1a:nd1b, nd2a:nd2b, 1:6)
      real phi(nd1a:nd1b, nd2a:nd2b)
c
      integer i1, i2
c
      if( solution_order .eq.4 ) then
        ! 4th order
        do i2 = n2a, n2b
        do i1 = n1a, n1b
          E(i1, i2, 1) = (phi(i1-2, i2)-8.0*phi(i1-1, i2)+
     *         8.0*phi(i1+1, i2)-phi(i1+2, i2))/(12.0*dx(1))
          E(i1, i2, 2) = (phi(i1, i2-2)-8.0*phi(i1, i2-1)+
     *         8.0*phi(i1, i2+1)-phi(i1, i2+2))/(12.0*dx(2))
        end do
        end do
      else
        ! 6th order
        do i2 = n2a,n2b
        do i1 = n1a,n1b
          E(i1,i2,1) = (
     *         -1.0  *phi(i1-3,i2)
     *         +9.0  *phi(i1-2,i2)
     *         -45.0 *phi(i1-1,i2)
     *         +45.0 *phi(i1+1,i2)
     *         -9.0  *phi(i1+2,i2)
     *         +1.0  *phi(i1+3,i2))/(60.0*dx(1))
          E(i1,i2,2) = (
     *         -1.0  *phi(i1,i2-3)
     *         +9.0  *phi(i1,i2-2)
     *         -45.0 *phi(i1,i2-1)
     *         +45.0 *phi(i1,i2+1)
     *         -9.0  *phi(i1,i2+2)
     *         +1.0  *phi(i1,i2+3))/(60.0*dx(2))
        end do
        end do
      end if
c
      return
      end
