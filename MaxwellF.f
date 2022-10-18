c
c Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
c See the top-level LICENSE file for details.
c Produced at the Lawrence Livermore National Laboratory
c
c SPDX-License-Identifier: MIT
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
     &     y,
     &     b,
     &     nd1a, nd1b, nd2a, nd2b,
     &     n1a, n1b, n2a, n2b,
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
     &     m1a, m1b, m2a, m2b,
     &     xlo, xhi, dx,
     &     c, avWeak, avStrong, solution_order,
     &     supergrid_lo, supergrid_hi,
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
      integer solution_order, xPeriodic, yPeriodic
      real supergrid_lo(1:2), supergrid_hi(1:2)
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
      real xmin, xmax, ymin, ymax
      real nu, nux, nuy
      real SGxa,SGxb
      real SGya,SGyb
      real x, y
c
      csquared = c**2
c
      xmin = xlo(1)
      ymin = xlo(2)
      xmax = xhi(1)
      ymax = xhi(2)
c
      SGxa = supergrid_lo(1)
      SGya = supergrid_lo(2)
      SGxb = supergrid_hi(1)
      SGyb = supergrid_hi(2)
c
      if( solution_order .eq. 4 ) then
        ! 4th order
        do i2 = m2a, m2b
        do i1 = m1a, m1b
          x = xmin+(0.5+i1)*dx(1)
          y = ymin+(0.5+i2)*dx(2)
          call SGMetricFunction( nux,x,xmin,xmax,SGxa,SGxb,
     *                           solution_order )
          call SGMetricFunction( nuy,y,ymin,ymax,SGya,SGyb,
     *                           solution_order )
          Exdy = nuy*(
     *              EMvars(i1, i2-2, 1)
     *         -8.0*EMvars(i1, i2-1, 1)
     *         +8.0*EMvars(i1, i2+1, 1)
     *             -EMvars(i1, i2+2, 1))/(12.0*dx(2))
          Eydx = nux*(
     *              EMvars(i1-2, i2, 2)
     *         -8.0*EMvars(i1-1, i2, 2)
     *         +8.0*EMvars(i1+1, i2, 2)
     *             -EMvars(i1+2, i2, 2))/(12.0*dx(1))
          Ezdx = nux*(
     *              EMvars(i1-2, i2, 3)
     *         -8.0*EMvars(i1-1, i2, 3)
     *         +8.0*EMvars(i1+1, i2, 3)
     *             -EMvars(i1+2, i2, 3))/(12.0*dx(1))
          Ezdy = nuy*(
     *              EMvars(i1, i2-2, 3)
     *         -8.0*EMvars(i1, i2-1, 3)
     *         +8.0*EMvars(i1, i2+1, 3)
     *             -EMvars(i1, i2+2, 3))/(12.0*dx(2))

          Bxdy = nuy*(
     *              EMvars(i1, i2-2, 4)
     *         -8.0*EMvars(i1, i2-1, 4)
     *         +8.0*EMvars(i1, i2+1, 4)
     *             -EMvars(i1, i2+2, 4))/(12.0*dx(2))
          Bydx = nux*(
     *              EMvars(i1-2, i2, 5)
     *         -8.0*EMvars(i1-1, i2, 5)
     *         +8.0*EMvars(i1+1, i2, 5)
     *             -EMvars(i1+2, i2, 5))/(12.0*dx(1))
          Bzdx = nux*(
     *              EMvars(i1-2, i2, 6)
     *         -8.0*EMvars(i1-1, i2, 6)
     *         +8.0*EMvars(i1+1, i2, 6)
     *             -EMvars(i1+2, i2, 6))/(12.0*dx(1))
          Bzdy = nuy*(
     *              EMvars(i1, i2-2, 6)
     *         -8.0*EMvars(i1, i2-1, 6)
     *         +8.0*EMvars(i1, i2+1, 6)
     *             -EMvars(i1, i2+2, 6))/(12.0*dx(2))
          
          dEMvars(i1, i2, 1) =  csquared*(Bzdy)-Jx(i1, i2)
          dEMvars(i1, i2, 2) = -csquared*(Bzdx)-Jy(i1, i2)
          dEMvars(i1, i2, 3) =  csquared*(Bydx-Bxdy)-Jz(i1, i2)

c          dEMvars(i1, i2, 1) = 0.0
c          dEMvars(i1, i2, 2) = 0.0
c          dEMvars(i1, i2, 3) = 0.0
          
          dEMvars(i1, i2, 4) = -Ezdy
          dEMvars(i1, i2, 5) =  Ezdx
          dEMvars(i1, i2, 6) =  Exdy-Eydx
          
        end do
        end do
      else
        ! 6th order
        do i2 = m2a, m2b
        do i1 = m1a, m1b
          x = xmin+(0.5+i1)*dx(1)
          y = ymin+(0.5+i2)*dx(2)
          call SGMetricFunction( nux,x,xmin,xmax,SGxa,SGxb,
     *                           solution_order )
          call SGMetricFunction( nuy,y,ymin,ymax,SGya,SGyb,
     *                           solution_order )
c          nux = 1.0
c          nuy = 1.0
          Exdy = nuy*(
     *         -1.0 *EMvars(i1, i2-3, 1)
     *         +9.0 *EMvars(i1, i2-2, 1)
     *         -45.0*EMvars(i1, i2-1, 1)
     *         +45.0*EMvars(i1, i2+1, 1)
     *         -9.0 *EMvars(i1, i2+2, 1)
     *         +1.0 *EMvars(i1, i2+3, 1))/(60.0*dx(2))
          Eydx = nux*(
     *         -1.0 *EMvars(i1-3, i2, 2)
     *         +9.0 *EMvars(i1-2, i2, 2)
     *         -45.0*EMvars(i1-1, i2, 2)
     *         +45.0*EMvars(i1+1, i2, 2)
     *         -9.0 *EMvars(i1+2, i2, 2)
     *         +1.0 *EMvars(i1+3, i2, 2))/(60.0*dx(1))
          Ezdx = nux*(
     *         -1.0 *EMvars(i1-3, i2, 3)
     *         +9.0 *EMvars(i1-2, i2, 3)
     *         -45.0*EMvars(i1-1, i2, 3)
     *         +45.0*EMvars(i1+1, i2, 3)
     *         -9.0 *EMvars(i1+2, i2, 3)
     *         +1.0 *EMvars(i1+3, i2, 3))/(60.0*dx(1))

          Ezdy = nuy*(
     *         -1.0 *EMvars(i1, i2-3, 3)
     *         +9.0 *EMvars(i1, i2-2, 3)
     *         -45.0*EMvars(i1, i2-1, 3)
     *         +45.0*EMvars(i1, i2+1, 3)
     *         -9.0 *EMvars(i1, i2+2, 3)
     *         +1.0 *EMvars(i1, i2+3, 3))/(60.0*dx(2))

          Bxdy = nuy*(
     *         -1.0 *EMvars(i1, i2-3, 4)
     *         +9.0 *EMvars(i1, i2-2, 4)
     *         -45.0*EMvars(i1, i2-1, 4)
     *         +45.0*EMvars(i1, i2+1, 4)
     *         -9.0 *EMvars(i1, i2+2, 4)
     *         +1.0 *EMvars(i1, i2+3, 4))/(60.0*dx(2))
          Bydx = nux*(
     *         -1.0 *EMvars(i1-3, i2, 5)
     *         +9.0 *EMvars(i1-2, i2, 5)
     *         -45.0*EMvars(i1-1, i2, 5)
     *         +45.0*EMvars(i1+1, i2, 5)
     *         -9.0 *EMvars(i1+2, i2, 5)
     *         +1.0 *EMvars(i1+3, i2, 5))/(60.0*dx(1))
          Bzdx = nux*(
     *         -1.0 *EMvars(i1-3, i2, 6)
     *         +9.0 *EMvars(i1-2, i2, 6)
     *         -45.0*EMvars(i1-1, i2, 6)
     *         +45.0*EMvars(i1+1, i2, 6)
     *         -9.0 *EMvars(i1+2, i2, 6)
     *         +1.0 *EMvars(i1+3, i2, 6))/(60.0*dx(1))
          Bzdy = nuy*(
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
     *           -(avWeak*c*dx(1)**4+avStrong*c*dx(1)**3)/16.0*uxxxx
     *           -(avWeak*c*dx(2)**4+avStrong*c*dx(2)**3)/16.0*uyyyy
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
     *           +(avWeak*c*dx(1)**6+avStrong*c*dx(1)**5)/64.0*uxxxxxx
     *           +(avWeak*c*dx(2)**6+avStrong*c*dx(2)**5)/64.0*uyyyyyy
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
      subroutine maxwelladdantennasource(
     &     md1a, md1b, md2a, md2b,
     &     m1a, m1b, m2a,  m2b,
     &     xlo, xhi, dx,
     &     antenna_source,
     &     dEMvars)
c
c.. compute rhs of Maxwell's equations
      implicit none
c
c.. declaration of incoming variables
      integer md1a, md1b, md2a, md2b
      integer m1a, m1b, m2a, m2b
      real xlo(1:4), xhi(1:4), dx(1:4)
      real antenna_source(md1a:md1b, md2a:md2b, 1:6)
      real dEMvars(md1a:md1b, md2a:md2b, 1:6)
c
c.. declaration of local variables
      integer i1, i2, comp
c
      do comp = 1, 6
        do i2 = m2a, m2b
          do i1 = m1a, m1b
            dEMvars(i1, i2, comp) =
     *        dEMvars(i1, i2, comp) - antenna_source(i1, i2, comp)
          end do
        end do
      end do
c
      return
      end
c
c ++++++++++++++
c
      subroutine SGMetricFunction( nu,x,xmin,xmax,SGxa,SGxb,
     &                             solution_order )
c
c .. coopute window function for supergrid layer
      implicit none
c
c.. declaration of incoming variables
      real nu
      real x
      real xmin,xmax
      real SGxa,SGxb
      integer solution_order
c
c.. declaration of local variables
      real xa,xb
      real xi
c
      if( x .lt. SGxa .and. x.gt.xmin ) then
        xa = SGxa
        xb = xmin
        xi = (x-xa)/(xb-xa)
      else if( x .gt. SGxb .and. x.lt.xmax ) then
        xa = SGxb
        xb = xmax
        xi = (x-xa)/(xb-xa)
      else
        xi = 0.0
      end if
      if (solution_order .eq. 4) then
        nu = 1.0+xi**4*(
     *       +20.0*xi**3
     *       -70.0*xi**2
     *       +84.0*xi
     *       -35.0)
      else
        nu = 1.0 + xi**6*(
     *       +252.0 *xi**5
     *       -1386.0*xi**4
     *       +3080.0*xi**3
     *       -3465.0*xi**2
     *       +1980.0*xi
     *       -462.0)
      end if
c
      return
      end
c
c ++++++++++++++
c
      subroutine maxwellevalvzrhs(
     &     md1a, md1b, md2a, md2b,
     &     m1a, m1b, m2a, m2b,
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
c ++++++++++++++
c
      subroutine maxwellsetembcs(
     &     md1a, md1b, md2a, md2b,
     &     m1a, m1b, m2a, m2b,
     &     EMvars,
     &     nx, ny,
     &     xPeriodic, yPeriodic,
     &     solution_order,
     &     c)
c
c..   set boundary conditions on Maxwell EM fields
      implicit none
c
c..   declaration of incoming variables
      integer md1a, md1b, md2a, md2b
      integer m1a, m1b, m2a, m2b
      real EMvars(md1a:md1b, md2a:md2b, 1:6)
      integer nx, ny
      integer xPeriodic, yPeriodic, solution_order
      real c
c
c..   declaration of local variables
      integer i1, i2, i3, i4, nghosts
      real u1,u2
      real w1,w2
c
      if (solution_order .eq. 4) then
        nghosts = 2
      else
        nghosts = 3
      end if

c     x direction boundary condition terms
      if (xPeriodic .eq. 0) then
c       Left edge
        if (m1a .eq. 0) then
          i1 = m1a
          do i2 = m2a, m2b
            do i4 = 1, nghosts
              do i3 = 1, 6
                ! first extrapolate
                EMvars(i1-i4, i2, i3) =
     *               +3.0*EMvars(i1-i4+1,i2,i3)
     *               -3.0*EMvars(i1-i4+2,i2,i3)
     *               +1.0*EMvars(i1-i4+3,i2,i3)
c                EMvars(i1-i4, i2, i3) =
c     *               +1.0*EMvars(i1-i4+1,i2,i3)
              end do
              u1 = EMvars(i1-i4,i2,2) ! Ey
              u2 = EMvars(i1-i4,i2,6) ! Bz
              w1 = +u1/(2.*c) + u2/2. ! right going
              w2 = -u1/(2.*c) + u2/2. ! left  going
              w1 = 0.0
              u1 = c*(w1-w2)
              u2 =    w1+w2
              EMvars(i1-i4,i2,2) = u1
              EMvars(i1-i4,i2,6) = u2
c
              u1 =-EMvars(i1-i4,i2,3) !-Ez
              u2 = EMvars(i1-i4,i2,5) ! By
              w1 = +u1/(2.*c) + u2/2. ! right going
              w2 = -u1/(2.*c) + u2/2. ! left  going
              w1 = 0.0
              u1 = c*(w1-w2)
              u2 =    w1+w2
              EMvars(i1-i4,i2,3) =-u1
              EMvars(i1-i4,i2,5) = u2
            end do
          end do
        end if
c       Right edge
        if (m1b .eq. nx-1) then
          i1 = m1b
          do i2 = m2a, m2b
            do i4 = 1, nghosts
              do i3 = 1, 6
                ! first extrapolate
                EMvars(i1+i4, i2, i3) =
     *               +3.0*EMvars(i1+i4-1,i2,i3)
     *               -3.0*EMvars(i1+i4-2,i2,i3)
     *               +1.0*EMvars(i1+i4-3,i2,i3)
c                EMvars(i1+i4, i2, i3) =
c     *               +1.0*EMvars(i1+i4-1,i2,i3)
              end do
              u1 = EMvars(i1+i4,i2,2) ! Ey
              u2 = EMvars(i1+i4,i2,6) ! Bz
              w1 = +u1/(2.*c) + u2/2. ! right going
              w2 = -u1/(2.*c) + u2/2. ! left  going
              w2 = 0.0
              u1 = c*(w1-w2)
              u2 =    w1+w2
              EMvars(i1+i4,i2,2) = u1
              EMvars(i1+i4,i2,6) = u2
c
              u1 =-EMvars(i1+i4,i2,3) !-Ez
              u2 = EMvars(i1+i4,i2,5) ! By
              w1 = +u1/(2.*c) + u2/2. ! right going
              w2 = -u1/(2.*c) + u2/2. ! left  going
              w2 = 0.0
              u1 = c*(w1-w2)
              u2 =    w1+w2
              EMvars(i1+i4,i2,3) =-u1
              EMvars(i1+i4,i2,5) = u2
            end do
          end do
        end if
      end if

c     y direction boundary condition terms
      if (yPeriodic .eq. 0) then
c       Bottom edge
        if (m2a .eq. 0) then
          i2 = m2a
          do i1 = m1a, m1b
            do i4 = 1, nghosts
              do i3 = 1, 6
                ! first extrapolate
                EMvars(i1, i2-i4, i3) =
     *               +3.0*EMvars(i1,i2-i4+1,i3)
     *               -3.0*EMvars(i1,i2-i4+2,i3)
     *               +1.0*EMvars(i1,i2-i4+3,i3)
c                EMvars(i1, i2-i4, i3) =
c     *               +1.0*EMvars(i1,i2-i4+1,i3)
              end do
              u1 =-EMvars(i1,i2-i4,1) !-Ex
              u2 = EMvars(i1,i2-i4,6) ! Bz
              w1 = +u1/(2.*c) + u2/2. ! right going
              w2 = -u1/(2.*c) + u2/2. ! left  going
              w1 = 0.0
              u1 = c*(w1-w2)
              u2 =    w1+w2
              EMvars(i1,i2-i4,1) =-u1
              EMvars(i1,i2-i4,6) = u2
c
              u1 = EMvars(i1,i2-i4,3) ! Ez
              u2 = EMvars(i1,i2-i4,4) ! Bx
              w1 = +u1/(2.*c) + u2/2. ! right going
              w2 = -u1/(2.*c) + u2/2. ! left  going
              w1 = 0.0
              u1 = c*(w1-w2)
              u2 =    w1+w2
              EMvars(i1,i2-i4,3) = u1
              EMvars(i1,i2-i4,4) = u2
            end do
          end do
        end if
c       Top edge
        if (m2b .eq. ny-1) then
          i2 = m2b
          do i1 = m1a, m1b
            do i4 = 1, nghosts
              do i3 = 1, 6
                ! first extrapolate
                EMvars(i1, i2+i4, i3) =
     *               +3.0*EMvars(i1,i2+i4-1,i3)
     *               -3.0*EMvars(i1,i2+i4-2,i3)
     *               +1.0*EMvars(i1,i2+i4-3,i3)
c                EMvars(i1, i2+i4, i3) =
c     *               +1.0*EMvars(i1,i2+i4-1,i3)
              end do
              u1 =-EMvars(i1,i2+i4,1) !-Ex
              u2 = EMvars(i1,i2+i4,6) ! Bz
              w1 = +u1/(2.*c) + u2/2. ! right going
              w2 = -u1/(2.*c) + u2/2. ! left  going
              w2 = 0.0
              u1 = c*(w1-w2)
              u2 =    w1+w2
              EMvars(i1,i2+i4,1) =-u1
              EMvars(i1,i2+i4,6) = u2
c
              u1 = EMvars(i1,i2+i4,3) ! Ez
              u2 = EMvars(i1,i2+i4,4) ! Bx
              w1 = +u1/(2.*c) + u2/2. ! right going
              w2 = -u1/(2.*c) + u2/2. ! left  going
              w2 = 0.0
              u1 = c*(w1-w2)
              u2 =    w1+w2
              EMvars(i1,i2+i4,3) = u1
              EMvars(i1,i2+i4,4) = u2
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
      subroutine maxwellsetvzbcs(
     &     md1a, md1b, md2a, md2b,
     &     m1a, m1b, m2a, m2b,
     &     vz,
     &     nx, ny,
     &     xPeriodic, yPeriodic,
     &     solution_order)
c
c..   set boundary conditions on Maxwell VZ fields
      implicit none
c
c..   declaration of incoming variables
      integer md1a, md1b, md2a, md2b
      integer m1a, m1b, m2a, m2b
      real vz(md1a:md1b, md2a:md2b)
      integer nx, ny, xPeriodic, yPeriodic, solution_order
c
c..   declaration of local variables
      integer i1, i2, i3, nghosts
c
      if (solution_order .eq. 4) then
        nghosts = 2
      else
        nghosts = 3
      end if
c     x direction boundary condition terms
      if (xPeriodic .eq. 0) then
c       Left edge
        if (m1a .eq. 0) then
          i1 = m1a
          do i2 = md2a, md2b
            do i3 = 1, nghosts
              vz(i1-i3, i2) = vz(i1+i3, i2)
            end do
          end do
        end if
c       Right edge
        if (m1b .eq. nx-1) then
          i1 = m1b
          do i2 = md2a, md2b
            do i3 = 1, nghosts
              vz(i1+i3, i2) = vz(i1-i3, i2)
            end do
          end do
        end if
      end if

c     y direction boundary condition terms
      if (yPeriodic .eq. 0) then
c       Bottom edge
        if (m2a .eq. 0) then
          i2 = m2a
          do i1 = md1a, md1b
            do i3 = 1, nghosts
              vz(i1, i2-i3) = vz(i1, i2+i3)
            end do
          end do
        end if
c       Top edge
        if (m2b .eq. ny-1) then
          i2 = m2b
          do i1 = md1a, md1b
            do i3 = 1, nghosts
              vz(i1, i2+i3) = vz(i1, i2-i3)
            end do
          end do
        end if
      end if
c
      return
      end
