c
c Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
c See the top-level LICENSE file for details.
c Produced at the Lawrence Livermore National Laboratory
c
c SPDX-License-Identifier: MIT
c
c Fortran functions called by ShapedRampedCosineElectricPotentialDriver.
c
      subroutine jeval_ghat_ext_pot(
     *   ghat, x, t, e_ext_param, pi)
c
c.. function to evaluate x-coordinate part of spatial external potential
      implicit none
c
      real ghat, x, t
      real e_ext_param(*)
      real pi
c
      real xwidth, omega, t0, x_shape, lwidth, x0
c
      xwidth  = e_ext_param(1)
      omega   = e_ext_param(4)
      t0      = e_ext_param(6)
      x_shape = e_ext_param(10)
      lwidth  = e_ext_param(11)
      x0      = e_ext_param(12)

      if (.false.) then
        !! temporary for JWB
c        write(6,*)'xwidth  ', pi/xwidth
        ghat = cos(pi*x/xwidth-omega*(t-t0))
      else
        if (abs(x-x0) .lt. 0.5*lwidth) then
          ghat = 1.0-x_shape*(sin(pi*(x-x0)/lwidth))**2
        else
          ghat = 1.0-x_shape
        end if
        ghat = ghat*cos(pi*x/xwidth-omega*(t-t0))
      end if

      return
      end
c
c+++++++++++
c
      subroutine jeval_h_ext_pot(
     *   h, y, e_ext_param, pi)
c
c.. function to evaluate y-coordinate part of spatial external potential
      implicit none
c
      real h, y
      real e_ext_param(*)
      real pi
c
      real ywidth, shape
c
      ywidth  = e_ext_param(2)
      shape   = e_ext_param(3)

      if (.false.) then
        !! temporary for JWB
c        write(6,*)'ywidth  ', pi/ywidth
        h = cos(pi*y/ywidth)
c        h = 1.0
      else
        if (abs(y) .lt. 0.5*ywidth) then
          h = 1.0-shape*(sin(pi*y/ywidth))**2
        else
          h = 1.0-shape
        end if
      end if

      return
      end
c
c+++++++++++
c
      subroutine evaluateShapedRampedPotentialDriver(
     &     phi,
     &     nd1a, nd1b, nd2a, nd2b,
     &     xlo, xhi, dx,
     &     t,
     &     e_ext_param)
c
c.. function to compute external potential
      implicit none
c
      integer nd1a, nd1b, nd2a, nd2b
      real xlo(2), xhi(2), dx(2), t
      real phi(nd1a:nd1b, nd2a:nd2b)
      real e_ext_param(*)
c
      integer i1, i2
c
      real phi_0, t0, t_rampup, t_hold, t_rampdown
      real envel, xcoord, ycoord, phi_ext
      real g
      real h
      real pi
c
      real one, four
c
      one = 1.0
      four = 4.0
c
      pi = four*atan(one)
c
      phi_0      = e_ext_param(5)
      t0         = e_ext_param(6)
      t_rampup   = e_ext_param(7)
      t_hold     = e_ext_param(8)
      t_rampdown = e_ext_param(9)
c
        ! define Gaussian quadrature weights

      if ((t .lt. (t0+t_rampup+t_hold+t_rampdown)) .and. t .ge. t0) then
        if (t .lt. t0+t_rampup) then
          envel = 0.5+0.5*tanh(4.0*(2.0*(t-t0)/t_rampup-1.0))
        else if (t .lt. t0+t_rampup+t_hold) then
          envel = 0.5+0.5*tanh(4.0)
        else
          envel = 0.5-0.5*
     *      tanh(4.0*(2.0*(t-t0-t_rampup-t_hold)/t_rampdown-1.0))
        end if
        do i2 = nd2a, nd2b
        do i1 = nd1a, nd1b
          xcoord = xlo(1)+dx(1)*(0.5+i1)
          ycoord = xlo(2)+dx(2)*(0.5+i2)
c
          call jeval_ghat_ext_pot(g, xcoord, t, e_ext_param, pi)
c
          call jeval_h_ext_pot(h, ycoord, e_ext_param, pi)
c
          phi_ext = phi_0*h*g
c
          phi(i1, i2) = phi(i1, i2)+envel*phi_ext
        end do
        end do
      end if
c
      return
      end
