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
      x_shape = e_ext_param(9)
      lwidth  = e_ext_param(10)
      x0      = e_ext_param(11)

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
      real phi_0, t0, t_ramp, t_off
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
      phi_0   = e_ext_param(5)
      t0      = e_ext_param(6)
      t_ramp  = e_ext_param(7)
      t_off   = e_ext_param(8)
c
        ! define Gaussian quadrature weights

      if ((t .lt. (t0+t_ramp+t_off)) .and. t .ge. t0) then
        if (t .lt. t0+t_ramp) then
          envel = 0.5+0.5*tanh(4.0*(2.0*(t-t0)/t_ramp-1.0))
        else
          envel = 0.5-0.5*tanh(4.0*(2.0*(t-t0-t_ramp)/t_off-1.0))
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
