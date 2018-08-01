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
c Fortran functions called by ShapedRampedCosineCurrentDriver.
c
      subroutine evaluateShapedRampedCurrentDriver( 
     &     Jx,
     &     Jy,
     &     Jz,
     &     nd1a, nd1b, nd2a, nd2b,
     &     n1a, n1b, n2a, n2b,
     &     xlo, xhi, dx,
     &     t,
     &     e_ext_param)
c
c.. function to compute applied current
      implicit none
c
      integer nd1a, nd1b, nd2a, nd2b
      integer n1a, n1b, n2a, n2b
      real xlo(2), xhi(2), dx(2), t
      real Jx(nd1a:nd1b, nd2a:nd2b)
      real Jy(nd1a:nd1b, nd2a:nd2b)
      real Jz(nd1a:nd1b, nd2a:nd2b)
      real e_ext_param(*)
c
      integer i1_apply, i2
c
      real width, apply_dir, shape, omega, J_0, t0, t_ramp, t_off, x0
      real envel, ycoord, J_ext, pi, i1_exact
c
      real one, four
c
      one = 1.0
      four = 4.0
c
      pi = four*atan(one)
c
      width     = e_ext_param(1)
      apply_dir = e_ext_param(2)
      shape     = e_ext_param(3)
      omega     = e_ext_param(4)
      J_0       = e_ext_param(5)
      t0        = e_ext_param(6)
      t_ramp    = e_ext_param(7)
      t_off     = e_ext_param(8)
      x0        = e_ext_param(9)
c
      if ((t .lt. (t0+t_ramp+t_off)) .and. t .ge. t0) then
        i1_exact = (x0 - xlo(1)) / dx(1)
        if (i1_exact .lt. 0) then
           i1_exact = i1_exact - 0.5
        else
           i1_exact = i1_exact + 0.5
        end if
        i1_apply = i1_exact
        if (t .lt. t0+t_ramp) then
          envel = 0.5+0.5*tanh(4.0*(2.0*(t-t0)/t_ramp-1.0))
        else
          envel = 0.5-0.5*tanh(4.0*(2.0*(t-t0-t_ramp)/t_off-1.0))
        end if
        do i2 = n2a, n2b
c
          ycoord = xlo(2)+dx(2)*i2
          if (abs(ycoord) .lt. 0.5*width) then
            J_ext = J_0*envel*cos(omega*(t-t0))*
     *            (1.0-shape*(sin(pi*ycoord/width))**2)
          else
            J_ext = J_0*envel*cos(omega*(t-t0))*(1.0-shape)
          endif
c
          if (apply_dir .eq. 1.0) then
             Jx(i1_apply, i2) = Jx(i1_apply, i2) + J_ext
          else if (apply_dir .eq. 2.0) then
             Jy(i1_apply, i2) = Jy(i1_apply, i2) + J_ext
          else
             Jz(i1_apply, i2) = Jz(i1_apply, i2) + J_ext
          end if
        end do
      end if
c
      return
      end
