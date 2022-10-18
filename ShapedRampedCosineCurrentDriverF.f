c
c Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
c See the top-level LICENSE file for details.
c Produced at the Lawrence Livermore National Laboratory
c
c SPDX-License-Identifier: MIT
c
c Fortran functions called by ShapedRampedCosineCurrentDriver.
c
      subroutine evaluateShapedRampedCurrentDriver( 
     &     antenna_source,
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
      real antenna_source(nd1a:nd1b, nd2a:nd2b, 1:6)
      real e_ext_param(*)
c
      integer i1_apply, i2_apply, i1, i2
c
      real width, apply_dir, shape, omega, J_0, t0, t_rampup, t_hold
      real t_rampdown, x0, plane, envel, xcoord, ycoord, J_ext, pi
      real i1_exact, i2_exact
c
      real one, four
c
      one = 1.0
      four = 4.0
c
      pi = four*atan(one)
c
      width      = e_ext_param(1)
      apply_dir  = e_ext_param(2)
      shape      = e_ext_param(3)
      omega      = e_ext_param(4)
      J_0        = e_ext_param(5)
      t0         = e_ext_param(6)
      t_rampup   = e_ext_param(7)
      t_hold     = e_ext_param(8)
      t_rampdown = e_ext_param(9)
      x0         = e_ext_param(10)
      plane      = e_ext_param(11)
c
      if ((t .lt. (t0+t_rampup+t_hold+t_rampdown)) .and. t .ge. t0) then
        if (plane .eq. 0.0) then
          i1_exact = (x0 - xlo(1)) / dx(1)
          if (i1_exact .lt. 0) then
            i1_exact = i1_exact - 0.5
          else
            i1_exact = i1_exact + 0.5
          end if
          i1_apply = i1_exact
          if (i1_apply .ge. nd1a .and. i1_apply .le. nd1b) then
            if (t .lt. t0+t_rampup) then
              envel = 0.5+0.5*tanh(4.0*(2.0*(t-t0)/t_rampup-1.0))
            else if (t .lt. t0+t_rampup+t_hold) then
              envel = 0.5+0.5*tanh(4.0)
            else
              envel = 0.5-0.5*
     *          tanh(4.0*(2.0*(t-t0-t_rampup-t_hold)/t_rampdown-1.0))
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
              if (apply_dir .eq. 2.0) then
                antenna_source(i1_apply, i2, 2) =
     *            antenna_source(i1_apply, i2, 2) + J_ext
              else
                antenna_source(i1_apply, i2, 3) =
     *             antenna_source(i1_apply, i2, 3) + J_ext
              end if
            end do
          end if
        else
          i2_exact = (x0 - xlo(2)) / dx(2)
          if (i2_exact .lt. 0) then
            i2_exact = i2_exact - 0.5
          else
            i2_exact = i2_exact + 0.5
          end if
          i2_apply = i2_exact
          if (i2_apply .ge. nd2a .and. i2_apply .le. nd2b) then
            if (t .lt. t0+t_rampup) then
              envel = 0.5+0.5*tanh(4.0*(2.0*(t-t0)/t_rampup-1.0))
            else if (t .lt. t0+t_rampup+t_hold) then
              envel = 0.5+0.5*tanh(4.0)
            else
              envel = 0.5-0.5*
     *          tanh(4.0*(2.0*(t-t0-t_rampup-t_hold)/t_rampdown-1.0))
            end if
            do i1 = n1a, n1b
c
              xcoord = xlo(1)+dx(1)*i1
              if (abs(xcoord) .lt. 0.5*width) then
                J_ext = J_0*envel*cos(omega*(t-t0))*
     *            (1.0-shape*(sin(pi*xcoord/width))**2)
              else
                J_ext = J_0*envel*cos(omega*(t-t0))*(1.0-shape)
              endif
c
              if (apply_dir .eq. 1.0) then
                antenna_source(i1, i2_apply, 1) =
     *            antenna_source(i1, i2_apply, 1) + J_ext
              else
                antenna_source(i1, i2_apply, 3) =
     *            antenna_source(i1, i2_apply, 3) + J_ext
              end if
            end do
          end if
        end if
      end if
c
      return
      end
