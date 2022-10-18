c
c Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
c See the top-level LICENSE file for details.
c Produced at the Lawrence Livermore National Laboratory
c
c SPDX-License-Identifier: MIT
c
c Fortran functions called by GaussianCurrentDriver.
c
      subroutine evaluateGaussianCurrentDriver( 
     &     antenna_source,
     &     nd1a, nd1b, nd2a, nd2b,
     &     n1a, n1b, n2a, n2b,
     &     xlo, xhi, dx,
     &     t,
     &     shape_lo, shape_hi,
     &     e_ext_param,
     &     shaping,
     &     phase,
     &     antenna_loc_idx)
c
c.. function to compute applied current
      implicit none
c
      integer nd1a, nd1b, nd2a, nd2b
      integer n1a, n1b, n2a, n2b
      integer shape_lo, shape_hi, antenna_loc_idx
      real xlo(2), xhi(2), dx(2), t
      real antenna_source(nd1a:nd1b, nd2a:nd2b, 1:6)
      real e_ext_param(*)
      real shaping(shape_lo:shape_hi)
      real phase(shape_lo:shape_hi)
c
      integer i1, i2
c
      real plane, omega, J_0, t0, t_rampup, t_hold, t_rampdown
      real envel, xcoord, ycoord, J_amp, pi
c
      real one, four
c
      one = 1.0
      four = 4.0
c
      pi = four*atan(one)
c
      plane      = e_ext_param(1)
      omega      = e_ext_param(6)
      J_0        = e_ext_param(8)
      t0         = e_ext_param(9)
      t_rampup   = e_ext_param(10)
      t_hold     = e_ext_param(11)
      t_rampdown = e_ext_param(12)
c
      if ((t .lt. (t0+t_rampup+t_hold+t_rampdown)) .and. t .ge. t0) then
        if (t .lt. t0+t_rampup) then
          envel = 0.5+0.5*tanh(4.0*(2.0*(t-t0)/t_rampup-1.0))
        else if (t .lt. t0+t_rampup+t_hold) then
          envel = 0.5+0.5*tanh(4.0)
        else
          envel = 0.5-0.5*
     *      tanh(4.0*(2.0*(t-t0-t_rampup-t_hold)/t_rampdown-1.0))
        end if
        J_amp = J_0*envel
        if (plane .eq. 0.0) then
          do i2 = n2a, n2b
            antenna_source(antenna_loc_idx, i2, 3) =
     *        antenna_source(antenna_loc_idx, i2, 3) +
     *        J_amp*shaping(i2)*cos(omega*(t-t0)-phase(i2))
          end do
        else
          do i1 = n1a, n1b
            antenna_source(i1, antenna_loc_idx, 3) =
     *        antenna_source(i1, antenna_loc_idx, 3) +
     *        J_amp*shaping(i1)*cos(omega*(t-t0)-phase(i1))
          end do
        end if
      end if
c
      return
      end
