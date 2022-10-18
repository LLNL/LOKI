c
c Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
c See the top-level LICENSE file for details.
c Produced at the Lawrence Livermore National Laboratory
c
c SPDX-License-Identifier: MIT
c
c Fortran functions called by UnidirectionalCurrentDriver.
c
      subroutine evaluateUnidirectionalCurrentDriver( 
     &     antenna_source,
     &     status,
     &     omega_eff2,
     &     nd1a, nd1b, nd2a, nd2b,
     &     n1a, n1b, n2a, n2b,
     &     xlo, xhi, dx,
     &     t,
     &     e_ext_param, R)
c
c.. function to compute applied current
      implicit none
c
c.. declarations of incoming variables
      integer status
      integer nd1a, nd1b, nd2a, nd2b
      integer n1a, n1b, n2a, n2b
      real xlo(2), xhi(2), dx(2), t
      real antenna_source(nd1a:nd1b, nd2a:nd2b, 1:6)
      real omega_eff2(nd1a:nd1b, nd2a:nd2b)
      real e_ext_param(*)
      real R(1:6, 1:2)
c
c.. declarations of local variables
      integer i1, i2, comp
      real omega, J_0, t0, t_rampup, t_hold, t_rampdown, polarAngle
      real x0x, x0y, x1x, x1y, beta, light_speed, shape, width
      real J_0_inPlane, J_0_outPlane, L, phi
      real one, four, pi, norm_fact, omega2, this_omega_eff2, tau
      real K(1:2), envel
      real xcoord, ycoord, coord_dx, coord_dy, alpha, Rad, delta, D
      real perp_falloff, par_shape
c
      status = 0
      one = 1.0
      four = 4.0
c
      pi = four*atan(one)
c
      omega        = e_ext_param(1)
      J_0          = e_ext_param(2)
      t0           = e_ext_param(3)
      t_rampup     = e_ext_param(4)
      t_hold       = e_ext_param(5)
      t_rampdown   = e_ext_param(6)
      polarAngle   = e_ext_param(7)
      x0x          = e_ext_param(8)
      x0y          = e_ext_param(9)
      x1x          = e_ext_param(10)
      x1y          = e_ext_param(11)
      beta         = e_ext_param(12)
      light_speed  = e_ext_param(13)
      shape        = e_ext_param(14)
      width        = e_ext_param(15)
      J_0_inplane  = e_ext_param(16)
      J_0_outplane = e_ext_param(17)
      L            = e_ext_param(18)
      phi          = e_ext_param(19)
      norm_fact = beta/sqrt(2*pi)
      omega2 = omega**2
c
c     Only do anything if the driver is on.
      if ((t .lt. (t0+t_rampup+t_hold+t_rampdown)) .and. t .ge. t0) then
c       Compute the appropriate modulating envelope.
        if (t .lt. t0+t_rampup) then
          envel = 0.5+0.5*tanh(4.0*(2.0*(t-t0)/t_rampup-1.0))
        else if (t .lt. t0+t_rampup+t_hold) then
          envel = 0.5+0.5*tanh(4.0)
        else
          envel = 0.5-0.5*
     *      tanh(4.0*(2.0*(t-t0-t_rampup-t_hold)/t_rampdown-1.0))
        end if
        do i2 = n2a, n2b
          ycoord = xlo(2)+dx(2)*(i2+0.5)
          do i1 = n1a, n1b
            xcoord = xlo(1)+dx(1)*(i1+0.5)
c           Compute the angle of the vector from the "from_pt" to this
c           grid point.  Then compute the angle between this vector
c           and the vector from the "from_pt" to the "to_pt".
            coord_dx = xcoord-x0x
            coord_dy = ycoord-x0y
            alpha = phi-atan2(coord_dy, coord_dx)
c           Now compute the perpendicular distance of this grid point
c           from the center of the emitter.  Also determine the overlap
c           with the emitter.
            Rad = sqrt(coord_dx**2+coord_dy**2)
            delta = Rad*sin(alpha)
            D = Rad*cos(alpha)
c           If this point overlaps the emitter then it contributes to
c           the antenna source term,
            if (D .gt. 0.0 .and. D .lt. L) then
c             Compute the current falloff as a Gaussian of the
c             perpendicular distance of the grid point from the
c             emitter.
              ! avoid denormalized results from the exponential
              if (abs(beta*delta) .le. 37.0) then
c               Account for the dispersion relation.
                this_omega_eff2 = omega_eff2(i1,i2)
                if (omega2 .le. this_omega_eff2) then
                  status = 1
                  return
                else
                  tau =
     *              delta*sqrt(1.0-this_omega_eff2/omega2)/light_speed
                endif
                perp_falloff = cos(omega*(t-t0-tau))*
     *            norm_fact*exp(-0.5*(beta*delta)**2)
              else
                perp_falloff = 0.0
              end if
c             Compute any current shaping.
              if (shape .ne. 0.0) then
                if (D .gt. (L-width)/2.0 .and.
     *              D .lt. (L+width)/2.0) then
                  par_shape = 1-shape*(sin(pi*(D-L/2.0)/width))**2
                else
                  par_shape = 1-shape
                endif
              else
                par_shape = 1.0
              end if
c             Now compute the actual antenna source at this grid point.
              K(1) = J_0_inPlane*perp_falloff*par_shape
              K(2) = J_0_outPlane*perp_falloff*par_shape
              do comp = 1, 6
                antenna_source(i1,i2,comp) = antenna_source(i1,i2,comp)+
     *            envel*(R(comp,1)*K(1)+R(comp,2)*K(2))
              end do
            end if
          end do
        end do
      end if
c
      return
      end
