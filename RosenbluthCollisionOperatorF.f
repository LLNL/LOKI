c
c Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
c See the top-level LICENSE file for details.
c Produced at the Lawrence Livermore National Laboratory
c
c SPDX-License-Identifier: MIT
c
c Fortran functions called by RosenbluthCollisionOperator.
c

      subroutine getRamp(A, v, vra, vrb, vmin, vmax, solution_order)
c
c.. evaluate ramp function
      implicit none
c
c.. declarations of incoming variables 
      real A, v, vra, vrb, vmin, vmax
      integer solution_order
c
c.. declarations of local variables
      real va, vb, xi
c
      if (v .lt. vra .and. v .ge. vmin) then
        va = vra
        vb = vmin
        xi = (v-va)/(vb-va)
      else if (v .gt. vrb .and. v .le. vmax) then
        va = vrb
        vb = vmax
        xi = (v-va)/(vb-va)
      else if (v .lt. vmin .or. v .gt. vmax) then
        xi = 1.0
      else
        xi = 0.0
      end if

      if (solution_order .eq. 4) then
         A = 1.0 + xi**4*(
     *     +20.0*xi**3
     *     -70.0*xi**2
     *     +84.0*xi
     *     -35.0)
      else
         A = 1.0 + xi**6*(
     *     +252.0 *xi**5
     *     -1386.0*xi**4
     *     +3080.0*xi**3
     *     -3465.0*xi**2
     *     +1980.0*xi
     *     -462.0)
      end if
c
      return
      end
c
c +++++++++++++
c
      subroutine getD(d1, d2, d3, vTh, C, vx, vy, vxgrid, vygrid,
     *     range_lo, range_hi, vxmin, vxmax, vymin, vymax,
     *     solution_order, sqrt2divpi)
c
c.. compute diffusion tensor
      implicit none
c
c.. declarations of incoming variables
      real d1, d2, d3, vTh, C, vx, vy, vxgrid, vygrid
      real vxmin, vxmax, vymin, vymax, sqrt2divpi
      real range_lo(1:2), range_hi(1:2)
      integer solution_order
c
c.. declarations of local variables
      real vx2, vy2, v, v2, v3, v4, v6, v8
      real f, g, Kh, Hh, Ax, Ay, factor, eps
      real theta, ctheta, ctheta2, sctheta
c
      eps = 1.e-10
      call getRamp(Ax, vxgrid, range_lo(1), range_hi(1), vxmin, vxmax,
     *             solution_order)
      call getRamp(Ay, vygrid, range_lo(2), range_hi(2), vymin, vymax,
     *             solution_order)
      vx2 = vx*vx
      vy2 = vy*vy
      v  = sqrt(vx2+vy2)
      v = v/vTh
c     For small velocity (<.01) compute D with a Taylor series expansion
C     through v**8 which will be accurate to machine precission,  For all
c     other velocities use the explicit expression for D.
      if (v .lt. 0.01) then
        factor = Ax*Ay*C*vTh*vTh*sqrt2divpi
c     Handle the cases where vx and/or vy are essentially 0 and the
c     conversion to polar coordinates is ill defined.
        if (abs(vx) .le. eps .and. abs(vy) .le. eps) then
          theta = 0.0
          v = 0.0
        else if (abs(vx) .le. eps) then
          theta = 2.0*atan(1.0)
        else if (abs(vy) .le. eps) then
          theta = 0.0
        else
          theta = atan(vy/vx)
        end if
        v2 = v*v
        v4 = v2*v2
        v6 = v4*v2
        v8 = v6*v2
        ctheta = cos(theta)
        ctheta2 = ctheta*ctheta
        sctheta = sin(theta)*ctheta
        d1 = factor*
     *    (2.0/3.0-
     *     (2.0*ctheta2+1)*v2/15.0+
     *     (4.0*ctheta2+1)*v4/140.0-
     *     (6.0*ctheta2+1)*v6/1512.0+
     *     (8.0*ctheta2+1)*v8/19008.0)
        d2 = factor*
     *    (-2.0*sctheta*v2/15.0+
     *     sctheta*v4/35.0-
     *     sctheta*v6/252.0+
     *     sctheta*v8/2376.0)
        d3 = factor*
     *    (2.0/3.0+
     *     (2.0*ctheta2-3.0)*v2/15.0-
     *     (4.0*ctheta2-5.0)*v4/140.0+
     *     (6.0*ctheta2-7.0)*v6/1512.0-
     *     (8.0*ctheta2-9.0)*v8/19008.0)
      else
        factor = Ax*Ay*C
        v2 = v*v
        v3 = v*v2
        g = derf(v/sqrt(2.0))
        f = sqrt2divpi*v*exp(-v2/2.0)
        Kh = 1.0/v3*((v2-1.0)*g+f)
        Hh = 1.0/v3*(g-f)
        d1 = factor*(2.0*Hh*vx2+Kh*vy2)/v2
        d2 = factor*vx*vy*(2.0*Hh-Kh)/v2
        d3 = factor*(2.0*Hh*vy2+Kh*vx2)/v2
      end if
      return
      end
c
c +++++++++++++
c
      subroutine computeRosenbluthDiffusionTensor(
     *     d, reLam, mass_ratio, rVx0, rVy0, IvTh, N,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,      
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     xlo, xhi, dx, velocities, BR, compute_re_lam,
     *     range_lo, range_hi, dparams, iparams)
c
c.. compute constants needed by collision operator
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
c
      real dparams(*)
      integer iparams(*)
c
      real d(n1a:n1b, n2a:n2b, n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2), 1:3)
      real reLam(n1a:n1b, n2a:n2b, n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2))
      real mass_ratio
      real rVx0(nd1a:nd1b, nd2a:nd2b)
      real rVy0(nd1a:nd1b, nd2a:nd2b)
      real IvTh(nd1a:nd1b, nd2a:nd2b)
      real N(nd1a:nd1b, nd2a:nd2b)
c
      real xlo(1:4), xhi(1:4), dx(1:4)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      integer BR, compute_re_lam
      real range_lo(1:2), range_hi(1:2)
c
c.. declarations of local variables
      integer i1, i2, i3, i4, solution_order, modulate_nu, do_relativity
      real vrolloff, C, pi, sqrt2divpi
      real vxmin, vxmax, vymin, vymax, vxgrid, vygrid, vx, vy
      real vTh, nu_base, nu
      real d1, d2, d3
      real vxlo, vxhi, vylo, vyhi, dvx, dvy
c
      nu_base = dparams(1)
      solution_order = iparams(1)
      modulate_nu = iparams(4)
      do_relativity = iparams(5)
      if (solution_order .eq. 4) then
        vrolloff = 3.0
      else
        vrolloff = 4.0
      end if
      vxlo = xlo(3)
      vxhi = xhi(3)
      vylo = xlo(4)
      vyhi = xhi(4)
      dvx = dx(3)
      dvy = dx(4)
      pi = 4.0*atan(1.0)
      sqrt2divpi = sqrt(2.0/pi)
      vxmin = vxlo + vrolloff*dvx
      vxmax = vxhi - vrolloff*dvx
      vymin = vylo + vrolloff*dvy
      vymax = vyhi - vrolloff*dvy

      do i4 = n4a-iparams(2), n4b+iparams(2)
        do i3 = n3a-iparams(2), n3b+iparams(2)
          vxgrid = velocities(i3, i4, 0)
          vygrid = velocities(i3, i4, 1)

          do i2 = n2a, n2b
            do i1 = n1a, n1b
              if (BR .eq. 1) then
                vTh = IvTh(i1, i2) / sqrt(mass_ratio)
                C = 1.0/(vTh**2)
              else
                vTh = IvTh(i1, i2)
                C = N(i1, i2) / ((vTh**3) * mass_ratio)
              end if
              vx = vxgrid - rVx0(i1, i2)
              vy = vygrid - rVy0(i1, i2)
              ! For relativity the grid is momentum based
              if (do_relativity .eq. 1) then
                vxgrid = vxlo + (i3+0.5)*dvx
                vygrid = vylo + (i4+0.5)*dvy
              end if
              call getD(d1, d2, d3, vTh, C, vx, vy, vxgrid, vygrid,
     *                  range_lo, range_hi, vxmin, vxmax, vymin, vymax,
     *                  solution_order, sqrt2divpi)
              d(i1,i2,i3,i4,1) = d1
              d(i1,i2,i3,i4,2) = d2
              d(i1,i2,i3,i4,3) = d3
              if (BR .eq. 0 .and. compute_re_lam .eq. 1) then
                if (modulate_nu .eq. 1) then
                  nu = nu_base*N(i1,i2)/IvTh(i1,i2)
                else
                  nu = nu_base
                end if
                reLam(i1,i2,i3,i4) = reLam(i1,i2,i3,i4) + nu*max(d1, d2)
              end if
            end do
          end do
        end do
      end do
      
      return
      end
c
c +++++++++++++
c
      subroutine addToDoubleArray4D(
     *     a, b, ca, cb,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b)
c
c.. compute a = ca*a + cb*b
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
c
      real a(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real b(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real ca, cb
c
c.. declarations of local variables
      integer i1, i2, i3, i4
c
      do i4 = nd4a, nd4b
        do i3 = nd3a, nd3b
          do i2 = nd2a, nd2b
             do i1 = nd1a, nd1b
                a(i1, i2, i3, i4) = ca * a(i1, i2, i3, i4) +
     *               cb * b(i1, i2, i3, i4)
            end do
          end do
        end do
      end do

      return
      end
c
c +++++++++++++
c
      subroutine computeRosenbluthTemps(
     *     uMomx, uMomy, uKE, IVx, IVy, Vth, n,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     velocities)
c
c.. compute temps needed for back reaction terms
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
c
      real uMomx(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real uMomy(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real uKE(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real IVx(nd1a:nd1b, nd2a:nd2b)
      real IVy(nd1a:nd1b, nd2a:nd2b)
      real Vth(nd1a:nd1b, nd2a:nd2b)
      real n(nd1a:nd1b, nd2a:nd2b)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
c     
c.. declarations of local variables
      integer i1, i2, i3, i4
      real alpha, pi, vx, vy, vxt, vyt, fMTemp
c
      pi = 4.0*atan(1.0)

      do i4 = nd4a, nd4b
        do i3 = nd3a, nd3b
          vxt = velocities(i3, i4, 0)
          vyt = velocities(i3, i4, 1)
          do i2 = nd2a, nd2b
            do i1 = nd1a, nd1b
               vx = vxt - IVx(i1,i2)
               vy = vyt - IVy(i1,i2)
               alpha = 1.0 / (Vth(i1,i2) * Vth(i1,i2))
               fMTemp = n(i1,i2)*alpha/(2.0*pi)*
     *           exp(-0.5*alpha*(vx**2+vy**2))
               uMomx(i1, i2, i3, i4) = vx*fMTemp
               uMomy(i1, i2, i3, i4) = vy*fMTemp
               uKE(i1, i2, i3, i4) = (vx**2+vy**2)*fMTemp
            end do
          end do
        end do
      end do

      return
      end
c
c +++++++++++++
c
      subroutine computeRosenbluthDenoms(
     *     rMomx, rMomy, rKE,
     *     cMomx, cMomy, cKE,
     *     IVx, IVy,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     velocities)
c
c.. compute back reaction denominator terms that will be reduced
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
c
      real rMomx(nd1a:nd1b, nd2a:nd2b)
      real rMomy(nd1a:nd1b, nd2a:nd2b)
      real rKE(nd1a:nd1b, nd2a:nd2b)
      real cMomx(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real cMomy(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real cKE(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real IVx(nd1a:nd1b, nd2a:nd2b)
      real IVy(nd1a:nd1b, nd2a:nd2b)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      real vx, vy, vxt, vyt
c
      do i4 = n4a, n4b
        do i3 = n3a, n3b
          vxt = velocities(i3, i4, 0)
          vyt = velocities(i3, i4, 1)
          do i2 = nd2a, nd2b
             do i1 = nd1a, nd1b
                vx = vxt - IVx(i1,i2)
                vy = vyt - IVy(i1,i2)

                rMomx(i1, i2) = rMomx(i1, i2) + vx*cMomx(i1, i2, i3, i4)
                rMomy(i1, i2) = rMomy(i1, i2) + vy*cMomy(i1, i2, i3, i4)
                rKE(i1, i2) =
     *               rKE(i1, i2) + (vx**2+vy**2)*cKE(i1, i2, i3, i4)
             end do
          end do
       end do
      end do

      return
      end
c
c +++++++++++++
c
      subroutine computeRosenbluthNumerators(
     *     rMomx, rMomy, rKE,
     *     IVx, IVy, dfa,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     velocities)
c
c.. compute back reaction denominator terms that will be reduced
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
c
      real rMomx(nd1a:nd1b, nd2a:nd2b)
      real rMomy(nd1a:nd1b, nd2a:nd2b)
      real rKE(nd1a:nd1b, nd2a:nd2b)
      real IVx(nd1a:nd1b, nd2a:nd2b)
      real IVy(nd1a:nd1b, nd2a:nd2b)
      real dfa(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
c
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      real vx, vy, vxt, vyt
c
      do i4 = n4a, n4b
        do i3 = n3a, n3b
          vxt = velocities(i3, i4, 0)
          vyt = velocities(i3, i4, 1)
          do i2 = nd2a, nd2b
             do i1 = nd1a, nd1b
                vx = vxt - IVx(i1,i2)
                vy = vyt - IVy(i1,i2)

              rMomx(i1, i2) = rMomx(i1, i2)+ vx * dfa(i1,i2,i3,i4)
              rMomy(i1, i2) = rMomy(i1, i2)+ vy * dfa(i1,i2,i3,i4)
              rKE(i1, i2) = rKE(i1, i2)+(vx**2+vy**2)*dfa(i1,i2,i3,i4)
            end do
          end do
        end do
      end do

      return
      end
c
c +++++++++++++
c
      subroutine computeRosenbluthEntropyChange(
     *     rEnt, u, dudt,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b)
c
c.. compute interspecies moments
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
c
      real rEnt(nd1a:nd1b, nd2a:nd2b)
      real u(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real dudt(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      real t1
c
      rEnt = 0.0
      
      do i4 = n4a, n4b
        do i3 = n3a, n3b
          do i2 = nd2a, nd2b
             do i1 = nd1a, nd1b
                t1 = log(abs(u(i1, i2, i3, i4)) + 1.0e-15) + 1.0
                rEnt(i1, i2) = rEnt(i1, i2) - dudt(i1, i2, i3, i4) * t1
            end do
          end do
        end do
      end do

      return
      end
c
c +++++++++++++
c
      subroutine computeRosenbluthSpeciesMoments(
     *     rN, rGammax, rGammay, rKE, u, diagnostics,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     velocities)
c
c.. compute interspecies moments
      implicit none
c
c.. declarations of incoming variables
      integer diagnostics
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
c
      real rN(nd1a:nd1b, nd2a:nd2b)
      real rGammax(nd1a:nd1b, nd2a:nd2b)
      real rGammay(nd1a:nd1b, nd2a:nd2b)
      real rKE(nd1a:nd1b, nd2a:nd2b)
      real u(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      real vx, vy, uval, eps
c
      eps = 1.0e-10
      rN = 0.0
      rGammax = 0.0
      rGammay = 0.0
      rKE = 0.0
      
      do i4 = n4a, n4b
        do i3 = n3a, n3b
          vx = velocities(i3, i4, 0)
          vy = velocities(i3, i4, 1)
          do i2 = nd2a, nd2b
            do i1 = nd1a, nd1b
              if (diagnostics .eq. 0) then
                uval = max(abs(u(i1, i2, i3, i4)), eps)
              else
                uval = u(i1, i2, i3, i4)
              end if
              if (diagnostics .eq. 0) then
                rN(i1, i2) = rN(i1, i2) + uval
              else
                rKE(i1, i2) = rKE(i1, i2) + (vx**2+vy**2)*uval
              end if
              rGammax(i1, i2) = rGammax(i1, i2) + vx*uval
              rGammay(i1, i2) = rGammay(i1, i2) + vy*uval
            end do
          end do
        end do
      end do

      return
      end
c
c +++++++++++++
c
      subroutine computeRosenbluthSpeciesKEC(
     *     rKEC, rVx0, rVy0, u,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     velocities)
c
c.. compute interspecies moments
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
c
      real rKEC(nd1a:nd1b, nd2a:nd2b)
      real rVx0(nd1a:nd1b, nd2a:nd2b)
      real rVy0(nd1a:nd1b, nd2a:nd2b)
      real u(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      real vx, vy, uval, eps
c
      eps = 1.0e-10
      do i4 = n4a, n4b
        do i3 = n3a, n3b
          vx = velocities(i3, i4, 0)
          vy = velocities(i3, i4, 1)
          do i2 = nd2a, nd2b
            do i1 = nd1a, nd1b
               uval = max(abs(u(i1, i2, i3, i4)), eps)
               rKEC(i1, i2) = rKEC(i1, i2) + ((vx - rVx0(i1, i2))**2
     *           +(vy - rVy0(i1, i2))**2)*uval
            end do
          end do
        end do
      end do

      return
      end
c
c +++++++++++++
c
      subroutine computeRosenbluthSpeciesReducedFields(
     *     vx, vy, n, Gammax, Gammay,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b)
c
c.. compute interspecies moments
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
c
      real vx(nd1a:nd1b, nd2a:nd2b)
      real vy(nd1a:nd1b, nd2a:nd2b)
      real n(nd1a:nd1b, nd2a:nd2b)
      real Gammax(nd1a:nd1b, nd2a:nd2b)
      real Gammay(nd1a:nd1b, nd2a:nd2b)
c
c.. declarations of local variables
      integer i1, i2
c
      do i2 = nd2a, nd2b
        do i1 = nd1a, nd1b
          vx(i1, i2) = Gammax(i1, i2) / n(i1, i2)
          vy(i1, i2) = Gammay(i1, i2) / n(i1, i2)
        end do
      end do

      return
      end
c
c +++++++++++++
c
      subroutine getCommunicatedSpeciesMoments(
     *     data, n, vx, vy, vthsq,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b, nf)
c
c.. compute interspecies moments
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b, nf
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
c
      real data(n1a:n1b, n2a:n2b, 1:nf)
      real n(nd1a:nd1b, nd2a:nd2b)
      real vx(nd1a:nd1b, nd2a:nd2b)
      real vy(nd1a:nd1b, nd2a:nd2b)
      real vthsq(nd1a:nd1b, nd2a:nd2b)
c
c.. declarations of local variables
      integer i1, i2
c
      do i2 = n2a, n2b
        do i1 = n1a, n1b
          data(i1, i2, 1) = n(i1, i2)
          data(i1, i2, 2) = vx(i1, i2)
          data(i1, i2, 3) = vy(i1, i2)
          data(i1, i2, 4) = vthsq(i1, i2)
        end do
      end do

      return
      end      
c
c +++++++++++++
c
      subroutine getCommunicatedBRSpeciesMoments(
     *     data, IMomxN, IMomyN, IKEN,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,nf,other_species_loc)
c
c.. compute interspecies moments
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
      integer nf, other_species_loc
c
      real data(n1a:n1b, n2a:n2b, 1:nf)
      real IMomxN(nd1a:nd1b, nd2a:nd2b)
      real IMomyN(nd1a:nd1b, nd2a:nd2b)
      real IKEN(nd1a:nd1b, nd2a:nd2b)
c
c.. declarations of local variables
      integer i1, i2
c
      do i2 = n2a, n2b
        do i1 = n1a, n1b
          data(i1, i2, 3*other_species_loc+1) = IKEN(i1, i2)
          data(i1, i2, 3*other_species_loc+2) = IMomxN(i1, i2)
          data(i1, i2, 3*other_species_loc+3) = IMomyN(i1, i2)
        end do
      end do

      return
      end      
c
c +++++++++++++
c
      subroutine setCommunicatedSpeciesMoments(
     *     n, vx, vy, vthsq, data,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,nf)
c
c.. compute interspecies moments
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b, nf
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
c
      real n(nd1a:nd1b, nd2a:nd2b)
      real vx(nd1a:nd1b, nd2a:nd2b)
      real vy(nd1a:nd1b, nd2a:nd2b)
      real vthsq(nd1a:nd1b, nd2a:nd2b)
      real data(n1a:n1b, n2a:n2b, 1:nf)
c
c.. declarations of local variables
      integer i1, i2
c
      do i2 = n2a, n2b
        do i1 = n1a, n1b
          n(i1, i2) = data(i1, i2, 1)
          vx(i1, i2) = data(i1, i2, 2)
          vy(i1, i2) = data(i1, i2, 3)
          vthsq(i1, i2) = data(i1, i2, 4)
        end do
      end do

      return
      end      
c
c +++++++++++++
c
      subroutine setCommunicatedSpeciesBRMoments(
     *     IMomxN, IMomyN, IKEN, data,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,nf)
c
c.. compute interspecies moments
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
      integer nf
c
      real IMomxN(nd1a:nd1b, nd2a:nd2b)
      real IMomyN(nd1a:nd1b, nd2a:nd2b)
      real IKEN(nd1a:nd1b, nd2a:nd2b)
      real data(n1a:n1b, n2a:n2b, 1:nf)
c
c.. declarations of local variables
      integer i1, i2
c
      do i2 = n2a, n2b
        do i1 = n1a, n1b
          IKEN(i1, i2) = data(i1, i2, 1)
          IMomxN(i1, i2) = data(i1, i2, 2)
          IMomyN(i1, i2) = data(i1, i2, 3)
        end do
      end do

      return
      end      
c
c +++++++++++++
c
      subroutine computeRosenbluthSpeciesVthermal(
     *     vthsq, KEC, n,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b)
c
c.. compute interspecies moments
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
c
      real vthsq(nd1a:nd1b, nd2a:nd2b)
      real KEC(nd1a:nd1b, nd2a:nd2b)
      real n(nd1a:nd1b, nd2a:nd2b)
c
c.. declarations of local variables
      integer i1, i2
c
      do i2 = nd2a, nd2b
        do i1 = nd1a, nd1b
          vthsq(i1, i2) = sqrt(0.5 * KEC(i1, i2) / n(i1, i2))
        end do
      end do

      return
      end
c
c +++++++++++++
c
      subroutine appendRosenbluthBackReaction(
     *     rhs,
     *     cMomx, cMomy, cKE,
     *     IMomxN, IMomyN, IKEN,
     *     IMomxD, IMomyD, IKED,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     massR)
c
c.. compute back reaction denominator terms that will be reduced
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
c
      real rhs(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real cMomx(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real cMomy(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real cKE(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real IMomxN(nd1a:nd1b, nd2a:nd2b)
      real IMomyN(nd1a:nd1b, nd2a:nd2b)
      real IKEN(nd1a:nd1b, nd2a:nd2b)
      real IMomxD(nd1a:nd1b, nd2a:nd2b)
      real IMomyD(nd1a:nd1b, nd2a:nd2b)
      real IKED(nd1a:nd1b, nd2a:nd2b)
c
      real massR
c
c.. declarations of local variables
      integer i1, i2, i3, i4
c
      do i4 = n4a, n4b
        do i3 = n3a, n3b
          do i2 = n2a, n2b
            do i1 = n1a, n1b
              rhs(i1, i2, i3, i4) = rhs(i1, i2, i3, i4) - massR*(
     *          IMomxN(i1, i2)/IMomxD(i1, i2)*cMomx(i1, i2, i3, i4) +
     *          IMomyN(i1, i2)/IMomyD(i1, i2)*cMomy(i1, i2, i3, i4) +
     *          IKEN(i1, i2)/IKED(i1, i2)*cKE(i1, i2, i3, i4))
            end do
          end do
        end do
      end do

      return
      end
c
c +++++++++++++
c
      function safelog(arg)
c
c.. compute the log of max(arg, tol)
      implicit none
c
c.. declarations of incoming variables
      real safelog, arg
c
c.. declarations of local variables
      real tol
c
      tol = 1.0e-16
      safelog = log(max(arg, tol))
      return
      end function
c
c +++++++++++++
c
      subroutine getmaxrelam(
     *     maxreLam, reLam,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     iparams)
c
c.. find the maximum of reLam
      implicit none
c
c.. declarations of incoming variables
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
c
      integer iparams(*)
c
      real maxreLam
      real reLam(n1a:n1b, n2a:n2b,
     #n3a-iparams(2):n3b+iparams(2), n4a-iparams(2):n4b+iparams(2))
c
c..   declarations of local variables
      integer i1, i2, i3, i4
      real thisreLam
c
      maxreLam = 0.0
      do i4 = n4a, n4b
      do i3 = n3a, n2b
      do i2 = n2a, n2b
      do i1 = n1a, n1b
        thisreLam = reLam(i1,i2,i3,i4)
        if (thisreLam > maxreLam) then
          maxreLam = thisreLam
        end if
      end do
      end do
      end do
      end do
      return
      end
c
c +++++++++++++
c
      subroutine appendRosenbluthCollision(
     *     rhs,f,velocities,IVx,IVy,d,Vth,N,massR,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     dx,
     *     dparams,
     *     iparams)
c
c.. compute and add the Rosenbluth collision operator to the RHS ...
      implicit none
c
c.. declarations of incoming variables 
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
c
      real dparams(*)
      integer iparams(*)
c
      real rhs(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real f(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real d(n1a:n1b, n2a:n2b, n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2), 1:3)
      real IVx(nd1a:nd1b, nd2a:nd2b)
      real IVy(nd1a:nd1b, nd2a:nd2b)
      real Vth(nd1a:nd1b, nd2a:nd2b)
      real N(nd1a:nd1b, nd2a:nd2b)
      real massR
c
      real dx(1:4)
c
c.. declarations of local variables
      integer solution_order, kernel_alg
c
      solution_order = iparams(1)
      kernel_alg = iparams(3)
c
      if( solution_order .eq. 4 ) then
        ! 4th order
        if (kernel_alg .eq. 0) then
          call appendRosenbluthCollision_4thMaplePrimitive(
     *       rhs,f,velocities,IVx,IVy,d,Vth,N,massR,
     *       nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *       n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *       dx,
     *       dparams,
     *       iparams)
        else if (kernel_alg .eq. 1) then
          call appendRosenbluthCollision_4thMapleLog(
     *       rhs,f,velocities,IVx,IVy,d,Vth,N,massR,
     *       nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *       n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *       dx,
     *       dparams,
     *       iparams)
        else if (kernel_alg .eq. 2) then
          call appendRosenbluthCollision_4thMapleOrig(
     *       rhs,f,velocities,IVx,IVy,d,Vth,N,massR,
     *       nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *       n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *       dx,
     *       dparams,
     *       iparams)
        else
          call appendRosenbluthCollision_4thMapleAsinh(
     *       rhs,f,velocities,IVx,IVy,d,Vth,N,massR,
     *       nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *       n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *       dx,
     *       dparams,
     *       iparams)
        end if
      else
        ! 6th order
        if (kernel_alg .eq. 0) then
          call appendRosenbluthCollision_6thMaplePrimitive(
     *       rhs,f,velocities,IVx,IVy,d,Vth,N,massR,
     *       nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *       n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *       dx,
     *       dparams,
     *       iparams)
        else if (kernel_alg .eq. 1) then
          call appendRosenbluthCollision_6thMapleLog(
     *       rhs,f,velocities,IVx,IVy,d,Vth,N,massR,
     *       nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *       n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *       dx,
     *       dparams,
     *       iparams)
        else if (kernel_alg .eq. 2) then
          call appendRosenbluthCollision_6thMapleOrig(
     *       rhs,f,velocities,IVx,IVy,d,Vth,N,massR,
     *       nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *       n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *       dx,
     *       dparams,
     *       iparams)
        else
          call appendRosenbluthCollision_6thMapleAsinh(
     *       rhs,f,velocities,IVx,IVy,d,Vth,N,massR,
     *       nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *       n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *       dx,
     *       dparams,
     *       iparams)
        end if
      end if
c
c
      return
      end
c
c +++++++++++++
c
      subroutine appendRosenbluthCollision_4thMaplePrimitive(
     *     rhs,fin,velocities,IVx,IVy,din,Vth,N,massR,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     dx, 
     *     dparams,
     *     iparams)
c
c.. compute and add the Rosenbluth collision operator to the RHS ...
      implicit none
c
c.. declarations of incoming variables 
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
c
      real dparams(*)
      integer iparams(*)
c
      real rhs(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real fin(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real din(n1a:n1b, n2a:n2b, n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2), 1:3)
      real IVx(nd1a:nd1b, nd2a:nd2b)
      real IVy(nd1a:nd1b, nd2a:nd2b)
      real Vth(nd1a:nd1b, nd2a:nd2b)
      real N(nd1a:nd1b, nd2a:nd2b)
      real massR
c     
      real dx(1:4)
c
c.. declarations of local variables
      real nu_base, nu
      integer i1, i2, i3, i4, i, j, modulate_nu
      real cg0
c
      real dvx, dvy, alpha
      real vx(-2:2,-2:2), vy(-2:2,-2:2)
      real f(-2:2,-2:2)
      real DD(-2:2,-2:2,1:3)
c
      real t1
      real t101
      real t106
      real t11
      real t12
      real t137
      real t139
      real t14
      real t148
      real t15
      real t158
      real t16
      real t177
      real t19
      real t196
      real t198
      real t199
      real t2
      real t200
      real t202
      real t203
      real t204
      real t205
      real t206
      real t208
      real t209
      real t21
      real t211
      real t212
      real t213
      real t216
      real t218
      real t220
      real t222
      real t226
      real t228
      real t23
      real t231
      real t233
      real t239
      real t240
      real t244
      real t245
      real t25
      real t253
      real t257
      real t263
      real t266
      real t279
      real t285
      real t287
      real t289
      real t29
      real t294
      real t297
      real t3
      real t302
      real t31
      real t333
      real t335
      real t34
      real t344
      real t355
      real t36
      real t365
      real t367
      real t368
      real t370
      real t375
      real t376
      real t378
      real t383
      real t384
      real t386
      real t39
      real t390
      real t392
      real t43
      real t436
      real t44
      real t48
      real t49
      real t5
      real t57
      real t6
      real t62
      real t67
      real t7
      real t70
      real t8
      real t83
      real t89
      real t9
      real t91
      real t93
      real t98
c
      modulate_nu = iparams(4)
      nu_base = dparams(1)
c
      dvx = dx(3)
      dvy = dx(4)
c      
      do i4 = n4a, n4b
      do i3 = n3a, n3b
      do i2 = n2a, n2b
      do i1 = n1a, n1b
        if (modulate_nu .eq. 1) then
          nu = nu_base*N(i1,i2)/Vth(i1,i2)
        else
          nu = nu_base
        end if
        alpha = 1.0 / (massR * (Vth(i1,i2)**2))
        do i = -2, 2
          do j = -2, 2
            vx(i,j) = velocities(i3+i,i4+j, 0) - IVx(i1,i2)
            vy(i,j) = velocities(i3+i,i4+j, 1) - IVy(i1,i2)
            f(i,j) = fin(i1,i2,i3+i,i4+j)
            DD(i,j,1) = din(i1,i2,i3+i,i4+j,1)
            DD(i,j,2) = din(i1,i2,i3+i,i4+j,2)
            DD(i,j,3) = din(i1,i2,i3+i,i4+j,3)
c            alpha(i,j) = 1.0 / (dparams(5) * (Vth(i1,i2)**2))
          end do
        end do
        t1 = vy(0,-2)
        t2 = 0.81E2 / 0.874E3 * t1
        t3 = vy(0,-1)
        t5 = vy(0,0)
        t6 = vy(0,2)
        t7 = 0.9E1 / 0.874E3 * t6
        t8 = -t2 + 0.392E3 / 0.437E3 * t3 + t5 - t7
        t9 = DD(0,-1,3)
        t11 = 0.9E1 / 0.874E3 * t1
        t12 = vy(0,1)
        t14 = 0.81E2 / 0.874E3 * t6
        t15 = t11 - t5 - 0.392E3 / 0.437E3 * t12 + t14
        t16 = DD(0,1,3)
        t19 = DD(0,0,3)
        t21 = DD(0,-2,3)
        t23 = DD(0,2,3)
        t25 = -0.81E2 / 0.874E3 * t21 - 0.9E1 / 0.874E3 * t23
        t29 = 0.9E1 / 0.874E3 * t21 + 0.81E2 / 0.874E3 * t23
        t31 = -t21 + t23
        t34 = t1 * t21
        t36 = t6 * t23
        t39 = f(0,0)
        t43 = 0.392E3 / 0.437E3 * t5
        t44 = 0.45E2 / 0.437E3 * t12
        t48 = 0.45E2 / 0.437E3 * t3
        t49 = -t7 + t11 - t48 + t44
        t57 = -0.9E1 / 0.874E3 * t31
        t62 = DD(0,-1,2)
        t67 = f(0,-1)
        t70 = t1 / 0.874E3
        t83 = DD(0,1,2)
        t89 = f(0,1)
        t91 = f(0,2)
        t93 = f(0,-2)
        t98 = -0.9E1 / 0.874E3 * t91 + 0.9E1 / 0.874E3 * t93
        t101 = 0.81E2 / 0.874E3 * t5
        t106 = t91 * (t5 - t6 / 0.9E1)
        t137 = t93 * t21
        t139 = t91 * t23
        t148 = DD(0,-2,2)
        t158 = DD(0,2,2)
        t177 = 0.4352E4 / 0.1311E4 * t19
        t196 = dvx ** 2
        t198 = vx(-2,0)
        t199 = 0.81E2 / 0.874E3 * t198
        t200 = vx(-1,0)
        t202 = vx(0,0)
        t203 = vx(2,0)
        t204 = 0.9E1 / 0.874E3 * t203
        t205 = t199 - 0.392E3 / 0.437E3 * t200 - t202 + t204
        t206 = DD(-1,0,1)
        t208 = 0.9E1 / 0.874E3 * t198
        t209 = vx(1,0)
        t211 = 0.81E2 / 0.874E3 * t203
        t212 = -t208 + t202 + 0.392E3 / 0.437E3 * t209 - t211
        t213 = DD(1,0,1)
        t216 = DD(0,0,1)
        t218 = DD(-2,0,1)
        t220 = DD(2,0,1)
        t222 = 0.81E2 / 0.874E3 * t218 + 0.9E1 / 0.874E3 * t220
        t226 = -0.9E1 / 0.874E3 * t218 - 0.81E2 / 0.874E3 * t220
        t228 = t218 - t220
        t231 = t218 * t198
        t233 = t220 * t203
        t239 = 0.392E3 / 0.437E3 * t202
        t240 = 0.45E2 / 0.437E3 * t209
        t244 = 0.45E2 / 0.437E3 * t200
        t245 = -t208 + t244 - t240 + t204
        t253 = -0.9E1 / 0.874E3 * t228
        t257 = DD(-1,0,2)
        t263 = f(-1,0)
        t266 = t198 / 0.874E3
        t279 = DD(1,0,2)
        t285 = f(1,0)
        t287 = f(-2,0)
        t289 = f(2,0)
        t294 = -0.9E1 / 0.874E3 * t287 + 0.9E1 / 0.874E3 * t289
        t297 = 0.81E2 / 0.874E3 * t202
        t302 = t289 * (t202 - t203 / 0.9E1)
        t333 = t287 * t218
        t335 = t289 * t220
        t344 = DD(-2,0,2)
        t355 = DD(2,0,2)
        t365 = 0.2560E4 / 0.1311E4 * f(-1,-1)
        t367 = 0.2560E4 / 0.1311E4 * f(-1,1)
        t368 = f(-1,-2)
        t370 = f(-1,2)
        t375 = 0.2560E4 / 0.1311E4 * f(1,-1)
        t376 = f(-2,-1)
        t378 = f(2,-1)
        t383 = 0.2560E4 / 0.1311E4 * f(1,1)
        t384 = f(-2,1)
        t386 = f(2,1)
        t390 = f(1,-2)
        t392 = f(1,2)
        t436 = dvy ** 2
        cg0 = -0.437E3 / 0.2048E4 / t196 / t436 * nu * (t196 * (dvy * (t
     #39 * (t9 * t8 + t16 * t15 + t19 * (-t2 + t3 - t12 + t14) + t3 * t2
     #5 + t12 * t29 + 0.81E2 / 0.874E3 * t5 * t31 + 0.9E1 / 0.874E3 * t3
     #4 - 0.9E1 / 0.874E3 * t36) + t67 * (t9 * (0.61E2 / 0.69E2 * t3 - 0
     #.307E3 / 0.2622E4 * t1 + t43 - t44 + t6 / 0.874E3) + t16 * t49 + t
     #19 * t8 + t3 * (t23 / 0.874E3 - 0.307E3 / 0.2622E4 * t21) + t12 * 
     #t57 + t5 * t25 + t36 / 0.874E3 - 0.37E2 / 0.2622E4 * t34 + 0.4096E
     #4 / 0.1311E4 * vx(0,-1) * t62) + t89 * (t9 * t49 + t16 * (-t70 + t
     #48 - t43 - 0.61E2 / 0.69E2 * t12 + 0.307E3 / 0.2622E4 * t6) + t19 
     #* t15 + t3 * t57 + t12 * (-t21 / 0.874E3 + 0.307E3 / 0.2622E4 * t2
     #3) + t5 * t29 - 0.4096E4 / 0.1311E4 * t83 * vx(0,1) - t34 / 0.874E
     #3 + 0.37E2 / 0.2622E4 * t36) + t9 * (t3 * (t91 / 0.874E3 - 0.307E3
     # / 0.2622E4 * t93) + t12 * t98 + t93 * (-0.37E2 / 0.2622E4 * t1 - 
     #t101) - 0.9E1 / 0.874E3 * t106) + t16 * (t3 * t98 + t12 * (0.307E3
     # / 0.2622E4 * t91 - t93 / 0.874E3) + t93 * (-t70 + 0.9E1 / 0.874E3
     # * t5) + 0.81E2 / 0.874E3 * (t5 + 0.37E2 / 0.243E3 * t6) * t91) + 
     #t19 * (t3 * (-0.81E2 / 0.874E3 * t93 - 0.9E1 / 0.874E3 * t91) + t1
     #2 * (0.9E1 / 0.874E3 * t93 + 0.81E2 / 0.874E3 * t91) + t93 * (t11 
     #- t101) + 0.81E2 / 0.874E3 * t106) + t3 * (-0.37E2 / 0.2622E4 * t1
     #37 + t139 / 0.874E3) + t12 * (-t137 / 0.874E3 + 0.37E2 / 0.2622E4 
     #* t139) + t93 * (-0.512E3 / 0.1311E4 * t148 * vx(0,-2) + 0.9E1 / 0
     #.874E3 * t5 * t21 - 0.67E2 / 0.2622E4 * t34) - 0.9E1 / 0.874E3 * t
     #91 * (t5 * t23 - 0.1024E4 / 0.27E2 * t158 * vx(0,2) - 0.67E2 / 0.2
     #7E2 * t36)) * alpha + t39 * (-0.128E3 / 0.437E3 * t23 - 0.128E3 / 
     #0.437E3 * t21 + 0.1280E4 / 0.437E3 * t9 + 0.2816E4 / 0.437E3 * t19
     # + 0.1280E4 / 0.437E3 * t16) + t67 * (0.256E3 / 0.1311E4 * t21 - 0
     #.4352E4 / 0.1311E4 * t9 - t177 + 0.256E3 / 0.1311E4 * t16) + t89 *
     # (0.256E3 / 0.1311E4 * t23 + 0.256E3 / 0.1311E4 * t9 - t177 - 0.43
     #52E4 / 0.1311E4 * t16) + 0.256E3 / 0.1311E4 * t93 * t9 + 0.256E3 /
     # 0.1311E4 * t91 * t16 + t19 * (0.128E3 / 0.1311E4 * t93 + 0.128E3 
     #/ 0.1311E4 * t91) + 0.128E3 / 0.1311E4 * t137 + 0.128E3 / 0.1311E4
     # * t139) - dvx * (dvy * (t39 * (t206 * t205 + t213 * t212 + t216 *
     # (t199 - t200 + t209 - t211) + t200 * t222 + t209 * t226 + 0.81E2 
     #/ 0.874E3 * t202 * t228 - 0.9E1 / 0.874E3 * t231 + 0.9E1 / 0.874E3
     # * t233) + t263 * (t206 * (0.307E3 / 0.2622E4 * t198 - 0.61E2 / 0.
     #69E2 * t200 - t239 + t240 - t203 / 0.874E3) + t213 * t245 + t216 *
     # t205 + t200 * (0.307E3 / 0.2622E4 * t218 - t220 / 0.874E3) + t209
     # * t253 + t202 * t222 - 0.4096E4 / 0.1311E4 * t257 * vy(-1,0) + 0.
     #37E2 / 0.2622E4 * t231 - t233 / 0.874E3) + t285 * (t206 * t245 + t
     #213 * (t266 - t244 + t239 + 0.61E2 / 0.69E2 * t209 - 0.307E3 / 0.2
     #622E4 * t203) + t216 * t212 + t200 * t253 + t209 * (t218 / 0.874E3
     # - 0.307E3 / 0.2622E4 * t220) + t202 * t226 + 0.4096E4 / 0.1311E4 
     #* t279 * vy(1,0) + t231 / 0.874E3 - 0.37E2 / 0.2622E4 * t233) + t2
     #06 * (t200 * (0.307E3 / 0.2622E4 * t287 - t289 / 0.874E3) + t209 *
     # t294 + t287 * (0.37E2 / 0.2622E4 * t198 + t297) + 0.9E1 / 0.874E3
     # * t302) + t213 * (t200 * t294 + t209 * (t287 / 0.874E3 - 0.307E3 
     #/ 0.2622E4 * t289) + t287 * (t266 - 0.9E1 / 0.874E3 * t202) - 0.81
     #E2 / 0.874E3 * (t202 + 0.37E2 / 0.243E3 * t203) * t289) + t216 * (
     #t200 * (0.9E1 / 0.874E3 * t289 + 0.81E2 / 0.874E3 * t287) + t209 *
     # (-0.81E2 / 0.874E3 * t289 - 0.9E1 / 0.874E3 * t287) + t287 * (-t2
     #08 + t297) - 0.81E2 / 0.874E3 * t302) + t200 * (0.37E2 / 0.2622E4 
     #* t333 - t335 / 0.874E3) + t209 * (t333 / 0.874E3 - 0.37E2 / 0.262
     #2E4 * t335) + t287 * (0.512E3 / 0.1311E4 * t344 * vy(-2,0) + 0.67E
     #2 / 0.2622E4 * t231 - 0.9E1 / 0.874E3 * t218 * t202) + 0.9E1 / 0.8
     #74E3 * t289 * (t220 * t202 - 0.67E2 / 0.27E2 * t233 - 0.1024E4 / 0
     #.27E2 * t355 * vy(2,0))) * alpha + t257 * (t365 - t367 - 0.256E3 /
     # 0.1311E4 * t368 + 0.256E3 / 0.1311E4 * t370) + t62 * (t365 - t375
     # - 0.256E3 / 0.1311E4 * t376 + 0.256E3 / 0.1311E4 * t378) + t83 * 
     #(-t367 + t383 + 0.256E3 / 0.1311E4 * t384 - 0.256E3 / 0.1311E4 * t
     #386) + t279 * (0.256E3 / 0.1311E4 * t390 - 0.256E3 / 0.1311E4 * t3
     #92 - t375 + t383) + t344 * (-0.256E3 / 0.1311E4 * t376 + 0.256E3 /
     # 0.1311E4 * t384) + t148 * (-0.256E3 / 0.1311E4 * t368 + 0.256E3 /
     # 0.1311E4 * t390) + t158 * (0.256E3 / 0.1311E4 * t370 - 0.256E3 / 
     #0.1311E4 * t392) + 0.256E3 / 0.1311E4 * (t378 - t386) * t355) * dv
     #y - 0.4352E4 / 0.1311E4 * t436 * (t39 * (-0.33E2 / 0.17E2 * t216 +
     # 0.3E1 / 0.34E2 * t218 - 0.15E2 / 0.17E2 * t206 - 0.15E2 / 0.17E2 
     #* t213 + 0.3E1 / 0.34E2 * t220) + t263 * (t216 - t218 / 0.17E2 + t
     #206 - t213 / 0.17E2) + t285 * (t216 - t206 / 0.17E2 + t213 - t220 
     #/ 0.17E2) - t287 * t206 / 0.17E2 - t289 * t213 / 0.17E2 + t216 * (
     #-t289 / 0.34E2 - t287 / 0.34E2) - t335 / 0.34E2 - t333 / 0.34E2))
        rhs(i1,i2,i3,i4) = rhs(i1,i2,i3,i4)+cg0
      end do
      end do
      end do
      end do
      return
      end
c
c +++++++++++++
c
      subroutine appendRosenbluthCollision_4thMapleLog(
     *     rhs,fin,velocities,IVx,IVy,din,Vth,N,massR,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     dx, 
     *     dparams,
     *     iparams)
c
c.. compute and add the Rosenbluth collision operator to the RHS ...
      implicit none
c
c.. declarations of incoming variables 
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
c
      real dparams(*)
      integer iparams(*)
c
      real rhs(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real fin(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real din(n1a:n1b, n2a:n2b, n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2), 1:3)
      real IVx(nd1a:nd1b, nd2a:nd2b)
      real IVy(nd1a:nd1b, nd2a:nd2b)
      real Vth(nd1a:nd1b, nd2a:nd2b)
      real N(nd1a:nd1b, nd2a:nd2b)
      real massR
c     
      real dx(1:4)
c
c.. declarations of local variables
      real nu_base, nu
      integer i1, i2, i3, i4, i, j, modulate_nu
      real cg1
c
      real dvx, dvy, alpha
      real vx(-2:2,-2:2), vy(-2:2,-2:2)
      real f(-2:2,-2:2)
      real D(-2:2,-2:2,1:3)
c
      real t1
      real t10
      real t100
      real t107
      real t108
      real t11
      real t111
      real t112
      real t113
      real t114
      real t115
      real t118
      real t119
      real t123
      real t124
      real t127
      real t128
      real t13
      real t132
      real t133
      real t134
      real t135
      real t136
      real t139
      real t14
      real t140
      real t146
      real t147
      real t15
      real t152
      real t153
      real t156
      real t157
      real t160
      real t161
      real t165
      real t166
      real t169
      real t17
      real t170
      real t176
      real t177
      real t18
      real t182
      real t183
      real t188
      real t189
      real t19
      real t194
      real t197
      real t198
      real t2
      real t204
      real t205
      real t210
      real t212
      real t214
      real t217
      real t22
      real t220
      real t223
      real t226
      real t229
      real t232
      real t235
      real t24
      real t25
      real t26
      real t27
      real t278
      real t28
      real t280
      real t281
      real t283
      real t284
      real t286
      real t288
      real t29
      real t290
      real t291
      real t293
      real t294
      real t296
      real t3
      real t301
      real t303
      real t304
      real t306
      real t308
      real t31
      real t310
      real t311
      real t313
      real t318
      real t32
      real t320
      real t321
      real t323
      real t325
      real t327
      real t328
      real t330
      real t335
      real t337
      real t339
      real t34
      real t341
      real t35
      real t36
      real t361
      real t364
      real t367
      real t370
      real t373
      real t376
      real t379
      real t38
      real t382
      real t39
      real t40
      real t42
      real t427
      real t46
      real t47
      real t5
      real t53
      real t54
      real t6
      real t61
      real t62
      real t69
      real t7
      real t70
      real t77
      real t78
      real t84
      real t85
      real t88
      real t9
      real t92
      real t93
      real t99
c
      modulate_nu = iparams(4)
      nu_base = dparams(1)
c
      dvx = dx(3)
      dvy = dx(4)
c      
      do i4 = n4a, n4b
      do i3 = n3a, n3b
      do i2 = n2a, n2b
      do i1 = n1a, n1b
        if (modulate_nu .eq. 1) then
          nu = nu_base*N(i1,i2)/Vth(i1,i2)
        else
          nu = nu_base
        end if
        alpha = 1.0 / (massR * (Vth(i1,i2)**2))
        do i = -2, 2
          do j = -2, 2
            vx(i,j) = velocities(i3+i,i4+j, 0) - IVx(i1,i2)
            vy(i,j) = velocities(i3+i,i4+j, 1) - IVy(i1,i2)
            f(i,j) = fin(i1,i2,i3+i,i4+j)
            D(i,j,1) = din(i1,i2,i3+i,i4+j,1)
            D(i,j,2) = din(i1,i2,i3+i,i4+j,2)
            D(i,j,3) = din(i1,i2,i3+i,i4+j,3)
c            alpha(i,j) = 1.0 / (dparams(5) * (Vth(i1,i2)**2))
          end do
        end do
        t1 = D(0,-2,3)
        t2 = f(0,-2)
        t3 = t2 * t1
        t5 = D(0,-1,3)
        t6 = f(0,-1)
        t7 = t6 * t5
        t9 = f(0,0)
        t10 = D(0,0,3)
        t11 = t10 * t9
        t13 = f(0,1)
        t14 = D(0,1,3)
        t15 = t14 * t13
        t17 = D(0,2,3)
        t18 = f(0,2)
        t19 = t18 * t17
        t22 = dvx ** 2
        t24 = f(1,0)
        t25 = D(1,0,1)
        t26 = t25 * t24
        t27 = f(2,0)
        t28 = D(2,0,1)
        t29 = t28 * t27
        t31 = D(0,0,1)
        t32 = t31 * t9
        t34 = D(-2,0,1)
        t35 = f(-2,0)
        t36 = t35 * t34
        t38 = f(-1,0)
        t39 = D(-1,0,1)
        t40 = t39 * t38
        t42 = dvy ** 2
        t46 = t9 ** 2
        t47 = log(t46)
        t53 = t38 ** 2
        t54 = log(t53)
        t61 = t6 ** 2
        t62 = log(t61)
        t69 = t13 ** 2
        t70 = log(t69)
        t77 = t24 ** 2
        t78 = log(t77)
        t84 = t35 ** 2
        t85 = log(t84)
        t88 = t11 / 0.2E1
        t92 = t2 ** 2
        t93 = log(t92)
        t99 = t18 ** 2
        t100 = log(t99)
        t107 = t27 ** 2
        t108 = log(t107)
        t111 = dvy * dvx
        t112 = D(-2,0,2)
        t113 = t35 * t112
        t114 = D(0,-1,2)
        t115 = t6 * t114
        t118 = f(-2,-1) ** 2
        t119 = log(t118)
        t123 = D(0,1,2)
        t124 = t13 * t123
        t127 = f(-2,1) ** 2
        t128 = log(t127)
        t132 = t47 * (t22 * (-0.3E1 / 0.20E2 * t3 + 0.3E1 / 0.2E1 * t7 +
     # 0.33E2 / 0.10E2 * t11 + 0.3E1 / 0.2E1 * t15 - 0.3E1 / 0.20E2 * t1
     #9) + 0.3E1 / 0.2E1 * t42 * (t26 - t29 / 0.10E2 + 0.11E2 / 0.5E1 * 
     #t32 - t36 / 0.10E2 + t40)) + t54 * (t36 - 0.17E2 * t40 - 0.17E2 * 
     #t32 + t26) * t42 / 0.10E2 - 0.17E2 / 0.10E2 * t62 * t22 * (t7 - t1
     #5 / 0.17E2 + t11 - t3 / 0.17E2) + t70 * (t7 - 0.17E2 * t11 - 0.17E
     #2 * t15 + t19) * t22 / 0.10E2 - 0.17E2 / 0.10E2 * t78 * t42 * (t26
     # - t29 / 0.17E2 + t32 - t40 / 0.17E2) + t85 * (t36 + 0.2E1 * t40 +
     # t32) * t42 / 0.20E2 + t93 * t22 * (t7 + t88 + t3 / 0.2E1) / 0.10E
     #2 + t100 * t22 * (t15 + t19 / 0.2E1 + t88) / 0.10E2 + t108 * (t26 
     #+ t29 / 0.2E1 + t32 / 0.2E1) * t42 / 0.10E2 + t119 * (t113 + t115)
     # * t111 / 0.10E2 - t128 * (t113 + t124) * t111 / 0.10E2
        t133 = D(-1,0,2)
        t134 = t38 * t133
        t135 = D(0,-2,2)
        t136 = t2 * t135
        t139 = f(-1,-2) ** 2
        t140 = log(t139)
        t146 = f(-1,-1) ** 2
        t147 = log(t146)
        t152 = f(-1,1) ** 2
        t153 = log(t152)
        t156 = D(0,2,2)
        t157 = t18 * t156
        t160 = f(-1,2) ** 2
        t161 = log(t160)
        t165 = D(1,0,2)
        t166 = t24 * t165
        t169 = f(1,-2) ** 2
        t170 = log(t169)
        t176 = f(1,-1) ** 2
        t177 = log(t176)
        t182 = f(1,1) ** 2
        t183 = log(t182)
        t188 = f(1,2) ** 2
        t189 = log(t188)
        t194 = t27 * D(2,0,2)
        t197 = f(2,-1) ** 2
        t198 = log(t197)
        t204 = f(2,1) ** 2
        t205 = log(t204)
        t210 = vx(0,0) ** 2
        t212 = vy(0,0) ** 2
        t214 = vx(0,-2) ** 2
        t217 = vx(0,-1) ** 2
        t220 = vx(0,1) ** 2
        t223 = vx(0,2) ** 2
        t226 = vy(0,-2) ** 2
        t229 = vy(0,-1) ** 2
        t232 = vy(0,1) ** 2
        t235 = vy(0,2) ** 2
        t278 = vx(-1,-2) ** 2
        t280 = vx(-1,-1) ** 2
        t281 = 0.10E2 * t280
        t283 = vx(-1,1) ** 2
        t284 = 0.10E2 * t283
        t286 = vx(-1,2) ** 2
        t288 = vy(-1,-2) ** 2
        t290 = vy(-1,-1) ** 2
        t291 = 0.10E2 * t290
        t293 = vy(-1,1) ** 2
        t294 = 0.10E2 * t293
        t296 = vy(-1,2) ** 2
        t301 = vx(-2,-1) ** 2
        t303 = vx(1,-1) ** 2
        t304 = 0.10E2 * t303
        t306 = vx(2,-1) ** 2
        t308 = vy(-2,-1) ** 2
        t310 = vy(1,-1) ** 2
        t311 = 0.10E2 * t310
        t313 = vy(2,-1) ** 2
        t318 = vx(-2,1) ** 2
        t320 = vx(1,1) ** 2
        t321 = 0.10E2 * t320
        t323 = vx(2,1) ** 2
        t325 = vy(-2,1) ** 2
        t327 = vy(1,1) ** 2
        t328 = 0.10E2 * t327
        t330 = vy(2,1) ** 2
        t335 = vx(1,-2) ** 2
        t337 = vx(1,2) ** 2
        t339 = vy(1,-2) ** 2
        t341 = vy(1,2) ** 2
        t361 = vx(-2,0) ** 2
        t364 = vx(-1,0) ** 2
        t367 = vx(1,0) ** 2
        t370 = vx(2,0) ** 2
        t373 = vy(-2,0) ** 2
        t376 = vy(-1,0) ** 2
        t379 = vy(1,0) ** 2
        t382 = vy(2,0) ** 2
        t427 = t140 * (t134 + t136) * t111 / 0.10E2 - t147 * (t134 + t11
     #5) * t111 + t153 * (t134 + t124) * t111 - t161 * (t134 + t157) * t
     #111 / 0.10E2 - t170 * (t136 + t166) * t111 / 0.10E2 + t177 * (t115
     # + t166) * t111 - t183 * (t124 + t166) * t111 + t189 * (t157 + t16
     #6) * t111 / 0.10E2 - t198 * (t115 + t194) * t111 / 0.10E2 + t205 *
     # (t124 + t194) * t111 / 0.10E2 + 0.3E1 / 0.2E1 * (t22 * (0.11E2 / 
     #0.5E1 * t9 * (t210 + t212 + t214 / 0.66E2 - 0.17E2 / 0.33E2 * t217
     # - 0.17E2 / 0.33E2 * t220 + t223 / 0.66E2 + t226 / 0.66E2 - 0.17E2
     # / 0.33E2 * t229 - 0.17E2 / 0.33E2 * t232 + t235 / 0.66E2) * t10 +
     # t6 * (-0.17E2 / 0.15E2 * t217 + t220 / 0.15E2 - 0.17E2 / 0.15E2 *
     # t229 + t212 + t232 / 0.15E2 + t214 / 0.15E2 + t226 / 0.15E2 + t21
     #0) * t5 + t13 * (t217 / 0.15E2 + t210 - 0.17E2 / 0.15E2 * t220 + t
     #229 / 0.15E2 + t212 - 0.17E2 / 0.15E2 * t232 + t223 / 0.15E2 + t23
     #5 / 0.15E2) * t14 - t2 * t1 * (t210 + t212 - t214 / 0.3E1 - 0.2E1 
     #/ 0.3E1 * t217 - t226 / 0.3E1 - 0.2E1 / 0.3E1 * t229) / 0.10E2 - t
     #17 * t18 * (t210 + t212 - 0.2E1 / 0.3E1 * t220 - t223 / 0.3E1 - 0.
     #2E1 / 0.3E1 * t232 - t235 / 0.3E1) / 0.10E2) + dvx * (t38 * (t278 
     #- t281 + t284 - t286 + t288 - t291 + t294 - t296) * t133 + t6 * (t
     #301 - t281 + t304 - t306 + t308 - t291 + t311 - t313) * t114 - t13
     # * (t318 - t284 + t321 - t323 + t325 - t294 + t328 - t330) * t123 
     #- t24 * (t335 - t304 + t321 - t337 + t339 - t311 + t328 - t341) * 
     #t165 + t35 * (t301 - t318 + t308 - t325) * t112 + t2 * (t278 - t33
     #5 + t288 - t339) * t135 - t18 * (t286 - t337 + t296 - t341) * t156
     # - (t306 - t323 + t313 - t330) * t194) * dvy / 0.15E2 + (0.11E2 / 
     #0.5E1 * t9 * (t210 + t212 + t361 / 0.66E2 - 0.17E2 / 0.33E2 * t364
     # - 0.17E2 / 0.33E2 * t367 + t370 / 0.66E2 + t373 / 0.66E2 - 0.17E2
     # / 0.33E2 * t376 - 0.17E2 / 0.33E2 * t379 + t382 / 0.66E2) * t31 +
     # t38 * (t212 + t210 - 0.17E2 / 0.15E2 * t364 + t367 / 0.15E2 - 0.1
     #7E2 / 0.15E2 * t376 + t379 / 0.15E2 + t361 / 0.15E2 + t373 / 0.15E
     #2) * t39 + t24 * (t210 + t212 + t364 / 0.15E2 - 0.17E2 / 0.15E2 * 
     #t367 + t376 / 0.15E2 - 0.17E2 / 0.15E2 * t379 + t370 / 0.15E2 + t3
     #82 / 0.15E2) * t25 - t35 * t34 * (t210 + t212 - t361 / 0.3E1 - 0.2
     #E1 / 0.3E1 * t364 - t373 / 0.3E1 - 0.2E1 / 0.3E1 * t376) / 0.10E2 
     #- t28 * (t210 + t212 - 0.2E1 / 0.3E1 * t367 - t370 / 0.3E1 - 0.2E1
     # / 0.3E1 * t379 - t382 / 0.3E1) * t27 / 0.10E2) * t42) * alpha
        cg1 = -0.5E1 / 0.24E2 / t22 / t42 * (t132 + t427) * nu
ccccc
ccccc
        rhs(i1,i2,i3,i4) = rhs(i1,i2,i3,i4)+cg1
      end do
      end do
      end do
      end do
c
      return
      end
c
c +++++++++++++
c
      subroutine appendRosenbluthCollision_4thMapleOrig(
     *     rhs,fin,velocities,IVx,IVy,din,Vth,N,massR,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     dx, 
     *     dparams,
     *     iparams)
c
c.. compute and add the Rosenbluth collision operator to the RHS ...
      implicit none
c
c.. declarations of incoming variables 
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
c
      real dparams(*)
      integer iparams(*)
c
      real rhs(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real fin(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real din(n1a:n1b, n2a:n2b, n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2), 1:3)
      real IVx(nd1a:nd1b, nd2a:nd2b)
      real IVy(nd1a:nd1b, nd2a:nd2b)
      real Vth(nd1a:nd1b, nd2a:nd2b)
      real N(nd1a:nd1b, nd2a:nd2b)
      real massR
c     
      real dx(1:4)
c
c.. declarations of local variables
      real nu_base, nu
      integer i1, i2, i3, i4, i, j, modulate_nu
      real cg0
c
      real dvx, dvy, alpha
      real vx(-2:2,-2:2), vy(-2:2,-2:2)
      real f(-2:2,-2:2)
      real DD(-2:2,-2:2,1:3)
c
      real t1
      real t100
      real t105
      real t108
      real t110
      real t113
      real t114
      real t115
      real t117
      real t12
      real t120
      real t123
      real t127
      real t13
      real t130
      real t132
      real t138
      real t14
      real t141
      real t145
      real t148
      real t15
      real t150
      real t154
      real t156
      real t159
      real t16
      real t161
      real t164
      real t165
      real t166
      real t168
      real t17
      real t174
      real t177
      real t181
      real t184
      real t186
      real t191
      real t194
      real t196
      real t199
      real t2
      real t20
      real t200
      real t201
      real t203
      real t209
      real t212
      real t216
      real t219
      real t221
      real t226
      real t229
      real t23
      real t231
      real t234
      real t235
      real t236
      real t238
      real t241
      real t244
      real t248
      real t251
      real t253
      real t259
      real t262
      real t266
      real t269
      real t27
      real t271
      real t276
      real t279
      real t281
      real t284
      real t285
      real t286
      real t288
      real t291
      real t294
      real t298
      real t30
      real t301
      real t303
      real t309
      real t312
      real t316
      real t319
      real t32
      real t321
      real t328
      real t338
      real t340
      real t341
      real t349
      real t350
      real t358
      real t361
      real t363
      real t366
      real t367
      real t368
      real t370
      real t377
      real t378
      real t379
      real t38
      real t380
      real t381
      real t385
      real t386
      real t390
      real t391
      real t395
      real t397
      real t398
      real t399
      real t4
      real t403
      real t407
      real t408
      real t41
      real t412
      real t413
      real t417
      real t419
      real t420
      real t424
      real t425
      real t429
      real t430
      real t445
      real t448
      real t45
      real t451
      real t454
      real t461
      real t462
      real t465
      real t466
      real t467
      real t471
      real t48
      real t481
      real t485
      real t486
      real t492
      real t495
      real t496
      real t50
      real t503
      real t518
      real t520
      real t533
      real t534
      real t537
      real t539
      real t541
      real t55
      real t550
      real t551
      real t58
      real t60
      real t600
      real t63
      real t64
      real t65
      real t66
      real t67
      real t7
      real t70
      real t73
      real t77
      real t80
      real t82
      real t88
      real t9
      real t91
      real t95
      real t98
c
      modulate_nu = iparams(4)
      nu_base = dparams(1)
c
      dvx = dx(3)
      dvy = dx(4)
c      
      do i4 = n4a, n4b
      do i3 = n3a, n3b
      do i2 = n2a, n2b
      do i1 = n1a, n1b
        if (modulate_nu .eq. 1) then
          nu = nu_base*N(i1,i2)/Vth(i1,i2)
        else
          nu = nu_base
        end if
        alpha = 1.0 / (massR * (Vth(i1,i2)**2))
        do i = -2, 2
          do j = -2, 2
            vx(i,j) = velocities(i3+i,i4+j, 0) - IVx(i1,i2)
            vy(i,j) = velocities(i3+i,i4+j, 1) - IVy(i1,i2)
            f(i,j) = fin(i1,i2,i3+i,i4+j)
            DD(i,j,1) = din(i1,i2,i3+i,i4+j,1)
            DD(i,j,2) = din(i1,i2,i3+i,i4+j,2)
            DD(i,j,3) = din(i1,i2,i3+i,i4+j,3)
c            alpha(i,j) = 1.0 / (dparams(5) * (Vth(i1,i2)**2))
          end do
        end do
        t1 = 0.1E1 / dvy
        t2 = nu * t1
        t4 = vx(0,1) ** 2
        t7 = exp(-t4 * alpha / 0.2E1)
        t9 = vy(0,1) ** 2
        t12 = exp(-t9 * alpha / 0.2E1)
        t13 = t12 * t7
        t14 = t13 * t2
        t15 = DD(0,1,2)
        t16 = 0.1E1 / dvx
        t17 = t16 * t15
        t20 = vx(1,1) ** 2
        t23 = exp(-t20 * alpha / 0.2E1)
        t27 = vy(1,1) ** 2
        t30 = exp(-t27 * alpha / 0.2E1)
        t32 = 0.1E1 / t30 / t23 * f(1,1)
        t38 = vx(-1,1) ** 2
        t41 = exp(-t38 * alpha / 0.2E1)
        t45 = vy(-1,1) ** 2
        t48 = exp(-t45 * alpha / 0.2E1)
        t50 = 0.1E1 / t48 / t41 * f(-1,1)
        t55 = vx(0,-1) ** 2
        t58 = exp(-t55 * alpha / 0.2E1)
        t60 = vy(0,-1) ** 2
        t63 = exp(-t60 * alpha / 0.2E1)
        t64 = t63 * t58
        t65 = t64 * t2
        t66 = DD(0,-1,2)
        t67 = t16 * t66
        t70 = vx(1,-1) ** 2
        t73 = exp(-t70 * alpha / 0.2E1)
        t77 = vy(1,-1) ** 2
        t80 = exp(-t77 * alpha / 0.2E1)
        t82 = 0.1E1 / t80 / t73 * f(1,-1)
        t88 = vx(-1,-1) ** 2
        t91 = exp(-t88 * alpha / 0.2E1)
        t95 = vy(-1,-1) ** 2
        t98 = exp(-t95 * alpha / 0.2E1)
        t100 = 0.1E1 / t98 / t91 * f(-1,-1)
        t105 = vx(0,2) ** 2
        t108 = exp(-t105 * alpha / 0.2E1)
        t110 = vy(0,2) ** 2
        t113 = exp(-t110 * alpha / 0.2E1)
        t114 = t113 * t108
        t115 = t114 * t2
        t117 = t16 * DD(0,2,2)
        t120 = vx(1,2) ** 2
        t123 = exp(-t120 * alpha / 0.2E1)
        t127 = vy(1,2) ** 2
        t130 = exp(-t127 * alpha / 0.2E1)
        t132 = 0.1E1 / t130 / t123 * f(1,2)
        t138 = vx(-1,2) ** 2
        t141 = exp(-t138 * alpha / 0.2E1)
        t145 = vy(-1,2) ** 2
        t148 = exp(-t145 * alpha / 0.2E1)
        t150 = 0.1E1 / t148 / t141 * f(-1,2)
        t154 = nu * t16
        t156 = vx(1,0) ** 2
        t159 = exp(-t156 * alpha / 0.2E1)
        t161 = vy(1,0) ** 2
        t164 = exp(-t161 * alpha / 0.2E1)
        t165 = t164 * t159
        t166 = t165 * t154
        t168 = t1 * DD(1,0,2)
        t174 = vx(1,-2) ** 2
        t177 = exp(-t174 * alpha / 0.2E1)
        t181 = vy(1,-2) ** 2
        t184 = exp(-t181 * alpha / 0.2E1)
        t186 = 0.1E1 / t184 / t177 * f(1,-2)
        t191 = vx(-1,0) ** 2
        t194 = exp(-t191 * alpha / 0.2E1)
        t196 = vy(-1,0) ** 2
        t199 = exp(-t196 * alpha / 0.2E1)
        t200 = t199 * t194
        t201 = t200 * t154
        t203 = t1 * DD(-1,0,2)
        t209 = vx(-1,-2) ** 2
        t212 = exp(-t209 * alpha / 0.2E1)
        t216 = vy(-1,-2) ** 2
        t219 = exp(-t216 * alpha / 0.2E1)
        t221 = 0.1E1 / t219 / t212 * f(-1,-2)
        t226 = vx(-2,0) ** 2
        t229 = exp(-t226 * alpha / 0.2E1)
        t231 = vy(-2,0) ** 2
        t234 = exp(-t231 * alpha / 0.2E1)
        t235 = t234 * t229
        t236 = t235 * t154
        t238 = t1 * DD(-2,0,2)
        t241 = vx(-2,1) ** 2
        t244 = exp(-t241 * alpha / 0.2E1)
        t248 = vy(-2,1) ** 2
        t251 = exp(-t248 * alpha / 0.2E1)
        t253 = 0.1E1 / t251 / t244 * f(-2,1)
        t259 = vx(-2,-1) ** 2
        t262 = exp(-t259 * alpha / 0.2E1)
        t266 = vy(-2,-1) ** 2
        t269 = exp(-t266 * alpha / 0.2E1)
        t271 = 0.1E1 / t269 / t262 * f(-2,-1)
        t276 = vx(2,0) ** 2
        t279 = exp(-t276 * alpha / 0.2E1)
        t281 = vy(2,0) ** 2
        t284 = exp(-t281 * alpha / 0.2E1)
        t285 = t284 * t279
        t286 = t285 * t154
        t288 = t1 * DD(2,0,2)
        t291 = vx(2,1) ** 2
        t294 = exp(-t291 * alpha / 0.2E1)
        t298 = vy(2,1) ** 2
        t301 = exp(-t298 * alpha / 0.2E1)
        t303 = 0.1E1 / t301 / t294 * f(2,1)
        t309 = vx(2,-1) ** 2
        t312 = exp(-t309 * alpha / 0.2E1)
        t316 = vy(2,-1) ** 2
        t319 = exp(-t316 * alpha / 0.2E1)
        t321 = 0.1E1 / t319 / t312 * f(2,-1)
        t328 = 0.5E1 / 0.12E2 * t32 * t17 * t14 - 0.5E1 / 0.12E2 * t50 *
     # t17 * t14 - 0.5E1 / 0.12E2 * t82 * t67 * t65 + 0.5E1 / 0.12E2 * t
     #100 * t67 * t65 - t132 * t117 * t115 / 0.24E2 + t150 * t117 * t115
     # / 0.24E2 - t132 * t168 * t166 / 0.24E2 + t186 * t168 * t166 / 0.2
     #4E2 + t150 * t203 * t201 / 0.24E2 - t221 * t203 * t201 / 0.24E2 + 
     #t253 * t238 * t236 / 0.24E2 - t271 * t238 * t236 / 0.24E2 - t303 *
     # t288 * t286 / 0.24E2 + t321 * t288 * t286 / 0.24E2 + 0.5E1 / 0.12
     #E2 * t32 * t168 * t166
        t338 = t1 * t16
        t340 = t7 * nu * t338
        t341 = t15 * t12
        t349 = t58 * nu * t338
        t350 = t66 * t63
        t358 = vx(0,-2) ** 2
        t361 = exp(-t358 * alpha / 0.2E1)
        t363 = vy(0,-2) ** 2
        t366 = exp(-t363 * alpha / 0.2E1)
        t367 = t366 * t361
        t368 = t367 * t2
        t370 = t16 * DD(0,-2,2)
        t377 = dvx ** 2
        t378 = 0.1E1 / t377
        t379 = nu * t378
        t380 = DD(1,0,1)
        t381 = f(1,0)
        t385 = DD(0,0,1)
        t386 = f(0,0)
        t390 = DD(-1,0,1)
        t391 = f(-1,0)
        t395 = dvy ** 2
        t397 = nu / t395
        t398 = DD(0,1,3)
        t399 = f(0,1)
        t403 = DD(0,0,3)
        t407 = DD(0,-1,3)
        t408 = f(0,-1)
        t412 = DD(0,2,3)
        t413 = f(0,2)
        t417 = -0.5E1 / 0.12E2 * t82 * t168 * t166 - 0.5E1 / 0.12E2 * t5
     #0 * t203 * t201 + 0.5E1 / 0.12E2 * t100 * t203 * t201 - t303 * t34
     #1 * t340 / 0.24E2 + t253 * t341 * t340 / 0.24E2 + t321 * t350 * t3
     #49 / 0.24E2 - t271 * t350 * t349 / 0.24E2 + t186 * t370 * t368 / 0
     #.24E2 - t221 * t370 * t368 / 0.24E2 + 0.17E2 / 0.24E2 * t381 * t38
     #0 * t379 - 0.11E2 / 0.8E1 * t386 * t385 * t379 + 0.17E2 / 0.24E2 *
     # t391 * t390 * t379 + 0.17E2 / 0.24E2 * t399 * t398 * t397 - 0.11E
     #2 / 0.8E1 * t386 * t403 * t397 + 0.17E2 / 0.24E2 * t408 * t407 * t
     #397 - t413 * t412 * t397 / 0.48E2
        t419 = DD(0,-2,3)
        t420 = f(0,-2)
        t424 = DD(2,0,1)
        t425 = f(2,0)
        t429 = DD(-2,0,1)
        t430 = f(-2,0)
        t445 = vx(0,0) ** 2
        t448 = exp(-t445 * alpha / 0.2E1)
        t451 = vy(0,0) ** 2
        t454 = exp(-t451 * alpha / 0.2E1)
        t466 = t385 * t454 * t448 * nu
        t471 = 0.1E1 / t234 / t229 * t430 * t378
        t496 = t454 * t448
        t461 = 0.1E1 / t63 / t58 * t408
        t462 = t398 * t13 * t397
        t465 = 0.1E1 / t454 / t448 * t386
        t467 = t419 * t367 * t397
        t481 = t380 * t165 * t379
        t485 = 0.1E1 / t164 / t159 * t381
        t486 = t424 * t285 * t379
        t492 = t385 * t496 * t379
        t495 = t390 * t200 * t379
        t503 = 0.1E1 / t199 / t194 * t391
        t518 = -t420 * t419 * t397 / 0.48E2 - t425 * t424 * t379 / 0.48E
     #2 - t430 * t429 * t379 / 0.48E2 - t462 * t461 / 0.24E2 + t467 * t4
     #65 / 0.16E2 - t467 * t461 / 0.24E2 - t471 * t466 / 0.48E2 - t471 *
     # t390 * t199 * t194 * nu / 0.24E2 - 0.5E1 / 0.8E1 * t481 * t465 - 
     #t486 * t485 / 0.24E2 + t486 * t465 / 0.16E2 + 0.17E2 / 0.24E2 * t4
     #92 * t485 - t495 * t485 / 0.24E2 - 0.5E1 / 0.8E1 * t495 * t465 + 0
     #.17E2 / 0.24E2 * t492 * t503
        t539 = 0.1E1 / t284 / t279 * t425 * t378
        t520 = t429 * t235 * t379
        t533 = 0.1E1 / t113 / t108 * t413
        t534 = t403 * t496 * t397
        t537 = 0.1E1 / t366 / t361 * t420
        t541 = t407 * t64 * t397
        t550 = 0.1E1 / t12 / t7 * t399
        t551 = t412 * t114 * t397
        t600 = -t481 * t503 / 0.24E2 + t520 * t465 / 0.16E2 - t520 * t50
     #3 / 0.24E2 - t539 * t380 * t164 * t159 * nu / 0.24E2 - t539 * t466
     # / 0.48E2 - t534 * t533 / 0.48E2 - t534 * t537 / 0.48E2 - t541 * t
     #537 / 0.24E2 - t462 * t533 / 0.24E2 - 0.5E1 / 0.8E1 * t462 * t465 
     #- t551 * t550 / 0.24E2 + t551 * t465 / 0.16E2 + 0.17E2 / 0.24E2 * 
     #t534 * t550 - t541 * t550 / 0.24E2 - 0.5E1 / 0.8E1 * t541 * t465 +
     # 0.17E2 / 0.24E2 * t534 * t461
        cg0 = t328 + t417 + t518 + t600
        rhs(i1,i2,i3,i4) = rhs(i1,i2,i3,i4)+cg0
      end do
      end do
      end do
      end do
      return
      end
c
c +++++++++++++
c
      subroutine appendRosenbluthCollision_4thMapleAsinh(
     *     rhs,fin,velocities,IVx,IVy,din,Vth,N,massR,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     dx, 
     *     dparams,
     *     iparams)
c
c.. compute and add the Rosenbluth collision operator to the RHS ...
      implicit none
c
c.. declarations of incoming variables 
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
c
      real dparams(*)
      integer iparams(*)
c
      real rhs(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real fin(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real din(n1a:n1b, n2a:n2b, n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2), 1:3)
      real IVx(nd1a:nd1b, nd2a:nd2b)
      real IVy(nd1a:nd1b, nd2a:nd2b)
      real Vth(nd1a:nd1b, nd2a:nd2b)
      real N(nd1a:nd1b, nd2a:nd2b)
      real massR
c     
      real dx(1:4)
c
c.. declarations of local variables
      real nu_base, nu, safelog
      integer i1, i2, i3, i4, i, j, modulate_nu
      real cg0
c
      real dvx, dvy, alpha, beta
      real vx(-2:2,-2:2), vy(-2:2,-2:2)
      real f(-2:2,-2:2)
      real DD(-2:2,-2:2,1:3)
c
      real t103
      real t106
      real t109
      real t11
      real t112
      real t115
      real t118
      real t121
      real t124
      real t127
      real t131
      real t133
      real t137
      real t14
      real t140
      real t142
      real t144
      real t146
      real t150
      real t153
      real t155
      real t157
      real t159
      real t163
      real t166
      real t168
      real t17
      real t170
      real t172
      real t176
      real t179
      real t181
      real t184
      real t192
      real t195
      real t196
      real t198
      real t199
      real t20
      real t201
      real t204
      real t207
      real t208
      real t210
      real t211
      real t213
      real t217
      real t219
      real t223
      real t226
      real t228
      real t229
      real t23
      real t231
      real t235
      real t238
      real t240
      real t241
      real t242
      real t244
      real t248
      real t251
      real t253
      real t254
      real t255
      real t257
      real t26
      real t261
      real t264
      real t266
      real t29
      real t3
      real t306
      real t309
      real t310
      real t312
      real t315
      real t318
      real t319
      real t32
      real t321
      real t325
      real t327
      real t331
      real t334
      real t336
      real t337
      real t339
      real t343
      real t346
      real t348
      real t349
      real t350
      real t352
      real t356
      real t359
      real t36
      real t361
      real t38
      real t387
      real t39
      real t390
      real t391
      real t393
      real t396
      real t399
      real t400
      real t402
      real t406
      real t408
      real t412
      real t415
      real t417
      real t418
      real t420
      real t424
      real t427
      real t429
      real t43
      real t430
      real t431
      real t433
      real t437
      real t440
      real t442
      real t453
      real t456
      real t459
      real t46
      real t462
      real t466
      real t468
      real t472
      real t475
      real t477
      real t478
      real t48
      real t480
      real t484
      real t487
      real t489
      real t49
      real t5
      real t51
      real t55
      real t58
      real t60
      real t62
      real t64
      real t68
      real t7
      real t71
      real t73
      real t75
      real t77
      real t81
      real t84
      real t86
      real t88
      real t9
      real t90
      real t94
      real t97
      real t99
c
      modulate_nu = iparams(4)
      nu_base = dparams(1)
c
      dvx = dx(3)
      dvy = dx(4)
c      
      do i4 = n4a, n4b
      do i3 = n3a, n3b
      do i2 = n2a, n2b
      do i1 = n1a, n1b
        if (modulate_nu .eq. 1) then
          nu = nu_base*N(i1,i2)/Vth(i1,i2)
        else
          nu = nu_base
        end if
        alpha = 1.0 / (massR * (Vth(i1,i2)**2))
        beta = 2.0*alpha
        do i = -2, 2
          do j = -2, 2
            vx(i,j) = velocities(i3+i,i4+j, 0) - IVx(i1,i2)
            vy(i,j) = velocities(i3+i,i4+j, 1) - IVy(i1,i2)
            f(i,j) = fin(i1,i2,i3+i,i4+j)
            DD(i,j,1) = din(i1,i2,i3+i,i4+j,1)
            DD(i,j,2) = din(i1,i2,i3+i,i4+j,2)
            DD(i,j,3) = din(i1,i2,i3+i,i4+j,3)
c            alpha(i,j) = 1.0 / (dparams(5) * (Vth(i1,i2)**2))
          end do
        end do
        t3 = vx(0,0) ** 2
        t5 = t3 * alpha / 0.2E1
        t7 = vy(0,0) ** 2
        t9 = t7 * alpha / 0.2E1
        t11 = vx(0,-1) ** 2
        t14 = vx(0,1) ** 2
        t17 = vy(0,-1) ** 2
        t20 = vy(0,1) ** 2
        t23 = vx(0,-2) ** 2
        t26 = vx(0,2) ** 2
        t29 = vy(0,-2) ** 2
        t32 = vy(0,2) ** 2
        t36 = f(0,0)
        t38 = beta ** 2
        t39 = t36 ** 2
        t43 = exp(-(t3 + t7) * alpha)
        t46 = sqrt(t38 * t39 + 0.4E1 * t43)
        t48 = safelog(beta * t36 + t46)
        t49 = f(0,-1)
        t51 = t49 ** 2
        t55 = exp(-(t11 + t17) * alpha)
        t58 = sqrt(t38 * t51 + 0.4E1 * t55)
        t60 = safelog(beta * t49 + t58)
        t62 = f(0,1)
        t64 = t62 ** 2
        t68 = exp(-(t14 + t20) * alpha)
        t71 = sqrt(t38 * t64 + 0.4E1 * t68)
        t73 = safelog(beta * t62 + t71)
        t75 = f(0,-2)
        t77 = t75 ** 2
        t81 = exp(-(t23 + t29) * alpha)
        t84 = sqrt(t38 * t77 + 0.4E1 * t81)
        t86 = safelog(beta * t75 + t84)
        t88 = f(0,2)
        t90 = t88 ** 2
        t94 = exp(-(t26 + t32) * alpha)
        t97 = sqrt(t38 * t90 + 0.4E1 * t94)
        t99 = safelog(beta * t88 + t97)
        t103 = dvx ** 2
        t106 = vx(-1,0) ** 2
        t109 = vx(1,0) ** 2
        t112 = vy(-1,0) ** 2
        t115 = vy(1,0) ** 2
        t118 = vx(-2,0) ** 2
        t121 = vx(2,0) ** 2
        t124 = vy(-2,0) ** 2
        t127 = vy(2,0) ** 2
        t131 = f(-1,0)
        t133 = t131 ** 2
        t137 = exp(-(t106 + t112) * alpha)
        t140 = sqrt(t133 * t38 + 0.4E1 * t137)
        t142 = safelog(beta * t131 + t140)
        t144 = f(1,0)
        t146 = t144 ** 2
        t150 = exp(-(t109 + t115) * alpha)
        t153 = sqrt(t146 * t38 + 0.4E1 * t150)
        t155 = safelog(beta * t144 + t153)
        t157 = f(-2,0)
        t159 = t157 ** 2
        t163 = exp(-(t118 + t124) * alpha)
        t166 = sqrt(t159 * t38 + 0.4E1 * t163)
        t168 = safelog(beta * t157 + t166)
        t170 = f(2,0)
        t172 = t170 ** 2
        t176 = exp(-(t121 + t127) * alpha)
        t179 = sqrt(t172 * t38 + 0.4E1 * t176)
        t181 = safelog(beta * t170 + t179)
        t184 = dvy ** 2
        t192 = vx(-1,-2) ** 2
        t195 = vx(-1,-1) ** 2
        t196 = 0.5E1 * t195
        t198 = vx(-1,1) ** 2
        t199 = 0.5E1 * t198
        t201 = vx(-1,2) ** 2
        t204 = vy(-1,-2) ** 2
        t207 = vy(-1,-1) ** 2
        t208 = 0.5E1 * t207
        t210 = vy(-1,1) ** 2
        t211 = 0.5E1 * t210
        t213 = vy(-1,2) ** 2
        t217 = f(-1,-2)
        t219 = t217 ** 2
        t223 = exp(-(t192 + t204) * alpha)
        t226 = sqrt(t219 * t38 + 0.4E1 * t223)
        t228 = safelog(beta * t217 + t226)
        t229 = f(-1,-1)
        t231 = t229 ** 2
        t235 = exp(-(t195 + t207) * alpha)
        t238 = sqrt(t231 * t38 + 0.4E1 * t235)
        t240 = safelog(beta * t229 + t238)
        t241 = 0.10E2 * t240
        t242 = f(-1,1)
        t244 = t242 ** 2
        t248 = exp(-(t198 + t210) * alpha)
        t251 = sqrt(t244 * t38 + 0.4E1 * t248)
        t253 = safelog(beta * t242 + t251)
        t254 = 0.10E2 * t253
        t255 = f(-1,2)
        t257 = t255 ** 2
        t261 = exp(-(t201 + t213) * alpha)
        t264 = sqrt(t257 * t38 + 0.4E1 * t261)
        t266 = safelog(beta * t255 + t264)
        t306 = vx(-2,-1) ** 2
        t309 = vx(1,-1) ** 2
        t310 = 0.5E1 * t309
        t312 = vx(2,-1) ** 2
        t315 = vy(-2,-1) ** 2
        t318 = vy(1,-1) ** 2
        t319 = 0.5E1 * t318
        t321 = vy(2,-1) ** 2
        t325 = f(-2,-1)
        t327 = t325 ** 2
        t331 = exp(-(t306 + t315) * alpha)
        t334 = sqrt(t327 * t38 + 0.4E1 * t331)
        t336 = safelog(beta * t325 + t334)
        t337 = f(1,-1)
        t339 = t337 ** 2
        t343 = exp(-(t309 + t318) * alpha)
        t346 = sqrt(t339 * t38 + 0.4E1 * t343)
        t348 = safelog(beta * t337 + t346)
        t349 = 0.10E2 * t348
        t350 = f(2,-1)
        t352 = t350 ** 2
        t356 = exp(-(t312 + t321) * alpha)
        t359 = sqrt(t352 * t38 + 0.4E1 * t356)
        t361 = safelog(beta * t350 + t359)
        t387 = vx(-2,1) ** 2
        t390 = vx(1,1) ** 2
        t391 = 0.5E1 * t390
        t393 = vx(2,1) ** 2
        t396 = vy(-2,1) ** 2
        t399 = vy(1,1) ** 2
        t400 = 0.5E1 * t399
        t402 = vy(2,1) ** 2
        t406 = f(-2,1)
        t408 = t406 ** 2
        t412 = exp(-(t387 + t396) * alpha)
        t415 = sqrt(t38 * t408 + 0.4E1 * t412)
        t417 = safelog(beta * t406 + t415)
        t418 = f(1,1)
        t420 = t418 ** 2
        t424 = exp(-(t390 + t399) * alpha)
        t427 = sqrt(t38 * t420 + 0.4E1 * t424)
        t429 = safelog(beta * t418 + t427)
        t430 = 0.10E2 * t429
        t431 = f(2,1)
        t433 = t431 ** 2
        t437 = exp(-(t393 + t402) * alpha)
        t440 = sqrt(t38 * t433 + 0.4E1 * t437)
        t442 = safelog(beta * t431 + t440)
        t453 = vx(1,-2) ** 2
        t456 = vx(1,2) ** 2
        t459 = vy(1,-2) ** 2
        t462 = vy(1,2) ** 2
        t466 = f(1,-2)
        t468 = t466 ** 2
        t472 = exp(-(t453 + t459) * alpha)
        t475 = sqrt(t38 * t468 + 0.4E1 * t472)
        t477 = safelog(beta * t466 + t475)
        t478 = f(1,2)
        t480 = t478 ** 2
        t484 = exp(-(t456 + t462) * alpha)
        t487 = sqrt(t38 * t480 + 0.4E1 * t484)
        t489 = safelog(beta * t478 + t487)
        cg0 = -0.11E2 / 0.8E1 / beta / t103 / t184 * nu * (t46 * (t103 *
     # (t5 + t9 + alpha * (-0.17E2 / 0.66E2 * t11 - 0.17E2 / 0.66E2 * t1
     #4 - 0.17E2 / 0.66E2 * t17 - 0.17E2 / 0.66E2 * t20 + t23 / 0.132E3 
     #+ t26 / 0.132E3 + t29 / 0.132E3 + t32 / 0.132E3) + t48 - 0.17E2 / 
     #0.33E2 * t60 - 0.17E2 / 0.33E2 * t73 + t86 / 0.66E2 + t99 / 0.66E2
     #) * DD(0,0,3) + DD(0,0,1) * t184 * (t5 + t9 + alpha * (-0.17E2 / 0
     #.66E2 * t106 - 0.17E2 / 0.66E2 * t109 - 0.17E2 / 0.66E2 * t112 - 0
     #.17E2 / 0.66E2 * t115 + t118 / 0.132E3 + t121 / 0.132E3 + t124 / 0
     #.132E3 + t127 / 0.132E3) - 0.17E2 / 0.33E2 * t142 - 0.17E2 / 0.33E
     #2 * t155 + t168 / 0.66E2 + t181 / 0.66E2 + t48)) + 0.5E1 / 0.11E2 
     #* t140 * dvy * (dvx * (alpha * (t192 / 0.2E1 - t196 + t199 - t201 
     #/ 0.2E1 + t204 / 0.2E1 - t208 + t211 - t213 / 0.2E1) + t228 - t241
     # + t254 - t266) * DD(-1,0,2) / 0.15E2 + dvy * (t5 + t9 + alpha * (
     #t118 / 0.30E2 - 0.17E2 / 0.30E2 * t106 + t109 / 0.30E2 + t124 / 0.
     #30E2 - 0.17E2 / 0.30E2 * t112 + t115 / 0.30E2) + t48 - 0.17E2 / 0.
     #15E2 * t142 + t155 / 0.15E2 + t168 / 0.15E2) * DD(-1,0,1)) + 0.5E1
     # / 0.11E2 * t58 * (dvx * DD(0,-1,3) * (t5 + t9 + alpha * (t23 / 0.
     #30E2 - 0.17E2 / 0.30E2 * t11 + t14 / 0.30E2 + t29 / 0.30E2 - 0.17E
     #2 / 0.30E2 * t17 + t20 / 0.30E2) + t48 - 0.17E2 / 0.15E2 * t60 + t
     #73 / 0.15E2 + t86 / 0.15E2) + DD(0,-1,2) * dvy * (alpha * (t306 / 
     #0.2E1 - t196 + t310 - t312 / 0.2E1 + t315 / 0.2E1 - t208 + t319 - 
     #t321 / 0.2E1) + t336 - t241 + t349 - t361) / 0.15E2) * dvx + 0.5E1
     # / 0.11E2 * t71 * (dvx * (t5 + t9 + alpha * (t11 / 0.30E2 - 0.17E2
     # / 0.30E2 * t14 + t26 / 0.30E2 + t17 / 0.30E2 - 0.17E2 / 0.30E2 * 
     #t20 + t32 / 0.30E2) + t48 + t60 / 0.15E2 - 0.17E2 / 0.15E2 * t73 +
     # t99 / 0.15E2) * DD(0,1,3) - DD(0,1,2) * dvy * (alpha * (t387 / 0.
     #2E1 - t199 + t391 - t393 / 0.2E1 + t396 / 0.2E1 - t211 + t400 - t4
     #02 / 0.2E1) + t417 - t254 + t430 - t442) / 0.15E2) * dvx + 0.5E1 /
     # 0.11E2 * t153 * dvy * (-dvx * DD(1,0,2) * (alpha * (t453 / 0.2E1 
     #- t310 + t391 - t456 / 0.2E1 + t459 / 0.2E1 - t319 + t400 - t462 /
     # 0.2E1) + t477 - t349 + t430 - t489) / 0.15E2 + dvy * DD(1,0,1) * 
     #(t5 + t9 + alpha * (t106 / 0.30E2 - 0.17E2 / 0.30E2 * t109 + t121 
     #/ 0.30E2 + t112 / 0.30E2 - 0.17E2 / 0.30E2 * t115 + t127 / 0.30E2)
     # + t48 + t142 / 0.15E2 - 0.17E2 / 0.15E2 * t155 + t181 / 0.15E2)) 
     #- t166 * dvy * (-0.2E1 / 0.3E1 * dvx * (alpha * (t306 / 0.2E1 - t3
     #87 / 0.2E1 + t315 / 0.2E1 - t396 / 0.2E1) + t336 - t417) * DD(-2,0
     #,2) + dvy * (t5 + t9 + alpha * (-t118 / 0.6E1 - t106 / 0.3E1 - t12
     #4 / 0.6E1 - t112 / 0.3E1) + t48 - 0.2E1 / 0.3E1 * t142 - t168 / 0.
     #3E1) * DD(-2,0,1)) / 0.22E2 - t84 * (dvx * DD(0,-2,3) * (t5 + t9 +
     # alpha * (-t23 / 0.6E1 - t11 / 0.3E1 - t29 / 0.6E1 - t17 / 0.3E1) 
     #+ t48 - 0.2E1 / 0.3E1 * t60 - t86 / 0.3E1) - 0.2E1 / 0.3E1 * (alph
     #a * (t192 / 0.2E1 - t453 / 0.2E1 + t204 / 0.2E1 - t459 / 0.2E1) + 
     #t228 - t477) * dvy * DD(0,-2,2)) * dvx / 0.22E2 - t97 * dvx * (dvx
     # * (t5 + t9 + alpha * (-t14 / 0.3E1 - t26 / 0.6E1 - t20 / 0.3E1 - 
     #t32 / 0.6E1) + t48 - 0.2E1 / 0.3E1 * t73 - t99 / 0.3E1) * DD(0,2,3
     #) + 0.2E1 / 0.3E1 * dvy * DD(0,2,2) * (alpha * (t201 / 0.2E1 - t45
     #6 / 0.2E1 + t213 / 0.2E1 - t462 / 0.2E1) + t266 - t489)) / 0.22E2 
     #- (0.2E1 / 0.3E1 * dvx * DD(2,0,2) * (alpha * (t312 / 0.2E1 - t393
     # / 0.2E1 + t321 / 0.2E1 - t402 / 0.2E1) + t361 - t442) + (t5 + t9 
     #+ alpha * (-t109 / 0.3E1 - t121 / 0.6E1 - t115 / 0.3E1 - t127 / 0.
     #6E1) + t48 - 0.2E1 / 0.3E1 * t155 - t181 / 0.3E1) * dvy * DD(2,0,1
     #)) * t179 * dvy / 0.22E2)
        rhs(i1,i2,i3,i4) = rhs(i1,i2,i3,i4)+cg0
      end do
      end do
      end do
      end do
      return
      end
c
c +++++++++++++
c      
      subroutine appendRosenbluthCollision_6thMaplePrimitive(
     *     rhs,fin,velocities,IVx,IVy,din,Vth,N,massR,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     dx,
     *     dparams,
     *     iparams)
c
c.. compute and add the Rosenbluth collision operator to the RHS ...
      implicit none
c
c.. declarations of incoming variables 
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
c
      real dparams(*)
      integer iparams(*)
c
      real rhs(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real fin(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real din(n1a:n1b, n2a:n2b, n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2), 1:3)
      real IVx(nd1a:nd1b, nd2a:nd2b)
      real IVy(nd1a:nd1b, nd2a:nd2b)
      real Vth(nd1a:nd1b, nd2a:nd2b)
      real N(nd1a:nd1b, nd2a:nd2b)
      real massR
c     
      real dx(1:4)
c
c.. declarations of local variables
      real nu_base, nu
      integer i1, i2, i3, i4, i, j, modulate_nu
      real cg0
c
      real dvx, dvy, alpha
      real vx(-3:3,-3:3), vy(-3:3,-3:3)
      real f(-3:3,-3:3)
      real DD(-3:3,-3:3,1:3)
c
      real t1
      real t101
      real t103
      real t108
      real t11
      real t111
      real t114
      real t115
      real t117
      real t12
      real t120
      real t126
      real t127
      real t13
      real t130
      real t138
      real t14
      real t140
      real t142
      real t146
      real t149
      real t150
      real t158
      real t16
      real t161
      real t167
      real t168
      real t170
      real t172
      real t174
      real t176
      real t178
      real t18
      real t181
      real t184
      real t185
      real t187
      real t188
      real t192
      real t194
      real t196
      real t202
      real t203
      real t204
      real t212
      real t216
      real t22
      real t222
      real t226
      real t23
      real t231
      real t234
      real t243
      real t248
      real t252
      real t273
      real t277
      real t28
      real t285
      real t288
      real t3
      real t30
      real t302
      real t31
      real t318
      real t32
      real t320
      real t339
      real t34
      real t346
      real t352
      real t353
      real t368
      real t39
      real t392
      real t40
      real t408
      real t409
      real t411
      real t413
      real t415
      real t417
      real t418
      real t420
      real t421
      real t422
      real t423
      real t425
      real t427
      real t431
      real t432
      real t435
      real t438
      real t440
      real t441
      real t443
      real t448
      real t449
      real t456
      real t458
      real t459
      real t460
      real t462
      real t466
      real t47
      real t470
      real t473
      real t474
      real t476
      real t479
      real t481
      real t483
      real t485
      real t489
      real t49
      real t491
      real t494
      real t497
      real t499
      real t5
      real t500
      real t501
      real t504
      real t506
      real t508
      real t51
      real t510
      real t512
      real t517
      real t519
      real t52
      real t523
      real t524
      real t526
      real t53
      real t530
      real t535
      real t536
      real t539
      real t548
      real t549
      real t551
      real t555
      real t558
      real t559
      real t567
      real t57
      real t571
      real t576
      real t577
      real t579
      real t581
      real t583
      real t585
      real t587
      real t590
      real t593
      real t594
      real t596
      real t597
      real t601
      real t603
      real t605
      real t61
      real t611
      real t612
      real t613
      real t621
      real t625
      real t63
      real t631
      real t636
      real t640
      real t643
      real t65
      real t652
      real t656
      real t661
      real t67
      real t681
      real t686
      real t69
      real t694
      real t697
      real t7
      real t71
      real t711
      real t728
      real t73
      real t730
      real t74
      real t750
      real t757
      real t76
      real t763
      real t766
      real t768
      real t771
      real t773
      real t775
      real t777
      real t780
      real t782
      real t785
      real t787
      real t789
      real t792
      real t794
      real t797
      real t799
      real t8
      real t80
      real t801
      real t804
      real t806
      real t809
      real t811
      real t815
      real t817
      real t82
      real t821
      real t825
      real t842
      real t846
      real t85
      real t88
      real t880
      real t896
      real t9
      real t90
      real t94
      real t95
      real t96
      real t97
      real t99
c
      modulate_nu = iparams(4)
      nu_base = dparams(1)
c     
      dvx = dx(3)
      dvy = dx(4)
c     
      do i4 = n4a, n4b
      do i3 = n3a, n3b
      do i2 = n2a, n2b
      do i1 = n1a, n1b
        if (modulate_nu .eq. 1) then
          nu = nu_base*N(i1,i2)/Vth(i1,i2)
        else
          nu = nu_base
        end if
        alpha = 1.0 / (massR * (Vth(i1,i2)**2))
        do i = -3, 3
          do j = -3, 3
            vx(i,j) = velocities(i3+i,i4+j, 0) - IVx(i1,i2)
            vy(i,j) = velocities(i3+i,i4+j, 1) - IVy(i1,i2)
            f(i,j) = fin(i1,i2,i3+i,i4+j)
            DD(i,j,1) = din(i1,i2,i3+i,i4+j,1)
            DD(i,j,2) = din(i1,i2,i3+i,i4+j,2)
            DD(i,j,3) = din(i1,i2,i3+i,i4+j,3)
          end do
        end do
ccccc
cc generated code 
ccccc
        t1 = vy(0,-3)
        t3 = vy(0,-2)
        t5 = vy(0,-1)
        t7 = vy(0,0)
        t8 = 0.3336741E7 / 0.50753E5 * t7
        t9 = vy(0,2)
        t11 = vy(0,3)
        t12 = 0.16875E5 / 0.101506E6 * t11
        t13 = 0.49473E5 / 0.50753E5 * t1 - 0.868347E6 / 0.101506E6 * t3 
     #+ 0.5643027E7 / 0.101506E6 * t5 + t8 - 0.46287E5 / 0.101506E6 * t9
     # + t12
        t14 = DD(0,-1,3)
        t16 = 0.16875E5 / 0.101506E6 * t1
        t18 = vy(0,1)
        t22 = -t16 + 0.46287E5 / 0.101506E6 * t3 - t8 - 0.5643027E7 / 0.
     #101506E6 * t18 + 0.868347E6 / 0.101506E6 * t9 - 0.49473E5 / 0.5075
     #3E5 * t11
        t23 = DD(0,1,3)
        t28 = 0.504756E6 / 0.50753E5 * t7
        t30 = 0.2025E4 / 0.101506E6 * t11
        t31 = -0.19179E5 / 0.101506E6 * t1 + 0.83124E5 / 0.50753E5 * t3 
     #- 0.868347E6 / 0.101506E6 * t5 - t28 + 0.46287E5 / 0.101506E6 * t1
     #8 - t30
        t32 = DD(0,-2,3)
        t34 = 0.2025E4 / 0.101506E6 * t1
        t39 = t34 - 0.46287E5 / 0.101506E6 * t5 + t28 + 0.868347E6 / 0.1
     #01506E6 * t18 - 0.83124E5 / 0.50753E5 * t9 + 0.19179E5 / 0.101506E
     #6 * t11
        t40 = DD(0,2,3)
        t47 = DD(0,0,3)
        t49 = DD(0,-3,3)
        t51 = DD(0,3,3)
        t52 = 0.2025E4 / 0.101506E6 * t51
        t53 = -0.19179E5 / 0.101506E6 * t49 - t52
        t57 = 0.49473E5 / 0.50753E5 * t49 + 0.16875E5 / 0.101506E6 * t51
        t61 = -0.16875E5 / 0.101506E6 * t49 - 0.49473E5 / 0.50753E5 * t5
     #1
        t63 = 0.2025E4 / 0.101506E6 * t49
        t65 = t63 + 0.19179E5 / 0.101506E6 * t51
        t67 = t49 - t51
        t69 = t1 * t49
        t71 = t11 * t51
        t73 = t14 * t13 + t23 * t22 + t32 * t31 + t40 * t39 + t47 * (t1 
     #- 0.504756E6 / 0.50753E5 * t3 + 0.3336741E7 / 0.50753E5 * t5 - 0.3
     #336741E7 / 0.50753E5 * t18 + 0.504756E6 / 0.50753E5 * t9 - t11) + 
     #t3 * t53 + t5 * t57 + t18 * t61 + t9 * t65 + t7 * t67 + 0.2281E4 /
     # 0.101506E6 * t69 - 0.2281E4 / 0.101506E6 * t71
        t74 = f(0,0)
        t76 = 0.1053495E7 / 0.101506E6 * t18
        t80 = 0.5643027E7 / 0.101506E6 * t7
        t82 = 0.5625E4 / 0.203012E6 * t11
        t85 = 0.1053495E7 / 0.101506E6 * t5
        t88 = -t85 - t16 + 0.82206E5 / 0.50753E5 * t3 + t76 - 0.82206E5 
     #/ 0.50753E5 * t9 + t12
        t90 = 0.39375E5 / 0.203012E6 * t9
        t94 = 0.868347E6 / 0.101506E6 * t7
        t95 = 0.82206E5 / 0.50753E5 * t18
        t96 = 0.675E3 / 0.203012E6 * t11
        t97 = -t90 + 0.3861E4 / 0.101506E6 * t1 - 0.224955E6 / 0.203012E
     #6 * t3 - 0.2338245E7 / 0.203012E6 * t5 - t94 + t95 + t96
        t99 = 0.39375E5 / 0.203012E6 * t3
        t101 = 0.46287E5 / 0.101506E6 * t7
        t103 = -t99 + t34 + 0.250911E6 / 0.203012E6 * t5 - t101 - t95 + 
     #0.52461E5 / 0.203012E6 * t9 - t82
        t108 = 0.3861E4 / 0.101506E6 * t49 + 0.675E3 / 0.203012E6 * t51
        t111 = 0.5625E4 / 0.203012E6 * t51
        t114 = -t67
        t115 = 0.16875E5 / 0.101506E6 * t114
        t117 = t63 - t111
        t120 = DD(0,-1,2)
        t126 = t14 * (-t76 + 0.60993E5 / 0.50753E5 * t1 - 0.2338245E7 / 
     #0.203012E6 * t3 + 0.11067363E8 / 0.203012E6 * t5 + t80 + 0.250911E
     #6 / 0.203012E6 * t9 - t82) + t23 * t88 + t32 * t97 + t40 * t103 + 
     #t47 * t13 + t3 * t108 + t5 * (0.60993E5 / 0.50753E5 * t49 - t111) 
     #+ t18 * t115 + t9 * t117 + t7 * t57 + 0.9437184E7 / 0.50753E5 * vx
     #(0,-1) * t120 - 0.279E3 / 0.101506E6 * t69 + 0.675E3 / 0.203012E6 
     #* t71
        t127 = f(0,-1)
        t130 = 0.5625E4 / 0.203012E6 * t1
        t138 = 0.82206E5 / 0.50753E5 * t5
        t140 = t90 + t130 - 0.52461E5 / 0.203012E6 * t3 + t138 + t101 - 
     #0.250911E6 / 0.203012E6 * t18 - t30
        t142 = 0.675E3 / 0.203012E6 * t1
        t146 = t99 - t142 - t138 + t94 + 0.2338245E7 / 0.203012E6 * t18 
     #+ 0.224955E6 / 0.203012E6 * t9 - 0.3861E4 / 0.101506E6 * t11
        t149 = 0.5625E4 / 0.203012E6 * t49
        t150 = t149 - t52
        t158 = -0.675E3 / 0.203012E6 * t49 - 0.3861E4 / 0.101506E6 * t51
        t161 = DD(0,1,2)
        t167 = t14 * t88 + t23 * (t85 + t130 - 0.250911E6 / 0.203012E6 *
     # t3 - t80 - 0.11067363E8 / 0.203012E6 * t18 + 0.2338245E7 / 0.2030
     #12E6 * t9 - 0.60993E5 / 0.50753E5 * t11) + t32 * t140 + t40 * t146
     # + t47 * t22 + t3 * t150 + t5 * t115 + t18 * (t149 - 0.60993E5 / 0
     #.50753E5 * t51) + t9 * t158 + t7 * t61 - 0.9437184E7 / 0.50753E5 *
     # vx(0,1) * t161 - 0.675E3 / 0.203012E6 * t69 + 0.279E3 / 0.101506E
     #6 * t71
        t168 = f(0,1)
        t170 = f(0,-2)
        t172 = f(0,2)
        t174 = f(0,-3)
        t176 = f(0,3)
        t178 = 0.3861E4 / 0.101506E6 * t174 + 0.675E3 / 0.203012E6 * t17
     #6
        t181 = 0.5625E4 / 0.203012E6 * t176
        t184 = -t174 + t176
        t185 = 0.16875E5 / 0.101506E6 * t184
        t187 = 0.2025E4 / 0.101506E6 * t174
        t188 = t187 - t181
        t192 = 0.49473E5 / 0.50753E5 * t174 + 0.16875E5 / 0.101506E6 * t
     #176
        t194 = t176 * t11
        t196 = t1 * t174
        t202 = 0.5625E4 / 0.203012E6 * t174
        t203 = 0.2025E4 / 0.101506E6 * t176
        t204 = t202 - t203
        t212 = -0.675E3 / 0.203012E6 * t174 - 0.3861E4 / 0.101506E6 * t1
     #76
        t216 = -0.16875E5 / 0.101506E6 * t174 - 0.49473E5 / 0.50753E5 * 
     #t176
        t222 = 0.1575E4 / 0.50753E5 * t9
        t226 = 0.83124E5 / 0.50753E5 * t7
        t231 = 0.1575E4 / 0.50753E5 * t3
        t234 = t231 - t142 - 0.39375E5 / 0.203012E6 * t5 + 0.39375E5 / 0
     #.203012E6 * t18 - t222 + t96
        t243 = 0.675E3 / 0.203012E6 * t114
        t248 = DD(0,-2,2)
        t252 = t32 * (t222 + 0.382941E6 / 0.1015060E7 * t1 - 0.674631E6 
     #/ 0.253765E6 * t3 - 0.224955E6 / 0.203012E6 * t5 + t226 - 0.52461E
     #5 / 0.203012E6 * t18 - 0.81E2 / 0.203012E6 * t11) + t40 * t234 + t
     #47 * t31 + t3 * (0.382941E6 / 0.1015060E7 * t49 - 0.81E2 / 0.20301
     #2E6 * t51) + t5 * t108 + t18 * t150 + t9 * t243 + t7 * t53 + 0.121
     #041E6 / 0.1015060E7 * t69 - 0.81E2 / 0.203012E6 * t71 - 0.9437184E
     #7 / 0.253765E6 * vx(0,-2) * t248
        t273 = DD(0,2,2)
        t277 = t32 * t234 + t40 * (-t231 + 0.81E2 / 0.203012E6 * t1 + 0.
     #52461E5 / 0.203012E6 * t5 - t226 + 0.224955E6 / 0.203012E6 * t18 +
     # 0.674631E6 / 0.253765E6 * t9 - 0.382941E6 / 0.1015060E7 * t11) + 
     #t47 * t39 + t3 * t243 + t5 * t117 + t18 * t158 + t9 * (0.81E2 / 0.
     #203012E6 * t49 - 0.382941E6 / 0.1015060E7 * t51) + t7 * t65 + 0.81
     #E2 / 0.203012E6 * t69 - 0.121041E6 / 0.1015060E7 * t71 + 0.9437184
     #E7 / 0.253765E6 * vx(0,2) * t273
        t285 = 0.675E3 / 0.203012E6 * t184
        t288 = -0.19179E5 / 0.101506E6 * t174 - t203
        t302 = t187 + 0.19179E5 / 0.101506E6 * t176
        t318 = t174 * t49
        t320 = t176 * t51
        t339 = DD(0,-3,2)
        t346 = DD(0,3,2)
        t353 = t74 * t73 + t127 * t126 + t168 * t167 + t14 * (t170 * t97
     # + t172 * t103 + t3 * t178 + t5 * (0.60993E5 / 0.50753E5 * t174 - 
     #t181) + t18 * t185 + t9 * t188 + t7 * t192 + 0.675E3 / 0.203012E6 
     #* t194 - 0.279E3 / 0.101506E6 * t196) + t23 * (t170 * t140 + t172 
     #* t146 + t3 * t204 + t5 * t185 + t18 * (t202 - 0.60993E5 / 0.50753
     #E5 * t176) + t9 * t212 + t7 * t216 + 0.279E3 / 0.101506E6 * t194 -
     # 0.675E3 / 0.203012E6 * t196) + t170 * t252 + t172 * t277 + t32 * 
     #(t3 * (0.382941E6 / 0.1015060E7 * t174 - 0.81E2 / 0.203012E6 * t17
     #6) + t5 * t178 + t18 * t204 + t9 * t285 + t7 * t288 + 0.121041E6 /
     # 0.1015060E7 * t196 - 0.81E2 / 0.203012E6 * t194) + t40 * (t3 * t2
     #85 + t5 * t188 + t18 * t212 + t9 * (0.81E2 / 0.203012E6 * t174 - 0
     #.382941E6 / 0.1015060E7 * t176) + t7 * t302 + 0.81E2 / 0.203012E6 
     #* t196 - 0.121041E6 / 0.1015060E7 * t194) + t47 * (t3 * t288 + t5 
     #* t192 + t18 * t216 + t9 * t302 - t7 * t184 + 0.2281E4 / 0.101506E
     #6 * t196 - 0.2281E4 / 0.101506E6 * t194) + t3 * (0.121041E6 / 0.10
     #15060E7 * t318 - 0.81E2 / 0.203012E6 * t320) + t5 * (-0.279E3 / 0.
     #101506E6 * t318 + 0.675E3 / 0.203012E6 * t320) + t18 * (-0.675E3 /
     # 0.203012E6 * t318 + 0.279E3 / 0.101506E6 * t320) + t9 * (0.81E2 /
     # 0.203012E6 * t318 - 0.121041E6 / 0.1015060E7 * t320) + t7 * (0.22
     #81E4 / 0.101506E6 * t318 - 0.2281E4 / 0.101506E6 * t320) + t174 * 
     #(0.1048576E7 / 0.253765E6 * vx(0,-3) * t339 + 0.150421E6 / 0.10150
     #60E7 * t69) - 0.1048576E7 / 0.253765E6 * (vx(0,3) * t346 + 0.15042
     #1E6 / 0.4194304E7 * t71) * t176
        t368 = 0.10272768E8 / 0.50753E5 * t47
        t392 = 0.2359296E7 / 0.253765E6 * t47
        t352 = alpha * dvy
        t408 = t353 * t352 + t74 * (0.147456E6 / 0.50753E5 * t49 - 0.380
     #1088E7 / 0.152259E6 * t32 + 0.25247744E8 / 0.152259E6 * t14 + 0.19
     #660800E8 / 0.50753E5 * t47 + 0.25247744E8 / 0.152259E6 * t23 - 0.3
     #801088E7 / 0.152259E6 * t40 + 0.147456E6 / 0.50753E5 * t51) + t127
     # * (-0.114688E6 / 0.50753E5 * t49 + 0.950272E6 / 0.50753E5 * t32 -
     # 0.10272768E8 / 0.50753E5 * t14 - t368 + 0.950272E6 / 0.50753E5 * 
     #t23 - 0.114688E6 / 0.50753E5 * t40) + t168 * (-0.114688E6 / 0.5075
     #3E5 * t32 + 0.950272E6 / 0.50753E5 * t14 - t368 - 0.10272768E8 / 0
     #.50753E5 * t23 + 0.950272E6 / 0.50753E5 * t40 - 0.114688E6 / 0.507
     #53E5 * t51) + t14 * (0.950272E6 / 0.50753E5 * t170 - 0.16384E5 / 0
     #.253765E6 * t172 - 0.606208E6 / 0.761295E6 * t174) + t23 * (-0.163
     #84E5 / 0.253765E6 * t170 + 0.950272E6 / 0.50753E5 * t172 - 0.60620
     #8E6 / 0.761295E6 * t176) + t170 * (-0.16384E5 / 0.253765E6 * t49 +
     # 0.2359296E7 / 0.253765E6 * t32 + t392) + t172 * (t392 + 0.2359296
     #E7 / 0.253765E6 * t40 - 0.16384E5 / 0.253765E6 * t51) - 0.606208E6
     # / 0.761295E6 * t174 * t32 - 0.606208E6 / 0.761295E6 * t176 * t40 
     #+ t47 * (-0.147456E6 / 0.253765E6 * t174 - 0.147456E6 / 0.253765E6
     # * t176) - 0.147456E6 / 0.253765E6 * t318 - 0.147456E6 / 0.253765E
     #6 * t320
        t409 = dvx ** 2
        t411 = vx(-3,0)
        t413 = vx(-2,0)
        t415 = vx(-1,0)
        t417 = vx(0,0)
        t418 = vx(2,0)
        t420 = vx(3,0)
        t421 = 0.625E3 / 0.247166E6 * t420
        t422 = -0.5497E4 / 0.370749E6 * t411 + 0.32161E5 / 0.247166E6 * 
     #t413 - 0.209001E6 / 0.247166E6 * t415 - t417 + 0.5143E4 / 0.741498
     #E6 * t418 - t421
        t423 = DD(-1,0,1)
        t425 = 0.625E3 / 0.247166E6 * t411
        t427 = vx(1,0)
        t431 = t425 - 0.5143E4 / 0.741498E6 * t413 + t417 + 0.209001E6 /
     # 0.247166E6 * t427 - 0.32161E5 / 0.247166E6 * t418 + 0.5497E4 / 0.
     #370749E6 * t420
        t432 = DD(1,0,1)
        t435 = 0.75E2 / 0.247166E6 * t420
        t438 = 0.56084E5 / 0.370749E6 * t417
        t440 = -0.5143E4 / 0.741498E6 * t427 + t435 - 0.9236E4 / 0.37074
     #9E6 * t413 + 0.32161E5 / 0.247166E6 * t415 + t438 + 0.2131E4 / 0.7
     #41498E6 * t411
        t441 = DD(-2,0,1)
        t443 = 0.75E2 / 0.247166E6 * t411
        t448 = -t443 + 0.5143E4 / 0.741498E6 * t415 - t438 - 0.32161E5 /
     # 0.247166E6 * t427 + 0.9236E4 / 0.370749E6 * t418 - 0.2131E4 / 0.7
     #41498E6 * t420
        t449 = DD(2,0,1)
        t456 = DD(0,0,1)
        t458 = DD(3,0,1)
        t459 = 0.75E2 / 0.247166E6 * t458
        t460 = DD(-3,0,1)
        t462 = t459 + 0.2131E4 / 0.741498E6 * t460
        t466 = -0.625E3 / 0.247166E6 * t458 - 0.5497E4 / 0.370749E6 * t4
     #60
        t470 = 0.5497E4 / 0.370749E6 * t458 + 0.625E3 / 0.247166E6 * t46
     #0
        t473 = 0.75E2 / 0.247166E6 * t460
        t474 = -0.2131E4 / 0.741498E6 * t458 - t473
        t476 = t458 - t460
        t479 = t458 * t420
        t481 = t411 * t460
        t483 = t423 * t422 + t432 * t431 + t441 * t440 + t449 * t448 + t
     #456 * (0.50753E5 / 0.3336741E7 * t420 - 0.50753E5 / 0.3336741E7 * 
     #t411 + 0.56084E5 / 0.370749E6 * t413 - t415 + t427 - 0.56084E5 / 0
     #.370749E6 * t418) + t413 * t462 + t415 * t466 + t427 * t470 + t418
     # * t474 + 0.50753E5 / 0.3336741E7 * t417 * t476 + 0.2281E4 / 0.667
     #3482E7 * t479 - 0.2281E4 / 0.6673482E7 * t481
        t485 = 0.117055E6 / 0.741498E6 * t427
        t489 = 0.209001E6 / 0.247166E6 * t417
        t491 = 0.625E3 / 0.1482996E7 * t420
        t494 = 0.117055E6 / 0.741498E6 * t415
        t497 = t494 + t425 - 0.9134E4 / 0.370749E6 * t413 - t485 + 0.913
     #4E4 / 0.370749E6 * t418 - t421
        t499 = 0.4375E4 / 0.1482996E7 * t418
        t500 = 0.9134E4 / 0.370749E6 * t427
        t501 = 0.25E2 / 0.494332E6 * t420
        t504 = 0.32161E5 / 0.247166E6 * t417
        t506 = t499 - t500 - t501 + 0.24995E5 / 0.1482996E7 * t413 + 0.2
     #59805E6 / 0.1482996E7 * t415 + t504 - 0.143E3 / 0.247166E6 * t411
        t508 = 0.4375E4 / 0.1482996E7 * t413
        t510 = 0.5143E4 / 0.741498E6 * t417
        t512 = t508 - t443 - 0.9293E4 / 0.494332E6 * t415 + t510 + t500 
     #- 0.1943E4 / 0.494332E6 * t418 + t491
        t517 = -0.25E2 / 0.494332E6 * t458 - 0.143E3 / 0.247166E6 * t460
        t519 = 0.625E3 / 0.1482996E7 * t458
        t523 = -t476
        t524 = 0.625E3 / 0.247166E6 * t523
        t526 = t519 - t473
        t530 = DD(-1,0,2)
        t535 = t423 * (t485 - 0.2259E4 / 0.123583E6 * t411 + 0.259805E6 
     #/ 0.1482996E7 * t413 - 0.1229707E7 / 0.1482996E7 * t415 - t489 - 0
     #.9293E4 / 0.494332E6 * t418 + t491) + t432 * t497 + t441 * t506 + 
     #t449 * t512 + t456 * t422 + t413 * t517 + t415 * (t519 - 0.2259E4 
     #/ 0.123583E6 * t460) + t427 * t524 + t418 * t526 + t417 * t466 - 0
     #.1048576E7 / 0.370749E6 * t530 * vy(-1,0) - 0.25E2 / 0.494332E6 * 
     #t479 + 0.31E2 / 0.741498E6 * t481
        t536 = f(-1,0)
        t539 = 0.625E3 / 0.1482996E7 * t411
        t548 = 0.9134E4 / 0.370749E6 * t415
        t549 = -t499 + 0.9293E4 / 0.494332E6 * t427 + t435 + 0.1943E4 / 
     #0.494332E6 * t413 - t548 - t510 - t539
        t551 = 0.25E2 / 0.494332E6 * t411
        t555 = -t508 + t551 + t548 - t504 - 0.259805E6 / 0.1482996E7 * t
     #427 - 0.24995E5 / 0.1482996E7 * t418 + 0.143E3 / 0.247166E6 * t420
        t558 = 0.625E3 / 0.1482996E7 * t460
        t559 = t459 - t558
        t567 = 0.143E3 / 0.247166E6 * t458 + 0.25E2 / 0.494332E6 * t460
        t571 = DD(1,0,2)
        t576 = t423 * t497 + t432 * (-t494 - t539 + 0.9293E4 / 0.494332E
     #6 * t413 + t489 + 0.1229707E7 / 0.1482996E7 * t427 - 0.259805E6 / 
     #0.1482996E7 * t418 + 0.2259E4 / 0.123583E6 * t420) + t441 * t549 +
     # t449 * t555 + t456 * t431 + t413 * t559 + t415 * t524 + t427 * (0
     #.2259E4 / 0.123583E6 * t458 - t558) + t418 * t567 + t417 * t470 + 
     #0.1048576E7 / 0.370749E6 * t571 * vy(1,0) - 0.31E2 / 0.741498E6 * 
     #t479 + 0.25E2 / 0.494332E6 * t481
        t577 = f(1,0)
        t579 = f(-2,0)
        t581 = f(2,0)
        t583 = f(-3,0)
        t585 = f(3,0)
        t587 = -0.143E3 / 0.247166E6 * t583 - 0.25E2 / 0.494332E6 * t585
        t590 = 0.625E3 / 0.1482996E7 * t585
        t593 = t583 - t585
        t594 = 0.625E3 / 0.247166E6 * t593
        t596 = 0.75E2 / 0.247166E6 * t583
        t597 = -t596 + t590
        t601 = -0.5497E4 / 0.370749E6 * t583 - 0.625E3 / 0.247166E6 * t5
     #85
        t603 = t583 * t411
        t605 = t585 * t420
        t611 = 0.625E3 / 0.1482996E7 * t583
        t612 = 0.75E2 / 0.247166E6 * t585
        t613 = -t611 + t612
        t621 = 0.25E2 / 0.494332E6 * t583 + 0.143E3 / 0.247166E6 * t585
        t625 = 0.625E3 / 0.247166E6 * t583 + 0.5497E4 / 0.370749E6 * t58
     #5
        t631 = 0.175E3 / 0.370749E6 * t418
        t636 = 0.9236E4 / 0.370749E6 * t417
        t640 = 0.175E3 / 0.370749E6 * t413
        t643 = -t640 + t551 + 0.4375E4 / 0.1482996E7 * t415 - 0.4375E4 /
     # 0.1482996E7 * t427 + t631 - t501
        t652 = 0.25E2 / 0.494332E6 * t523
        t656 = DD(-2,0,2)
        t661 = t441 * (-t631 + 0.1943E4 / 0.494332E6 * t427 + 0.3E1 / 0.
     #494332E6 * t420 + 0.74959E5 / 0.1853745E7 * t413 + 0.24995E5 / 0.1
     #482996E7 * t415 - t636 - 0.14183E5 / 0.2471660E7 * t411) + t449 * 
     #t643 + t456 * t440 + t413 * (0.3E1 / 0.494332E6 * t458 - 0.14183E5
     # / 0.2471660E7 * t460) + t415 * t517 + t427 * t559 + t418 * t652 +
     # t417 * t462 + 0.1048576E7 / 0.1853745E7 * t656 * vy(-2,0) + 0.3E1
     # / 0.494332E6 * t479 - 0.4483E4 / 0.2471660E7 * t481
        t681 = DD(2,0,2)
        t686 = t441 * t643 + t449 * (t640 - 0.3E1 / 0.494332E6 * t411 - 
     #0.1943E4 / 0.494332E6 * t415 + t636 - 0.24995E5 / 0.1482996E7 * t4
     #27 - 0.74959E5 / 0.1853745E7 * t418 + 0.14183E5 / 0.2471660E7 * t4
     #20) + t456 * t448 + t413 * t652 + t415 * t526 + t427 * t567 + t418
     # * (0.14183E5 / 0.2471660E7 * t458 - 0.3E1 / 0.494332E6 * t460) + 
     #t417 * t474 - 0.1048576E7 / 0.1853745E7 * t681 * vy(2,0) + 0.4483E
     #4 / 0.2471660E7 * t479 - 0.3E1 / 0.494332E6 * t481
        t694 = 0.25E2 / 0.494332E6 * t593
        t697 = 0.2131E4 / 0.741498E6 * t583 + t612
        t711 = -t596 - 0.2131E4 / 0.741498E6 * t585
        t728 = t460 * t583
        t730 = t458 * t585
        t750 = DD(-3,0,2)
        t757 = DD(3,0,2)
        t763 = t74 * t483 + t536 * t535 + t577 * t576 + t423 * (t579 * t
     #506 + t581 * t512 + t413 * t587 + t415 * (-0.2259E4 / 0.123583E6 *
     # t583 + t590) + t427 * t594 + t418 * t597 + t417 * t601 + 0.31E2 /
     # 0.741498E6 * t603 - 0.25E2 / 0.494332E6 * t605) + t432 * (t579 * 
     #t549 + t581 * t555 + t413 * t613 + t415 * t594 + t427 * (-t611 + 0
     #.2259E4 / 0.123583E6 * t585) + t418 * t621 + t417 * t625 + 0.25E2 
     #/ 0.494332E6 * t603 - 0.31E2 / 0.741498E6 * t605) + t579 * t661 + 
     #t581 * t686 + t441 * (t413 * (-0.14183E5 / 0.2471660E7 * t583 + 0.
     #3E1 / 0.494332E6 * t585) + t415 * t587 + t427 * t613 + t418 * t694
     # + t417 * t697 - 0.4483E4 / 0.2471660E7 * t603 + 0.3E1 / 0.494332E
     #6 * t605) + t449 * (t413 * t694 + t415 * t597 + t427 * t621 + t418
     # * (-0.3E1 / 0.494332E6 * t583 + 0.14183E5 / 0.2471660E7 * t585) +
     # t417 * t711 - 0.3E1 / 0.494332E6 * t603 + 0.4483E4 / 0.2471660E7 
     #* t605) + t456 * (t413 * t697 + t415 * t601 + t427 * t625 + t418 *
     # t711 - 0.50753E5 / 0.3336741E7 * t417 * t593 - 0.2281E4 / 0.66734
     #82E7 * t603 + 0.2281E4 / 0.6673482E7 * t605) + t413 * (-0.4483E4 /
     # 0.2471660E7 * t728 + 0.3E1 / 0.494332E6 * t730) + t415 * (0.31E2 
     #/ 0.741498E6 * t728 - 0.25E2 / 0.494332E6 * t730) + t427 * (0.25E2
     # / 0.494332E6 * t728 - 0.31E2 / 0.741498E6 * t730) + t418 * (-0.3E
     #1 / 0.494332E6 * t728 + 0.4483E4 / 0.2471660E7 * t730) + t417 * (-
     #0.2281E4 / 0.6673482E7 * t728 + 0.2281E4 / 0.6673482E7 * t730) + t
     #583 * (-0.1048576E7 / 0.16683705E8 * t750 * vy(-3,0) - 0.150421E6 
     #/ 0.66734820E8 * t481) + 0.150421E6 / 0.66734820E8 * t585 * (t479 
     #+ 0.4194304E7 / 0.150421E6 * t757 * vy(3,0))
        t766 = f(-1,3)
        t768 = f(-1,-3)
        t771 = 0.16777216E8 / 0.50051115E8 * f(-1,2)
        t773 = 0.16777216E8 / 0.50051115E8 * f(-1,-2)
        t775 = 0.19922944E8 / 0.10010223E8 * f(-1,-1)
        t777 = 0.19922944E8 / 0.10010223E8 * f(-1,1)
        t780 = f(3,-1)
        t782 = f(-3,-1)
        t785 = 0.19922944E8 / 0.10010223E8 * f(1,-1)
        t787 = 0.16777216E8 / 0.50051115E8 * f(2,-1)
        t789 = 0.16777216E8 / 0.50051115E8 * f(-2,-1)
        t792 = f(3,1)
        t794 = f(-3,1)
        t797 = 0.19922944E8 / 0.10010223E8 * f(1,1)
        t799 = 0.16777216E8 / 0.50051115E8 * f(2,1)
        t801 = 0.16777216E8 / 0.50051115E8 * f(-2,1)
        t804 = f(1,3)
        t806 = f(1,-3)
        t809 = 0.16777216E8 / 0.50051115E8 * f(1,2)
        t811 = 0.16777216E8 / 0.50051115E8 * f(1,-2)
        t815 = 0.262144E6 / 0.10010223E8 * f(-2,2)
        t817 = 0.262144E6 / 0.10010223E8 * f(-2,-2)
        t821 = 0.262144E6 / 0.10010223E8 * f(2,-2)
        t825 = 0.262144E6 / 0.10010223E8 * f(2,2)
        t842 = t763 * t352 + t530 * (-0.524288E6 / 0.16683705E8 * t766 +
     # 0.524288E6 / 0.16683705E8 * t768 + t771 - t773 + t775 - t777) + t
     #120 * (-0.524288E6 / 0.16683705E8 * t780 + 0.524288E6 / 0.16683705
     #E8 * t782 + t775 - t785 + t787 - t789) + t161 * (0.524288E6 / 0.16
     #683705E8 * t792 - 0.524288E6 / 0.16683705E8 * t794 - t777 + t797 -
     # t799 + t801) + t571 * (0.524288E6 / 0.16683705E8 * t804 - 0.52428
     #8E6 / 0.16683705E8 * t806 - t809 + t811 - t785 + t797) + t656 * (-
     #t815 + t817 - t789 + t801) + t248 * (-t821 + t817 - t773 + t811) +
     # t273 * (t825 - t815 + t771 - t809) + t681 * (t825 - t821 + t787 -
     # t799) + t750 * (0.524288E6 / 0.16683705E8 * t782 - 0.524288E6 / 0
     #.16683705E8 * t794) + t339 * (0.524288E6 / 0.16683705E8 * t768 - 0
     #.524288E6 / 0.16683705E8 * t806) + t346 * (-0.524288E6 / 0.1668370
     #5E8 * t766 + 0.524288E6 / 0.16683705E8 * t804) - 0.524288E6 / 0.16
     #683705E8 * (t780 - t792) * t757
        t846 = dvy ** 2
        t880 = 0.48E2 / 0.1045E4 * t456
        t896 = t74 * (-0.3E1 / 0.209E3 * t460 + 0.232E3 / 0.1881E4 * t44
     #1 - 0.1541E4 / 0.1881E4 * t423 - 0.400E3 / 0.209E3 * t456 - 0.1541
     #E4 / 0.1881E4 * t432 + 0.232E3 / 0.1881E4 * t449 - 0.3E1 / 0.209E3
     # * t458) + t536 * (0.7E1 / 0.627E3 * t460 - 0.58E2 / 0.627E3 * t44
     #1 + t423 + t456 - 0.58E2 / 0.627E3 * t432 + 0.7E1 / 0.627E3 * t449
     #) + t577 * (0.7E1 / 0.627E3 * t441 - 0.58E2 / 0.627E3 * t423 + t45
     #6 + t432 - 0.58E2 / 0.627E3 * t449 + 0.7E1 / 0.627E3 * t458) + t42
     #3 * (-0.58E2 / 0.627E3 * t579 + t581 / 0.3135E4 + 0.37E2 / 0.9405E
     #4 * t583) + t432 * (t579 / 0.3135E4 - 0.58E2 / 0.627E3 * t581 + 0.
     #37E2 / 0.9405E4 * t585) + t579 * (t460 / 0.3135E4 - 0.48E2 / 0.104
     #5E4 * t441 - t880) + t581 * (-t880 - 0.48E2 / 0.1045E4 * t449 + t4
     #58 / 0.3135E4) + 0.37E2 / 0.9405E4 * t441 * t583 + 0.37E2 / 0.9405
     #E4 * t449 * t585 + t456 * (0.3E1 / 0.1045E4 * t583 + 0.3E1 / 0.104
     #5E4 * t585) + 0.3E1 / 0.1045E4 * t728 + 0.3E1 / 0.1045E4 * t730
        cg0 = -0.50753E5 / 0.12582912E8 / t409 / t846 * nu * (t409 * t40
     #8 - 0.3336741E7 / 0.50753E5 * dvx * dvy * t842 - 0.10272768E8 / 0.
     #50753E5 * t896 * t846)
ccccc
ccccc
        rhs(i1,i2,i3,i4) = rhs(i1,i2,i3,i4)+cg0
      end do
      end do
      end do
      end do
c
      return
      end
c
c +++++++++++++
c      
      subroutine appendRosenbluthCollision_6thMapleLog(
     *     rhs,fin,velocities,IVx,IVy,din,Vth,N,massR,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     dx,
     *     dparams,
     *     iparams)
c
c.. compute and add the Rosenbluth collision operator to the RHS ...
      implicit none
c
c.. declarations of incoming variables 
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
c
      real dparams(*)
      integer iparams(*)
c
      real rhs(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real fin(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real din(n1a:n1b, n2a:n2b, n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2), 1:3)
      real IVx(nd1a:nd1b, nd2a:nd2b)
      real IVy(nd1a:nd1b, nd2a:nd2b)
      real Vth(nd1a:nd1b, nd2a:nd2b)
      real N(nd1a:nd1b, nd2a:nd2b)
      real massR
c     
      real dx(1:4)
c
c.. declarations of local variables
      real nu_base, nu
      integer i1, i2, i3, i4, i, j, modulate_nu
      real cg
c
      real dvx, dvy, alpha
      real vx(-3:3,-3:3), vy(-3:3,-3:3)
      real f(-3:3,-3:3)
      real D(-3:3,-3:3,1:3)
c
      real t1
      real t10
      real t100
      real t101
      real t102
      real t105
      real t106
      real t11
      real t112
      real t113
      real t119
      real t120
      real t126
      real t127
      real t132
      real t133
      real t136
      real t137
      real t142
      real t143
      real t146
      real t147
      real t148
      real t15
      real t151
      real t152
      real t157
      real t158
      real t16
      real t161
      real t162
      real t163
      real t164
      real t167
      real t168
      real t17
      real t174
      real t175
      real t181
      real t182
      real t188
      real t189
      real t195
      real t196
      real t199
      real t2
      real t20
      real t200
      real t201
      real t202
      real t203
      real t204
      real t205
      real t206
      real t208
      real t209
      real t21
      real t213
      real t214
      real t217
      real t218
      real t219
      real t220
      real t222
      real t223
      real t229
      real t230
      real t237
      real t238
      real t241
      real t242
      real t243
      real t244
      real t245
      real t246
      real t247
      real t249
      real t25
      real t250
      real t254
      real t255
      real t258
      real t259
      real t26
      real t260
      real t261
      real t263
      real t264
      real t268
      real t269
      real t27
      real t272
      real t277
      real t278
      real t285
      real t286
      real t294
      real t295
      real t3
      real t30
      real t304
      real t305
      real t31
      real t311
      real t314
      real t315
      real t323
      real t324
      real t333
      real t334
      real t35
      real t355
      real t356
      real t359
      real t36
      real t361
      real t363
      real t366
      real t369
      real t37
      real t372
      real t375
      real t378
      real t381
      real t384
      real t387
      real t390
      real t393
      real t396
      real t398
      real t4
      real t40
      real t41
      real t412
      real t426
      real t45
      real t46
      real t47
      real t477
      real t479
      real t481
      real t484
      real t487
      real t490
      real t493
      real t496
      real t499
      real t5
      real t50
      real t502
      real t505
      real t507
      real t508
      real t51
      real t512
      real t514
      real t517
      real t520
      real t523
      real t526
      real t528
      real t531
      real t534
      real t537
      real t539
      real t54
      real t544
      real t546
      real t549
      real t55
      real t552
      real t555
      real t558
      real t56
      real t560
      real t563
      real t566
      real t569
      real t57
      real t571
      real t576
      real t578
      real t58
      real t583
      real t586
      real t588
      real t59
      real t590
      real t595
      real t598
      real t599
      real t6
      real t603
      real t605
      real t609
      real t611
      real t619
      real t62
      real t623
      real t63
      real t633
      real t637
      real t661
      real t665
      real t668
      real t671
      real t674
      real t677
      real t68
      real t680
      real t683
      real t686
      real t689
      real t69
      real t692
      real t695
      real t698
      real t7
      real t700
      real t714
      real t727
      real t75
      real t76
      real t779
      real t82
      real t83
      real t89
      real t90
      real t96
      real t97
c
      modulate_nu = iparams(4)
      nu_base = dparams(1)
c     
      dvx = dx(3)
      dvy = dx(4)
c     
      do i4 = n4a, n4b
      do i3 = n3a, n3b
      do i2 = n2a, n2b
      do i1 = n1a, n1b
        if (modulate_nu .eq. 1) then
          nu = nu_base*N(i1,i2)/Vth(i1,i2)
        else
          nu = nu_base
        end if
        alpha = 1.0 / (massR * (Vth(i1,i2)**2))
        do i = -3, 3
          do j = -3, 3
            vx(i,j) = velocities(i3+i,i4+j, 0) - IVx(i1,i2)
            vy(i,j) = velocities(i3+i,i4+j, 1) - IVy(i1,i2)
            f(i,j) = fin(i1,i2,i3+i,i4+j)
            D(i,j,1) = din(i1,i2,i3+i,i4+j,1)
            D(i,j,2) = din(i1,i2,i3+i,i4+j,2)
            D(i,j,3) = din(i1,i2,i3+i,i4+j,3)
          end do
        end do
ccccc
cc generated code 
ccccc
        t1 = dvy * dvx
        t2 = f(-1,0)
        t3 = D(-1,0,2)
        t4 = t3 * t2
        t5 = f(0,-2)
        t6 = D(0,-2,2)
        t7 = t6 * t5
        t10 = f(-1,-2) ** 2
        t11 = log(t10)
        t15 = f(0,-1)
        t16 = D(0,-1,2)
        t17 = t16 * t15
        t20 = f(-1,-1) ** 2
        t21 = log(t20)
        t25 = f(0,1)
        t26 = D(0,1,2)
        t27 = t26 * t25
        t30 = f(-1,1) ** 2
        t31 = log(t30)
        t35 = f(0,2)
        t36 = D(0,2,2)
        t37 = t36 * t35
        t40 = f(-1,2) ** 2
        t41 = log(t40)
        t45 = f(0,3)
        t46 = D(0,3,2)
        t47 = t46 * t45
        t50 = f(-1,3) ** 2
        t51 = log(t50)
        t54 = f(0,-3)
        t55 = D(0,-3,2)
        t56 = t55 * t54
        t57 = f(1,0)
        t58 = D(1,0,2)
        t59 = t58 * t57
        t62 = f(1,-3) ** 2
        t63 = log(t62)
        t68 = f(1,-2) ** 2
        t69 = log(t68)
        t75 = f(1,-1) ** 2
        t76 = log(t75)
        t82 = f(1,1) ** 2
        t83 = log(t82)
        t89 = f(1,2) ** 2
        t90 = log(t89)
        t96 = f(1,3) ** 2
        t97 = log(t96)
        t100 = D(2,0,2)
        t101 = f(2,0)
        t102 = t101 * t100
        t105 = f(2,-2) ** 2
        t106 = log(t105)
        t112 = f(2,-1) ** 2
        t113 = log(t112)
        t119 = f(2,1) ** 2
        t120 = log(t119)
        t126 = f(2,2) ** 2
        t127 = log(t126)
        t132 = f(3,0)
        t133 = t132 * D(3,0,2)
        t136 = f(3,-1) ** 2
        t137 = log(t136)
        t142 = f(3,1) ** 2
        t143 = log(t142)
        t146 = f(-3,0)
        t147 = D(-3,0,2)
        t148 = t147 * t146
        t151 = f(-3,-1) ** 2
        t152 = log(t151)
        t157 = f(-3,1) ** 2
        t158 = log(t157)
        t161 = -0.32E2 / 0.3E1 * t11 * (t4 + t7) * t1 + 0.190E3 / 0.3E1 
     #* t21 * (t4 + t17) * t1 - 0.190E3 / 0.3E1 * t31 * (t4 + t27) * t1 
     #+ 0.32E2 / 0.3E1 * t41 * (t4 + t37) * t1 - t51 * (t4 + t47) * t1 -
     # t63 * (t56 + t59) * t1 + 0.32E2 / 0.3E1 * t69 * (t7 + t59) * t1 -
     # 0.190E3 / 0.3E1 * t76 * (t17 + t59) * t1 + 0.190E3 / 0.3E1 * t83 
     #* (t27 + t59) * t1 - 0.32E2 / 0.3E1 * t90 * (t37 + t59) * t1 + t97
     # * (t47 + t59) * t1 - 0.5E1 / 0.6E1 * t106 * (t7 + t102) * t1 + 0.
     #32E2 / 0.3E1 * t113 * (t17 + t102) * t1 - 0.32E2 / 0.3E1 * t120 * 
     #(t27 + t102) * t1 + 0.5E1 / 0.6E1 * t127 * (t37 + t102) * t1 - t13
     #7 * (t17 + t133) * t1 + t143 * (t27 + t133) * t1 + t152 * (t148 + 
     #t17) * t1 - t158 * (t148 + t27) * t1
        t162 = f(-2,0)
        t163 = D(-2,0,2)
        t164 = t163 * t162
        t167 = f(-2,-2) ** 2
        t168 = log(t167)
        t174 = f(-2,-1) ** 2
        t175 = log(t174)
        t181 = f(-2,1) ** 2
        t182 = log(t181)
        t188 = f(-2,2) ** 2
        t189 = log(t188)
        t195 = f(-1,-3) ** 2
        t196 = log(t195)
        t199 = dvy ** 2
        t200 = D(0,0,1)
        t201 = f(0,0)
        t202 = t201 * t200
        t203 = D(-3,0,1)
        t204 = t146 * t203
        t205 = D(-2,0,1)
        t206 = t162 * t205
        t208 = D(-1,0,1)
        t209 = t2 * t208
        t213 = t146 ** 2
        t214 = log(t213)
        t217 = D(1,0,1)
        t218 = t57 * t217
        t219 = D(2,0,1)
        t220 = t101 * t219
        t222 = D(3,0,1)
        t223 = t132 * t222
        t229 = t101 ** 2
        t230 = log(t229)
        t237 = t132 ** 2
        t238 = log(t237)
        t241 = dvx ** 2
        t242 = D(0,1,3)
        t243 = t25 * t242
        t244 = D(0,2,3)
        t245 = t35 * t244
        t246 = D(0,3,3)
        t247 = t45 * t246
        t249 = D(0,0,3)
        t250 = t201 * t249
        t254 = t45 ** 2
        t255 = log(t254)
        t258 = D(0,-3,3)
        t259 = t54 * t258
        t260 = D(0,-2,3)
        t261 = t5 * t260
        t263 = D(0,-1,3)
        t264 = t15 * t263
        t268 = t54 ** 2
        t269 = log(t268)
        t272 = 0.144E3 * t250
        t277 = t35 ** 2
        t278 = log(t277)
        t285 = t5 ** 2
        t286 = log(t285)
        t294 = t162 ** 2
        t295 = log(t294)
        t304 = t57 ** 2
        t305 = log(t304)
        t311 = 0.627E3 / 0.7E1 * t250
        t314 = t25 ** 2
        t315 = log(t314)
        t323 = t15 ** 2
        t324 = log(t323)
        t333 = t2 ** 2
        t334 = log(t333)
        t355 = t201 ** 2
        t356 = log(t355)
        t359 = vx(0,0) ** 2
        t361 = vy(0,0) ** 2
        t363 = vx(0,-3) ** 2
        t366 = vx(0,-2) ** 2
        t369 = vx(0,-1) ** 2
        t372 = vx(0,1) ** 2
        t375 = vx(0,2) ** 2
        t378 = vx(0,3) ** 2
        t381 = vy(0,-3) ** 2
        t384 = vy(0,-2) ** 2
        t387 = vy(0,-1) ** 2
        t390 = vy(0,1) ** 2
        t393 = vy(0,2) ** 2
        t396 = vy(0,3) ** 2
        t398 = t359 + t361 - 0.3E1 / 0.2000E4 * t363 + 0.3E1 / 0.125E3 *
     # t366 - 0.209E3 / 0.400E3 * t369 - 0.209E3 / 0.400E3 * t372 + 0.3E
     #1 / 0.125E3 * t375 - 0.3E1 / 0.2000E4 * t378 - 0.3E1 / 0.2000E4 * 
     #t381 + 0.3E1 / 0.125E3 * t384 - 0.209E3 / 0.400E3 * t387 - 0.209E3
     # / 0.400E3 * t390 + 0.3E1 / 0.125E3 * t393 - 0.3E1 / 0.2000E4 * t3
     #96
        t412 = t359 + t361 - 0.37E2 / 0.7705E4 * t363 + 0.174E3 / 0.1541
     #E4 * t366 - 0.1881E4 / 0.1541E4 * t369 + 0.174E3 / 0.1541E4 * t372
     # - 0.3E1 / 0.7705E4 * t375 - 0.37E2 / 0.7705E4 * t381 + 0.174E3 / 
     #0.1541E4 * t384 - 0.1881E4 / 0.1541E4 * t387 + 0.174E3 / 0.1541E4 
     #* t390 - 0.3E1 / 0.7705E4 * t393
        t426 = t359 + t361 - 0.3E1 / 0.7705E4 * t366 + 0.174E3 / 0.1541E
     #4 * t369 - 0.1881E4 / 0.1541E4 * t372 + 0.174E3 / 0.1541E4 * t375 
     #- 0.37E2 / 0.7705E4 * t378 - 0.3E1 / 0.7705E4 * t384 + 0.174E3 / 0
     #.1541E4 * t387 - 0.1881E4 / 0.1541E4 * t390 + 0.174E3 / 0.1541E4 *
     # t393 - 0.37E2 / 0.7705E4 * t396
        t477 = vx(-1,-3) ** 2
        t479 = vy(-1,-3) ** 2
        t481 = vx(-1,-1) ** 2
        t484 = vx(-1,1) ** 2
        t487 = vy(-1,-1) ** 2
        t490 = vy(-1,1) ** 2
        t493 = vx(-1,-2) ** 2
        t496 = vx(-1,2) ** 2
        t499 = vy(-1,-2) ** 2
        t502 = vy(-1,2) ** 2
        t505 = vy(-1,3) ** 2
        t507 = vx(-1,3) ** 2
        t508 = t477 + t479 + 0.190E3 / 0.3E1 * t481 - 0.190E3 / 0.3E1 * 
     #t484 + 0.190E3 / 0.3E1 * t487 - 0.190E3 / 0.3E1 * t490 - 0.32E2 / 
     #0.3E1 * t493 + 0.32E2 / 0.3E1 * t496 - 0.32E2 / 0.3E1 * t499 + 0.3
     #2E2 / 0.3E1 * t502 - t505 - t507
        t512 = vx(1,-1) ** 2
        t514 = vx(2,-1) ** 2
        t517 = vx(3,-1) ** 2
        t520 = vy(-3,-1) ** 2
        t523 = vy(-2,-1) ** 2
        t526 = vy(1,-1) ** 2
        t528 = vy(2,-1) ** 2
        t531 = vy(3,-1) ** 2
        t534 = vx(-3,-1) ** 2
        t537 = vx(-2,-1) ** 2
        t539 = t481 - t512 + 0.16E2 / 0.95E2 * t514 - 0.3E1 / 0.190E3 * 
     #t517 + 0.3E1 / 0.190E3 * t520 - 0.16E2 / 0.95E2 * t523 + t487 - t5
     #26 + 0.16E2 / 0.95E2 * t528 - 0.3E1 / 0.190E3 * t531 + 0.3E1 / 0.1
     #90E3 * t534 - 0.16E2 / 0.95E2 * t537
        t544 = vx(1,1) ** 2
        t546 = vx(2,1) ** 2
        t549 = vx(3,1) ** 2
        t552 = vy(-3,1) ** 2
        t555 = vy(-2,1) ** 2
        t558 = vy(1,1) ** 2
        t560 = vy(2,1) ** 2
        t563 = vy(3,1) ** 2
        t566 = vx(-3,1) ** 2
        t569 = vx(-2,1) ** 2
        t571 = t484 - t544 + 0.16E2 / 0.95E2 * t546 - 0.3E1 / 0.190E3 * 
     #t549 + 0.3E1 / 0.190E3 * t552 - 0.16E2 / 0.95E2 * t555 + t490 - t5
     #58 + 0.16E2 / 0.95E2 * t560 - 0.3E1 / 0.190E3 * t563 + 0.3E1 / 0.1
     #90E3 * t566 - 0.16E2 / 0.95E2 * t569
        t576 = vx(1,-3) ** 2
        t578 = vx(1,-2) ** 2
        t583 = vx(1,2) ** 2
        t586 = vx(1,3) ** 2
        t588 = vy(1,-3) ** 2
        t590 = vy(1,-2) ** 2
        t595 = vy(1,2) ** 2
        t598 = vy(1,3) ** 2
        t599 = t576 - 0.32E2 / 0.3E1 * t578 + 0.190E3 / 0.3E1 * t512 - 0
     #.190E3 / 0.3E1 * t544 + 0.32E2 / 0.3E1 * t583 - t586 + t588 - 0.32
     #E2 / 0.3E1 * t590 + 0.190E3 / 0.3E1 * t526 - 0.190E3 / 0.3E1 * t55
     #8 + 0.32E2 / 0.3E1 * t595 - t598
        t603 = vx(-2,2) ** 2
        t605 = vy(-2,-2) ** 2
        t609 = vy(-2,2) ** 2
        t611 = vx(-2,-2) ** 2
        t619 = vx(2,-2) ** 2
        t623 = vy(2,-2) ** 2
        t633 = vx(2,2) ** 2
        t637 = vy(2,2) ** 2
        t661 = t2 * t508 * t3 + 0.190E3 / 0.3E1 * t15 * t16 * t539 - 0.1
     #90E3 / 0.3E1 * t25 * t571 * t26 - t57 * t58 * t599 - 0.5E1 / 0.6E1
     # * t162 * t163 * (t603 - t605 + 0.64E2 / 0.5E1 * t523 - 0.64E2 / 0
     #.5E1 * t555 + t609 - t611 + 0.64E2 / 0.5E1 * t537 - 0.64E2 / 0.5E1
     # * t569) - 0.32E2 / 0.3E1 * t5 * t6 * (t493 - t578 + 0.5E1 / 0.64E
     #2 * t619 - 0.5E1 / 0.64E2 * t605 + t499 - t590 + 0.5E1 / 0.64E2 * 
     #t623 - 0.5E1 / 0.64E2 * t611) - 0.5E1 / 0.6E1 * t35 * (t603 - 0.64
     #E2 / 0.5E1 * t496 + 0.64E2 / 0.5E1 * t583 - t633 + t609 - 0.64E2 /
     # 0.5E1 * t502 + 0.64E2 / 0.5E1 * t595 - t637) * t36 - 0.5E1 / 0.6E
     #1 * t101 * t100 * (t619 - 0.64E2 / 0.5E1 * t514 + 0.64E2 / 0.5E1 *
     # t546 - t633 + t623 - 0.64E2 / 0.5E1 * t528 + 0.64E2 / 0.5E1 * t56
     #0 - t637) + t146 * (t534 - t566 + t520 - t552) * t147 + t54 * (t47
     #7 - t576 + t479 - t588) * t55 - t45 * (t507 - t586 + t505 - t598) 
     #* t46 - (t517 - t549 + t531 - t563) * t133
        t665 = vx(-3,0) ** 2
        t668 = vx(-2,0) ** 2
        t671 = vx(-1,0) ** 2
        t674 = vx(1,0) ** 2
        t677 = vx(2,0) ** 2
        t680 = vx(3,0) ** 2
        t683 = vy(-3,0) ** 2
        t686 = vy(-2,0) ** 2
        t689 = vy(-1,0) ** 2
        t692 = vy(1,0) ** 2
        t695 = vy(2,0) ** 2
        t698 = vy(3,0) ** 2
        t700 = t359 + t361 - 0.3E1 / 0.2000E4 * t665 + 0.3E1 / 0.125E3 *
     # t668 - 0.209E3 / 0.400E3 * t671 - 0.209E3 / 0.400E3 * t674 + 0.3E
     #1 / 0.125E3 * t677 - 0.3E1 / 0.2000E4 * t680 - 0.3E1 / 0.2000E4 * 
     #t683 + 0.3E1 / 0.125E3 * t686 - 0.209E3 / 0.400E3 * t689 - 0.209E3
     # / 0.400E3 * t692 + 0.3E1 / 0.125E3 * t695 - 0.3E1 / 0.2000E4 * t6
     #98
        t714 = t359 + t361 - 0.37E2 / 0.7705E4 * t683 + 0.174E3 / 0.1541
     #E4 * t668 - 0.1881E4 / 0.1541E4 * t671 + 0.174E3 / 0.1541E4 * t674
     # - 0.3E1 / 0.7705E4 * t677 + 0.174E3 / 0.1541E4 * t686 - 0.1881E4 
     #/ 0.1541E4 * t689 + 0.174E3 / 0.1541E4 * t692 - 0.3E1 / 0.7705E4 *
     # t695 - 0.37E2 / 0.7705E4 * t665
        t727 = t359 + t361 - 0.3E1 / 0.7705E4 * t668 + 0.174E3 / 0.1541E
     #4 * t671 - 0.1881E4 / 0.1541E4 * t674 + 0.174E3 / 0.1541E4 * t677 
     #- 0.3E1 / 0.7705E4 * t686 + 0.174E3 / 0.1541E4 * t689 - 0.1881E4 /
     # 0.1541E4 * t692 + 0.174E3 / 0.1541E4 * t695 - 0.37E2 / 0.7705E4 *
     # t680 - 0.37E2 / 0.7705E4 * t698
        t779 = 0.5E1 / 0.6E1 * t168 * (t164 + t7) * t1 - 0.32E2 / 0.3E1 
     #* t175 * (t164 + t17) * t1 + 0.32E2 / 0.3E1 * t182 * (t164 + t27) 
     #* t1 - 0.5E1 / 0.6E1 * t189 * (t164 + t37) * t1 + t196 * (t4 + t56
     #) * t1 + 0.9E1 / 0.32E2 * t214 * (t202 + t204 + 0.37E2 / 0.27E2 * 
     #t206 + 0.37E2 / 0.27E2 * t209) * t199 - 0.145E3 / 0.16E2 * t230 * 
     #(t218 + 0.72E2 / 0.145E3 * t220 - t223 / 0.290E3 + 0.72E2 / 0.145E
     #3 * t202 - t209 / 0.290E3) * t199 + 0.37E2 / 0.96E2 * t238 * (t218
     # + t220 + 0.27E2 / 0.37E2 * t223 + 0.27E2 / 0.37E2 * t202) * t199 
     #+ 0.37E2 / 0.96E2 * t255 * (t243 + t245 + 0.27E2 / 0.37E2 * t247 +
     # 0.27E2 / 0.37E2 * t250) * t241 + 0.9E1 / 0.32E2 * t269 * (t259 + 
     #0.37E2 / 0.27E2 * t261 + 0.37E2 / 0.27E2 * t264 + t250) * t241 + t
     #278 * (t264 - t272 - 0.290E3 * t243 - 0.144E3 * t245 + t247) * t24
     #1 / 0.32E2 + t286 * (t259 - 0.144E3 * t261 - 0.290E3 * t264 - t272
     # + t243) * t241 / 0.32E2 + t295 * (t204 - 0.144E3 * t206 - 0.290E3
     # * t209 - 0.144E3 * t202 + t218) * t199 / 0.32E2 + 0.3135E4 / 0.32
     #E2 * t305 * (t218 - 0.58E2 / 0.627E3 * t220 + 0.7E1 / 0.627E3 * t2
     #23 + t202 + 0.7E1 / 0.627E3 * t206 - 0.58E2 / 0.627E3 * t209) * t1
     #99 + 0.35E2 / 0.32E2 * t315 * t241 * (t261 - 0.58E2 / 0.7E1 * t264
     # + 0.627E3 / 0.7E1 * t243 - 0.58E2 / 0.7E1 * t245 + t247 + t311) +
     # 0.35E2 / 0.32E2 * t324 * (t259 - 0.58E2 / 0.7E1 * t261 + 0.627E3 
     #/ 0.7E1 * t264 - 0.58E2 / 0.7E1 * t243 + t245 + t311) * t241 - 0.1
     #45E3 / 0.16E2 * t334 * (t218 - 0.7E1 / 0.58E2 * t220 - 0.627E3 / 0
     #.58E2 * t202 - 0.7E1 / 0.58E2 * t204 + t206 - 0.627E3 / 0.58E2 * t
     #209) * t199 + t356 * (t241 * (-0.45E2 / 0.32E2 * t259 - 0.45E2 / 0
     #.32E2 * t247 + 0.145E3 / 0.12E2 * t261 - 0.7705E4 / 0.96E2 * t264 
     #- 0.375E3 / 0.2E1 * t250 - 0.7705E4 / 0.96E2 * t243 + 0.145E3 / 0.
     #12E2 * t245) - 0.7705E4 / 0.96E2 * (t218 - 0.232E3 / 0.1541E4 * t2
     #20 + 0.27E2 / 0.1541E4 * t223 + 0.3600E4 / 0.1541E4 * t202 + 0.27E
     #2 / 0.1541E4 * t204 - 0.232E3 / 0.1541E4 * t206 + t209) * t199) + 
     #(t241 * (-0.375E3 / 0.2E1 * t201 * t249 * t398 - 0.7705E4 / 0.96E2
     # * t15 * t263 * t412 - 0.7705E4 / 0.96E2 * t25 * t426 * t242 + 0.1
     #45E3 / 0.12E2 * t5 * (t359 + t361 + 0.37E2 / 0.1160E4 * t363 - 0.5
     #4E2 / 0.145E3 * t366 - 0.3E1 / 0.4E1 * t369 + 0.21E2 / 0.232E3 * t
     #372 + 0.37E2 / 0.1160E4 * t381 - 0.54E2 / 0.145E3 * t384 - 0.3E1 /
     # 0.4E1 * t387 + 0.21E2 / 0.232E3 * t390) * t260 + 0.145E3 / 0.12E2
     # * t35 * t244 * (t359 + t361 + 0.21E2 / 0.232E3 * t369 - 0.3E1 / 0
     #.4E1 * t372 - 0.54E2 / 0.145E3 * t375 + 0.37E2 / 0.1160E4 * t378 +
     # 0.21E2 / 0.232E3 * t387 - 0.3E1 / 0.4E1 * t390 - 0.54E2 / 0.145E3
     # * t393 + 0.37E2 / 0.1160E4 * t396) - 0.45E2 / 0.32E2 * t54 * (t35
     #9 + t361 - t363 / 0.5E1 - t366 / 0.45E2 - 0.7E1 / 0.9E1 * t369 - t
     #381 / 0.5E1 - t384 / 0.45E2 - 0.7E1 / 0.9E1 * t387) * t258 - 0.45E
     #2 / 0.32E2 * t246 * t45 * (t359 + t361 - 0.7E1 / 0.9E1 * t372 - t3
     #75 / 0.45E2 - t378 / 0.5E1 - 0.7E1 / 0.9E1 * t390 - t393 / 0.45E2 
     #- t396 / 0.5E1)) + dvx * dvy * t661 - 0.7705E4 / 0.96E2 * (0.3600E
     #4 / 0.1541E4 * t201 * t200 * t700 + t2 * t714 * t208 + t57 * t727 
     #* t217 - 0.232E3 / 0.1541E4 * t162 * (t359 + t361 + 0.37E2 / 0.116
     #0E4 * t665 - 0.54E2 / 0.145E3 * t668 - 0.3E1 / 0.4E1 * t671 + 0.21
     #E2 / 0.232E3 * t674 + 0.37E2 / 0.1160E4 * t683 - 0.54E2 / 0.145E3 
     #* t686 - 0.3E1 / 0.4E1 * t689 + 0.21E2 / 0.232E3 * t692) * t205 - 
     #0.232E3 / 0.1541E4 * t101 * t219 * (t359 + t361 + 0.21E2 / 0.232E3
     # * t671 - 0.3E1 / 0.4E1 * t674 - 0.54E2 / 0.145E3 * t677 + 0.37E2 
     #/ 0.1160E4 * t680 + 0.21E2 / 0.232E3 * t689 - 0.3E1 / 0.4E1 * t692
     # - 0.54E2 / 0.145E3 * t695 + 0.37E2 / 0.1160E4 * t698) + 0.27E2 / 
     #0.1541E4 * t146 * (t359 + t361 - t665 / 0.5E1 - t668 / 0.45E2 - 0.
     #7E1 / 0.9E1 * t671 - t683 / 0.5E1 - t686 / 0.45E2 - 0.7E1 / 0.9E1 
     #* t689) * t203 + 0.27E2 / 0.1541E4 * t222 * t132 * (t359 + t361 - 
     #0.7E1 / 0.9E1 * t674 - t677 / 0.45E2 - t680 / 0.5E1 - 0.7E1 / 0.9E
     #1 * t692 - t695 / 0.45E2 - t698 / 0.5E1)) * t199) * alpha
        cg = 0.1E1 / t241 / t199 * nu * (t161 + t779) / 0.240E3
ccccc
ccccc
        rhs(i1,i2,i3,i4) = rhs(i1,i2,i3,i4)+cg
      end do
      end do
      end do
      end do
c
      return
      end
c
c +++++++++++++
c      
      subroutine appendRosenbluthCollision_6thMapleOrig(
     *     rhs,fin,velocities,IVx,IVy,din,Vth,N,massR,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     dx,
     *     dparams,
     *     iparams)
c
c.. compute and add the Rosenbluth collision operator to the RHS ...
      implicit none
c
c.. declarations of incoming variables 
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
c
      real dparams(*)
      integer iparams(*)
c
      real rhs(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real fin(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real din(n1a:n1b, n2a:n2b, n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2), 1:3)
      real IVx(nd1a:nd1b, nd2a:nd2b)
      real IVy(nd1a:nd1b, nd2a:nd2b)
      real Vth(nd1a:nd1b, nd2a:nd2b)
      real N(nd1a:nd1b, nd2a:nd2b)
      real massR
c     
      real dx(1:4)
c
c.. declarations of local variables
      real nu_base, nu
      integer i1, i2, i3, i4, i, j, modulate_nu
      real cg0
c
      real dvx, dvy, alpha
      real vx(-3:3,-3:3), vy(-3:3,-3:3)
      real f(-3:3,-3:3)
      real DD(-3:3,-3:3,1:3)
c
      real t1
      real t1000
      real t101
      real t1012
      real t1032
      real t105
      real t107
      real t108
      real t1083
      real t111
      real t114
      real t118
      real t12
      real t121
      real t123
      real t128
      real t13
      real t131
      real t133
      real t135
      real t138
      real t139
      real t14
      real t140
      real t143
      real t146
      real t15
      real t150
      real t153
      real t155
      real t16
      real t161
      real t164
      real t168
      real t17
      real t171
      real t173
      real t178
      real t181
      real t183
      real t185
      real t188
      real t189
      real t190
      real t193
      real t196
      real t2
      real t20
      real t200
      real t203
      real t205
      real t211
      real t214
      real t218
      real t221
      real t223
      real t228
      real t23
      real t231
      real t233
      real t235
      real t238
      real t239
      real t240
      real t243
      real t246
      real t250
      real t253
      real t255
      real t261
      real t264
      real t268
      real t27
      real t271
      real t273
      real t279
      real t282
      real t286
      real t289
      real t291
      real t296
      real t299
      real t30
      real t301
      real t304
      real t305
      real t306
      real t308
      real t311
      real t314
      real t318
      real t32
      real t321
      real t323
      real t329
      real t332
      real t336
      real t339
      real t341
      real t346
      real t349
      real t351
      real t354
      real t355
      real t356
      real t358
      real t361
      real t364
      real t368
      real t371
      real t373
      real t377
      real t38
      real t380
      real t383
      real t387
      real t390
      real t392
      real t4
      real t401
      real t404
      real t408
      real t41
      real t411
      real t413
      real t422
      real t425
      real t429
      real t432
      real t434
      real t438
      real t439
      real t440
      real t443
      real t446
      real t45
      real t450
      real t453
      real t455
      real t461
      real t464
      real t468
      real t471
      real t473
      real t477
      real t478
      real t479
      real t48
      real t482
      real t485
      real t489
      real t492
      real t494
      real t50
      real t500
      real t503
      real t507
      real t510
      real t512
      real t516
      real t517
      real t518
      real t521
      real t524
      real t528
      real t531
      real t533
      real t539
      real t54
      real t542
      real t546
      real t549
      real t551
      real t556
      real t559
      real t56
      real t561
      real t564
      real t565
      real t566
      real t568
      real t573
      real t576
      real t578
      real t581
      real t582
      real t583
      real t585
      real t59
      real t593
      real t596
      real t598
      real t601
      real t602
      real t603
      real t605
      real t609
      real t61
      real t624
      real t627
      real t629
      real t632
      real t633
      real t634
      real t636
      real t64
      real t644
      real t647
      real t649
      real t65
      real t652
      real t653
      real t654
      real t656
      real t66
      real t68
      real t684
      real t694
      real t696
      real t698
      real t699
      real t7
      real t708
      real t709
      real t71
      real t718
      real t724
      real t725
      real t727
      real t728
      real t732
      real t733
      real t734
      real t737
      real t738
      real t739
      real t74
      real t742
      real t743
      real t745
      real t746
      real t747
      real t751
      real t752
      real t755
      real t757
      real t758
      real t761
      real t763
      real t770
      real t774
      real t776
      real t78
      real t786
      real t792
      real t798
      real t800
      real t81
      real t813
      real t814
      real t824
      real t83
      real t841
      real t842
      real t843
      real t845
      real t846
      real t851
      real t854
      real t859
      real t864
      real t866
      real t871
      real t874
      real t875
      real t877
      real t879
      real t880
      real t882
      real t884
      real t89
      real t895
      real t9
      real t901
      real t910
      real t92
      real t937
      real t946
      real t952
      real t96
      real t961
      real t964
      real t975
      real t976
      real t977
      real t980
      real t981
      real t986
      real t99
      real t995
      real t996
      real t999
c
      modulate_nu = iparams(4)
      nu_base = dparams(1)
c     
      dvx = dx(3)
      dvy = dx(4)
c     
      do i4 = n4a, n4b
      do i3 = n3a, n3b
      do i2 = n2a, n2b
      do i1 = n1a, n1b
        if (modulate_nu .eq. 1) then
          nu = nu_base*N(i1,i2)/Vth(i1,i2)
        else
          nu = nu_base
        end if
        alpha = 1.0 / (massR * (Vth(i1,i2)**2))
        do i = -3, 3
          do j = -3, 3
            vx(i,j) = velocities(i3+i,i4+j, 0) - IVx(i1,i2)
            vy(i,j) = velocities(i3+i,i4+j, 1) - IVy(i1,i2)
            f(i,j) = fin(i1,i2,i3+i,i4+j)
            DD(i,j,1) = din(i1,i2,i3+i,i4+j,1)
            DD(i,j,2) = din(i1,i2,i3+i,i4+j,2)
            DD(i,j,3) = din(i1,i2,i3+i,i4+j,3)
          end do
        end do
ccccc
cc generated code 
ccccc
        t1 = 0.1E1 / dvy
        t2 = nu * t1
        t4 = vx(0,2) ** 2
        t7 = exp(-t4 * alpha / 0.2E1)
        t9 = vy(0,2) ** 2
        t12 = exp(-t9 * alpha / 0.2E1)
        t13 = t12 * t7
        t14 = t13 * t2
        t15 = DD(0,2,2)
        t16 = 0.1E1 / dvx
        t17 = t16 * t15
        t20 = vx(1,2) ** 2
        t23 = exp(-t20 * alpha / 0.2E1)
        t27 = vy(1,2) ** 2
        t30 = exp(-t27 * alpha / 0.2E1)
        t32 = 0.1E1 / t30 / t23 * f(1,2)
        t38 = vx(-1,2) ** 2
        t41 = exp(-t38 * alpha / 0.2E1)
        t45 = vy(-1,2) ** 2
        t48 = exp(-t45 * alpha / 0.2E1)
        t50 = 0.1E1 / t48 / t41 * f(-1,2)
        t54 = nu * t16
        t56 = vx(2,0) ** 2
        t59 = exp(-t56 * alpha / 0.2E1)
        t61 = vy(2,0) ** 2
        t64 = exp(-t61 * alpha / 0.2E1)
        t65 = t64 * t59
        t66 = t65 * t54
        t68 = t1 * DD(2,0,2)
        t71 = vx(2,1) ** 2
        t74 = exp(-t71 * alpha / 0.2E1)
        t78 = vy(2,1) ** 2
        t81 = exp(-t78 * alpha / 0.2E1)
        t83 = 0.1E1 / t81 / t74 * f(2,1)
        t89 = vx(2,-1) ** 2
        t92 = exp(-t89 * alpha / 0.2E1)
        t96 = vy(2,-1) ** 2
        t99 = exp(-t96 * alpha / 0.2E1)
        t101 = 0.1E1 / t99 / t92 * f(2,-1)
        t105 = t1 * t16
        t107 = t7 * nu * t105
        t108 = t15 * t12
        t111 = vx(-2,2) ** 2
        t114 = exp(-t111 * alpha / 0.2E1)
        t118 = vy(-2,2) ** 2
        t121 = exp(-t118 * alpha / 0.2E1)
        t123 = 0.1E1 / t121 / t114 * f(-2,2)
        t128 = vx(0,-1) ** 2
        t131 = exp(-t128 * alpha / 0.2E1)
        t133 = t131 * nu * t105
        t135 = vy(0,-1) ** 2
        t138 = exp(-t135 * alpha / 0.2E1)
        t139 = DD(0,-1,2)
        t140 = t139 * t138
        t143 = vx(3,-1) ** 2
        t146 = exp(-t143 * alpha / 0.2E1)
        t150 = vy(3,-1) ** 2
        t153 = exp(-t150 * alpha / 0.2E1)
        t155 = 0.1E1 / t153 / t146 * f(3,-1)
        t161 = vx(-3,-1) ** 2
        t164 = exp(-t161 * alpha / 0.2E1)
        t168 = vy(-3,-1) ** 2
        t171 = exp(-t168 * alpha / 0.2E1)
        t173 = 0.1E1 / t171 / t164 * f(-3,-1)
        t178 = vx(0,-2) ** 2
        t181 = exp(-t178 * alpha / 0.2E1)
        t183 = t181 * nu * t105
        t185 = vy(0,-2) ** 2
        t188 = exp(-t185 * alpha / 0.2E1)
        t189 = DD(0,-2,2)
        t190 = t189 * t188
        t193 = vx(2,-2) ** 2
        t196 = exp(-t193 * alpha / 0.2E1)
        t200 = vy(2,-2) ** 2
        t203 = exp(-t200 * alpha / 0.2E1)
        t205 = 0.1E1 / t203 / t196 * f(2,-2)
        t211 = vx(-2,-2) ** 2
        t214 = exp(-t211 * alpha / 0.2E1)
        t218 = vy(-2,-2) ** 2
        t221 = exp(-t218 * alpha / 0.2E1)
        t223 = 0.1E1 / t221 / t214 * f(-2,-2)
        t228 = vx(0,1) ** 2
        t231 = exp(-t228 * alpha / 0.2E1)
        t233 = t231 * nu * t105
        t235 = vy(0,1) ** 2
        t238 = exp(-t235 * alpha / 0.2E1)
        t239 = DD(0,1,2)
        t240 = t239 * t238
        t243 = vx(3,1) ** 2
        t246 = exp(-t243 * alpha / 0.2E1)
        t250 = vy(3,1) ** 2
        t253 = exp(-t250 * alpha / 0.2E1)
        t255 = 0.1E1 / t253 / t246 * f(3,1)
        t261 = vx(-3,1) ** 2
        t264 = exp(-t261 * alpha / 0.2E1)
        t268 = vy(-3,1) ** 2
        t271 = exp(-t268 * alpha / 0.2E1)
        t273 = 0.1E1 / t271 / t264 * f(-3,1)
        t279 = vx(2,2) ** 2
        t282 = exp(-t279 * alpha / 0.2E1)
        t286 = vy(2,2) ** 2
        t289 = exp(-t286 * alpha / 0.2E1)
        t291 = 0.1E1 / t289 / t282 * f(2,2)
        t296 = vx(0,-3) ** 2
        t299 = exp(-t296 * alpha / 0.2E1)
        t301 = vy(0,-3) ** 2
        t304 = exp(-t301 * alpha / 0.2E1)
        t305 = t304 * t299
        t306 = t305 * t2
        t308 = t16 * DD(0,-3,2)
        t311 = vx(1,-3) ** 2
        t314 = exp(-t311 * alpha / 0.2E1)
        t318 = vy(1,-3) ** 2
        t321 = exp(-t318 * alpha / 0.2E1)
        t323 = 0.1E1 / t321 / t314 * f(1,-3)
        t329 = vx(-1,-3) ** 2
        t332 = exp(-t329 * alpha / 0.2E1)
        t336 = vy(-1,-3) ** 2
        t339 = exp(-t336 * alpha / 0.2E1)
        t341 = 0.1E1 / t339 / t332 * f(-1,-3)
        t346 = vx(0,3) ** 2
        t349 = exp(-t346 * alpha / 0.2E1)
        t351 = vy(0,3) ** 2
        t354 = exp(-t351 * alpha / 0.2E1)
        t355 = t354 * t349
        t356 = t355 * t2
        t358 = t16 * DD(0,3,2)
        t361 = vx(1,3) ** 2
        t364 = exp(-t361 * alpha / 0.2E1)
        t368 = vy(1,3) ** 2
        t371 = exp(-t368 * alpha / 0.2E1)
        t373 = 0.1E1 / t371 / t364 * f(1,3)
        t377 = -0.4E1 / 0.45E2 * t32 * t17 * t14 + 0.4E1 / 0.45E2 * t50 
     #* t17 * t14 - 0.4E1 / 0.45E2 * t83 * t68 * t66 + 0.4E1 / 0.45E2 * 
     #t101 * t68 * t66 - t123 * t108 * t107 / 0.144E3 - t155 * t140 * t1
     #33 / 0.120E3 + t173 * t140 * t133 / 0.120E3 - t205 * t190 * t183 /
     # 0.144E3 + t223 * t190 * t183 / 0.144E3 + t255 * t240 * t233 / 0.1
     #20E3 - t273 * t240 * t233 / 0.120E3 + t291 * t108 * t107 / 0.144E3
     # - t323 * t308 * t306 / 0.120E3 + t341 * t308 * t306 / 0.120E3 + t
     #373 * t358 * t356 / 0.120E3
        t380 = vx(-1,3) ** 2
        t383 = exp(-t380 * alpha / 0.2E1)
        t387 = vy(-1,3) ** 2
        t390 = exp(-t387 * alpha / 0.2E1)
        t392 = 0.1E1 / t390 / t383 * f(-1,3)
        t401 = vx(-2,1) ** 2
        t404 = exp(-t401 * alpha / 0.2E1)
        t408 = vy(-2,1) ** 2
        t411 = exp(-t408 * alpha / 0.2E1)
        t413 = 0.1E1 / t411 / t404 * f(-2,1)
        t422 = vx(-2,-1) ** 2
        t425 = exp(-t422 * alpha / 0.2E1)
        t429 = vy(-2,-1) ** 2
        t432 = exp(-t429 * alpha / 0.2E1)
        t434 = 0.1E1 / t432 / t425 * f(-2,-1)
        t438 = t188 * t181
        t439 = t438 * t2
        t440 = t16 * t189
        t443 = vx(1,-2) ** 2
        t446 = exp(-t443 * alpha / 0.2E1)
        t450 = vy(1,-2) ** 2
        t453 = exp(-t450 * alpha / 0.2E1)
        t455 = 0.1E1 / t453 / t446 * f(1,-2)
        t461 = vx(-1,-2) ** 2
        t464 = exp(-t461 * alpha / 0.2E1)
        t468 = vy(-1,-2) ** 2
        t471 = exp(-t468 * alpha / 0.2E1)
        t473 = 0.1E1 / t471 / t464 * f(-1,-2)
        t477 = t238 * t231
        t478 = t477 * t2
        t479 = t16 * t239
        t482 = vx(1,1) ** 2
        t485 = exp(-t482 * alpha / 0.2E1)
        t489 = vy(1,1) ** 2
        t492 = exp(-t489 * alpha / 0.2E1)
        t494 = 0.1E1 / t492 / t485 * f(1,1)
        t500 = vx(-1,1) ** 2
        t503 = exp(-t500 * alpha / 0.2E1)
        t507 = vy(-1,1) ** 2
        t510 = exp(-t507 * alpha / 0.2E1)
        t512 = 0.1E1 / t510 / t503 * f(-1,1)
        t516 = t138 * t131
        t517 = t516 * t2
        t518 = t16 * t139
        t521 = vx(1,-1) ** 2
        t524 = exp(-t521 * alpha / 0.2E1)
        t528 = vy(1,-1) ** 2
        t531 = exp(-t528 * alpha / 0.2E1)
        t533 = 0.1E1 / t531 / t524 * f(1,-1)
        t539 = vx(-1,-1) ** 2
        t542 = exp(-t539 * alpha / 0.2E1)
        t546 = vy(-1,-1) ** 2
        t549 = exp(-t546 * alpha / 0.2E1)
        t551 = 0.1E1 / t549 / t542 * f(-1,-1)
        t556 = vx(1,0) ** 2
        t559 = exp(-t556 * alpha / 0.2E1)
        t561 = vy(1,0) ** 2
        t564 = exp(-t561 * alpha / 0.2E1)
        t565 = t564 * t559
        t566 = t565 * t54
        t568 = t1 * DD(1,0,2)
        t573 = vx(-2,0) ** 2
        t576 = exp(-t573 * alpha / 0.2E1)
        t578 = vy(-2,0) ** 2
        t581 = exp(-t578 * alpha / 0.2E1)
        t582 = t581 * t576
        t583 = t582 * t54
        t585 = t1 * DD(-2,0,2)
        t593 = vx(-1,0) ** 2
        t596 = exp(-t593 * alpha / 0.2E1)
        t598 = vy(-1,0) ** 2
        t601 = exp(-t598 * alpha / 0.2E1)
        t602 = t601 * t596
        t603 = t602 * t54
        t605 = t1 * DD(-1,0,2)
        t609 = -t392 * t358 * t356 / 0.120E3 - 0.4E1 / 0.45E2 * t83 * t2
     #40 * t233 + 0.4E1 / 0.45E2 * t413 * t240 * t233 + 0.4E1 / 0.45E2 *
     # t101 * t140 * t133 - 0.4E1 / 0.45E2 * t434 * t140 * t133 + 0.4E1 
     #/ 0.45E2 * t455 * t440 * t439 - 0.4E1 / 0.45E2 * t473 * t440 * t43
     #9 + 0.19E2 / 0.36E2 * t494 * t479 * t478 - 0.19E2 / 0.36E2 * t512 
     #* t479 * t478 - 0.19E2 / 0.36E2 * t533 * t518 * t517 + 0.19E2 / 0.
     #36E2 * t551 * t518 * t517 - t323 * t568 * t566 / 0.120E3 - t123 * 
     #t585 * t583 / 0.144E3 + t223 * t585 * t583 / 0.144E3 - t392 * t605
     # * t603 / 0.120E3
        t624 = vx(-3,0) ** 2
        t627 = exp(-t624 * alpha / 0.2E1)
        t629 = vy(-3,0) ** 2
        t632 = exp(-t629 * alpha / 0.2E1)
        t633 = t632 * t627
        t634 = t633 * t54
        t636 = t1 * DD(-3,0,2)
        t644 = vx(3,0) ** 2
        t647 = exp(-t644 * alpha / 0.2E1)
        t649 = vy(3,0) ** 2
        t652 = exp(-t649 * alpha / 0.2E1)
        t653 = t652 * t647
        t654 = t653 * t54
        t656 = t1 * DD(3,0,2)
        t684 = t341 * t605 * t603 / 0.120E3 + t291 * t68 * t66 / 0.144E3
     # - t205 * t68 * t66 / 0.144E3 + t373 * t568 * t566 / 0.120E3 - t27
     #3 * t636 * t634 / 0.120E3 + t173 * t636 * t634 / 0.120E3 + t255 * 
     #t656 * t654 / 0.120E3 - t155 * t656 * t654 / 0.120E3 - 0.4E1 / 0.4
     #5E2 * t32 * t568 * t566 + 0.4E1 / 0.45E2 * t455 * t568 * t566 + 0.
     #4E1 / 0.45E2 * t50 * t605 * t603 - 0.4E1 / 0.45E2 * t473 * t605 * 
     #t603 + 0.4E1 / 0.45E2 * t413 * t585 * t583 - 0.4E1 / 0.45E2 * t434
     # * t585 * t583 + 0.19E2 / 0.36E2 * t494 * t568 * t566
        t694 = dvy ** 2
        t696 = nu / t694
        t698 = DD(0,1,3)
        t699 = f(0,3)
        t708 = DD(0,-2,3)
        t709 = f(0,-3)
        t718 = DD(0,2,3)
        t724 = DD(0,-1,3)
        t725 = f(0,-2)
        t734 = vx(0,0) ** 2
        t737 = exp(-t734 * alpha / 0.2E1)
        t739 = vy(0,0) ** 2
        t742 = exp(-t739 * alpha / 0.2E1)
        t743 = t742 * t737
        t745 = DD(0,0,3)
        t751 = DD(0,-3,3)
        t757 = DD(0,3,3)
        t758 = f(0,2)
        t774 = f(0,1)
        t786 = f(0,0)
        t727 = 0.1E1 / t354 / t349 * t699
        t728 = t698 * t477 * t696
        t732 = 0.1E1 / t304 / t299 * t709
        t733 = t708 * t438 * t696
        t738 = t718 * t13 * t696
        t746 = 0.1E1 / t188 / t181 * t725
        t747 = t724 * t516 * t696
        t752 = t745 * t743 * t696
        t755 = t751 * t305 * t696
        t761 = 0.1E1 / t12 / t7 * t758
        t763 = t757 * t355 * t696
        t770 = 0.1E1 / t238 / t231 * t774
        t776 = 0.1E1 / t742 / t737 * t786
        t798 = -0.19E2 / 0.36E2 * t533 * t568 * t566 - 0.19E2 / 0.36E2 *
     # t512 * t605 * t603 + 0.19E2 / 0.36E2 * t551 * t605 * t603 + 0.37E
     #2 / 0.11520E5 * t728 * t727 + 0.37E2 / 0.11520E5 * t733 * t732 + 0
     #.37E2 / 0.11520E5 * t738 * t727 - 0.29E2 / 0.384E3 * t747 * t746 -
     # 0.3E1 / 0.80E2 * t752 * t746 + t755 * t746 / 0.3840E4 + t763 * t7
     #61 / 0.3840E4 - 0.29E2 / 0.384E3 * t728 * t761 - 0.3E1 / 0.80E2 * 
     #t752 * t761 + 0.209E3 / 0.256E3 * t752 * t770 + 0.7E1 / 0.768E3 * 
     #t763 * t770 - 0.3E1 / 0.256E3 * t763 * t776 - 0.29E2 / 0.384E3 * t
     #747 * t770
        t813 = f(0,-1)
        t841 = dvx ** 2
        t842 = 0.1E1 / t841
        t843 = nu * t842
        t845 = DD(0,0,1)
        t846 = f(-3,0)
        t866 = f(3,0)
        t875 = DD(-1,0,1)
        t877 = t875 * t601 * t596 * nu
        t792 = t846 / t632 / t627
        t879 = t842 * t792
        t800 = 0.1E1 / t138 / t131 * t813
        t814 = t845 * t743 * t843
        t824 = 0.1E1 / t652 / t647 * t866
        t882 = -0.1541E4 / 0.2304E4 * t747 * t776 + 0.7E1 / 0.768E3 * t7
     #33 * t770 + 0.29E2 / 0.288E3 * t733 * t776 + 0.209E3 / 0.256E3 * t
     #752 * t800 - 0.29E2 / 0.384E3 * t728 * t800 + 0.7E1 / 0.768E3 * t7
     #38 * t800 - 0.29E2 / 0.384E3 * t733 * t800 - 0.3E1 / 0.256E3 * t75
     #5 * t776 + 0.7E1 / 0.768E3 * t755 * t800 + 0.3E1 / 0.1280E4 * t814
     # * t792 - 0.1541E4 / 0.2304E4 * t728 * t776 - 0.29E2 / 0.384E3 * t
     #738 * t770 + 0.29E2 / 0.288E3 * t738 * t776 + 0.3E1 / 0.1280E4 * t
     #814 * t824 + 0.37E2 / 0.11520E5 * t879 * t877
        t884 = DD(-2,0,1)
        t901 = f(1,0)
        t910 = DD(3,0,1)
        t937 = f(-1,0)
        t946 = DD(1,0,1)
        t952 = DD(2,0,1)
        t851 = 0.1E1 / t564 / t559 * t901
        t854 = t910 * t653 * t843
        t859 = t875 * t602 * t843
        t864 = t884 * t582 * t843
        t871 = 0.1E1 / t601 / t596 * t937
        t874 = t946 * t565 * t843
        t880 = t952 * t65 * t843
        t961 = 0.37E2 / 0.11520E5 * t879 * t884 * t581 * t576 * nu + 0.3
     #E1 / 0.1280E4 * t752 * t732 + 0.37E2 / 0.11520E5 * t747 * t732 + 0
     #.3E1 / 0.1280E4 * t752 * t727 + 0.209E3 / 0.256E3 * t814 * t851 + 
     #0.7E1 / 0.768E3 * t854 * t851 - 0.3E1 / 0.256E3 * t854 * t776 - 0.
     #29E2 / 0.384E3 * t859 * t851 - 0.1541E4 / 0.2304E4 * t859 * t776 +
     # 0.7E1 / 0.768E3 * t864 * t851 + 0.29E2 / 0.288E3 * t864 * t776 + 
     #0.209E3 / 0.256E3 * t814 * t871 - 0.29E2 / 0.384E3 * t874 * t871 +
     # 0.7E1 / 0.768E3 * t880 * t871 - 0.29E2 / 0.384E3 * t864 * t871
        t964 = DD(-3,0,1)
        t975 = t946 * t564 * t559 * nu
        t976 = f(2,0)
        t977 = t976 * t842
        t980 = 0.1E1 / t64 / t59
        t981 = t980 * t977
        t986 = t845 * t742 * t737 * nu
        t995 = f(-2,0)
        t996 = t995 * t842
        t999 = 0.1E1 / t581 / t576
        t1000 = t999 * t996
        t1012 = t842 * t824
        t895 = t964 * t633 * t843
        t1032 = -0.3E1 / 0.256E3 * t895 * t776 + 0.7E1 / 0.768E3 * t895 
     #* t871 - 0.29E2 / 0.384E3 * t981 * t975 - 0.3E1 / 0.80E2 * t981 * 
     #t986 + t981 * t877 / 0.3840E4 + t854 * t980 * t976 / 0.3840E4 - 0.
     #3E1 / 0.80E2 * t1000 * t986 + t1000 * t975 / 0.3840E4 - 0.29E2 / 0
     #.384E3 * t1000 * t877 + t895 * t999 * t995 / 0.3840E4 + 0.37E2 / 0
     #.11520E5 * t1012 * t975 + 0.37E2 / 0.11520E5 * t1012 * t952 * t64 
     #* t59 * nu + t728 * t746 / 0.3840E4 + t747 * t761 / 0.3840E4 - 0.1
     #541E4 / 0.2304E4 * t874 * t776
        t1083 = -0.29E2 / 0.384E3 * t880 * t851 + 0.29E2 / 0.288E3 * t88
     #0 * t776 - 0.3E1 / 0.80E2 * t725 * t708 * t696 + 0.3E1 / 0.1280E4 
     #* t699 * t757 * t696 + 0.3E1 / 0.1280E4 * t846 * t964 * t843 + 0.3
     #E1 / 0.1280E4 * t866 * t910 * t843 + 0.209E3 / 0.256E3 * t901 * t9
     #46 * t843 - 0.25E2 / 0.16E2 * t786 * t845 * t843 + 0.209E3 / 0.256
     #E3 * t937 * t875 * t843 - 0.3E1 / 0.80E2 * t996 * t884 * nu - 0.3E
     #1 / 0.80E2 * t977 * t952 * nu + 0.3E1 / 0.1280E4 * t709 * t751 * t
     #696 + 0.209E3 / 0.256E3 * t774 * t698 * t696 - 0.25E2 / 0.16E2 * t
     #786 * t745 * t696 + 0.209E3 / 0.256E3 * t813 * t724 * t696 - 0.3E1
     # / 0.80E2 * t758 * t718 * t696
        cg0 = t377 + t609 + t684 + t798 + t882 + t961 + t1032 + t1083
ccccc
ccccc
        rhs(i1,i2,i3,i4) = rhs(i1,i2,i3,i4)+cg0
      end do
      end do
      end do
      end do
c
      return
      end
c
c +++++++++++++
c      
      subroutine appendRosenbluthCollision_6thMapleAsinh(
     *     rhs,fin,velocities,IVx,IVy,din,Vth,N,massR,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     dx,
     *     dparams,
     *     iparams)
c
c.. compute and add the Rosenbluth collision operator to the RHS ...
      implicit none
c
c.. declarations of incoming variables 
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
c
      real dparams(*)
      integer iparams(*)
c
      real rhs(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real fin(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real din(n1a:n1b, n2a:n2b, n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2), 1:3)
      real IVx(nd1a:nd1b, nd2a:nd2b)
      real IVy(nd1a:nd1b, nd2a:nd2b)
      real Vth(nd1a:nd1b, nd2a:nd2b)
      real N(nd1a:nd1b, nd2a:nd2b)
      real massR
c     
      real dx(1:4)
c
c.. declarations of local variables
      real nu_base, nu, safelog
      integer i1, i2, i3, i4, i, j, modulate_nu
      real cg0
c
      real dvx, dvy, alpha, beta
      real vx(-3:3,-3:3), vy(-3:3,-3:3)
      real f(-3:3,-3:3)
      real DD(-3:3,-3:3,1:3)
c
      real t10
      real t101
      real t105
      real t108
      real t1086
      real t110
      real t112
      real t114
      real t118
      real t121
      real t123
      real t125
      real t127
      real t13
      real t131
      real t134
      real t136
      real t141
      real t145
      real t148
      real t151
      real t154
      real t157
      real t16
      real t160
      real t163
      real t166
      real t169
      real t172
      real t175
      real t178
      real t180
      real t182
      real t184
      real t188
      real t19
      real t191
      real t193
      real t195
      real t197
      real t2
      real t201
      real t204
      real t206
      real t208
      real t210
      real t214
      real t217
      real t219
      real t22
      real t221
      real t223
      real t227
      real t230
      real t232
      real t234
      real t236
      real t240
      real t243
      real t245
      real t247
      real t249
      real t25
      real t253
      real t256
      real t258
      real t262
      real t269
      real t272
      real t273
      real t275
      real t276
      real t278
      real t279
      real t28
      real t281
      real t284
      real t287
      real t290
      real t291
      real t293
      real t294
      real t296
      real t297
      real t299
      real t302
      real t304
      real t306
      real t308
      real t31
      real t312
      real t315
      real t317
      real t318
      real t320
      real t324
      real t327
      real t329
      real t331
      real t333
      real t337
      real t34
      real t340
      real t342
      real t344
      real t346
      real t350
      real t353
      real t355
      real t356
      real t358
      real t362
      real t365
      real t367
      real t368
      real t369
      real t37
      real t371
      real t375
      real t378
      real t380
      real t381
      real t4
      real t40
      real t43
      real t434
      real t437
      real t440
      real t441
      real t443
      real t446
      real t449
      real t45
      real t452
      real t455
      real t456
      real t458
      real t461
      real t463
      real t465
      real t467
      real t47
      real t471
      real t474
      real t476
      real t477
      real t479
      real t483
      real t486
      real t488
      real t49
      real t490
      real t492
      real t496
      real t499
      real t50
      real t501
      real t503
      real t505
      real t509
      real t512
      real t514
      real t515
      real t517
      real t521
      real t524
      real t526
      real t54
      real t558
      real t561
      real t564
      real t565
      real t567
      real t57
      real t570
      real t573
      real t576
      real t579
      real t580
      real t582
      real t585
      real t587
      real t589
      real t59
      real t591
      real t595
      real t598
      real t6
      real t60
      real t600
      real t601
      real t603
      real t607
      real t610
      real t612
      real t614
      real t616
      real t62
      real t620
      real t623
      real t625
      real t626
      real t628
      real t632
      real t635
      real t637
      real t639
      real t641
      real t645
      real t648
      real t650
      real t66
      real t662
      real t665
      real t670
      real t673
      real t676
      real t679
      real t684
      real t687
      real t689
      real t69
      real t691
      real t693
      real t697
      real t700
      real t702
      real t703
      real t705
      real t709
      real t71
      real t712
      real t714
      real t718
      real t720
      real t724
      real t727
      real t729
      real t73
      real t731
      real t733
      real t737
      real t740
      real t742
      real t75
      real t776
      real t777
      real t779
      real t784
      real t785
      real t787
      real t79
      real t791
      real t793
      real t797
      real t8
      real t800
      real t802
      real t805
      real t807
      real t811
      real t814
      real t816
      real t82
      real t84
      real t86
      real t864
      real t865
      real t869
      real t870
      real t875
      real t877
      real t88
      real t881
      real t884
      real t886
      real t887
      real t92
      real t920
      real t925
      real t929
      real t931
      real t935
      real t938
      real t940
      real t95
      real t97
      real t99
c
      modulate_nu = iparams(4)
      nu_base = dparams(1)
c     
      dvx = dx(3)
      dvy = dx(4)
c     
      do i4 = n4a, n4b
      do i3 = n3a, n3b
      do i2 = n2a, n2b
      do i1 = n1a, n1b
        if (modulate_nu .eq. 1) then
          nu = nu_base*N(i1,i2)/Vth(i1,i2)
        else
          nu = nu_base
        end if
        alpha = 1.0 / (massR * (Vth(i1,i2)**2))
        beta = 2.0*alpha
        do i = -3, 3
          do j = -3, 3
            vx(i,j) = velocities(i3+i,i4+j, 0) - IVx(i1,i2)
            vy(i,j) = velocities(i3+i,i4+j, 1) - IVy(i1,i2)
            f(i,j) = fin(i1,i2,i3+i,i4+j)
            DD(i,j,1) = din(i1,i2,i3+i,i4+j,1)
            DD(i,j,2) = din(i1,i2,i3+i,i4+j,2)
            DD(i,j,3) = din(i1,i2,i3+i,i4+j,3)
          end do
        end do
ccccc
cc generated code 
ccccc
        t2 = vx(0,0) ** 2
        t4 = t2 * alpha / 0.2E1
        t6 = vy(0,0) ** 2
        t8 = t6 * alpha / 0.2E1
        t10 = vx(0,-3) ** 2
        t13 = vx(0,-2) ** 2
        t16 = vx(0,-1) ** 2
        t19 = vx(0,1) ** 2
        t22 = vx(0,2) ** 2
        t25 = vx(0,3) ** 2
        t28 = vy(0,-3) ** 2
        t31 = vy(0,-2) ** 2
        t34 = vy(0,-1) ** 2
        t37 = vy(0,1) ** 2
        t40 = vy(0,2) ** 2
        t43 = vy(0,3) ** 2
        t45 = -0.3E1 / 0.4000E4 * t10 + 0.3E1 / 0.250E3 * t13 - 0.209E3 
     #/ 0.800E3 * t16 - 0.209E3 / 0.800E3 * t19 + 0.3E1 / 0.250E3 * t22 
     #- 0.3E1 / 0.4000E4 * t25 - 0.3E1 / 0.4000E4 * t28 + 0.3E1 / 0.250E
     #3 * t31 - 0.209E3 / 0.800E3 * t34 - 0.209E3 / 0.800E3 * t37 + 0.3E
     #1 / 0.250E3 * t40 - 0.3E1 / 0.4000E4 * t43
        t47 = f(0,0)
        t49 = beta ** 2
        t50 = t47 ** 2
        t54 = exp(-(t2 + t6) * alpha)
        t57 = sqrt(t49 * t50 + 0.4E1 * t54)
        t59 = safelog(beta * t47 + t57)
        t60 = f(0,3)
        t62 = t60 ** 2
        t66 = exp(-(t25 + t43) * alpha)
        t69 = sqrt(t49 * t62 + 0.4E1 * t66)
        t71 = safelog(beta * t60 + t69)
        t73 = f(0,-3)
        t75 = t73 ** 2
        t79 = exp(-(t10 + t28) * alpha)
        t82 = sqrt(t49 * t75 + 0.4E1 * t79)
        t84 = safelog(beta * t73 + t82)
        t86 = f(0,2)
        t88 = t86 ** 2
        t92 = exp(-(t22 + t40) * alpha)
        t95 = sqrt(t49 * t88 + 0.4E1 * t92)
        t97 = safelog(beta * t86 + t95)
        t99 = f(0,-2)
        t101 = t99 ** 2
        t105 = exp(-(t13 + t31) * alpha)
        t108 = sqrt(t101 * t49 + 0.4E1 * t105)
        t110 = safelog(beta * t99 + t108)
        t112 = f(0,1)
        t114 = t112 ** 2
        t118 = exp(-(t19 + t37) * alpha)
        t121 = sqrt(t114 * t49 + 0.4E1 * t118)
        t123 = safelog(beta * t112 + t121)
        t125 = f(0,-1)
        t127 = t125 ** 2
        t131 = exp(-(t16 + t34) * alpha)
        t134 = sqrt(t127 * t49 + 0.4E1 * t131)
        t136 = safelog(beta * t125 + t134)
        t141 = dvx ** 2
        t145 = vx(-3,0) ** 2
        t148 = vx(-2,0) ** 2
        t151 = vx(-1,0) ** 2
        t154 = vx(1,0) ** 2
        t157 = vx(2,0) ** 2
        t160 = vx(3,0) ** 2
        t163 = vy(-3,0) ** 2
        t166 = vy(-2,0) ** 2
        t169 = vy(-1,0) ** 2
        t172 = vy(1,0) ** 2
        t175 = vy(2,0) ** 2
        t178 = vy(3,0) ** 2
        t180 = -0.3E1 / 0.4000E4 * t145 + 0.3E1 / 0.250E3 * t148 - 0.209
     #E3 / 0.800E3 * t151 - 0.209E3 / 0.800E3 * t154 + 0.3E1 / 0.250E3 *
     # t157 - 0.3E1 / 0.4000E4 * t160 - 0.3E1 / 0.4000E4 * t163 + 0.3E1 
     #/ 0.250E3 * t166 - 0.209E3 / 0.800E3 * t169 - 0.209E3 / 0.800E3 * 
     #t172 + 0.3E1 / 0.250E3 * t175 - 0.3E1 / 0.4000E4 * t178
        t182 = f(3,0)
        t184 = t182 ** 2
        t188 = exp(-(t160 + t178) * alpha)
        t191 = sqrt(t184 * t49 + 0.4E1 * t188)
        t193 = safelog(beta * t182 + t191)
        t195 = f(-3,0)
        t197 = t195 ** 2
        t201 = exp(-(t145 + t163) * alpha)
        t204 = sqrt(t197 * t49 + 0.4E1 * t201)
        t206 = safelog(beta * t195 + t204)
        t208 = f(2,0)
        t210 = t208 ** 2
        t214 = exp(-(t157 + t175) * alpha)
        t217 = sqrt(t210 * t49 + 0.4E1 * t214)
        t219 = safelog(beta * t208 + t217)
        t221 = f(-2,0)
        t223 = t221 ** 2
        t227 = exp(-(t148 + t166) * alpha)
        t230 = sqrt(t223 * t49 + 0.4E1 * t227)
        t232 = safelog(beta * t221 + t230)
        t234 = f(1,0)
        t236 = t234 ** 2
        t240 = exp(-(t154 + t172) * alpha)
        t243 = sqrt(t236 * t49 + 0.4E1 * t240)
        t245 = safelog(beta * t234 + t243)
        t247 = f(-1,0)
        t249 = t247 ** 2
        t253 = exp(-(t151 + t169) * alpha)
        t256 = sqrt(t249 * t49 + 0.4E1 * t253)
        t258 = safelog(beta * t247 + t256)
        t262 = dvy ** 2
        t269 = vy(-1,-3) ** 2
        t272 = vy(-1,-2) ** 2
        t273 = t272 / 0.2E1
        t275 = vy(-1,-1) ** 2
        t276 = 0.95E2 / 0.32E2 * t275
        t278 = vy(-1,1) ** 2
        t279 = 0.95E2 / 0.32E2 * t278
        t281 = vy(-1,2) ** 2
        t284 = vy(-1,3) ** 2
        t287 = vx(-1,-3) ** 2
        t290 = vx(-1,-2) ** 2
        t291 = t290 / 0.2E1
        t293 = vx(-1,-1) ** 2
        t294 = 0.95E2 / 0.32E2 * t293
        t296 = vx(-1,1) ** 2
        t297 = 0.95E2 / 0.32E2 * t296
        t299 = vx(-1,2) ** 2
        t302 = vx(-1,3) ** 2
        t304 = -0.3E1 / 0.64E2 * t269 + t273 - t276 + t279 - t281 / 0.2E
     #1 + 0.3E1 / 0.64E2 * t284 - 0.3E1 / 0.64E2 * t287 + t291 - t294 + 
     #t297 - t299 / 0.2E1 + 0.3E1 / 0.64E2 * t302
        t306 = f(-1,-2)
        t308 = t306 ** 2
        t312 = exp(-(t290 + t272) * alpha)
        t315 = sqrt(t308 * t49 + 0.4E1 * t312)
        t317 = safelog(beta * t306 + t315)
        t318 = f(-1,-3)
        t320 = t318 ** 2
        t324 = exp(-(t287 + t269) * alpha)
        t327 = sqrt(t320 * t49 + 0.4E1 * t324)
        t329 = safelog(beta * t318 + t327)
        t331 = f(-1,3)
        t333 = t331 ** 2
        t337 = exp(-(t302 + t284) * alpha)
        t340 = sqrt(t333 * t49 + 0.4E1 * t337)
        t342 = safelog(beta * t331 + t340)
        t344 = f(-1,2)
        t346 = t344 ** 2
        t350 = exp(-(t299 + t281) * alpha)
        t353 = sqrt(t346 * t49 + 0.4E1 * t350)
        t355 = safelog(beta * t344 + t353)
        t356 = f(-1,1)
        t358 = t356 ** 2
        t362 = exp(-(t296 + t278) * alpha)
        t365 = sqrt(t358 * t49 + 0.4E1 * t362)
        t367 = safelog(beta * t356 + t365)
        t368 = 0.95E2 / 0.16E2 * t367
        t369 = f(-1,-1)
        t371 = t369 ** 2
        t375 = exp(-(t293 + t275) * alpha)
        t378 = sqrt(t371 * t49 + 0.4E1 * t375)
        t380 = safelog(beta * t369 + t378)
        t381 = 0.95E2 / 0.16E2 * t380
        t434 = vy(-2,-1) ** 2
        t437 = vy(1,-1) ** 2
        t440 = vy(2,-1) ** 2
        t441 = t440 / 0.2E1
        t443 = vy(3,-1) ** 2
        t446 = vx(-3,-1) ** 2
        t449 = vx(-2,-1) ** 2
        t452 = vx(1,-1) ** 2
        t455 = vx(2,-1) ** 2
        t456 = t455 / 0.2E1
        t458 = vx(3,-1) ** 2
        t461 = vy(-3,-1) ** 2
        t463 = t434 / 0.2E1 - t276 + 0.95E2 / 0.32E2 * t437 - t441 + 0.3
     #E1 / 0.64E2 * t443 - 0.3E1 / 0.64E2 * t446 + t449 / 0.2E1 - t294 +
     # 0.95E2 / 0.32E2 * t452 - t456 + 0.3E1 / 0.64E2 * t458 - 0.3E1 / 0
     #.64E2 * t461
        t465 = f(-2,-1)
        t467 = t465 ** 2
        t471 = exp(-(t449 + t434) * alpha)
        t474 = sqrt(t467 * t49 + 0.4E1 * t471)
        t476 = safelog(beta * t465 + t474)
        t477 = f(-3,-1)
        t479 = t477 ** 2
        t483 = exp(-(t446 + t461) * alpha)
        t486 = sqrt(t479 * t49 + 0.4E1 * t483)
        t488 = safelog(beta * t477 + t486)
        t490 = f(3,-1)
        t492 = t490 ** 2
        t496 = exp(-(t458 + t443) * alpha)
        t499 = sqrt(t49 * t492 + 0.4E1 * t496)
        t501 = safelog(beta * t490 + t499)
        t503 = f(2,-1)
        t505 = t503 ** 2
        t509 = exp(-(t455 + t440) * alpha)
        t512 = sqrt(t49 * t505 + 0.4E1 * t509)
        t514 = safelog(beta * t503 + t512)
        t515 = f(1,-1)
        t517 = t515 ** 2
        t521 = exp(-(t452 + t437) * alpha)
        t524 = sqrt(t49 * t517 + 0.4E1 * t521)
        t526 = safelog(beta * t515 + t524)
        t558 = vy(-2,1) ** 2
        t561 = vy(1,1) ** 2
        t564 = vy(2,1) ** 2
        t565 = t564 / 0.2E1
        t567 = vy(3,1) ** 2
        t570 = vx(-3,1) ** 2
        t573 = vx(-2,1) ** 2
        t576 = vx(1,1) ** 2
        t579 = vx(2,1) ** 2
        t580 = t579 / 0.2E1
        t582 = vx(3,1) ** 2
        t585 = vy(-3,1) ** 2
        t587 = t558 / 0.2E1 - t279 + 0.95E2 / 0.32E2 * t561 - t565 + 0.3
     #E1 / 0.64E2 * t567 - 0.3E1 / 0.64E2 * t570 + t573 / 0.2E1 - t297 +
     # 0.95E2 / 0.32E2 * t576 - t580 + 0.3E1 / 0.64E2 * t582 - 0.3E1 / 0
     #.64E2 * t585
        t589 = f(-2,1)
        t591 = t589 ** 2
        t595 = exp(-(t573 + t558) * alpha)
        t598 = sqrt(t49 * t591 + 0.4E1 * t595)
        t600 = safelog(beta * t589 + t598)
        t601 = f(-3,1)
        t603 = t601 ** 2
        t607 = exp(-(t570 + t585) * alpha)
        t610 = sqrt(t49 * t603 + 0.4E1 * t607)
        t612 = safelog(beta * t601 + t610)
        t614 = f(2,1)
        t616 = t614 ** 2
        t620 = exp(-(t579 + t564) * alpha)
        t623 = sqrt(t49 * t616 + 0.4E1 * t620)
        t625 = safelog(beta * t614 + t623)
        t626 = f(3,1)
        t628 = t626 ** 2
        t632 = exp(-(t582 + t567) * alpha)
        t635 = sqrt(t49 * t628 + 0.4E1 * t632)
        t637 = safelog(beta * t626 + t635)
        t639 = f(1,1)
        t641 = t639 ** 2
        t645 = exp(-(t576 + t561) * alpha)
        t648 = sqrt(t49 * t641 + 0.4E1 * t645)
        t650 = safelog(beta * t639 + t648)
        t662 = vy(1,-3) ** 2
        t665 = vy(1,-2) ** 2
        t670 = vy(1,2) ** 2
        t673 = vy(1,3) ** 2
        t676 = vx(1,-3) ** 2
        t679 = vx(1,-2) ** 2
        t684 = vx(1,2) ** 2
        t687 = vx(1,3) ** 2
        t689 = -t662 / 0.2E1 + 0.16E2 / 0.3E1 * t665 - 0.95E2 / 0.3E1 * 
     #t437 + 0.95E2 / 0.3E1 * t561 - 0.16E2 / 0.3E1 * t670 + t673 / 0.2E
     #1 - t676 / 0.2E1 + 0.16E2 / 0.3E1 * t679 - 0.95E2 / 0.3E1 * t452 +
     # 0.95E2 / 0.3E1 * t576 - 0.16E2 / 0.3E1 * t684 + t687 / 0.2E1
        t691 = f(1,3)
        t693 = t691 ** 2
        t697 = exp(-(t687 + t673) * alpha)
        t700 = sqrt(t49 * t693 + 0.4E1 * t697)
        t702 = safelog(beta * t691 + t700)
        t703 = f(1,2)
        t705 = t703 ** 2
        t709 = exp(-(t684 + t670) * alpha)
        t712 = sqrt(t49 * t705 + 0.4E1 * t709)
        t714 = safelog(beta * t703 + t712)
        t718 = f(1,-2)
        t720 = t718 ** 2
        t724 = exp(-(t679 + t665) * alpha)
        t727 = sqrt(t49 * t720 + 0.4E1 * t724)
        t729 = safelog(beta * t718 + t727)
        t731 = f(1,-3)
        t733 = t731 ** 2
        t737 = exp(-(t676 + t662) * alpha)
        t740 = sqrt(t49 * t733 + 0.4E1 * t737)
        t742 = safelog(beta * t731 + t740)
        t776 = vy(-2,2) ** 2
        t777 = t776 / 0.2E1
        t779 = vx(-2,-2) ** 2
        t784 = vx(-2,2) ** 2
        t785 = t784 / 0.2E1
        t787 = vy(-2,-2) ** 2
        t791 = f(-2,2)
        t793 = t791 ** 2
        t797 = exp(-(t784 + t776) * alpha)
        t800 = sqrt(t49 * t793 + 0.4E1 * t797)
        t802 = safelog(beta * t791 + t800)
        t805 = f(-2,-2)
        t807 = t805 ** 2
        t811 = exp(-(t779 + t787) * alpha)
        t814 = sqrt(t49 * t807 + 0.4E1 * t811)
        t816 = safelog(beta * t805 + t814)
        t864 = vy(2,-2) ** 2
        t865 = 0.5E1 / 0.128E3 * t864
        t869 = vx(2,-2) ** 2
        t870 = 0.5E1 / 0.128E3 * t869
        t875 = f(2,-2)
        t877 = t875 ** 2
        t881 = exp(-(t869 + t864) * alpha)
        t884 = sqrt(t49 * t877 + 0.4E1 * t881)
        t886 = safelog(beta * t875 + t884)
        t887 = 0.5E1 / 0.64E2 * t886
        t920 = vy(2,2) ** 2
        t925 = vx(2,2) ** 2
        t929 = f(2,2)
        t931 = t929 ** 2
        t935 = exp(-(t925 + t920) * alpha)
        t938 = sqrt(t49 * t931 + 0.4E1 * t935)
        t940 = safelog(beta * t929 + t938)
        t1086 = t57 * (0.400E3 / 0.3E1 * t141 * DD(0,0,3) * (t4 + t8 + a
     #lpha * t45 + t59 - 0.3E1 / 0.2000E4 * t71 - 0.3E1 / 0.2000E4 * t84
     # + 0.3E1 / 0.125E3 * t97 + 0.3E1 / 0.125E3 * t110 - 0.209E3 / 0.40
     #0E3 * t123 - 0.209E3 / 0.400E3 * t136) + 0.400E3 / 0.3E1 * t262 * 
     #(t4 + t8 + alpha * t180 + t59 - 0.3E1 / 0.2000E4 * t193 - 0.3E1 / 
     #0.2000E4 * t206 + 0.3E1 / 0.125E3 * t219 + 0.3E1 / 0.125E3 * t232 
     #- 0.209E3 / 0.400E3 * t245 - 0.209E3 / 0.400E3 * t258) * DD(0,0,1)
     #) + 0.1541E4 / 0.27E2 * t256 * (0.1024E4 / 0.7705E4 * dvx * (alpha
     # * t304 + t317 - 0.3E1 / 0.32E2 * t329 + 0.3E1 / 0.32E2 * t342 - t
     #355 + t368 - t381) * DD(-1,0,2) + (t4 + t8 + alpha * (-0.37E2 / 0.
     #15410E5 * t145 - 0.37E2 / 0.15410E5 * t163 + 0.87E2 / 0.1541E4 * t
     #148 - 0.1881E4 / 0.3082E4 * t151 + 0.87E2 / 0.1541E4 * t154 - 0.3E
     #1 / 0.15410E5 * t157 + 0.87E2 / 0.1541E4 * t166 - 0.1881E4 / 0.308
     #2E4 * t169 + 0.87E2 / 0.1541E4 * t172 - 0.3E1 / 0.15410E5 * t175) 
     #- 0.3E1 / 0.7705E4 * t219 + 0.174E3 / 0.1541E4 * t232 + t59 - 0.18
     #81E4 / 0.1541E4 * t258 - 0.37E2 / 0.7705E4 * t206 + 0.174E3 / 0.15
     #41E4 * t245) * dvy * DD(-1,0,1)) * dvy + 0.1541E4 / 0.27E2 * t134 
     #* (dvx * (t4 + t8 + alpha * (0.87E2 / 0.1541E4 * t13 - 0.1881E4 / 
     #0.3082E4 * t16 + 0.87E2 / 0.1541E4 * t19 - 0.3E1 / 0.15410E5 * t22
     # + 0.87E2 / 0.1541E4 * t31 - 0.1881E4 / 0.3082E4 * t34 + 0.87E2 / 
     #0.1541E4 * t37 - 0.3E1 / 0.15410E5 * t40 - 0.37E2 / 0.15410E5 * t1
     #0 - 0.37E2 / 0.15410E5 * t28) - 0.1881E4 / 0.1541E4 * t136 - 0.37E
     #2 / 0.7705E4 * t84 - 0.3E1 / 0.7705E4 * t97 + 0.174E3 / 0.1541E4 *
     # t110 + t59 + 0.174E3 / 0.1541E4 * t123) * DD(0,-1,3) + 0.1024E4 /
     # 0.7705E4 * dvy * (alpha * t463 + t476 - 0.3E1 / 0.32E2 * t488 + 0
     #.3E1 / 0.32E2 * t501 - t514 + 0.95E2 / 0.16E2 * t526 - t381) * DD(
     #0,-1,2)) * dvx + 0.1541E4 / 0.27E2 * t121 * (dvx * DD(0,1,3) * (t4
     # + t8 + alpha * (-0.37E2 / 0.15410E5 * t43 + 0.87E2 / 0.1541E4 * t
     #22 - 0.3E1 / 0.15410E5 * t31 + 0.87E2 / 0.1541E4 * t34 - 0.1881E4 
     #/ 0.3082E4 * t37 + 0.87E2 / 0.1541E4 * t40 - 0.3E1 / 0.15410E5 * t
     #13 + 0.87E2 / 0.1541E4 * t16 - 0.1881E4 / 0.3082E4 * t19 - 0.37E2 
     #/ 0.15410E5 * t25) + 0.174E3 / 0.1541E4 * t136 - 0.37E2 / 0.7705E4
     # * t71 + 0.174E3 / 0.1541E4 * t97 - 0.3E1 / 0.7705E4 * t110 - 0.18
     #81E4 / 0.1541E4 * t123 + t59) - 0.1024E4 / 0.7705E4 * dvy * DD(0,1
     #,2) * (alpha * t587 + t600 - 0.3E1 / 0.32E2 * t612 - t625 + 0.3E1 
     #/ 0.32E2 * t637 + 0.95E2 / 0.16E2 * t650 - t368)) * dvx + 0.1541E4
     # / 0.27E2 * t243 * (-0.96E2 / 0.7705E4 * dvx * DD(1,0,2) * (alpha 
     #* t689 + t702 - 0.32E2 / 0.3E1 * t714 + 0.190E3 / 0.3E1 * t650 - 0
     #.190E3 / 0.3E1 * t526 + 0.32E2 / 0.3E1 * t729 - t742) + (t4 + t8 +
     # alpha * (-0.1881E4 / 0.3082E4 * t154 + 0.87E2 / 0.1541E4 * t157 -
     # 0.3E1 / 0.15410E5 * t166 + 0.87E2 / 0.1541E4 * t169 - 0.37E2 / 0.
     #15410E5 * t160 - 0.1881E4 / 0.3082E4 * t172 + 0.87E2 / 0.1541E4 * 
     #t175 - 0.37E2 / 0.15410E5 * t178 - 0.3E1 / 0.15410E5 * t148 + 0.87
     #E2 / 0.1541E4 * t151) - 0.37E2 / 0.7705E4 * t193 + 0.174E3 / 0.154
     #1E4 * t219 - 0.3E1 / 0.7705E4 * t232 + t59 + 0.174E3 / 0.1541E4 * 
     #t258 - 0.1881E4 / 0.1541E4 * t245) * dvy * DD(1,0,1)) * dvy - 0.23
     #2E3 / 0.27E2 * t230 * dvy * (-0.2E1 / 0.29E2 * dvx * DD(-2,0,2) * 
     #(alpha * (0.32E2 / 0.5E1 * t434 - 0.32E2 / 0.5E1 * t558 + t777 - t
     #779 / 0.2E1 + 0.32E2 / 0.5E1 * t449 - 0.32E2 / 0.5E1 * t573 + t785
     # - t787 / 0.2E1) + t802 - 0.64E2 / 0.5E1 * t600 + 0.64E2 / 0.5E1 *
     # t476 - t816) + dvy * DD(-2,0,1) * (t4 + t8 + alpha * (0.37E2 / 0.
     #2320E4 * t145 + 0.37E2 / 0.2320E4 * t163 - 0.27E2 / 0.145E3 * t148
     # - 0.3E1 / 0.8E1 * t151 + 0.21E2 / 0.464E3 * t154 - 0.27E2 / 0.145
     #E3 * t166 - 0.3E1 / 0.8E1 * t169 + 0.21E2 / 0.464E3 * t172) - 0.54
     #E2 / 0.145E3 * t232 + t59 - 0.3E1 / 0.4E1 * t258 + 0.37E2 / 0.1160
     #E4 * t206 + 0.21E2 / 0.232E3 * t245)) - 0.232E3 / 0.27E2 * t108 * 
     #dvx * (dvx * DD(0,-2,3) * (t4 + t8 + alpha * (-0.27E2 / 0.145E3 * 
     #t13 - 0.3E1 / 0.8E1 * t16 + 0.21E2 / 0.464E3 * t19 - 0.27E2 / 0.14
     #5E3 * t31 - 0.3E1 / 0.8E1 * t34 + 0.21E2 / 0.464E3 * t37 + 0.37E2 
     #/ 0.2320E4 * t10 + 0.37E2 / 0.2320E4 * t28) - 0.3E1 / 0.4E1 * t136
     # + 0.37E2 / 0.1160E4 * t84 - 0.54E2 / 0.145E3 * t110 + t59 + 0.21E
     #2 / 0.232E3 * t123) - 0.128E3 / 0.145E3 * DD(0,-2,2) * dvy * (alph
     #a * (t273 - t665 / 0.2E1 + t865 - 0.5E1 / 0.128E3 * t779 + t291 - 
     #t679 / 0.2E1 + t870 - 0.5E1 / 0.128E3 * t787) + t317 - 0.5E1 / 0.6
     #4E2 * t816 + t887 - t729)) - 0.232E3 / 0.27E2 * t95 * dvx * (dvx *
     # (t4 + t8 + alpha * (0.21E2 / 0.464E3 * t16 - 0.3E1 / 0.8E1 * t19 
     #- 0.27E2 / 0.145E3 * t22 + 0.21E2 / 0.464E3 * t34 - 0.3E1 / 0.8E1 
     #* t37 - 0.27E2 / 0.145E3 * t40 + 0.37E2 / 0.2320E4 * t43 + 0.37E2 
     #/ 0.2320E4 * t25) + 0.21E2 / 0.232E3 * t136 + 0.37E2 / 0.1160E4 * 
     #t71 - 0.54E2 / 0.145E3 * t97 - 0.3E1 / 0.4E1 * t123 + t59) * DD(0,
     #2,3) - 0.2E1 / 0.29E2 * (alpha * (t777 - 0.32E2 / 0.5E1 * t281 + 0
     #.32E2 / 0.5E1 * t670 - t920 / 0.2E1 + t785 - 0.32E2 / 0.5E1 * t299
     # + 0.32E2 / 0.5E1 * t684 - t925 / 0.2E1) + t802 - t940 + 0.64E2 / 
     #0.5E1 * t714 - 0.64E2 / 0.5E1 * t355) * dvy * DD(0,2,2)) - 0.232E3
     # / 0.27E2 * t217 * (-0.128E3 / 0.145E3 * dvx * DD(2,0,2) * (alpha 
     #* (t865 - t441 + t565 - 0.5E1 / 0.128E3 * t920 + t870 - t456 + t58
     #0 - 0.5E1 / 0.128E3 * t925) + t625 - 0.5E1 / 0.64E2 * t940 - t514 
     #+ t887) + dvy * DD(2,0,1) * (t4 + t8 + alpha * (0.21E2 / 0.464E3 *
     # t151 - 0.3E1 / 0.8E1 * t154 - 0.27E2 / 0.145E3 * t157 + 0.21E2 / 
     #0.464E3 * t169 - 0.3E1 / 0.8E1 * t172 - 0.27E2 / 0.145E3 * t175 + 
     #0.37E2 / 0.2320E4 * t160 + 0.37E2 / 0.2320E4 * t178) + 0.37E2 / 0.
     #1160E4 * t193 - 0.54E2 / 0.145E3 * t219 + t59 + 0.21E2 / 0.232E3 *
     # t258 - 0.3E1 / 0.4E1 * t245)) * dvy + t204 * dvy * (0.32E2 / 0.45
     #E2 * dvx * (alpha * (-t446 / 0.2E1 + t570 / 0.2E1 - t461 / 0.2E1 +
     # t585 / 0.2E1) + t612 - t488) * DD(-3,0,2) + dvy * (t4 + t8 + alph
     #a * (-t145 / 0.10E2 - t163 / 0.10E2 - t148 / 0.90E2 - 0.7E1 / 0.18
     #E2 * t151 - t166 / 0.90E2 - 0.7E1 / 0.18E2 * t169) - t232 / 0.45E2
     # - 0.7E1 / 0.9E1 * t258 - t206 / 0.5E1 + t59) * DD(-3,0,1)) + t82 
     #* dvx * (dvx * DD(0,-3,3) * (t4 + t8 + alpha * (-t13 / 0.90E2 - 0.
     #7E1 / 0.18E2 * t16 - t31 / 0.90E2 - 0.7E1 / 0.18E2 * t34 - t10 / 0
     #.10E2 - t28 / 0.10E2) - t84 / 0.5E1 - t110 / 0.45E2 - 0.7E1 / 0.9E
     #1 * t136 + t59) - 0.32E2 / 0.45E2 * DD(0,-3,2) * dvy * (alpha * (t
     #287 / 0.2E1 - t676 / 0.2E1 + t269 / 0.2E1 - t662 / 0.2E1) + t329 -
     # t742)) + t69 * (dvx * (t4 + t8 + alpha * (-0.7E1 / 0.18E2 * t19 -
     # t22 / 0.90E2 - 0.7E1 / 0.18E2 * t37 - t40 / 0.90E2 - t25 / 0.10E2
     # - t43 / 0.10E2) - t71 / 0.5E1 - t97 / 0.45E2 - 0.7E1 / 0.9E1 * t1
     #23 + t59) * DD(0,3,3) - 0.32E2 / 0.45E2 * (alpha * (-t302 / 0.2E1 
     #+ t687 / 0.2E1 - t284 / 0.2E1 + t673 / 0.2E1) + t702 - t342) * dvy
     # * DD(0,3,2)) * dvx + dvy * (0.32E2 / 0.45E2 * dvx * (alpha * (t45
     #8 / 0.2E1 - t582 / 0.2E1 + t443 / 0.2E1 - t567 / 0.2E1) + t501 - t
     #637) * DD(3,0,2) + (t4 + t8 + alpha * (-t160 / 0.10E2 - t178 / 0.1
     #0E2 - 0.7E1 / 0.18E2 * t154 - t157 / 0.90E2 - 0.7E1 / 0.18E2 * t17
     #2 - t175 / 0.90E2) - t193 / 0.5E1 - t219 / 0.45E2 - 0.7E1 / 0.9E1 
     #* t245 + t59) * DD(3,0,1) * dvy) * t191
        cg0 = -0.3E1 / 0.256E3 / beta / t141 / t262 * nu * t1086
ccccc
ccccc
        rhs(i1,i2,i3,i4) = rhs(i1,i2,i3,i4)+cg0
      end do
      end do
      end do
      end do
c
      return
      end
