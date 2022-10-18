c
c Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
c See the top-level LICENSE file for details.
c Produced at the Lawrence Livermore National Laboratory
c
c SPDX-License-Identifier: MIT
c
c Fortran functions called by PitchAngleCollisionOperator.
c

      subroutine evaluateCollisionality(
     *     nuei, vx, vy, vxgrid, vygrid, range_lo, range_hi,
     *     vxmin, vxmax, vymin, vymax, vfloor, vthermal,
     *     nuCoefficient, solution_order)
c
c.. compute and add a simple collision operator to the RHS ...
      implicit none
c
c.. declarations of incoming variables 
      real nuei, vx, vy, vxgrid, vygrid, vxmin, vxmax, vymin, vymax
      real vfloor, vthermal, nuCoefficient
      real range_lo(1:2), range_hi(1:2)
      integer solution_order
c
c.. declarations of local variables
      real v, va, vb
      real alpha, vxra, vxrb, xi
      real beta, vyra, vyrb, eta
c
      vxra = range_lo(1)
      vxrb = range_hi(1)
      vyra = range_lo(2)
      vyrb = range_hi(2)
      v = max(sqrt(vx**2+vy**2), vfloor)
      nuei = nuCoefficient*(vthermal/v)**3
      if (vxgrid .lt. vxra .and. vxgrid .ge. vxmin) then
        va = vxra
        vb = vxmin
        xi = (vxgrid-va)/(vb-va)
      else if (vxgrid .gt. vxrb .and. vxgrid .le. vxmax) then
        va = vxrb
        vb = vxmax
        xi = (vxgrid-va)/(vb-va)
      else if (vxgrid .lt. vxmin .or. vxgrid .gt. vxmax) then
        xi = 1.0
      else
        xi = 0.0
      endif
      if (vygrid .lt. vyra .and. vygrid .ge. vymin) then
        va = vyra
        vb = vymin
        eta = (vygrid-va)/(vb-va)
      else if (vygrid .gt. vyrb .and. vygrid .le. vymax) then
        va = vyrb
        vb = vymax
        eta = (vygrid-va)/(vb-va)
      else if (vygrid .lt. vymin .or. vygrid .gt. vymax) then
        eta = 1.0
      else
        eta = 0.0
      endif

      if (solution_order .eq. 4) then
         alpha = 1.0 + xi**4*(
     *         +20.0*xi**3
     *         -70.0*xi**2
     *         +84.0*xi
     *         -35.0)
         beta = 1.0 + eta**4*(
     *        +20.0*eta**3
     *        -70.0*eta**2
     *        +84.0*eta
     *        -35.0)
      else
         alpha = 1.0 + xi**6*(
     *         +252.0 *xi**5
     *         -1386.0*xi**4
     *         +3080.0*xi**3
     *         -3465.0*xi**2
     *         +1980.0*xi
     *         -462.0)
         beta = 1.0 + eta**6*(
     *        +252.0 *eta**5
     *        -1386.0*eta**4
     *        +3080.0*eta**3
     *        -3465.0*eta**2
     *        +1980.0*eta
     *        -462.0)
      end if
      nuei = nuei*alpha*beta
c
      return
      end
c
c ++++++++++++++
c
      subroutine conservativePitchAngle_4th(
     *     rhs, f,
     *     nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b,
     *     n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b,
     *     xlo, xhi, dx,
     *     velocities, range_lo, range_hi,
     *     vfloor, Vth, nuCoefficient, IVx, IVy,
     *     do_relativity)
c
c.. compute and add a simple collision operator to the RHS ...
      implicit none
c
c.. declarations of incoming variables 
      integer nd1a, nd1b, nd2a, nd2b
      integer nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
      integer do_relativity
c
      real rhs(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real f(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real Vth(nd1a:nd1b, nd2a:nd2b)
      real IVx(nd1a:nd1b, nd2a:nd2b)
      real IVy(nd1a:nd1b, nd2a:nd2b)
      real range_lo(1:2), range_hi(1:2)
      real vfloor, nuCoefficient
      real xlo(1:4), xhi(1:4), dx(1:4)
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      integer gi3, gi4
      integer solution_order, vrolloff
      real vxgrid, vygrid, vxt, vyt
      real nuei
      real temp
      real vthermal, vflowx, vflowy
c
      real nu(-2:2, -2:2)
      real vx(-2:2, -2:2)
      real vy(-2:2, -2:2)
c
      real dvx, dvy, vxmin, vxmax, vymin, vymax, vxlo, vxhi, vylo, vyhi
c
c
c     generated declarations
        real t1
        real t10
        real t106
        real t107
        real t11
        real t112
        real t113
        real t116
        real t117
        integer t12
        real t121
        real t127
        real t13
        real t130
        real t133
        real t139
        real t143
        real t15
        real t155
        real t158
        integer t159
        real t160
        real t161
        real t164
        real t181
        real t182
        real t184
        real t185
        real t188
        real t189
        real t19
        real t191
        real t195
        real t196
        real t197
        real t199
        real t20
        real t202
        real t204
        real t208
        real t209
        real t21
        real t216
        real t22
        real t222
        real t223
        real t225
        real t228
        real t23
        real t230
        real t231
        real t24
        real t240
        real t244
        real t246
        real t248
        real t25
        real t252
        real t256
        real t26
        real t27
        real t271
        real t272
        real t273
        real t274
        real t276
        real t277
        real t278
        real t279
        real t28
        real t282
        real t283
        real t285
        real t289
        real t29
        real t290
        real t291
        real t292
        real t293
        real t295
        real t298
        real t3
        real t30
        real t303
        real t304
        real t305
        real t306
        real t318
        real t327
        real t331
        real t333
        real t335
        real t34
        real t354
        real t355
        real t357
        real t359
        real t36
        real t363
        real t367
        real t369
        real t37
        real t371
        real t372
        real t374
        real t376
        real t380
        real t382
        real t389
        real t39
        real t4
        real t41
        real t42
        real t426
        real t429
        integer t43
        real t437
        real t44
        real t45
        real t47
        real t48
        integer t5
        real t50
        real t54
        real t56
        real t58
        integer t6
        real t60
        integer t61
        real t62
        real t68
        real t7
        integer t72
        real t73
        real t8
        real t85
        real t88
        integer t89
        real t90
        real t91
        real t94
c
      solution_order = 4
      vrolloff = 3
      vxlo = xlo(3)
      vxhi = xhi(3)
      vylo = xlo(4)
      vyhi = xhi(4)
      dvx = dx(3)
      dvy = dx(4)
      vxmin = vxlo + vrolloff*dvx
      vxmax = vxhi - vrolloff*dvx
      vymin = vylo + vrolloff*dvy
      vymax = vyhi - vrolloff*dvy
c
c
      do i4 = n4a, n4b
      do i3 = n3a, n3b
      do i2 = n2a, n2b
      do i1 = n1a, n1b
        vthermal = Vth(i1,i2)
        vflowx = IVx(i1, i2)
        vflowy = IVy(i1, i2)

        ! setup variable coefficient arrays
        do gi4 = -2, 2
        do gi3 = -2, 2

            vxgrid = velocities(i3+gi3, i4+gi4, 0)
            vygrid = velocities(i3+gi3, i4+gi4, 1)
            vxt = vxgrid-vflowx
            vyt = vygrid-vflowy
            ! For relativity the grid is momentum based
            if (do_relativity .eq. 1) then
              vxgrid = vxlo + (i3+gi3+0.5)*dvx
              vygrid = vylo + (i4+gi4+0.5)*dvy
            end if

          call evaluateCollisionality(
     *       nuei, vxt, vyt, vxgrid, vygrid, range_lo, range_hi,
     *       vxmin, vxmax, vymin, vymax,
     *       vfloor, vthermal, nuCoefficient, solution_order)

          vx(gi3, gi4) = vxt
          vy(gi3, gi4) = vyt
          nu(gi3, gi4) = nuei
        end do
        end do
ccccc
cc generated code 
ccccc
        t1 = nu(1,0)
        t3 = t1 * vx(1,0)
        t4 = vy(1,0)
        t5 = i3 + 1
        t6 = i4 + 1
        t7 = f(i1,i2,t5,t6)
        t8 = f(i1,i2,t5,i4)
        t10 = 0.1E1 / dvy
        t11 = (t7 - t8) * t10
        t12 = i4 - 1
        t13 = f(i1,i2,t5,t12)
        t15 = (t8 - t13) * t10
        t19 = t3 * t4 * (t11 / 0.2E1 + t15 / 0.2E1)
        t20 = nu(0,0)
        t21 = vx(0,0)
        t22 = t20 * t21
        t23 = vy(0,0)
        t24 = f(i1,i2,i3,t6)
        t25 = f(i1,i2,i3,i4)
        t26 = t24 - t25
        t27 = t26 * t10
        t28 = f(i1,i2,i3,t12)
        t29 = t25 - t28
        t30 = t29 * t10
        t34 = t22 * t23 * (t27 / 0.2E1 + t30 / 0.2E1)
        t36 = 0.1E1 / dvx
        t37 = (t19 - t34) * t36
        t39 = nu(-1,0)
        t41 = t39 * vx(-1,0)
        t42 = vy(-1,0)
        t43 = i3 - 1
        t44 = f(i1,i2,t43,t6)
        t45 = f(i1,i2,t43,i4)
        t47 = (t44 - t45) * t10
        t48 = f(i1,i2,t43,t12)
        t50 = (t45 - t48) * t10
        t54 = t41 * t42 * (t47 / 0.2E1 + t50 / 0.2E1)
        t56 = (t34 - t54) * t36
        t58 = dvx ** 2
        t60 = dvy / t58
        t61 = i4 + 2
        t62 = f(i1,i2,t5,t61)
        t68 = t60 * (t11 - t15)
        t72 = i4 - 2
        t73 = f(i1,i2,t5,t72)
        t85 = nu(2,0)
        t88 = vy(2,0)
        t89 = i3 + 2
        t90 = f(i1,i2,t89,t6)
        t91 = f(i1,i2,t89,i4)
        t94 = f(i1,i2,t89,t12)
        t106 = f(i1,i2,i3,t61)
        t107 = t106 - t24
        t112 = t60 * (t27 - t30)
        t113 = -t60 * (t107 * t10 - t27) + t112
        t116 = f(i1,i2,i3,t72)
        t117 = t28 - t116
        t121 = -t112 + t60 * (t30 - t117 * t10)
        t127 = t22 * t23 * (t113 * t10 / 0.12E2 + t121 * t10 / 0.12E2)
        t130 = (t37 - t56) * t36 / 0.6E1
        t133 = f(i1,i2,t43,t61)
        t139 = t60 * (t47 - t50)
        t143 = f(i1,i2,t43,t72)
        t155 = nu(-2,0)
        t158 = vy(-2,0)
        t159 = i3 - 2
        t160 = f(i1,i2,t159,t6)
        t161 = f(i1,i2,t159,i4)
        t164 = f(i1,i2,t159,t12)
        t181 = t4 ** 2
        t182 = t1 * t181
        t184 = t88 ** 2
        t185 = t85 * t184
        t188 = t23 ** 2
        t189 = t20 * t188
        t191 = (t182 - t189) * t36
        t195 = t189 / 0.2E1
        t196 = t42 ** 2
        t197 = t39 * t196
        t199 = (t189 - t197) * t36
        t202 = dvx * (t191 - t199) / 0.16E2
        t204 = t8 - t25
        t208 = t158 ** 2
        t209 = t155 * t208
        t216 = t25 - t45
        t222 = t182 / 0.2E1 + t189 / 0.2E1
        t223 = t91 - t8
        t225 = t204 * t36
        t228 = t216 * t36
        t230 = (t225 - t228) * t36
        t231 = -(t223 * t36 - t225) * t36 + t230
        t240 = t222 * t204 * t36
        t244 = t189 / 0.2E1 + t197 / 0.2E1
        t246 = t244 * t216 * t36
        t248 = (t240 - t246) * t36
        t252 = t45 - t161
        t256 = -t230 + (t228 - t252 * t36) * t36
        t271 = nu(0,1)
        t272 = vx(0,1)
        t273 = t272 ** 2
        t274 = t271 * t273
        t276 = nu(0,2)
        t277 = vx(0,2)
        t278 = t277 ** 2
        t279 = t276 * t278
        t282 = t21 ** 2
        t283 = t20 * t282
        t285 = (t274 - t283) * t10
        t289 = t283 / 0.2E1
        t290 = nu(0,-1)
        t291 = vx(0,-1)
        t292 = t291 ** 2
        t293 = t290 * t292
        t295 = (t283 - t293) * t10
        t298 = dvy * (t285 - t295) / 0.16E2
        t303 = nu(0,-2)
        t304 = vx(0,-2)
        t305 = t304 ** 2
        t306 = t303 * t305
        t318 = t274 / 0.2E1 + t283 / 0.2E1
        t327 = t318 * t26 * t10
        t331 = t283 / 0.2E1 + t293 / 0.2E1
        t333 = t331 * t29 * t10
        t335 = (t327 - t333) * t10
        t354 = t271 * t272
        t355 = vy(0,1)
        t357 = (t7 - t24) * t36
        t359 = (t24 - t44) * t36
        t363 = t354 * t355 * (t357 / 0.2E1 + t359 / 0.2E1)
        t367 = t22 * t23 * (t225 / 0.2E1 + t228 / 0.2E1)
        t369 = (t363 - t367) * t10
        t371 = t290 * t291
        t372 = vy(0,-1)
        t374 = (t13 - t28) * t36
        t376 = (t28 - t48) * t36
        t380 = t371 * t372 * (t374 / 0.2E1 + t376 / 0.2E1)
        t382 = (t367 - t380) * t10
        t389 = (t357 - t359) * t36
        t426 = t22 * t23 * (t231 * t36 / 0.12E2 + t256 * t36 / 0.12E2)
        t429 = t60 * (t369 - t382) / 0.6E1
        t437 = (t374 - t376) * t36
        temp = 
     #-t37 / 0.2E1 - t56 / 0.2E1 + t58 * (-(t3 * t4 * ((-t60 * (
     #(t62 - t7) * t10 - t11) / 0.6E1 + t68 / 0.6E1) * t10 / 0.2E1 + (-t
     #68 / 0.6E1 + t60 * (t15 - (t13 - t73) * t10) / 0.6E1) * t10 / 0.2E
     #1) - ((t85 * vx(2,0) * t88 * ((t90 - t91) * t10 / 0.2E1 + (t91 - t
     #94) * t10 / 0.2E1) - t19) * t36 - t37) * t36 / 0.6E1 - t127 + t130
     #) * t36 / 0.2E1 - (t127 - t130 - t41 * t42 * ((-t60 * ((t133 - t44
     #) * t10 - t47) / 0.6E1 + t139 / 0.6E1) * t10 / 0.2E1 + (-t139 / 0.
     #6E1 + t60 * (t50 - (t48 - t143) * t10) / 0.6E1) * t10 / 0.2E1) + (
     #t56 - (t54 - t155 * vx(-2,0) * t158 * ((t160 - t161) * t10 / 0.2E1
     # + (t161 - t164) * t10 / 0.2E1)) * t36) * t36 / 0.6E1) * t36 / 0.2
     #E1) + ((t182 / 0.2E1 - dvx * ((t185 - t182) * t36 - t191) / 0.16E2
     # + t195 - t202) * t204 * t36 - (t195 - t202 + t197 / 0.2E1 - dvx *
     # (t199 - (t197 - t209) * t36) / 0.16E2) * t216 * t36) * t36 + dvx 
     #* (t222 * t231 * t36 / 0.24E2 - (((t185 / 0.2E1 + t182 / 0.2E1) * 
     #t223 * t36 - t240) * t36 - t248) * t36 / 0.24E2 - t244 * t256 * t3
     #6 / 0.24E2 + (t248 - (t246 - (t197 / 0.2E1 + t209 / 0.2E1) * t252 
     #* t36) * t36) * t36 / 0.24E2) + ((t274 / 0.2E1 - dvy * ((t279 - t2
     #74) * t10 - t285) / 0.16E2 + t289 - t298) * t26 * t10 - (t289 - t2
     #98 + t293 / 0.2E1 - dvy * (t295 - (t293 - t306) * t10) / 0.16E2) *
     # t29 * t10) * t10 + t58 * (t318 * t113 * t10 / 0.24E2 - t60 * (((t
     #279 / 0.2E1 + t274 / 0.2E1) * t107 * t10 - t327) * t10 - t335) / 0
     #.24E2 - t331 * t121 * t10 / 0.24E2 + t60 * (t335 - (t333 - (t293 /
     # 0.2E1 + t306 / 0.2E1) * t117 * t10) * t10) / 0.24E2) * t10 - t369
     # / 0.2E1 - t382 / 0.2E1 + t58 * (-(t354 * t355 * ((-((t90 - t7) * 
     #t36 - t357) * t36 / 0.6E1 + t389 / 0.6E1) * t36 / 0.2E1 + (-t389 /
     # 0.6E1 + (t359 - (t44 - t160) * t36) * t36 / 0.6E1) * t36 / 0.2E1)
     # - t60 * ((t276 * t277 * vy(0,2) * ((t62 - t106) * t36 / 0.2E1 + (
     #t106 - t133) * t36 / 0.2E1) - t363) * t10 - t369) / 0.6E1 - t426 +
     # t429) * t10 / 0.2E1 - (t426 - t429 - t371 * t372 * ((-((t94 - t13
     #) * t36 - t374) * t36 / 0.6E1 + t437 / 0.6E1) * t36 / 0.2E1 + (-t4
     #37 / 0.6E1 + (t376 - (t48 - t164) * t36) * t36 / 0.6E1) * t36 / 0.
     #2E1) + t60 * (t382 - (t380 - t303 * t304 * vy(0,-2) * ((t73 - t116
     #) * t36 / 0.2E1 + (t116 - t143) * t36 / 0.2E1)) * t10) / 0.6E1) * 
     #t10 / 0.2E1)

ccccc
ccccc
        rhs(i1, i2, i3, i4) = rhs(i1, i2, i3, i4)+temp
      end do
      end do
      end do
      end do
c
      return
      end
c
c ++++++++++++++
c
      subroutine conservativePitchAngle_6th(
     *     rhs, f,
     *     nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b,
     *     n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b,
     *     xlo, xhi, dx,
     *     velocities, range_lo, range_hi,
     *     vfloor, Vth, nuCoefficient, IVx, IVy,
     *     do_relativity)
c
c.. compute and add a simple collision operator to the RHS ...
      implicit none
c
c.. declarations of incoming variables 
      integer nd1a, nd1b, nd2a, nd2b
      integer nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
      integer do_relativity
c
      real rhs(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real f(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real Vth(nd1a:nd1b, nd2a:nd2b)
      real IVx(nd1a:nd1b, nd2a:nd2b)
      real IVy(nd1a:nd1b, nd2a:nd2b)
      real range_lo(1:2), range_hi(1:2)
      real vfloor, nuCoefficient
      real xlo(1:4), xhi(1:4), dx(1:4)
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      integer gi3, gi4
      integer solution_order, vrolloff
      real vxgrid, vygrid, vxt, vyt
      real nuei
      real temp
      real vthermal, vflowx, vflowy
c
      real nu(-3:3, -3:3)
      real vx(-3:3, -3:3)
      real vy(-3:3, -3:3)

      real dvx, dvy, vxmin, vxmax, vymin, vymax, vxlo, vxhi, vylo, vyhi
c
c
c     generated declarations
      real t1
        real t10
        real t100
        real t1000
        real t102
        real t1024
        real t104
        real t1042
        real t106
        real t1060
        real t1063
        real t1065
        real t107
        real t1076
        real t1078
        real t108
        real t1081
        real t1082
        real t1084
        real t1087
        real t109
        real t1097
        real t11
        real t110
        real t1101
        real t1103
        real t111
        real t112
        real t1127
        real t113
        real t116
        real t117
        real t118
        real t119
        integer t12
        real t120
        real t121
        real t127
        real t129
        real t13
        real t130
        real t133
        real t135
        real t136
        real t138
        real t139
        real t143
        real t145
        real t146
        real t15
        real t154
        real t155
        real t157
        real t158
        integer t159
        real t160
        real t161
        real t163
        real t164
        real t166
        real t170
        real t172
        real t174
        real t181
        real t182
        real t183
        real t185
        integer t186
        real t187
        real t19
        real t192
        real t195
        real t197
        real t20
        real t200
        real t202
        real t204
        integer t208
        real t209
        real t21
        real t22
        real t225
        real t23
        real t231
        real t235
        real t24
        real t25
        real t250
        real t254
        real t257
        integer t258
        real t259
        real t26
        real t260
        real t263
        real t27
        real t277
        real t28
        real t281
        real t282
        real t284
        real t286
        real t289
        real t29
        real t291
        real t294
        real t296
        real t298
        real t299
        real t3
        real t30
        real t302
        real t303
        real t305
        real t311
        real t317
        real t319
        real t322
        real t324
        real t327
        real t330
        real t335
        real t338
        real t34
        real t340
        real t343
        real t345
        real t347
        real t351
        real t36
        real t367
        real t37
        real t373
        real t377
        real t39
        real t394
        real t397
        integer t398
        real t399
        real t4
        real t400
        real t403
        real t41
        real t42
        real t424
        real t425
        real t426
        real t427
        real t428
        integer t43
        real t430
        real t431
        real t432
        real t434
        real t435
        real t437
        real t438
        real t439
        real t44
        real t440
        real t443
        real t445
        real t448
        real t449
        real t45
        real t451
        real t452
        real t453
        real t455
        real t459
        real t461
        real t462
        real t463
        real t465
        real t466
        real t467
        real t469
        real t47
        real t472
        real t474
        real t477
        real t479
        real t48
        real t480
        real t481
        real t484
        real t492
        real t497
        real t498
        real t499
        integer t5
        real t50
        real t500
        real t502
        real t503
        real t505
        real t506
        real t507
        real t517
        real t520
        real t522
        real t524
        real t528
        real t529
        real t531
        real t532
        real t533
        real t54
        real t550
        real t551
        real t555
        real t558
        real t56
        real t562
        real t564
        real t565
        real t570
        real t576
        real t58
        real t580
        real t582
        real t584
        real t593
        real t597
        real t599
        integer t6
        real t60
        real t603
        real t605
        real t607
        integer t61
        real t611
        real t613
        real t615
        real t617
        real t619
        real t62
        real t623
        real t627
        real t631
        real t64
        real t65
        real t659
        real t660
        real t661
        real t662
        real t663
        real t664
        real t665
        real t666
        real t667
        real t669
        real t67
        real t670
        real t671
        real t673
        real t674
        real t676
        real t677
        real t678
        real t679
        real t68
        real t680
        real t683
        real t685
        real t688
        real t689
        real t690
        real t691
        real t693
        real t694
        real t695
        real t697
        real t7
        real t701
        real t703
        real t704
        real t705
        real t706
        real t707
        real t709
        real t710
        real t711
        real t713
        real t716
        integer t72
        real t720
        real t722
        real t723
        real t724
        real t725
        real t726
        real t729
        real t73
        real t741
        real t742
        real t75
        real t752
        real t755
        real t757
        real t759
        real t76
        real t763
        real t781
        real t786
        real t793
        real t797
        real t799
        real t8
        real t801
        real t810
        real t814
        real t816
        real t820
        real t822
        real t824
        real t828
        real t830
        real t832
        real t834
        real t836
        real t84
        real t85
        real t869
        real t87
        real t870
        real t872
        real t874
        real t878
        real t88
        real t882
        real t884
        real t886
        real t887
        real t889
        integer t89
        real t891
        real t895
        real t897
        real t90
        real t900
        real t902
        real t904
        real t905
        real t909
        real t91
        real t911
        real t912
        real t918
        real t919
        real t920
        real t922
        real t924
        real t928
        real t93
        real t930
        real t931
        real t94
        real t941
        real t942
        real t944
        real t948
        real t950
        real t952
        real t953
        real t957
        real t959
        real t96
        real t960
        real t966
        real t967
        real t968
        real t970
        real t972
        real t976
        real t978
        real t979
        real t994
        real t998
c
      solution_order = 6
      vrolloff = 4
      vxlo = xlo(3)
      vxhi = xhi(3)
      vylo = xlo(4)
      vyhi = xhi(4)
      dvx = dx(3)
      dvy = dx(4)
      vxmin = vxlo + vrolloff*dvx
      vxmax = vxhi - vrolloff*dvx
      vymin = vylo + vrolloff*dvy
      vymax = vyhi - vrolloff*dvy
c
      do i4 = n4a, n4b
      do i3 = n3a, n3b
      do i2 = n2a, n2b
      do i1 = n1a, n1b
        vthermal = Vth(i1,i2)
        vflowx = IVx(i1, i2)
        vflowy = IVy(i1, i2)

        ! setup variable coefficient arrays
        do gi4 = -3, 3
        do gi3 = -3, 3

            vxgrid = velocities(i3+gi3, i4+gi4, 0)
            vygrid = velocities(i3+gi3, i4+gi4, 1)
            vxt = vxgrid-vflowx
            vyt = vygrid-vflowy
            ! For relativity the grid is momentum based
            if (do_relativity .eq. 1) then
              vxgrid = vxlo + (i3+gi3+0.5)*dvx
              vygrid = vylo + (i4+gi4+0.5)*dvy
            end if

          call evaluateCollisionality(
     *       nuei, vxt, vyt, vxgrid, vygrid, range_lo, range_hi,
     *       vxmin, vxmax, vymin, vymax,
     *       vfloor, vthermal, nuCoefficient, solution_order)

          vx(gi3, gi4) = vxt
          vy(gi3, gi4) = vyt
          nu(gi3, gi4) = nuei
        end do
        end do
ccccc
cc generated code 
ccccc
        t1 = nu(1,0)
        t3 = t1 * vx(1,0)
        t4 = vy(1,0)
        t5 = i3 + 1
        t6 = i4 + 1
        t7 = f(i1,i2,t5,t6)
        t8 = f(i1,i2,t5,i4)
        t10 = 0.1E1 / dvy
        t11 = (t7 - t8) * t10
        t12 = i4 - 1
        t13 = f(i1,i2,t5,t12)
        t15 = (t8 - t13) * t10
        t19 = t3 * t4 * (t11 / 0.2E1 + t15 / 0.2E1)
        t20 = nu(0,0)
        t21 = vx(0,0)
        t22 = t20 * t21
        t23 = vy(0,0)
        t24 = f(i1,i2,i3,t6)
        t25 = f(i1,i2,i3,i4)
        t26 = t24 - t25
        t27 = t26 * t10
        t28 = f(i1,i2,i3,t12)
        t29 = t25 - t28
        t30 = t29 * t10
        t34 = t22 * t23 * (t27 / 0.2E1 + t30 / 0.2E1)
        t36 = 0.1E1 / dvx
        t37 = (t19 - t34) * t36
        t39 = nu(-1,0)
        t41 = t39 * vx(-1,0)
        t42 = vy(-1,0)
        t43 = i3 - 1
        t44 = f(i1,i2,t43,t6)
        t45 = f(i1,i2,t43,i4)
        t47 = (t44 - t45) * t10
        t48 = f(i1,i2,t43,t12)
        t50 = (t45 - t48) * t10
        t54 = t41 * t42 * (t47 / 0.2E1 + t50 / 0.2E1)
        t56 = (t34 - t54) * t36
        t58 = dvx ** 2
        t60 = dvy / t58
        t61 = i4 + 2
        t62 = f(i1,i2,t5,t61)
        t64 = (t62 - t7) * t10
        t65 = t64 - t11
        t67 = t11 - t15
        t68 = t60 * t67
        t72 = i4 - 2
        t73 = f(i1,i2,t5,t72)
        t75 = (t13 - t73) * t10
        t76 = t15 - t75
        t84 = t3 * t4 * ((-t60 * t65 / 0.6E1 + t68 / 0.6E1) * t10 / 0.2E
     #1 + (t60 * t76 / 0.6E1 - t68 / 0.6E1) * t10 / 0.2E1)
        t85 = nu(2,0)
        t87 = t85 * vx(2,0)
        t88 = vy(2,0)
        t89 = i3 + 2
        t90 = f(i1,i2,t89,t6)
        t91 = f(i1,i2,t89,i4)
        t93 = (t90 - t91) * t10
        t94 = f(i1,i2,t89,t12)
        t96 = (t91 - t94) * t10
        t100 = t87 * t88 * (t93 / 0.2E1 + t96 / 0.2E1)
        t102 = (t100 - t19) * t36
        t104 = (t102 - t37) * t36
        t106 = f(i1,i2,i3,t61)
        t107 = t106 - t24
        t108 = t107 * t10
        t109 = t108 - t27
        t110 = t60 * t109
        t111 = t27 - t30
        t112 = t60 * t111
        t113 = -t110 + t112
        t116 = f(i1,i2,i3,t72)
        t117 = t28 - t116
        t118 = t117 * t10
        t119 = t30 - t118
        t120 = t60 * t119
        t121 = -t112 + t120
        t127 = t22 * t23 * (t113 * t10 / 0.12E2 + t121 * t10 / 0.12E2)
        t129 = (t37 - t56) * t36
        t130 = t129 / 0.6E1
        t133 = f(i1,i2,t43,t61)
        t135 = (t133 - t44) * t10
        t136 = t135 - t47
        t138 = t47 - t50
        t139 = t60 * t138
        t143 = f(i1,i2,t43,t72)
        t145 = (t48 - t143) * t10
        t146 = t50 - t145
        t154 = t41 * t42 * ((-t60 * t136 / 0.6E1 + t139 / 0.6E1) * t10 /
     # 0.2E1 + (t60 * t146 / 0.6E1 - t139 / 0.6E1) * t10 / 0.2E1)
        t155 = nu(-2,0)
        t157 = t155 * vx(-2,0)
        t158 = vy(-2,0)
        t159 = i3 - 2
        t160 = f(i1,i2,t159,t6)
        t161 = f(i1,i2,t159,i4)
        t163 = (t160 - t161) * t10
        t164 = f(i1,i2,t159,t12)
        t166 = (t161 - t164) * t10
        t170 = t157 * t158 * (t163 / 0.2E1 + t166 / 0.2E1)
        t172 = (t54 - t170) * t36
        t174 = (t56 - t172) * t36
        t181 = t58 ** 2
        t182 = dvy ** 2
        t183 = t182 * dvy
        t185 = t183 / t181
        t186 = i4 + 3
        t187 = f(i1,i2,t5,t186)
        t192 = t65 * t10
        t195 = t67 * t10
        t197 = (t192 - t195) * t10
        t200 = t76 * t10
        t202 = (t195 - t200) * t10
        t204 = t185 * (t197 - t202)
        t208 = i4 - 3
        t209 = f(i1,i2,t5,t208)
        t225 = f(i1,i2,t89,t61)
        t231 = t60 * (t93 - t96)
        t235 = f(i1,i2,t89,t72)
        t250 = (t84 - t127) * t36
        t254 = nu(3,0)
        t257 = vy(3,0)
        t258 = i3 + 3
        t259 = f(i1,i2,t258,t6)
        t260 = f(i1,i2,t258,i4)
        t263 = f(i1,i2,t258,t12)
        t277 = (t104 - t129) * t36
        t281 = f(i1,i2,i3,t186)
        t282 = t281 - t106
        t284 = t10 * t282 - t108
        t286 = t109 * t10
        t289 = t111 * t10
        t291 = (t286 - t289) * t10
        t294 = t119 * t10
        t296 = (t289 - t294) * t10
        t298 = t185 * (t291 - t296)
        t299 = t185 * ((t10 * t284 - t286) * t10 - t291) - t298
        t302 = f(i1,i2,i3,t208)
        t303 = t116 - t302
        t305 = -t10 * t303 + t118
        t311 = t298 - t185 * (t296 - (-t10 * t305 + t294) * t10)
        t317 = t22 * t23 * (t299 * t10 / 0.60E2 + t311 * t10 / 0.60E2)
        t319 = (t127 - t154) * t36
        t322 = (t250 - t319) * t36 / 0.6E1
        t324 = (t129 - t174) * t36
        t327 = (t277 - t324) * t36 / 0.30E2
        t330 = f(i1,i2,t43,t186)
        t335 = t136 * t10
        t338 = t138 * t10
        t340 = (t335 - t338) * t10
        t343 = t146 * t10
        t345 = (t338 - t343) * t10
        t347 = t185 * (t340 - t345)
        t351 = f(i1,i2,t43,t208)
        t367 = f(i1,i2,t159,t61)
        t373 = t60 * (t163 - t166)
        t377 = f(i1,i2,t159,t72)
        t394 = nu(-3,0)
        t397 = vy(-3,0)
        t398 = i3 - 3
        t399 = f(i1,i2,t398,t6)
        t400 = f(i1,i2,t398,i4)
        t403 = f(i1,i2,t398,t12)
        t424 = t4 ** 2
        t425 = t1 * t424
        t426 = t425 / 0.2E1
        t427 = t88 ** 2
        t428 = t85 * t427
        t430 = (-t425 + t428) * t36
        t431 = t23 ** 2
        t432 = t20 * t431
        t434 = (-t432 + t425) * t36
        t435 = t430 - t434
        t437 = dvx * t435 / 0.16E2
        t438 = t58 * dvx
        t439 = t257 ** 2
        t440 = t254 * t439
        t443 = (-t428 + t440) * t36 - t430
        t445 = t435 * t36
        t448 = t42 ** 2
        t449 = t39 * t448
        t451 = (-t449 + t432) * t36
        t452 = t434 - t451
        t453 = t452 * t36
        t455 = (t445 - t453) * t36
        t459 = t432 / 0.2E1
        t461 = dvx * t452 / 0.16E2
        t462 = t158 ** 2
        t463 = t155 * t462
        t465 = (-t463 + t449) * t36
        t466 = t451 - t465
        t467 = t466 * t36
        t469 = (t453 - t467) * t36
        t472 = 0.3E1 / 0.256E3 * t438 * (t455 - t469)
        t474 = t8 - t25
        t477 = t449 / 0.2E1
        t479 = dvx * t466 / 0.16E2
        t480 = t397 ** 2
        t481 = t394 * t480
        t484 = t465 - (-t481 + t463) * t36
        t492 = t25 - t45
        t497 = t426 - t437 + t459 - t461
        t498 = t91 - t8
        t499 = t498 * t36
        t500 = t474 * t36
        t502 = (t499 - t500) * t36
        t503 = t492 * t36
        t505 = (t500 - t503) * t36
        t506 = -t502 + t505
        t507 = t506 / 0.24E2
        t517 = t497 * t474 * t36
        t520 = t459 - t461 + t477 - t479
        t522 = t520 * t492 * t36
        t524 = (t517 - t522) * t36
        t528 = t45 - t161
        t529 = t528 * t36
        t531 = (t503 - t529) * t36
        t532 = -t505 + t531
        t533 = t532 / 0.24E2
        t550 = t432 / 0.2E1 + t425 / 0.2E1
        t551 = t260 - t91
        t555 = (t36 * t551 - t499) * t36 - t502
        t558 = -t506 * t36
        t562 = -t532 * t36
        t564 = (t558 - t562) * t36
        t565 = (t36 * t555 - t558) * t36 - t564
        t570 = t425 / 0.2E1 + t428 / 0.2E1
        t576 = t550 * t507 * t36
        t580 = t449 / 0.2E1 + t432 / 0.2E1
        t582 = t580 * t533 * t36
        t584 = (t576 - t582) * t36
        t593 = t570 * t498 * t36
        t597 = t550 * t474 * t36
        t599 = (t593 - t597) * t36
        t603 = t580 * t492 * t36
        t605 = (t597 - t603) * t36
        t607 = (t599 - t605) * t36
        t611 = t463 / 0.2E1 + t449 / 0.2E1
        t613 = t611 * t528 * t36
        t615 = (t603 - t613) * t36
        t617 = (t605 - t615) * t36
        t619 = (t607 - t617) * t36
        t623 = t161 - t400
        t627 = t531 - (-t36 * t623 + t529) * t36
        t631 = t564 - (-t36 * t627 + t562) * t36
        t659 = nu(0,1)
        t660 = vx(0,1)
        t661 = t660 ** 2
        t662 = t659 * t661
        t663 = t662 / 0.2E1
        t664 = nu(0,2)
        t665 = vx(0,2)
        t666 = t665 ** 2
        t667 = t664 * t666
        t669 = (-t662 + t667) * t10
        t670 = t21 ** 2
        t671 = t20 * t670
        t673 = (-t671 + t662) * t10
        t674 = t669 - t673
        t676 = dvy * t674 / 0.16E2
        t677 = nu(0,3)
        t678 = vx(0,3)
        t679 = t678 ** 2
        t680 = t677 * t679
        t683 = (-t667 + t680) * t10 - t669
        t685 = t674 * t10
        t688 = nu(0,-1)
        t689 = vx(0,-1)
        t690 = t689 ** 2
        t691 = t688 * t690
        t693 = (-t691 + t671) * t10
        t694 = t673 - t693
        t695 = t694 * t10
        t697 = (t685 - t695) * t10
        t701 = t671 / 0.2E1
        t703 = dvy * t694 / 0.16E2
        t704 = nu(0,-2)
        t705 = vx(0,-2)
        t706 = t705 ** 2
        t707 = t704 * t706
        t709 = (-t707 + t691) * t10
        t710 = t693 - t709
        t711 = t710 * t10
        t713 = (t695 - t711) * t10
        t716 = 0.3E1 / 0.256E3 * t183 * (t697 - t713)
        t720 = t691 / 0.2E1
        t722 = dvy * t710 / 0.16E2
        t723 = nu(0,-3)
        t724 = vx(0,-3)
        t725 = t724 ** 2
        t726 = t723 * t725
        t729 = t709 - (-t726 + t707) * t10
        t741 = t663 - t676 + t701 - t703
        t742 = t113 / 0.24E2
        t752 = t741 * t26 * t10
        t755 = t701 - t703 + t720 - t722
        t757 = t755 * t29 * t10
        t759 = (t752 - t757) * t10
        t763 = t121 / 0.24E2
        t781 = t671 / 0.2E1 + t662 / 0.2E1
        t786 = t662 / 0.2E1 + t667 / 0.2E1
        t793 = t781 * t742 * t10
        t797 = t691 / 0.2E1 + t671 / 0.2E1
        t799 = t797 * t763 * t10
        t801 = (t793 - t799) * t10
        t810 = t786 * t107 * t10
        t814 = t781 * t26 * t10
        t816 = (t810 - t814) * t10
        t820 = t797 * t29 * t10
        t822 = (t814 - t820) * t10
        t824 = (t816 - t822) * t10
        t828 = t707 / 0.2E1 + t691 / 0.2E1
        t830 = t828 * t117 * t10
        t832 = (t820 - t830) * t10
        t834 = (t822 - t832) * t10
        t836 = (t824 - t834) * t10
        t869 = t659 * t660
        t870 = vy(0,1)
        t872 = (t7 - t24) * t36
        t874 = (t24 - t44) * t36
        t878 = t869 * t870 * (t872 / 0.2E1 + t874 / 0.2E1)
        t882 = t22 * t23 * (t500 / 0.2E1 + t503 / 0.2E1)
        t884 = (t878 - t882) * t10
        t886 = t688 * t689
        t887 = vy(0,-1)
        t889 = (t13 - t28) * t36
        t891 = (t28 - t48) * t36
        t895 = t886 * t887 * (t889 / 0.2E1 + t891 / 0.2E1)
        t897 = (t882 - t895) * t10
        t900 = (t90 - t7) * t36
        t902 = (t900 - t872) * t36
        t904 = (t872 - t874) * t36
        t905 = -t902 + t904
        t909 = (t44 - t160) * t36
        t911 = (t874 - t909) * t36
        t912 = -t904 + t911
        t918 = t869 * t870 * (t905 * t36 / 0.12E2 + t912 * t36 / 0.12E2)
        t919 = t664 * t665
        t920 = vy(0,2)
        t922 = (t62 - t106) * t36
        t924 = (t106 - t133) * t36
        t928 = t919 * t920 * (t922 / 0.2E1 + t924 / 0.2E1)
        t930 = (t928 - t878) * t10
        t931 = t930 - t884
        t941 = t22 * t23 * (t506 * t36 / 0.12E2 + t532 * t36 / 0.12E2)
        t942 = t884 - t897
        t944 = t60 * t942 / 0.6E1
        t948 = (t94 - t13) * t36
        t950 = (t948 - t889) * t36
        t952 = (t889 - t891) * t36
        t953 = -t950 + t952
        t957 = (t48 - t164) * t36
        t959 = (t891 - t957) * t36
        t960 = -t952 + t959
        t966 = t886 * t887 * (t953 * t36 / 0.12E2 + t960 * t36 / 0.12E2)
        t967 = t704 * t705
        t968 = vy(0,-2)
        t970 = (t73 - t116) * t36
        t972 = (t116 - t143) * t36
        t976 = t967 * t968 * (t970 / 0.2E1 + t972 / 0.2E1)
        t978 = (t895 - t976) * t10
        t979 = t897 - t978
        t994 = -t905 * t36
        t998 = -t912 * t36
        t1000 = (t994 - t998) * t36
        t1024 = (t922 - t924) * t36
        t1042 = (t918 - t941) * t10
        t1060 = t931 * t10
        t1063 = t942 * t10
        t1065 = (t1060 - t1063) * t10
        t1076 = t22 * t23 * (t565 * t36 / 0.60E2 + t631 * t36 / 0.60E2)
        t1078 = (t941 - t966) * t10
        t1081 = t60 * (t1042 - t1078) / 0.6E1
        t1082 = t979 * t10
        t1084 = (t1063 - t1082) * t10
        t1087 = t185 * (t1065 - t1084) / 0.30E2
        t1097 = -t953 * t36
        t1101 = -t960 * t36
        t1103 = (t1097 - t1101) * t36
        t1127 = (t970 - t972) * t36
        temp = 
     #-t37 / 0.2E1 - t56 / 0.2E1 + t58 * (-(t84 - t104 / 0.6E1 -
     # t127 + t130) * t36 / 0.2E1 - (t127 - t130 - t154 + t174 / 0.6E1) 
     #* t36 / 0.2E1) + t181 * (-(t3 * t4 * ((t185 * ((((t187 - t62) * t1
     #0 - t64) * t10 - t192) * t10 - t197) / 0.30E2 - t204 / 0.30E2) * t
     #10 / 0.2E1 + (t204 / 0.30E2 - t185 * (t202 - (t200 - (t75 - (t73 -
     # t209) * t10) * t10) * t10) / 0.30E2) * t10 / 0.2E1) - ((t87 * t88
     # * ((-t60 * ((t225 - t90) * t10 - t93) / 0.6E1 + t231 / 0.6E1) * t
     #10 / 0.2E1 + (-t231 / 0.6E1 + t60 * (t96 - (t94 - t235) * t10) / 0
     #.6E1) * t10 / 0.2E1) - t84) * t36 - t250) * t36 / 0.6E1 + ((((t254
     # * vx(3,0) * t257 * ((t259 - t260) * t10 / 0.2E1 + (t260 - t263) *
     # t10 / 0.2E1) - t100) * t36 - t102) * t36 - t104) * t36 - t277) * 
     #t36 / 0.30E2 - t317 + t322 - t327) * t36 / 0.2E1 - (t317 - t322 + 
     #t327 - t41 * t42 * ((t185 * ((((t330 - t133) * t10 - t135) * t10 -
     # t335) * t10 - t340) / 0.30E2 - t347 / 0.30E2) * t10 / 0.2E1 + (t3
     #47 / 0.30E2 - t185 * (t345 - (t343 - (t145 - (t143 - t351) * t10) 
     #* t10) * t10) / 0.30E2) * t10 / 0.2E1) + (t319 - (t154 - t157 * t1
     #58 * ((-t60 * ((t367 - t160) * t10 - t163) / 0.6E1 + t373 / 0.6E1)
     # * t10 / 0.2E1 + (-t373 / 0.6E1 + t60 * (t166 - (t164 - t377) * t1
     #0) / 0.6E1) * t10 / 0.2E1)) * t36) * t36 / 0.6E1 - (t324 - (t174 -
     # (t172 - (t170 - t394 * vx(-3,0) * t397 * ((t399 - t400) * t10 / 0
     #.2E1 + (t400 - t403) * t10 / 0.2E1)) * t36) * t36) * t36) * t36 / 
     #0.30E2) * t36 / 0.2E1) + ((t426 - t437 + 0.3E1 / 0.256E3 * t438 * 
     #((t36 * t443 - t445) * t36 - t455) + t459 - t461 + t472) * t474 * 
     #t36 - (t459 - t461 + t472 + t477 - t479 + 0.3E1 / 0.256E3 * t438 *
     # (t469 - (-t36 * t484 + t467) * t36)) * t492 * t36) * t36 + dvx * 
     #(t497 * t507 * t36 - (((t428 / 0.2E1 - dvx * t443 / 0.16E2 + t426 
     #- t437) * t498 * t36 - t517) * t36 - t524) * t36 / 0.24E2 - t520 *
     # t533 * t36 + (t524 - (t522 - (t477 - t479 + t463 / 0.2E1 - dvx * 
     #t484 / 0.16E2) * t528 * t36) * t36) * t36 / 0.24E2) + t438 * (0.3E
     #1 / 0.640E3 * t550 * t565 * t36 - ((-t570 * t555 * t36 / 0.24E2 - 
     #t576) * t36 - t584) * t36 / 0.24E2 + 0.3E1 / 0.640E3 * (((((t428 /
     # 0.2E1 + t440 / 0.2E1) * t551 * t36 - t593) * t36 - t599) * t36 - 
     #t607) * t36 - t619) * t36 - 0.3E1 / 0.640E3 * t580 * t631 * t36 + 
     #(t584 - (t582 + t611 * t627 * t36 / 0.24E2) * t36) * t36 / 0.24E2 
     #- 0.3E1 / 0.640E3 * (t619 - (t617 - (t615 - (t613 - (t481 / 0.2E1 
     #+ t463 / 0.2E1) * t623 * t36) * t36) * t36) * t36) * t36) + ((t663
     # - t676 + 0.3E1 / 0.256E3 * t183 * ((t10 * t683 - t685) * t10 - t6
     #97) + t701 - t703 + t716) * t26 * t10 - (t701 - t703 + t716 + t720
     # - t722 + 0.3E1 / 0.256E3 * t183 * (t713 - (-t10 * t729 + t711) * 
     #t10)) * t29 * t10) * t10 + t58 * (t741 * t742 * t10 - t60 * (((t66
     #7 / 0.2E1 - dvy * t683 / 0.16E2 + t663 - t676) * t107 * t10 - t752
     #) * t10 - t759) / 0.24E2 - t755 * t763 * t10 + t60 * (t759 - (t757
     # - (t720 - t722 + t707 / 0.2E1 - dvy * t729 / 0.16E2) * t117 * t10
     #) * t10) / 0.24E2) * t10 + t181 * (0.3E1 / 0.640E3 * t781 * t299 *
     # t10 - t60 * ((t786 * (-t60 * t284 / 0.24E2 + t110 / 0.24E2) * t10
     # - t793) * t10 - t801) / 0.24E2 + 0.3E1 / 0.640E3 * t185 * (((((t6
     #67 / 0.2E1 + t680 / 0.2E1) * t282 * t10 - t810) * t10 - t816) * t1
     #0 - t824) * t10 - t836) - 0.3E1 / 0.640E3 * t797 * t311 * t10 + t6
     #0 * (t801 - (t799 - t828 * (t60 * t305 / 0.24E2 - t120 / 0.24E2) *
     # t10) * t10) / 0.24E2 - 0.3E1 / 0.640E3 * t185 * (t836 - (t834 - (
     #t832 - (t830 - (t726 / 0.2E1 + t707 / 0.2E1) * t303 * t10) * t10) 
     #* t10) * t10)) * t10 - t884 / 0.2E1 - t897 / 0.2E1 + t58 * (-(t918
     # - t60 * t931 / 0.6E1 - t941 + t944) * t10 / 0.2E1 - (t941 - t944 
     #- t966 + t60 * t979 / 0.6E1) * t10 / 0.2E1) + t181 * (-(t869 * t87
     #0 * ((((((t259 - t90) * t36 - t900) * t36 - t902) * t36 - t994) * 
     #t36 / 0.30E2 - t1000 / 0.30E2) * t36 / 0.2E1 + (t1000 / 0.30E2 - (
     #t998 - (t911 - (t909 - (t160 - t399) * t36) * t36) * t36) * t36 / 
     #0.30E2) * t36 / 0.2E1) - t60 * ((t919 * t920 * ((-((t225 - t62) * 
     #t36 - t922) * t36 / 0.6E1 + t1024 / 0.6E1) * t36 / 0.2E1 + (-t1024
     # / 0.6E1 + (t924 - (t133 - t367) * t36) * t36 / 0.6E1) * t36 / 0.2
     #E1) - t918) * t10 - t1042) / 0.6E1 + t185 * ((((t677 * t678 * vy(0
     #,3) * ((t187 - t281) * t36 / 0.2E1 + (t281 - t330) * t36 / 0.2E1) 
     #- t928) * t10 - t930) * t10 - t1060) * t10 - t1065) / 0.30E2 - t10
     #76 + t1081 - t1087) * t10 / 0.2E1 - (t1076 - t1081 + t1087 - t886 
     #* t887 * ((((((t263 - t94) * t36 - t948) * t36 - t950) * t36 - t10
     #97) * t36 / 0.30E2 - t1103 / 0.30E2) * t36 / 0.2E1 + (t1103 / 0.30
     #E2 - (t1101 - (t959 - (t957 - (t164 - t403) * t36) * t36) * t36) *
     # t36 / 0.30E2) * t36 / 0.2E1) + t60 * (t1078 - (t966 - t967 * t968
     # * ((-((t235 - t73) * t36 - t970) * t36 / 0.6E1 + t1127 / 0.6E1) *
     # t36 / 0.2E1 + (-t1127 / 0.6E1 + (t972 - (t143 - t377) * t36) * t3
     #6 / 0.6E1) * t36 / 0.2E1)) * t10) / 0.6E1 - t185 * (t1084 - (t1082
     # - (t978 - (t976 - t723 * t724 * vy(0,-3) * ((t209 - t302) * t36 /
     # 0.2E1 + (t302 - t351) * t36 / 0.2E1)) * t10) * t10) * t10) / 0.30
     #E2) * t10 / 0.2E1)

ccccc
ccccc
        rhs(i1, i2, i3, i4) = rhs(i1, i2, i3, i4)+temp
      end do
      end do
      end do
      end do
c
      return
      end
ccccc
ccccc
ccccc
ccccc
c
c ++++++++++++++
c
      subroutine nonConservativePitchAngle_4th(
     *     rhs, f,
     *     nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b,
     *     n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b,
     *     xlo, xhi, dx,
     *     velocities, range_lo, range_hi,
     *     vfloor, Vth, nuCoefficient, IVx, IVy,
     *     do_relativity)
c
c.. compute and add a simple collision operator to the RHS ...
      implicit none
c
c.. declarations of incoming variables 
      integer nd1a, nd1b, nd2a, nd2b
      integer nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
      integer do_relativity
c
      real rhs (nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real f   (nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real Vth(nd1a:nd1b, nd2a:nd2b)
      real IVx(nd1a:nd1b, nd2a:nd2b)
      real IVy(nd1a:nd1b, nd2a:nd2b)
      real range_lo(1:2), range_hi(1:2)
      real vfloor, nuCoefficient
      real xlo(1:4), xhi(1:4), dx(1:4)
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      integer solution_order, vrolloff
      real vxgrid, vygrid, vx, vy, dvx, dvy, vxmin, vxmax, vymin, vymax
      real fvxp2, fvxp1, fvxm1, fvxm2
      real fvxvx, fvyvy, fvxvy
      real fvx, fvy
      real nuei
      real temp
      real vthermal, vflowx, vflowy, vxlo, vxhi, vylo, vyhi
c
      solution_order = 4
      vrolloff = 3
      vxlo = xlo(3)
      vxhi = xhi(3)
      vylo = xlo(4)
      vyhi = xhi(4)
      dvx = dx(3)
      dvy = dx(4)
      vxmin = vxlo + vrolloff*dvx
      vxmax = vxhi - vrolloff*dvx
      vymin = vylo + vrolloff*dvy
      vymax = vyhi - vrolloff*dvy

      do i4 = n4a, n4b
      do i3 = n3a, n3b
      do i2 = n2a, n2b
      do i1 = n1a, n1b
      vthermal = Vth(i1,i2)
      vflowx = IVx(i1, i2)
      vflowy = IVy(i1, i2)
      vxgrid = velocities(i3, i4, 0)
      vygrid = velocities(i3, i4, 1)
      vx = vxgrid-vflowx
      vy = vygrid-vflowy
      ! For relativity the grid is momentum based
      if (do_relativity .eq. 1) then
        vxgrid = vxlo + (i3+0.5)*dvx
        vygrid = vylo + (i4+0.5)*dvy
      end if

      call evaluateCollisionality(
     *     nuei, vx, vy, vxgrid, vygrid, range_lo, range_hi,
     *     vxmin, vxmax, vymin, vymax,
     *     vfloor, vthermal, nuCoefficient, solution_order)

        fvxp2 = 
     *       (-1.0*f(i1, i2, i3+2, i4+2)
     *       +8.0 *f(i1, i2, i3+1, i4+2)
     *       -8.0 *f(i1, i2, i3-1, i4+2)
     *       +1.0 *f(i1, i2, i3-2, i4+2))/(12.0*dvx)

        fvxp1 = 
     *       (-1.0*f(i1, i2, i3+2, i4+1)
     *       +8.0 *f(i1, i2, i3+1, i4+1)
     *       -8.0 *f(i1, i2, i3-1, i4+1)
     *       +1.0 *f(i1, i2, i3-2, i4+1))/(12.0*dvx)

        fvxm1 = 
     *       (-1.0*f(i1, i2, i3+2, i4-1)
     *       +8.0 *f(i1, i2, i3+1, i4-1)
     *       -8.0 *f(i1, i2, i3-1, i4-1)
     *       +1.0 *f(i1, i2, i3-2, i4-1))/(12.0*dvx)

        fvxm2 = 
     *       (-1.0*f(i1, i2, i3+2, i4-2)
     *       +8.0 *f(i1, i2, i3+1, i4-2)
     *       -8.0 *f(i1, i2, i3-1, i4-2)
     *       +1.0 *f(i1, i2, i3-2, i4-2))/(12.0*dvx)

        fvxvy = 
     *       (-1.0*fvxp2
     *       +8.0 *fvxp1
     *       -8.0 *fvxm1
     *       +1.0 *fvxm2)/(12.0*dvy)

        fvxvx = 
     *       (-1.0*f(i1, i2, i3+2, i4)
     *       +16.0*f(i1, i2, i3+1, i4)
     *       -30.0*f(i1, i2, i3, i4)
     *       +16.0*f(i1, i2, i3-1, i4)
     *       -1.0 *f(i1, i2, i3-2, i4))/(12.0*(dvx)**2)
        
        fvyvy = 
     *       (-1.0*f(i1, i2, i3, i4+2)
     *       +16.0*f(i1, i2, i3, i4+1)
     *       -30.0*f(i1, i2, i3, i4)
     *       +16.0*f(i1, i2, i3, i4-1)
     *       -1.0 *f(i1, i2, i3, i4-2))/(12.0*(dvy)**2)
        
        fvx = 
     *       (-1.0*f(i1, i2, i3+2, i4)
     *        +8.0*f(i1, i2, i3+1, i4)
     *        -8.0*f(i1, i2, i3-1, i4)
     *        +1.0*f(i1, i2, i3-2, i4))/(12.0*dvx)

        fvy = 
     *       (-1.0*f(i1, i2, i3, i4+2)
     *        +8.0*f(i1, i2, i3, i4+1)
     *        -8.0*f(i1, i2, i3, i4-1)
     *        +1.0*f(i1, i2, i3, i4-2))/(12.0*dvy)
        
        temp = nuei*(
     *       vx*vx*fvyvy
     *       -2.0*vx*vy*fvxvy
     *       +vy*vy*fvxvx
     *       -vy*fvy
     *       -vx*fvx)

        rhs(i1, i2, i3, i4) = rhs(i1, i2, i3, i4)+temp
      end do
      end do
      end do
      end do
c
      return
      end
c
c ++++++++++++++
c
      subroutine appendPitchAngleCollision(
     *     rhs, f, velocities, IVx, IVy, Vth,
     *     nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b,
     *     n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b,
     *     xlo, xhi, dx,
     *     range_lo, range_hi,
     *     dparams,
     *     iparams)
c
c.. compute and add a simple pitch angle collision operator to the RHS ...
      implicit none
c
c.. declarations of incoming variables 
      integer nd1a, nd1b, nd2a, nd2b
      integer nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
c
      real rhs(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real f(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real IVx(nd1a:nd1b, nd2a:nd2b)
      real IVy(nd1a:nd1b, nd2a:nd2b)
      real Vth(nd1a:nd1b, nd2a:nd2b)
      real range_lo(1:2), range_hi(1:2)
c
      real dparams(*)
      integer iparams(*)
c
      real xlo(1:4), xhi(1:4), dx(1:4)
c
      real vfloor, nuCoefficient
      integer icons, solution_order, do_relativity
c
      vfloor = dparams(1)
      nuCoefficient = dparams(3)
      icons = iparams(1)
      solution_order = iparams(2)
      do_relativity = iparams(3)

      if (icons .eq. 1) then
        if( solution_order .eq. 4 ) then
          ! 4th order conservative implementation of collision operator
          call conservativePitchAngle_4th(
     *         rhs, f,
     *         nd1a, nd1b, nd2a, nd2b,
     *         nd3a, nd3b, nd4a, nd4b,
     *         n1a, n1b, n2a, n2b,
     *         n3a, n3b, n4a, n4b,
     *         xlo, xhi, dx,
     *         velocities,
     *         range_lo, range_hi,
     *         vfloor, Vth, nuCoefficient, IVx, IVy,
     *         do_relativity)
        else
          call conservativePitchAngle_6th(
     *         rhs, f,
     *         nd1a, nd1b, nd2a, nd2b,
     *         nd3a, nd3b, nd4a, nd4b,
     *         n1a, n1b, n2a, n2b,
     *         n3a, n3b, n4a, n4b,
     *         xlo, xhi, dx,
     *         velocities,
     *         range_lo, range_hi,
     *         vfloor, Vth, nuCoefficient, IVx, IVy,
     *         do_relativity)
        end if
      else
        if( solution_order .eq. 4 ) then
          !4th order nonconservative implementation of collision operator
          call nonConservativePitchAngle_4th(
     *         rhs, f,
     *         nd1a, nd1b, nd2a, nd2b,
     *         nd3a, nd3b, nd4a, nd4b,
     *         n1a, n1b, n2a, n2b,
     *         n3a, n3b, n4a, n4b,
     *         xlo, xhi, dx,
     *         velocities,
     *         range_lo, range_hi,
     *         vfloor, Vth, nuCoefficient, IVx, IVy,
     *         do_relativity)
        end if
      end if
c
      return
      end
c
c +++++++++++++
c
      subroutine computePitchAngleSpeciesMoments(
     *     rN, rGammax, rGammay, u,
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
      real rN(nd1a:nd1b, nd2a:nd2b)
      real rGammax(nd1a:nd1b, nd2a:nd2b)
      real rGammay(nd1a:nd1b, nd2a:nd2b)
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

      do i4 = n4a, n4b
        do i3 = n3a, n3b
          vx = velocities(i3, i4, 0)
          vy = velocities(i3, i4, 1)
          do i2 = nd2a, nd2b
            do i1 = nd1a, nd1b
              uval = max(abs(u(i1, i2, i3, i4)), eps)
              rN(i1, i2) = rN(i1, i2) + uval
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
      subroutine computePitchAngleSpeciesReducedFields(
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
      subroutine computePitchAngleSpeciesKEC(
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
      subroutine computePitchAngleSpeciesVthermal(
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
