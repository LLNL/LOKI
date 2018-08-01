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
c Fortran functions called by RosenbluthCollisionOperator.
c

      subroutine computeRosenbluthConstants(
     *     khat, hhat, vx, vy,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     nGhosts,
     *     xlo, xhi, dx,
     *     dparams, iparams)
c
c.. compute constants needed by collision operator
      implicit none
c
c.. declarations of incoming variables
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b,nGhosts
c
      real khat(n3a-nGhosts:n3b+nGhosts, n4a-nGhosts:n4b+nGhosts)
      real hhat(n3a-nGhosts:n3b+nGhosts, n4a-nGhosts:n4b+nGhosts)
      real vx(n3a-nGhosts:n3b+nGhosts, n4a-nGhosts:n4b+nGhosts)
      real vy(n3a-nGhosts:n3b+nGhosts, n4a-nGhosts:n4b+nGhosts)
c
      real dparams(*)
      integer iparams(*)
c
      real xlo(1:4), xhi(1:4), dx(1:4)
c
c.. declarations of local variables
      integer i3, i4, solution_order
      real vrolloff
      real erfwrap
      real vx0, vy0, dvx, dvy
      real vxt, vyt, vt, v, vmax, vceil, vfloor, vTh, tmp
      real vflowx, vflowy, Kh, Hh, A, pi
c
      solution_order = iparams(1)
      if (solution_order .eq. 4) then
        vrolloff = 3.0
      else
        vrolloff = 4.0
      end if
      pi = 4.0*atan(1.0)
      vceil = dparams(1)
      vfloor = dparams(2)
      vTh = dparams(3)
      vflowx = dparams(8)
      vflowy = dparams(9)
      vx0 = xlo(3)
      vy0 = xlo(4)
      dvx = dx(3)
      dvy = dx(4)
      vmax = min(max(abs(xhi(3)), abs(xlo(3)))-vrolloff*dvx,
     *           max(abs(xhi(4)), abs(xlo(4)))-vrolloff*dvy)

      do i4 = n4a-nGhosts, n4b+nGhosts
        do i3 = n3a-nGhosts, n3b+nGhosts
          vxt = vx0+(0.5+i3)*dvx-vflowx
          vyt = vy0+(0.5+i4)*dvy-vflowy
          vt  = sqrt(vxt**2+vyt**2)
          v   = max(vfloor, vt)
          call getRamp(A, v, vceil, vmax, pi)
          v   = v/vTh
          tmp = erfwrap(v/sqrt(2.0))
          Kh = 1.0/(v**3)*
     *      ((v**2-1.0)*tmp+sqrt(2.0/pi)*v*exp(-v**2/2.0))
          khat(i3, i4) = A*Kh/(v**2)
          Hh = 1.0/(v**3)*(tmp-sqrt(2.0/pi)*v*exp(-(v**2/2.0)))
          hhat(i3, i4) = A*Hh/(v**2)
          vx(i3, i4) = vxt
          vy(i3, i4) = vyt
        end do
      end do

      return
      end
c
c +++++++++++++
c
      subroutine computeRosenbluthTemps(
     *     fM, uMomx, uMomy, uKE,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     xlo, xhi, dx,
     *     dparams)
c
c.. compute temps needed for back reaction terms
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
c
      real fM(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real uMomx(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real uMomy(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real uKE(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
c
      real dparams(*)
c
      real xlo(1:4), xhi(1:4), dx(1:4)
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      real alpha, pi, vx0, vy0, dvx, dvy, vx, vy, fMTemp
      real vflowx, vflowy
c
      alpha = dparams(6)
      vflowx = dparams(8)
      vflowy = dparams(9)
      pi = 4.0*atan(1.0)
      vx0 = xlo(3)
      vy0 = xlo(4)
      dvx = dx(3)
      dvy = dx(4)

      do i4 = nd4a, nd4b
        vy = vy0 + (0.5+i4)*dvy-vflowy
        do i3 = nd3a, nd3b
          vx = vx0 + (0.5+i3)*dvx-vflowx
          do i2 = nd2a, nd2b
            do i1 = nd1a, nd1b
               fMTemp = alpha/(2.0*pi)*exp(-0.5*alpha*(vx**2+vy**2))
               fM(i1, i2, i3, i4) = fMTemp
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
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     xlo, xhi, dx,
     *     dparams)
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
c
      real xlo(1:4), xhi(1:4), dx(1:4)
c
      real dparams(*)
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      real vx0, vy0, dvx, dvy, vx, vy, vflowx, vflowy
c
      vx0 = xlo(3)
      vy0 = xlo(4)
      dvx = dx(3)
      dvy = dx(4)
      vflowx = dparams(8)
      vflowy = dparams(9)

      do i4 = n4a, n4b
        vy = vy0 + (0.5+i4)*dvy - vflowy
        do i3 = n3a, n3b
          vx = vx0 + (0.5+i3)*dvx - vflowx
          do i2 = nd2a, nd2b
            do i1 = nd1a, nd1b
              rMomx(i1, i2) = rMomx(i1, i2) + vx*cMomx(i1, i2, i3, i4)
              rMomy(i1, i2) = rMomy(i1, i2) + vy*cMomy(i1, i2, i3, i4)
              rKE(i1, i2) =
     *          rKE(i1, i2) + (vx**2+vy**2)*cKE(i1, i2, i3, i4)
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
     *     cMomx, cMomy, cKE, fM, u,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b)
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
      real fM(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real u(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      real temp
c
      do i4 = n4a, n4b
        do i3 = n3a, n3b
          do i2 = nd2a, nd2b
            do i1 = nd1a, nd1b
              temp = u(i1, i2, i3, i4) / fM(i1, i2, i3, i4)
              rMomx(i1, i2) =
     *          rMomx(i1, i2) + temp*cMomx(i1, i2, i3, i4)
              rMomy(i1, i2) =
     *          rMomy(i1, i2) + temp*cMomy(i1, i2, i3, i4)
              rKE(i1, i2) = rKE(i1, i2) + temp*cKE(i1, i2, i3, i4)
            end do
          end do
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
     *     dparams)
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
      real dparams(*)
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      real massR
c
      massR = dparams(7)
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
      subroutine appendRosenbluthCollision(
     *     rhs,f,vx,vy,Kh,Hh,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     xlo, xhi, dx,
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
      real vx(n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2))
      real vy(n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2))
      real Kh(n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2))
      real Hh(n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2))
c
      real xlo(1:4), xhi(1:4), dx(1:4)
c
c.. declarations of local variables
      integer solution_order
c
      solution_order = iparams(1)
c
      if( solution_order .eq. 4 ) then
        ! 4th order
        call appendRosenbluthCollision_4th(
     *     rhs,f,vx,vy,Kh,Hh,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     xlo, xhi, dx,
     *     dparams,
     *     iparams)
      else
        ! 6th order
        call appendRosenbluthCollision_6th(
     *     rhs,f,vx,vy,Kh,Hh,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     xlo, xhi, dx,
     *     dparams,
     *     iparams)
      end if
c
c
      return
      end
c
c +++++++++++++
c
      subroutine appendRosenbluthCollision_4th(
     *     rhs,f,vx,vy,Kh,Hh,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     xlo, xhi, dx,
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
      real vx(n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2))
      real vy(n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2))
      real Kh(n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2))
      real Hh(n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2))
c
      real xlo(1:4), xhi(1:4), dx(1:4)
c
c.. declarations of local variables
      real nu, alpha
      integer i1, i2, i3, i4
      real temp
c
      real dvx, dvy
c
      real t1
      real t10
      real t100
      real t104
      real t106
      real t107
      real t109
      real t110
      real t111
      real t112
      real t116
      real t117
      real t118
      real t122
      real t125
      integer t126
      real t129
      real t131
      real t135
      real t137
      real t14
      real t140
      real t142
      real t145
      real t147
      real t151
      real t153
      real t155
      real t159
      real t161
      real t163
      real t165
      integer t166
      real t169
      real t17
      real t171
      real t175
      real t177
      integer t18
      real t183
      integer t187
      integer t19
      real t190
      real t192
      real t196
      real t198
      real t2
      real t210
      real t211
      real t212
      real t213
      real t217
      real t218
      real t219
      real t22
      real t223
      integer t227
      real t230
      real t232
      real t236
      real t238
      real t24
      real t241
      real t246
      real t248
      real t252
      real t254
      real t267
      real t268
      real t269
      real t270
      real t274
      real t276
      real t277
      real t28
      real t282
      real t283
      real t287
      real t288
      real t289
      real t290
      real t294
      real t296
      real t297
      real t3
      real t30
      real t301
      real t307
      real t310
      real t315
      real t317
      real t321
      real t323
      real t329
      real t33
      real t335
      real t337
      real t341
      real t343
      real t35
      real t355
      real t356
      real t357
      real t358
      real t36
      real t362
      real t363
      real t364
      real t368
      integer t37
      integer t372
      real t375
      real t377
      real t381
      real t383
      real t386
      real t391
      real t393
      real t397
      real t399
      real t4
      real t40
      real t42
      real t421
      real t428
      real t436
      real t438
      real t442
      real t448
      real t450
      real t453
      real t455
      real t456
      real t46
      real t465
      real t472
      real t473
      real t48
      real t480
      real t481
      real t483
      real t486
      real t488
      real t489
      real t495
      real t499
      real t50
      real t504
      real t506
      real t508
      real t512
      real t516
      real t522
      real t532
      real t533
      real t536
      real t54
      real t540
      real t542
      real t543
      real t546
      real t55
      real t550
      real t558
      real t56
      real t560
      real t564
      real t565
      real t566
      real t569
      real t57
      real t573
      real t575
      real t578
      real t58
      real t580
      real t584
      real t585
      real t588
      real t592
      real t599
      real t606
      real t612
      real t616
      real t62
      real t621
      real t623
      real t625
      real t63
      real t634
      real t64
      real t650
      real t652
      real t654
      real t658
      real t662
      real t664
      real t671
      real t673
      real t675
      real t679
      real t68
      real t681
      real t688
      real t71
      real t729
      real t73
      real t732
      real t74
      real t740
      real t75
      real t76
      real t8
      real t80
      real t82
      real t85
      real t86
      real t87
      real t89
      real t9
      real t90
      real t91
      real t92
      real t96
      real t98
      real t99
c
      nu = dparams(5)
      alpha = dparams(6)
c
      dvx = dx(3)
      dvy = dx(4)
c
      do i4 = n4a, n4b
      do i3 = n3a, n3b
      do i2 = n2a, n2b
      do i1 = n1a, n1b
ccccc
cc generated code 
ccccc
        t1 = vx(i3+1,i4)
        t2 = t1 ** 2
        t3 = vy(i3+1,i4)
        t4 = t3 ** 2
        t8 = exp(-alpha * (t2 + t4) / 0.2E1)
        t9 = nu * t8
        t10 = Hh(i3+1,i4)
        t14 = Kh(i3+1,i4)
        t17 = 0.2E1 * t1 * t10 * t3 - t1 * t14 * t3
        t18 = i3 + 1
        t19 = i4 + 1
        t22 = vx(i3+1,i4+1) ** 2
        t24 = vy(i3+1,i4+1) ** 2
        t28 = exp(-alpha * (t22 + t24) / 0.2E1)
        t30 = f(i1,i2,t18,t19) / t28
        t33 = f(i1,i2,t18,i4) / t8
        t35 = 0.1E1 / dvy
        t36 = (t30 - t33) * t35
        t37 = i4 - 1
        t40 = vx(i3+1,i4-1) ** 2
        t42 = vy(i3+1,i4-1) ** 2
        t46 = exp(-alpha * (t40 + t42) / 0.2E1)
        t48 = f(i1,i2,t18,t37) / t46
        t50 = (t33 - t48) * t35
        t54 = t9 * t17 * (t36 / 0.2E1 + t50 / 0.2E1)
        t55 = vx(i3,i4)
        t56 = t55 ** 2
        t57 = vy(i3,i4)
        t58 = t57 ** 2
        t62 = exp(-alpha * (t56 + t58) / 0.2E1)
        t63 = nu * t62
        t64 = Hh(i3,i4)
        t68 = Kh(i3,i4)
        t71 = 0.2E1 * t55 * t57 * t64 - t55 * t57 * t68
        t73 = vx(i3,i4+1)
        t74 = t73 ** 2
        t75 = vy(i3,i4+1)
        t76 = t75 ** 2
        t80 = exp(-alpha * (t74 + t76) / 0.2E1)
        t82 = f(i1,i2,i3,t19) / t80
        t85 = f(i1,i2,i3,i4) / t62
        t86 = t82 - t85
        t87 = t86 * t35
        t89 = vx(i3,i4-1)
        t90 = t89 ** 2
        t91 = vy(i3,i4-1)
        t92 = t91 ** 2
        t96 = exp(-alpha * (t90 + t92) / 0.2E1)
        t98 = f(i1,i2,i3,t37) / t96
        t99 = t85 - t98
        t100 = t99 * t35
        t104 = t63 * t71 * (t87 / 0.2E1 + t100 / 0.2E1)
        t106 = 0.1E1 / dvx
        t107 = (t54 - t104) * t106
        t109 = vx(i3-1,i4)
        t110 = t109 ** 2
        t111 = vy(i3-1,i4)
        t112 = t111 ** 2
        t116 = exp(-alpha * (t110 + t112) / 0.2E1)
        t117 = nu * t116
        t118 = Hh(i3-1,i4)
        t122 = Kh(i3-1,i4)
        t125 = 0.2E1 * t109 * t111 * t118 - t109 * t111 * t122
        t126 = i3 - 1
        t129 = vx(i3-1,i4+1) ** 2
        t131 = vy(i3-1,i4+1) ** 2
        t135 = exp(-alpha * (t129 + t131) / 0.2E1)
        t137 = f(i1,i2,t126,t19) / t135
        t140 = f(i1,i2,t126,i4) / t116
        t142 = (t137 - t140) * t35
        t145 = vx(i3-1,i4-1) ** 2
        t147 = vy(i3-1,i4-1) ** 2
        t151 = exp(-alpha * (t145 + t147) / 0.2E1)
        t153 = f(i1,i2,t126,t37) / t151
        t155 = (t140 - t153) * t35
        t159 = t117 * t125 * (t142 / 0.2E1 + t155 / 0.2E1)
        t161 = (t104 - t159) * t106
        t163 = dvx ** 2
        t165 = dvy / t163
        t166 = i4 + 2
        t169 = vx(i3+1,i4+2) ** 2
        t171 = vy(i3+1,i4+2) ** 2
        t175 = exp(-alpha * (t169 + t171) / 0.2E1)
        t177 = f(i1,i2,t18,t166) / t175
        t183 = t165 * (t36 - t50)
        t187 = i4 - 2
        t190 = vx(i3+1,i4-2) ** 2
        t192 = vy(i3+1,i4-2) ** 2
        t196 = exp(-alpha * (t190 + t192) / 0.2E1)
        t198 = f(i1,i2,t18,t187) / t196
        t210 = vx(i3+2,i4)
        t211 = t210 ** 2
        t212 = vy(i3+2,i4)
        t213 = t212 ** 2
        t217 = exp(-alpha * (t211 + t213) / 0.2E1)
        t218 = nu * t217
        t219 = Hh(i3+2,i4)
        t223 = Kh(i3+2,i4)
        t227 = i3 + 2
        t230 = vx(i3+2,i4+1) ** 2
        t232 = vy(i3+2,i4+1) ** 2
        t236 = exp(-alpha * (t230 + t232) / 0.2E1)
        t238 = f(i1,i2,t227,t19) / t236
        t241 = f(i1,i2,t227,i4) / t217
        t246 = vx(i3+2,i4-1) ** 2
        t248 = vy(i3+2,i4-1) ** 2
        t252 = exp(-alpha * (t246 + t248) / 0.2E1)
        t254 = f(i1,i2,t227,t37) / t252
        t267 = vx(i3,i4+2)
        t268 = t267 ** 2
        t269 = vy(i3,i4+2)
        t270 = t269 ** 2
        t274 = exp(-alpha * (t268 + t270) / 0.2E1)
        t276 = f(i1,i2,i3,t166) / t274
        t277 = t276 - t82
        t282 = t165 * (t87 - t100)
        t283 = -t165 * (t277 * t35 - t87) + t282
        t287 = vx(i3,i4-2)
        t288 = t287 ** 2
        t289 = vy(i3,i4-2)
        t290 = t289 ** 2
        t294 = exp(-alpha * (t288 + t290) / 0.2E1)
        t296 = f(i1,i2,i3,t187) / t294
        t297 = t98 - t296
        t301 = -t282 + t165 * (-t297 * t35 + t100)
        t307 = t63 * t71 * (t283 * t35 / 0.12E2 + t301 * t35 / 0.12E2)
        t310 = (t107 - t161) * t106 / 0.6E1
        t315 = vx(i3-1,i4+2) ** 2
        t317 = vy(i3-1,i4+2) ** 2
        t321 = exp(-alpha * (t315 + t317) / 0.2E1)
        t323 = f(i1,i2,t126,t166) / t321
        t329 = t165 * (t142 - t155)
        t335 = vx(i3-1,i4-2) ** 2
        t337 = vy(i3-1,i4-2) ** 2
        t341 = exp(-alpha * (t335 + t337) / 0.2E1)
        t343 = f(i1,i2,t126,t187) / t341
        t355 = vx(i3-2,i4)
        t356 = t355 ** 2
        t357 = vy(i3-2,i4)
        t358 = t357 ** 2
        t362 = exp(-alpha * (t356 + t358) / 0.2E1)
        t363 = nu * t362
        t364 = Hh(i3-2,i4)
        t368 = Kh(i3-2,i4)
        t372 = i3 - 2
        t375 = vx(i3-2,i4+1) ** 2
        t377 = vy(i3-2,i4+1) ** 2
        t381 = exp(-alpha * (t375 + t377) / 0.2E1)
        t383 = f(i1,i2,t372,t19) / t381
        t386 = f(i1,i2,t372,i4) / t362
        t391 = vx(i3-2,i4-1) ** 2
        t393 = vy(i3-2,i4-1) ** 2
        t397 = exp(-alpha * (t391 + t393) / 0.2E1)
        t399 = f(i1,i2,t372,t37) / t397
        t421 = sqrt(t9 * (0.2E1 * t10 * t2 + t14 * t4))
        t428 = sqrt(t218 * (0.2E1 * t211 * t219 + t213 * t223))
        t436 = sqrt(t63 * (0.2E1 * t56 * t64 + t58 * t68))
        t438 = (t421 - t436) * t106
        t442 = t436 / 0.2E1
        t448 = sqrt(t117 * (0.2E1 * t110 * t118 + t112 * t122))
        t450 = (t436 - t448) * t106
        t453 = dvx * (t438 - t450) / 0.16E2
        t455 = (t421 / 0.2E1 - dvx * ((t428 - t421) * t106 - t438) / 0.1
     #6E2 + t442 - t453) ** 2
        t456 = t33 - t85
        t465 = sqrt(t363 * (0.2E1 * t356 * t364 + t358 * t368))
        t472 = (t442 - t453 + t448 / 0.2E1 - dvx * (t450 - (t448 - t465)
     # * t106) / 0.16E2) ** 2
        t473 = t85 - t140
        t480 = (t421 / 0.2E1 + t436 / 0.2E1) ** 2
        t481 = t241 - t33
        t483 = t456 * t106
        t486 = t473 * t106
        t488 = (t483 - t486) * t106
        t489 = -(t106 * t481 - t483) * t106 + t488
        t495 = (t428 / 0.2E1 + t421 / 0.2E1) ** 2
        t499 = t480 * t456 * t106
        t504 = (t436 / 0.2E1 + t448 / 0.2E1) ** 2
        t506 = t504 * t473 * t106
        t508 = (t499 - t506) * t106
        t512 = t140 - t386
        t516 = -t488 + (-t106 * t512 + t486) * t106
        t522 = (t448 / 0.2E1 + t465 / 0.2E1) ** 2
        t532 = nu * t80
        t533 = Hh(i3,i4+1)
        t536 = Kh(i3,i4+1)
        t540 = sqrt(t532 * (0.2E1 * t533 * t76 + t536 * t74))
        t542 = nu * t274
        t543 = Hh(i3,i4+2)
        t546 = Kh(i3,i4+2)
        t550 = sqrt(t542 * (t268 * t546 + 0.2E1 * t270 * t543))
        t558 = sqrt(t63 * (t56 * t68 + 0.2E1 * t58 * t64))
        t560 = (t540 - t558) * t35
        t564 = t558 / 0.2E1
        t565 = nu * t96
        t566 = Hh(i3,i4-1)
        t569 = Kh(i3,i4-1)
        t573 = sqrt(t565 * (0.2E1 * t566 * t92 + t569 * t90))
        t575 = (t558 - t573) * t35
        t578 = dvy * (t560 - t575) / 0.16E2
        t580 = (t540 / 0.2E1 - dvy * ((t550 - t540) * t35 - t560) / 0.16
     #E2 + t564 - t578) ** 2
        t584 = nu * t294
        t585 = Hh(i3,i4-2)
        t588 = Kh(i3,i4-2)
        t592 = sqrt(t584 * (t288 * t588 + 0.2E1 * t290 * t585))
        t599 = (t564 - t578 + t573 / 0.2E1 - dvy * (t575 - (t573 - t592)
     # * t35) / 0.16E2) ** 2
        t606 = (t540 / 0.2E1 + t558 / 0.2E1) ** 2
        t612 = (t550 / 0.2E1 + t540 / 0.2E1) ** 2
        t616 = t606 * t86 * t35
        t621 = (t558 / 0.2E1 + t573 / 0.2E1) ** 2
        t623 = t621 * t99 * t35
        t625 = (t616 - t623) * t35
        t634 = (t573 / 0.2E1 + t592 / 0.2E1) ** 2
        t650 = 0.2E1 * t533 * t73 * t75 - t536 * t73 * t75
        t652 = (t30 - t82) * t106
        t654 = (t82 - t137) * t106
        t658 = t532 * t650 * (t652 / 0.2E1 + t654 / 0.2E1)
        t662 = t63 * t71 * (t483 / 0.2E1 + t486 / 0.2E1)
        t664 = (t658 - t662) * t35
        t671 = 0.2E1 * t566 * t89 * t91 - t569 * t89 * t91
        t673 = (t48 - t98) * t106
        t675 = (t98 - t153) * t106
        t679 = t565 * t671 * (t673 / 0.2E1 + t675 / 0.2E1)
        t681 = (t662 - t679) * t35
        t688 = (t652 - t654) * t106
        t729 = t63 * t71 * (t489 * t106 / 0.12E2 + t516 * t106 / 0.12E2)
        t732 = t165 * (t664 - t681) / 0.6E1
        t740 = (t673 - t675) * t106
        temp = t107 / 0.2E1 + t161 / 0.2E1 + t163 * ((t9 * t17 * ((-t165
     # * ((t177 - t30) * t35 - t36) / 0.6E1 + t183 / 0.6E1) * t35 / 0.2E
     #1 + (-t183 / 0.6E1 + t165 * (t50 - (t48 - t198) * t35) / 0.6E1) * 
     #t35 / 0.2E1) - ((t218 * (0.2E1 * t210 * t212 * t219 - t210 * t212 
     #* t223) * ((t238 - t241) * t35 / 0.2E1 + (t241 - t254) * t35 / 0.2
     #E1) - t54) * t106 - t107) * t106 / 0.6E1 - t307 + t310) * t106 / 0
     #.2E1 + (t307 - t310 - t117 * t125 * ((-t165 * ((t323 - t137) * t35
     # - t142) / 0.6E1 + t329 / 0.6E1) * t35 / 0.2E1 + (-t329 / 0.6E1 + 
     #t165 * (t155 - (t153 - t343) * t35) / 0.6E1) * t35 / 0.2E1) + (t16
     #1 - (t159 - t363 * (0.2E1 * t355 * t357 * t364 - t355 * t357 * t36
     #8) * ((t383 - t386) * t35 / 0.2E1 + (t386 - t399) * t35 / 0.2E1)) 
     #* t106) * t106 / 0.6E1) * t106 / 0.2E1) + (t106 * t455 * t456 - t1
     #06 * t472 * t473) * t106 + dvx * (t480 * t489 * t106 / 0.24E2 - ((
     #t106 * t481 * t495 - t499) * t106 - t508) * t106 / 0.24E2 - t504 *
     # t516 * t106 / 0.24E2 + (t508 - (-t106 * t512 * t522 + t506) * t10
     #6) * t106 / 0.24E2) + (t35 * t580 * t86 - t35 * t599 * t99) * t35 
     #+ t163 * (t606 * t283 * t35 / 0.24E2 - t165 * ((t277 * t35 * t612 
     #- t616) * t35 - t625) / 0.24E2 - t621 * t301 * t35 / 0.24E2 + t165
     # * (t625 - (-t297 * t35 * t634 + t623) * t35) / 0.24E2) * t35 + t6
     #64 / 0.2E1 + t681 / 0.2E1 + t163 * ((t532 * t650 * ((-((t238 - t30
     #) * t106 - t652) * t106 / 0.6E1 + t688 / 0.6E1) * t106 / 0.2E1 + (
     #-t688 / 0.6E1 + (t654 - (t137 - t383) * t106) * t106 / 0.6E1) * t1
     #06 / 0.2E1) - t165 * ((t542 * (0.2E1 * t267 * t269 * t543 - t267 *
     # t269 * t546) * ((t177 - t276) * t106 / 0.2E1 + (t276 - t323) * t1
     #06 / 0.2E1) - t658) * t35 - t664) / 0.6E1 - t729 + t732) * t35 / 0
     #.2E1 + (t729 - t732 - t565 * t671 * ((-((t254 - t48) * t106 - t673
     #) * t106 / 0.6E1 + t740 / 0.6E1) * t106 / 0.2E1 + (-t740 / 0.6E1 +
     # (t675 - (t153 - t399) * t106) * t106 / 0.6E1) * t106 / 0.2E1) + t
     #165 * (t681 - (t679 - t584 * (0.2E1 * t287 * t289 * t585 - t287 * 
     #t289 * t588) * ((t198 - t296) * t106 / 0.2E1 + (t296 - t343) * t10
     #6 / 0.2E1)) * t35) / 0.6E1) * t35 / 0.2E1)
ccccc
ccccc
        rhs(i1,i2,i3,i4) = rhs(i1,i2,i3,i4)+temp
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
      subroutine appendRosenbluthCollision_6th(
     *     rhs,f,vx,vy,Kh,Hh,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     xlo, xhi, dx,
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
      real vx(n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2))
      real vy(n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2))
      real Kh(n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2))
      real Hh(n3a-iparams(2):n3b+iparams(2),
     #n4a-iparams(2):n4b+iparams(2))
c
      real xlo(1:4), xhi(1:4), dx(1:4)
c
c.. declarations of local variables
      real nu, alpha
      integer i1, i2, i3, i4
      real temp
c
      real dvx, dvy
c
c generated declarations
        real t1
        real t10
        real t100
        real t1002
        real t1006
        real t1008
        real t1010
        real t1019
        real t1023
        real t1025
        real t1029
        real t1031
        real t1033
        real t1037
        real t1039
        real t104
        real t1041
        real t1043
        real t1045
        real t1049
        real t1053
        real t1057
        real t106
        real t107
        real t1085
        real t1086
        real t109
        real t1092
        real t1093
        real t1094
        real t1095
        real t110
        real t1101
        real t1103
        real t1109
        real t111
        real t1111
        real t1112
        real t1114
        real t1115
        real t1116
        real t112
        real t1122
        real t1125
        real t1127
        real t1130
        real t1131
        real t1137
        real t1139
        real t1140
        real t1141
        real t1143
        real t1147
        real t1149
        real t1150
        real t1151
        real t1157
        real t1159
        real t116
        real t1160
        real t1161
        real t1163
        real t1166
        real t117
        real t1170
        real t1172
        real t1173
        real t1174
        real t118
        real t1180
        real t1183
        real t1195
        real t1196
        real t1206
        real t1209
        real t1211
        real t1213
        real t1217
        real t122
        real t1235
        real t1240
        real t1247
        real t125
        real t1251
        real t1253
        real t1255
        integer t126
        real t1264
        real t1268
        real t1270
        real t1274
        real t1276
        real t1278
        real t1282
        real t1284
        real t1286
        real t1288
        real t129
        real t1290
        real t131
        real t1329
        real t1331
        real t1333
        real t1337
        real t1341
        real t1343
        real t135
        real t1351
        real t1353
        real t1355
        real t1359
        real t1361
        real t1364
        real t1366
        real t1368
        real t1369
        real t137
        real t1373
        real t1375
        real t1376
        real t1382
        real t1389
        real t1391
        real t1393
        real t1397
        real t1399
        real t14
        real t140
        real t1400
        real t1410
        real t1411
        real t1413
        real t1417
        real t1419
        real t142
        real t1421
        real t1422
        real t1426
        real t1428
        real t1429
        real t1435
        real t1442
        real t1444
        real t1446
        real t145
        real t1450
        real t1452
        real t1453
        real t1468
        real t147
        real t1472
        real t1474
        real t1498
        real t151
        real t1516
        real t153
        real t1539
        real t1542
        real t1544
        real t155
        real t1555
        real t1557
        real t1560
        real t1561
        real t1563
        real t1566
        real t1576
        real t1580
        real t1582
        real t159
        real t1606
        real t161
        real t163
        real t165
        integer t166
        real t169
        real t17
        real t171
        real t175
        real t177
        real t179
        integer t18
        real t180
        real t182
        real t183
        integer t187
        integer t19
        real t190
        real t192
        real t196
        real t198
        real t2
        real t200
        real t201
        real t209
        real t210
        real t211
        real t212
        real t213
        real t217
        real t218
        real t219
        real t22
        real t223
        real t226
        integer t227
        real t230
        real t232
        real t236
        real t238
        real t24
        real t241
        real t243
        real t246
        real t248
        real t252
        real t254
        real t256
        real t260
        real t262
        real t264
        real t267
        real t268
        real t269
        real t270
        real t274
        real t276
        real t277
        real t278
        real t279
        real t28
        real t280
        real t281
        real t282
        real t283
        real t287
        real t288
        real t289
        real t290
        real t294
        real t296
        real t297
        real t298
        real t299
        real t3
        real t30
        real t300
        real t301
        real t307
        real t309
        real t310
        real t315
        real t317
        real t321
        real t323
        real t325
        real t326
        real t328
        real t329
        real t33
        real t335
        real t337
        real t341
        real t343
        real t345
        real t346
        real t35
        real t354
        real t355
        real t356
        real t357
        real t358
        real t36
        real t362
        real t363
        real t364
        real t368
        integer t37
        real t371
        integer t372
        real t375
        real t377
        real t381
        real t383
        real t386
        real t388
        real t391
        real t393
        real t397
        real t399
        real t4
        real t40
        real t401
        real t405
        real t407
        real t409
        real t416
        real t417
        real t418
        real t42
        real t420
        integer t421
        real t424
        real t426
        real t430
        real t432
        real t437
        real t440
        real t442
        real t445
        real t447
        real t449
        integer t453
        real t456
        real t458
        real t46
        real t462
        real t464
        real t48
        real t482
        real t484
        real t488
        real t490
        real t496
        real t50
        real t502
        real t504
        real t508
        real t510
        real t525
        real t529
        real t530
        real t531
        real t532
        real t536
        real t537
        real t538
        real t54
        real t542
        integer t546
        real t549
        real t55
        real t551
        real t555
        real t557
        real t56
        real t560
        real t565
        real t567
        real t57
        real t571
        real t573
        real t58
        real t587
        real t592
        real t593
        real t594
        real t595
        real t599
        real t601
        real t602
        real t604
        real t606
        real t609
        real t611
        real t614
        real t616
        real t618
        real t619
        real t62
        real t623
        real t624
        real t625
        real t626
        real t63
        real t630
        real t632
        real t633
        real t635
        real t64
        real t641
        real t647
        real t649
        real t652
        real t654
        real t657
        real t662
        real t664
        real t668
        real t670
        real t675
        real t678
        real t68
        real t680
        real t683
        real t685
        real t687
        real t693
        real t695
        real t699
        real t701
        real t71
        real t719
        real t721
        real t725
        real t727
        real t73
        real t733
        real t739
        real t74
        real t741
        real t745
        real t747
        real t75
        real t76
        real t764
        real t765
        real t766
        real t767
        real t771
        real t772
        real t773
        real t777
        integer t781
        real t784
        real t786
        real t790
        real t792
        real t795
        real t8
        real t80
        real t800
        real t802
        real t806
        real t808
        real t82
        real t833
        real t834
        real t839
        real t841
        real t846
        real t848
        real t849
        real t85
        real t851
        real t852
        real t857
        real t86
        real t860
        real t862
        real t869
        real t87
        real t871
        real t872
        real t873
        real t875
        real t879
        real t881
        real t886
        real t888
        real t889
        real t89
        real t890
        real t892
        real t895
        real t897
        real t9
        real t90
        real t900
        real t902
        real t907
        real t91
        real t910
        real t918
        real t92
        real t923
        real t924
        real t925
        real t926
        real t928
        real t929
        real t931
        real t932
        real t933
        real t943
        real t946
        real t948
        real t950
        real t954
        real t955
        real t957
        real t958
        real t959
        real t96
        real t976
        real t977
        real t98
        real t981
        real t984
        real t988
        real t99
        real t990
        real t991
        real t996
c
      nu = dparams(5)
      alpha = dparams(6)
c
      dvx = dx(3)
      dvy = dx(4)
c
      do i4 = n4a, n4b
      do i3 = n3a, n3b
      do i2 = n2a, n2b
      do i1 = n1a, n1b
ccccc
cc generated code 
ccccc
        t1 = vx(i3+1,i4)
        t2 = t1 ** 2
        t3 = vy(i3+1,i4)
        t4 = t3 ** 2
        t8 = exp(-alpha * (t2 + t4) / 0.2E1)
        t9 = nu * t8
        t10 = Hh(i3+1,i4)
        t14 = Kh(i3+1,i4)
        t17 = 0.2E1 * t1 * t10 * t3 - t1 * t14 * t3
        t18 = i3 + 1
        t19 = i4 + 1
        t22 = vx(i3+1,i4+1) ** 2
        t24 = vy(i3+1,i4+1) ** 2
        t28 = exp(-alpha * (t22 + t24) / 0.2E1)
        t30 = f(i1,i2,t18,t19) / t28
        t33 = f(i1,i2,t18,i4) / t8
        t35 = 0.1E1 / dvy
        t36 = (t30 - t33) * t35
        t37 = i4 - 1
        t40 = vx(i3+1,i4-1) ** 2
        t42 = vy(i3+1,i4-1) ** 2
        t46 = exp(-alpha * (t40 + t42) / 0.2E1)
        t48 = f(i1,i2,t18,t37) / t46
        t50 = (t33 - t48) * t35
        t54 = t9 * t17 * (t36 / 0.2E1 + t50 / 0.2E1)
        t55 = vx(i3,i4)
        t56 = t55 ** 2
        t57 = vy(i3,i4)
        t58 = t57 ** 2
        t62 = exp(-alpha * (t56 + t58) / 0.2E1)
        t63 = nu * t62
        t64 = Hh(i3,i4)
        t68 = Kh(i3,i4)
        t71 = 0.2E1 * t55 * t57 * t64 - t55 * t57 * t68
        t73 = vx(i3,i4+1)
        t74 = t73 ** 2
        t75 = vy(i3,i4+1)
        t76 = t75 ** 2
        t80 = exp(-alpha * (t74 + t76) / 0.2E1)
        t82 = f(i1,i2,i3,t19) / t80
        t85 = f(i1,i2,i3,i4) / t62
        t86 = t82 - t85
        t87 = t86 * t35
        t89 = vx(i3,i4-1)
        t90 = t89 ** 2
        t91 = vy(i3,i4-1)
        t92 = t91 ** 2
        t96 = exp(-alpha * (t90 + t92) / 0.2E1)
        t98 = f(i1,i2,i3,t37) / t96
        t99 = t85 - t98
        t100 = t99 * t35
        t104 = t63 * t71 * (t87 / 0.2E1 + t100 / 0.2E1)
        t106 = 0.1E1 / dvx
        t107 = (t54 - t104) * t106
        t109 = vx(i3-1,i4)
        t110 = t109 ** 2
        t111 = vy(i3-1,i4)
        t112 = t111 ** 2
        t116 = exp(-alpha * (t110 + t112) / 0.2E1)
        t117 = nu * t116
        t118 = Hh(i3-1,i4)
        t122 = Kh(i3-1,i4)
        t125 = 0.2E1 * t109 * t111 * t118 - t109 * t111 * t122
        t126 = i3 - 1
        t129 = vx(i3-1,i4+1) ** 2
        t131 = vy(i3-1,i4+1) ** 2
        t135 = exp(-alpha * (t129 + t131) / 0.2E1)
        t137 = f(i1,i2,t126,t19) / t135
        t140 = f(i1,i2,t126,i4) / t116
        t142 = (t137 - t140) * t35
        t145 = vx(i3-1,i4-1) ** 2
        t147 = vy(i3-1,i4-1) ** 2
        t151 = exp(-alpha * (t145 + t147) / 0.2E1)
        t153 = f(i1,i2,t126,t37) / t151
        t155 = (t140 - t153) * t35
        t159 = t117 * t125 * (t142 / 0.2E1 + t155 / 0.2E1)
        t161 = (t104 - t159) * t106
        t163 = dvx ** 2
        t165 = dvy / t163
        t166 = i4 + 2
        t169 = vx(i3+1,i4+2) ** 2
        t171 = vy(i3+1,i4+2) ** 2
        t175 = exp(-alpha * (t169 + t171) / 0.2E1)
        t177 = f(i1,i2,t18,t166) / t175
        t179 = (t177 - t30) * t35
        t180 = t179 - t36
        t182 = t36 - t50
        t183 = t165 * t182
        t187 = i4 - 2
        t190 = vx(i3+1,i4-2) ** 2
        t192 = vy(i3+1,i4-2) ** 2
        t196 = exp(-alpha * (t190 + t192) / 0.2E1)
        t198 = f(i1,i2,t18,t187) / t196
        t200 = (t48 - t198) * t35
        t201 = t50 - t200
        t209 = t9 * t17 * ((-t165 * t180 / 0.6E1 + t183 / 0.6E1) * t35 /
     # 0.2E1 + (t165 * t201 / 0.6E1 - t183 / 0.6E1) * t35 / 0.2E1)
        t210 = vx(i3+2,i4)
        t211 = t210 ** 2
        t212 = vy(i3+2,i4)
        t213 = t212 ** 2
        t217 = exp(-alpha * (t211 + t213) / 0.2E1)
        t218 = nu * t217
        t219 = Hh(i3+2,i4)
        t223 = Kh(i3+2,i4)
        t226 = 0.2E1 * t210 * t212 * t219 - t210 * t212 * t223
        t227 = i3 + 2
        t230 = vx(i3+2,i4+1) ** 2
        t232 = vy(i3+2,i4+1) ** 2
        t236 = exp(-alpha * (t230 + t232) / 0.2E1)
        t238 = f(i1,i2,t227,t19) / t236
        t241 = f(i1,i2,t227,i4) / t217
        t243 = (t238 - t241) * t35
        t246 = vx(i3+2,i4-1) ** 2
        t248 = vy(i3+2,i4-1) ** 2
        t252 = exp(-alpha * (t246 + t248) / 0.2E1)
        t254 = f(i1,i2,t227,t37) / t252
        t256 = (t241 - t254) * t35
        t260 = t218 * t226 * (t243 / 0.2E1 + t256 / 0.2E1)
        t262 = (t260 - t54) * t106
        t264 = (t262 - t107) * t106
        t267 = vx(i3,i4+2)
        t268 = t267 ** 2
        t269 = vy(i3,i4+2)
        t270 = t269 ** 2
        t274 = exp(-alpha * (t268 + t270) / 0.2E1)
        t276 = f(i1,i2,i3,t166) / t274
        t277 = t276 - t82
        t278 = t277 * t35
        t279 = t278 - t87
        t280 = t165 * t279
        t281 = t87 - t100
        t282 = t165 * t281
        t283 = -t280 + t282
        t287 = vx(i3,i4-2)
        t288 = t287 ** 2
        t289 = vy(i3,i4-2)
        t290 = t289 ** 2
        t294 = exp(-alpha * (t288 + t290) / 0.2E1)
        t296 = f(i1,i2,i3,t187) / t294
        t297 = t98 - t296
        t298 = t297 * t35
        t299 = t100 - t298
        t300 = t165 * t299
        t301 = -t282 + t300
        t307 = t63 * t71 * (t283 * t35 / 0.12E2 + t301 * t35 / 0.12E2)
        t309 = (t107 - t161) * t106
        t310 = t309 / 0.6E1
        t315 = vx(i3-1,i4+2) ** 2
        t317 = vy(i3-1,i4+2) ** 2
        t321 = exp(-alpha * (t315 + t317) / 0.2E1)
        t323 = f(i1,i2,t126,t166) / t321
        t325 = (t323 - t137) * t35
        t326 = t325 - t142
        t328 = t142 - t155
        t329 = t165 * t328
        t335 = vx(i3-1,i4-2) ** 2
        t337 = vy(i3-1,i4-2) ** 2
        t341 = exp(-alpha * (t335 + t337) / 0.2E1)
        t343 = f(i1,i2,t126,t187) / t341
        t345 = (t153 - t343) * t35
        t346 = t155 - t345
        t354 = t117 * t125 * ((-t165 * t326 / 0.6E1 + t329 / 0.6E1) * t3
     #5 / 0.2E1 + (t165 * t346 / 0.6E1 - t329 / 0.6E1) * t35 / 0.2E1)
        t355 = vx(i3-2,i4)
        t356 = t355 ** 2
        t357 = vy(i3-2,i4)
        t358 = t357 ** 2
        t362 = exp(-alpha * (t356 + t358) / 0.2E1)
        t363 = nu * t362
        t364 = Hh(i3-2,i4)
        t368 = Kh(i3-2,i4)
        t371 = 0.2E1 * t355 * t357 * t364 - t355 * t357 * t368
        t372 = i3 - 2
        t375 = vx(i3-2,i4+1) ** 2
        t377 = vy(i3-2,i4+1) ** 2
        t381 = exp(-alpha * (t375 + t377) / 0.2E1)
        t383 = f(i1,i2,t372,t19) / t381
        t386 = f(i1,i2,t372,i4) / t362
        t388 = (t383 - t386) * t35
        t391 = vx(i3-2,i4-1) ** 2
        t393 = vy(i3-2,i4-1) ** 2
        t397 = exp(-alpha * (t391 + t393) / 0.2E1)
        t399 = f(i1,i2,t372,t37) / t397
        t401 = (t386 - t399) * t35
        t405 = t363 * t371 * (t388 / 0.2E1 + t401 / 0.2E1)
        t407 = (t159 - t405) * t106
        t409 = (t161 - t407) * t106
        t416 = t163 ** 2
        t417 = dvy ** 2
        t418 = t417 * dvy
        t420 = t418 / t416
        t421 = i4 + 3
        t424 = vx(i3+1,i4+3) ** 2
        t426 = vy(i3+1,i4+3) ** 2
        t430 = exp(-alpha * (t424 + t426) / 0.2E1)
        t432 = f(i1,i2,t18,t421) / t430
        t437 = t180 * t35
        t440 = t182 * t35
        t442 = (t437 - t440) * t35
        t445 = t201 * t35
        t447 = (t440 - t445) * t35
        t449 = t420 * (t442 - t447)
        t453 = i4 - 3
        t456 = vx(i3+1,i4-3) ** 2
        t458 = vy(i3+1,i4-3) ** 2
        t462 = exp(-alpha * (t456 + t458) / 0.2E1)
        t464 = f(i1,i2,t18,t453) / t462
        t482 = vx(i3+2,i4+2) ** 2
        t484 = vy(i3+2,i4+2) ** 2
        t488 = exp(-alpha * (t482 + t484) / 0.2E1)
        t490 = f(i1,i2,t227,t166) / t488
        t496 = t165 * (t243 - t256)
        t502 = vx(i3+2,i4-2) ** 2
        t504 = vy(i3+2,i4-2) ** 2
        t508 = exp(-alpha * (t502 + t504) / 0.2E1)
        t510 = f(i1,i2,t227,t187) / t508
        t525 = (t209 - t307) * t106
        t529 = vx(i3+3,i4)
        t530 = t529 ** 2
        t531 = vy(i3+3,i4)
        t532 = t531 ** 2
        t536 = exp(-alpha * (t530 + t532) / 0.2E1)
        t537 = nu * t536
        t538 = Hh(i3+3,i4)
        t542 = Kh(i3+3,i4)
        t546 = i3 + 3
        t549 = vx(i3+3,i4+1) ** 2
        t551 = vy(i3+3,i4+1) ** 2
        t555 = exp(-alpha * (t549 + t551) / 0.2E1)
        t557 = f(i1,i2,t546,t19) / t555
        t560 = f(i1,i2,t546,i4) / t536
        t565 = vx(i3+3,i4-1) ** 2
        t567 = vy(i3+3,i4-1) ** 2
        t571 = exp(-alpha * (t565 + t567) / 0.2E1)
        t573 = f(i1,i2,t546,t37) / t571
        t587 = (t264 - t309) * t106
        t592 = vx(i3,i4+3)
        t593 = t592 ** 2
        t594 = vy(i3,i4+3)
        t595 = t594 ** 2
        t599 = exp(-alpha * (t593 + t595) / 0.2E1)
        t601 = f(i1,i2,i3,t421) / t599
        t602 = t601 - t276
        t604 = t35 * t602 - t278
        t606 = t279 * t35
        t609 = t281 * t35
        t611 = (t606 - t609) * t35
        t614 = t299 * t35
        t616 = (t609 - t614) * t35
        t618 = t420 * (t611 - t616)
        t619 = t420 * ((t35 * t604 - t606) * t35 - t611) - t618
        t623 = vx(i3,i4-3)
        t624 = t623 ** 2
        t625 = vy(i3,i4-3)
        t626 = t625 ** 2
        t630 = exp(-alpha * (t624 + t626) / 0.2E1)
        t632 = f(i1,i2,i3,t453) / t630
        t633 = t296 - t632
        t635 = -t35 * t633 + t298
        t641 = t618 - t420 * (t616 - (-t35 * t635 + t614) * t35)
        t647 = t63 * t71 * (t619 * t35 / 0.60E2 + t641 * t35 / 0.60E2)
        t649 = (t307 - t354) * t106
        t652 = (t525 - t649) * t106 / 0.6E1
        t654 = (t309 - t409) * t106
        t657 = (t587 - t654) * t106 / 0.30E2
        t662 = vx(i3-1,i4+3) ** 2
        t664 = vy(i3-1,i4+3) ** 2
        t668 = exp(-alpha * (t662 + t664) / 0.2E1)
        t670 = f(i1,i2,t126,t421) / t668
        t675 = t326 * t35
        t678 = t328 * t35
        t680 = (t675 - t678) * t35
        t683 = t346 * t35
        t685 = (t678 - t683) * t35
        t687 = t420 * (t680 - t685)
        t693 = vx(i3-1,i4-3) ** 2
        t695 = vy(i3-1,i4-3) ** 2
        t699 = exp(-alpha * (t693 + t695) / 0.2E1)
        t701 = f(i1,i2,t126,t453) / t699
        t719 = vx(i3-2,i4+2) ** 2
        t721 = vy(i3-2,i4+2) ** 2
        t725 = exp(-alpha * (t719 + t721) / 0.2E1)
        t727 = f(i1,i2,t372,t166) / t725
        t733 = t165 * (t388 - t401)
        t739 = vx(i3-2,i4-2) ** 2
        t741 = vy(i3-2,i4-2) ** 2
        t745 = exp(-alpha * (t739 + t741) / 0.2E1)
        t747 = f(i1,i2,t372,t187) / t745
        t764 = vx(i3-3,i4)
        t765 = t764 ** 2
        t766 = vy(i3-3,i4)
        t767 = t766 ** 2
        t771 = exp(-alpha * (t765 + t767) / 0.2E1)
        t772 = nu * t771
        t773 = Hh(i3-3,i4)
        t777 = Kh(i3-3,i4)
        t781 = i3 - 3
        t784 = vx(i3-3,i4+1) ** 2
        t786 = vy(i3-3,i4+1) ** 2
        t790 = exp(-alpha * (t784 + t786) / 0.2E1)
        t792 = f(i1,i2,t781,t19) / t790
        t795 = f(i1,i2,t781,i4) / t771
        t800 = vx(i3-3,i4-1) ** 2
        t802 = vy(i3-3,i4-1) ** 2
        t806 = exp(-alpha * (t800 + t802) / 0.2E1)
        t808 = f(i1,i2,t781,t37) / t806
        t833 = t9 * (0.2E1 * t10 * t2 + t14 * t4)
        t834 = t833 / 0.2E1
        t839 = t218 * (0.2E1 * t211 * t219 + t213 * t223)
        t841 = (t839 - t833) * t106
        t846 = t63 * (0.2E1 * t56 * t64 + t58 * t68)
        t848 = (t833 - t846) * t106
        t849 = t841 - t848
        t851 = dvx * t849 / 0.16E2
        t852 = t163 * dvx
        t857 = t537 * (0.2E1 * t530 * t538 + t532 * t542)
        t860 = (t857 - t839) * t106 - t841
        t862 = t849 * t106
        t869 = t117 * (0.2E1 * t110 * t118 + t112 * t122)
        t871 = (t846 - t869) * t106
        t872 = t848 - t871
        t873 = t872 * t106
        t875 = (t862 - t873) * t106
        t879 = t846 / 0.2E1
        t881 = dvx * t872 / 0.16E2
        t886 = t363 * (0.2E1 * t356 * t364 + t358 * t368)
        t888 = (t869 - t886) * t106
        t889 = t871 - t888
        t890 = t889 * t106
        t892 = (t873 - t890) * t106
        t895 = 0.3E1 / 0.256E3 * t852 * (t875 - t892)
        t897 = t33 - t85
        t900 = t869 / 0.2E1
        t902 = dvx * t889 / 0.16E2
        t907 = t772 * (0.2E1 * t765 * t773 + t767 * t777)
        t910 = t888 - (t886 - t907) * t106
        t918 = t85 - t140
        t923 = t834 - t851 + t879 - t881
        t924 = t241 - t33
        t925 = t924 * t106
        t926 = t897 * t106
        t928 = (t925 - t926) * t106
        t929 = t918 * t106
        t931 = (t926 - t929) * t106
        t932 = -t928 + t931
        t933 = t932 / 0.24E2
        t943 = t923 * t897 * t106
        t946 = t879 - t881 + t900 - t902
        t948 = t946 * t918 * t106
        t950 = (t943 - t948) * t106
        t954 = t140 - t386
        t955 = t954 * t106
        t957 = (t929 - t955) * t106
        t958 = -t931 + t957
        t959 = t958 / 0.24E2
        t976 = t833 / 0.2E1 + t846 / 0.2E1
        t977 = t560 - t241
        t981 = (t106 * t977 - t925) * t106 - t928
        t984 = -t932 * t106
        t988 = -t958 * t106
        t990 = (t984 - t988) * t106
        t991 = (t106 * t981 - t984) * t106 - t990
        t996 = t839 / 0.2E1 + t833 / 0.2E1
        t1002 = t976 * t933 * t106
        t1006 = t846 / 0.2E1 + t869 / 0.2E1
        t1008 = t1006 * t959 * t106
        t1010 = (t1002 - t1008) * t106
        t1019 = t996 * t924 * t106
        t1023 = t976 * t897 * t106
        t1025 = (t1019 - t1023) * t106
        t1029 = t1006 * t918 * t106
        t1031 = (t1023 - t1029) * t106
        t1033 = (t1025 - t1031) * t106
        t1037 = t869 / 0.2E1 + t886 / 0.2E1
        t1039 = t1037 * t954 * t106
        t1041 = (t1029 - t1039) * t106
        t1043 = (t1031 - t1041) * t106
        t1045 = (t1033 - t1043) * t106
        t1049 = t386 - t795
        t1053 = t957 - (-t1049 * t106 + t955) * t106
        t1057 = t990 - (-t1053 * t106 + t988) * t106
        t1085 = nu * t80
        t1086 = Hh(i3,i4+1)
        t1092 = t1085 * (0.2E1 * t1086 * t76 + t74 * Kh(i3,i4+1))
        t1093 = t1092 / 0.2E1
        t1094 = nu * t274
        t1095 = Hh(i3,i4+2)
        t1101 = t1094 * (0.2E1 * t1095 * t270 + t268 * Kh(i3,i4+2))
        t1103 = (t1101 - t1092) * t35
        t1109 = t63 * (t56 * Kh(i3,i4) + 0.2E1 * t58 * t64)
        t1111 = (t1092 - t1109) * t35
        t1112 = t1103 - t1111
        t1114 = dvy * t1112 / 0.16E2
        t1115 = nu * t599
        t1116 = Hh(i3,i4+3)
        t1122 = t1115 * (0.2E1 * t1116 * t595 + t593 * Kh(i3,i4+3))
        t1125 = (t1122 - t1101) * t35 - t1103
        t1127 = t1112 * t35
        t1130 = nu * t96
        t1131 = Hh(i3,i4-1)
        t1137 = t1130 * (0.2E1 * t1131 * t92 + t90 * Kh(i3,i4-1))
        t1139 = (t1109 - t1137) * t35
        t1140 = t1111 - t1139
        t1141 = t1140 * t35
        t1143 = (t1127 - t1141) * t35
        t1147 = t1109 / 0.2E1
        t1149 = dvy * t1140 / 0.16E2
        t1150 = nu * t294
        t1151 = Hh(i3,i4-2)
        t1157 = t1150 * (0.2E1 * t1151 * t290 + t288 * Kh(i3,i4-2))
        t1159 = (t1137 - t1157) * t35
        t1160 = t1139 - t1159
        t1161 = t1160 * t35
        t1163 = (t1141 - t1161) * t35
        t1166 = 0.3E1 / 0.256E3 * t418 * (t1143 - t1163)
        t1170 = t1137 / 0.2E1
        t1172 = dvy * t1160 / 0.16E2
        t1173 = nu * t630
        t1174 = Hh(i3,i4-3)
        t1180 = t1173 * (0.2E1 * t1174 * t626 + t624 * Kh(i3,i4-3))
        t1183 = t1159 - (t1157 - t1180) * t35
        t1195 = t1093 - t1114 + t1147 - t1149
        t1196 = t283 / 0.24E2
        t1206 = t1195 * t86 * t35
        t1209 = t1147 - t1149 + t1170 - t1172
        t1211 = t1209 * t99 * t35
        t1213 = (t1206 - t1211) * t35
        t1217 = t301 / 0.24E2
        t1235 = t1092 / 0.2E1 + t1109 / 0.2E1
        t1240 = t1101 / 0.2E1 + t1092 / 0.2E1
        t1247 = t1235 * t1196 * t35
        t1251 = t1109 / 0.2E1 + t1137 / 0.2E1
        t1253 = t1251 * t1217 * t35
        t1255 = (t1247 - t1253) * t35
        t1264 = t1240 * t277 * t35
        t1268 = t1235 * t86 * t35
        t1270 = (t1264 - t1268) * t35
        t1274 = t1251 * t99 * t35
        t1276 = (t1268 - t1274) * t35
        t1278 = (t1270 - t1276) * t35
        t1282 = t1137 / 0.2E1 + t1157 / 0.2E1
        t1284 = t1282 * t297 * t35
        t1286 = (t1274 - t1284) * t35
        t1288 = (t1276 - t1286) * t35
        t1290 = (t1278 - t1288) * t35
        t1329 = 0.2E1 * t1086 * t73 * t75 - t73 * t75 * Kh(i3,i4+1)
        t1331 = (t30 - t82) * t106
        t1333 = (t82 - t137) * t106
        t1337 = t1085 * t1329 * (t1331 / 0.2E1 + t1333 / 0.2E1)
        t1341 = t63 * t71 * (t926 / 0.2E1 + t929 / 0.2E1)
        t1343 = (t1337 - t1341) * t35
        t1351 = 0.2E1 * t1131 * t89 * t91 - t89 * t91 * Kh(i3,i4-1)
        t1353 = (t48 - t98) * t106
        t1355 = (t98 - t153) * t106
        t1359 = t1130 * t1351 * (t1353 / 0.2E1 + t1355 / 0.2E1)
        t1361 = (t1341 - t1359) * t35
        t1364 = (t238 - t30) * t106
        t1366 = (t1364 - t1331) * t106
        t1368 = (t1331 - t1333) * t106
        t1369 = -t1366 + t1368
        t1373 = (t137 - t383) * t106
        t1375 = (t1333 - t1373) * t106
        t1376 = -t1368 + t1375
        t1382 = t1085 * t1329 * (t1369 * t106 / 0.12E2 + t1376 * t106 / 
     #0.12E2)
        t1389 = 0.2E1 * t1095 * t267 * t269 - t267 * t269 * Kh(i3,i4+2)
        t1391 = (t177 - t276) * t106
        t1393 = (t276 - t323) * t106
        t1397 = t1094 * t1389 * (t1391 / 0.2E1 + t1393 / 0.2E1)
        t1399 = (t1397 - t1337) * t35
        t1400 = t1399 - t1343
        t1410 = t63 * t71 * (t932 * t106 / 0.12E2 + t958 * t106 / 0.12E2
     #)
        t1411 = t1343 - t1361
        t1413 = t165 * t1411 / 0.6E1
        t1417 = (t254 - t48) * t106
        t1419 = (t1417 - t1353) * t106
        t1421 = (t1353 - t1355) * t106
        t1422 = -t1419 + t1421
        t1426 = (t153 - t399) * t106
        t1428 = (t1355 - t1426) * t106
        t1429 = -t1421 + t1428
        t1435 = t1130 * t1351 * (t1422 * t106 / 0.12E2 + t1429 * t106 / 
     #0.12E2)
        t1442 = 0.2E1 * t1151 * t287 * t289 - t287 * t289 * Kh(i3,i4-2)
        t1444 = (t198 - t296) * t106
        t1446 = (t296 - t343) * t106
        t1450 = t1150 * t1442 * (t1444 / 0.2E1 + t1446 / 0.2E1)
        t1452 = (t1359 - t1450) * t35
        t1453 = t1361 - t1452
        t1468 = -t1369 * t106
        t1472 = -t1376 * t106
        t1474 = (t1468 - t1472) * t106
        t1498 = (t1391 - t1393) * t106
        t1516 = (t1382 - t1410) * t35
        t1539 = t1400 * t35
        t1542 = t1411 * t35
        t1544 = (t1539 - t1542) * t35
        t1555 = t63 * t71 * (t991 * t106 / 0.60E2 + t1057 * t106 / 0.60E
     #2)
        t1557 = (t1410 - t1435) * t35
        t1560 = t165 * (t1516 - t1557) / 0.6E1
        t1561 = t1453 * t35
        t1563 = (t1542 - t1561) * t35
        t1566 = t420 * (t1544 - t1563) / 0.30E2
        t1576 = -t1422 * t106
        t1580 = -t1429 * t106
        t1582 = (t1576 - t1580) * t106
        t1606 = (t1444 - t1446) * t106
        temp = 
     #t107 / 0.2E1 + t161 / 0.2E1 + t163 * ((t209 - t264 / 0.6E1
     # - t307 + t310) * t106 / 0.2E1 + (t307 - t310 - t354 + t409 / 0.6E
     #1) * t106 / 0.2E1) + t416 * ((t9 * t17 * ((t420 * ((((t432 - t177)
     # * t35 - t179) * t35 - t437) * t35 - t442) / 0.30E2 - t449 / 0.30E
     #2) * t35 / 0.2E1 + (t449 / 0.30E2 - t420 * (t447 - (t445 - (t200 -
     # (t198 - t464) * t35) * t35) * t35) / 0.30E2) * t35 / 0.2E1) - ((t
     #218 * t226 * ((-t165 * ((t490 - t238) * t35 - t243) / 0.6E1 + t496
     # / 0.6E1) * t35 / 0.2E1 + (-t496 / 0.6E1 + t165 * (t256 - (t254 - 
     #t510) * t35) / 0.6E1) * t35 / 0.2E1) - t209) * t106 - t525) * t106
     # / 0.6E1 + ((((t537 * (0.2E1 * t529 * t531 * t538 - t529 * t531 * 
     #t542) * ((t557 - t560) * t35 / 0.2E1 + (t560 - t573) * t35 / 0.2E1
     #) - t260) * t106 - t262) * t106 - t264) * t106 - t587) * t106 / 0.
     #30E2 - t647 + t652 - t657) * t106 / 0.2E1 + (t647 - t652 + t657 - 
     #t117 * t125 * ((t420 * ((((t670 - t323) * t35 - t325) * t35 - t675
     #) * t35 - t680) / 0.30E2 - t687 / 0.30E2) * t35 / 0.2E1 + (t687 / 
     #0.30E2 - t420 * (t685 - (t683 - (t345 - (t343 - t701) * t35) * t35
     #) * t35) / 0.30E2) * t35 / 0.2E1) + (t649 - (t354 - t363 * t371 * 
     #((-t165 * ((t727 - t383) * t35 - t388) / 0.6E1 + t733 / 0.6E1) * t
     #35 / 0.2E1 + (-t733 / 0.6E1 + t165 * (t401 - (t399 - t747) * t35) 
     #/ 0.6E1) * t35 / 0.2E1)) * t106) * t106 / 0.6E1 - (t654 - (t409 - 
     #(t407 - (t405 - t772 * (0.2E1 * t764 * t766 * t773 - t764 * t766 *
     # t777) * ((t792 - t795) * t35 / 0.2E1 + (t795 - t808) * t35 / 0.2E
     #1)) * t106) * t106) * t106) * t106 / 0.30E2) * t106 / 0.2E1) + ((t
     #834 - t851 + 0.3E1 / 0.256E3 * t852 * ((t106 * t860 - t862) * t106
     # - t875) + t879 - t881 + t895) * t897 * t106 - (t879 - t881 + t895
     # + t900 - t902 + 0.3E1 / 0.256E3 * t852 * (t892 - (-t106 * t910 + 
     #t890) * t106)) * t918 * t106) * t106 + dvx * (t923 * t933 * t106 -
     # (((t839 / 0.2E1 - dvx * t860 / 0.16E2 + t834 - t851) * t924 * t10
     #6 - t943) * t106 - t950) * t106 / 0.24E2 - t946 * t959 * t106 + (t
     #950 - (t948 - (t900 - t902 + t886 / 0.2E1 - dvx * t910 / 0.16E2) *
     # t954 * t106) * t106) * t106 / 0.24E2) + t852 * (0.3E1 / 0.640E3 *
     # t976 * t991 * t106 - ((-t996 * t981 * t106 / 0.24E2 - t1002) * t1
     #06 - t1010) * t106 / 0.24E2 + 0.3E1 / 0.640E3 * (((((t857 / 0.2E1 
     #+ t839 / 0.2E1) * t977 * t106 - t1019) * t106 - t1025) * t106 - t1
     #033) * t106 - t1045) * t106 - 0.3E1 / 0.640E3 * t1006 * t1057 * t1
     #06 + (t1010 - (t1008 + t1037 * t1053 * t106 / 0.24E2) * t106) * t1
     #06 / 0.24E2 - 0.3E1 / 0.640E3 * (t1045 - (t1043 - (t1041 - (t1039 
     #- (t886 / 0.2E1 + t907 / 0.2E1) * t1049 * t106) * t106) * t106) * 
     #t106) * t106) + ((t1093 - t1114 + 0.3E1 / 0.256E3 * t418 * ((t1125
     # * t35 - t1127) * t35 - t1143) + t1147 - t1149 + t1166) * t86 * t3
     #5 - (t1147 - t1149 + t1166 + t1170 - t1172 + 0.3E1 / 0.256E3 * t41
     #8 * (t1163 - (-t1183 * t35 + t1161) * t35)) * t99 * t35) * t35 + t
     #163 * (t1195 * t1196 * t35 - t165 * (((t1101 / 0.2E1 - dvy * t1125
     # / 0.16E2 + t1093 - t1114) * t277 * t35 - t1206) * t35 - t1213) / 
     #0.24E2 - t1209 * t1217 * t35 + t165 * (t1213 - (t1211 - (t1170 - t
     #1172 + t1157 / 0.2E1 - dvy * t1183 / 0.16E2) * t297 * t35) * t35) 
     #/ 0.24E2) * t35 + t416 * (0.3E1 / 0.640E3 * t1235 * t619 * t35 - t
     #165 * ((t1240 * (-t165 * t604 / 0.24E2 + t280 / 0.24E2) * t35 - t1
     #247) * t35 - t1255) / 0.24E2 + 0.3E1 / 0.640E3 * t420 * (((((t1122
     # / 0.2E1 + t1101 / 0.2E1) * t602 * t35 - t1264) * t35 - t1270) * t
     #35 - t1278) * t35 - t1290) - 0.3E1 / 0.640E3 * t1251 * t641 * t35 
     #+ t165 * (t1255 - (t1253 - t1282 * (t165 * t635 / 0.24E2 - t300 / 
     #0.24E2) * t35) * t35) / 0.24E2 - 0.3E1 / 0.640E3 * t420 * (t1290 -
     # (t1288 - (t1286 - (t1284 - (t1157 / 0.2E1 + t1180 / 0.2E1) * t633
     # * t35) * t35) * t35) * t35)) * t35 + t1343 / 0.2E1 + t1361 / 0.2E
     #1 + t163 * ((t1382 - t165 * t1400 / 0.6E1 - t1410 + t1413) * t35 /
     # 0.2E1 + (t1410 - t1413 - t1435 + t165 * t1453 / 0.6E1) * t35 / 0.
     #2E1) + t416 * ((t1085 * t1329 * ((((((t557 - t238) * t106 - t1364)
     # * t106 - t1366) * t106 - t1468) * t106 / 0.30E2 - t1474 / 0.30E2)
     # * t106 / 0.2E1 + (t1474 / 0.30E2 - (t1472 - (t1375 - (t1373 - (t3
     #83 - t792) * t106) * t106) * t106) * t106 / 0.30E2) * t106 / 0.2E1
     #) - t165 * ((t1094 * t1389 * ((-((t490 - t177) * t106 - t1391) * t
     #106 / 0.6E1 + t1498 / 0.6E1) * t106 / 0.2E1 + (-t1498 / 0.6E1 + (t
     #1393 - (t323 - t727) * t106) * t106 / 0.6E1) * t106 / 0.2E1) - t13
     #82) * t35 - t1516) / 0.6E1 + t420 * ((((t1115 * (0.2E1 * t1116 * t
     #592 * t594 - t592 * t594 * Kh(i3,i4+3)) * ((t432 - t601) * t106 / 
     #0.2E1 + (t601 - t670) * t106 / 0.2E1) - t1397) * t35 - t1399) * t3
     #5 - t1539) * t35 - t1544) / 0.30E2 - t1555 + t1560 - t1566) * t35 
     #/ 0.2E1 + (t1555 - t1560 + t1566 - t1130 * t1351 * ((((((t573 - t2
     #54) * t106 - t1417) * t106 - t1419) * t106 - t1576) * t106 / 0.30E
     #2 - t1582 / 0.30E2) * t106 / 0.2E1 + (t1582 / 0.30E2 - (t1580 - (t
     #1428 - (t1426 - (t399 - t808) * t106) * t106) * t106) * t106 / 0.3
     #0E2) * t106 / 0.2E1) + t165 * (t1557 - (t1435 - t1150 * t1442 * ((
     #-((t510 - t198) * t106 - t1444) * t106 / 0.6E1 + t1606 / 0.6E1) * 
     #t106 / 0.2E1 + (-t1606 / 0.6E1 + (t1446 - (t343 - t747) * t106) * 
     #t106 / 0.6E1) * t106 / 0.2E1)) * t35) / 0.6E1 - t420 * (t1563 - (t
     #1561 - (t1452 - (t1450 - t1173 * (0.2E1 * t1174 * t623 * t625 - t6
     #23 * t625 * Kh(i3,i4-3)) * ((t464 - t632) * t106 / 0.2E1 + (t632 -
     # t701) * t106 / 0.2E1)) * t35) * t35) * t35) / 0.30E2) * t35 / 0.2
     #E1)
ccccc
ccccc
        rhs(i1,i2,i3,i4) = rhs(i1,i2,i3,i4)+temp
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
      subroutine getRamp(A, v, vceil, vmax, pi)
c
c.. evaluate ramp function
      implicit none
c
c.. declarations of incoming variables 
      real A, v, vceil, vmax, pi
c
c.. declarations of local variables

c
      if (v .le. vceil) then
        A = 1.0
      elseif (v .gt. vmax) then
        A = 0.0
      else
        A = 1.0-(sin(pi/2.0*(v-vceil)/(vmax-vceil)))**2
      end if
c
      return
      end
