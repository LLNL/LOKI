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
c Fortran functions called by CollisionOperator.
c
      subroutine computeDistFuncAvg(
     &     nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b,
     &     n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b,
     &     xlo, xhi, deltax,
     &     vel_space_size,
     &     config_space_size,
     &     u,
     &     uxybar)
c
c.. function to compute average of distribution function over
c.. configuration space
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
      integer vel_space_size, config_space_size
      real xlo(1:4), xhi(1:4), deltax(1:4)
      real u(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real uxybar(1:vel_space_size)
c
c.. declarations of local variables
      integer idx, i1, i2, i3, i4
c
      idx = 1
      do i4 = n4a, n4b
        do i3 = n3a, n3b
          uxybar(idx) = 0.0
          do i2 = n2a, n2b
            do i1 = n1a, n1b
              uxybar(idx) = uxybar(idx) + u(i1, i2, i3, i4)
            end do
          end do
          uxybar(idx) = uxybar(idx) / config_space_size
          idx = idx + 1
        end do
      end do
      return
      end
c
c **************
c
      subroutine computeVThLocal(
     &     nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b,
     &     n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b,
     &     xlo, xhi, deltax,
     &     vflowx, vflowy,
     &     comm,
     &     u, vth)
c
c.. function to compute local kinetic species vthermal
      implicit none
      include 'mpif.h'
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
      real xlo(1:4), xhi(1:4), deltax(1:4)
      real vflowx, vflowy
      real u(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real vth
      integer comm
c
c.. declarations of local variables
      real vx0, dvx, vy0, dvy, vx, vy, vy2
      real uxybar, T, N, N_local, KE, KE_local
      integer NxNy, i1, i2, i3, i4, ierr
c
      N_local = 0.0
      KE_local = 0.0
      vx0 = xlo(3)
      dvx = deltax(3)
      vy0 = xlo(4)
      dvy = deltax(4)
      NxNy = (n1b-n1a+1)*(n2b-n2a+1)
      do i4 = n4a, n4b
        vy = vy0 + (0.5 + i4) * dvy - vflowy
        vy2 = vy * vy
        do i3 = n3a, n3b
          vx = vx0 + (0.5 + i3) * dvx - vflowx
          uxybar = 0.0
          do i2 = n2a, n2b
            do i1 = n1a, n1b
              uxybar = uxybar + u(i1, i2, i3, i4)
            end do
          end do
          uxybar = uxybar / NxNy
          N_local = N_local + uxybar
          KE_local = KE_local + (vx * vx + vy2) * uxybar
        end do
      end do
      call MPI_ALLREDUCE(KE_local, KE, 1, MPI_DOUBLE_PRECISION,
     *                   MPI_SUM, comm, ierr)
      call MPI_ALLREDUCE(N_local, N, 1, MPI_DOUBLE_PRECISION,
     *                   MPI_SUM, comm, ierr)
      T = 0.5 * KE / N
      vth = sqrt(T)
      return
      end
c
c **************
c
      subroutine computeVThGlobal(
     &     nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b,
     &     n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b,
     &     xlo, xhi, deltax,
     &     vel_space_size,
     &     vflowx, vflowy,
     &     uxybar, is_config_space_head_node, comm, vth)
c
c.. function to compute local kinetic species vthermal
      implicit none
      include 'mpif.h'
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
      integer vel_space_size
      real xlo(1:4), xhi(1:4), deltax(1:4)
      real vflowx, vflowy
      real uxybar(1:vel_space_size)
      real vth
      integer is_config_space_head_node, comm
c
c.. declarations of local variables
      real vx0, dvx, vy0, dvy, vx, vy, vy2
      real T, N, N_local, KE, KE_local
      integer idx, i1, i2, i3, i4, ierr
c
      N_local = 0.0
      KE_local = 0.0
      if (is_config_space_head_node .eq. 1) then
        vx0 = xlo(3)
        dvx = deltax(3)
        vy0 = xlo(4)
        dvy = deltax(4)
        idx = 1
        do i4 = n4a, n4b
          vy = vy0 + (0.5 + i4) * dvy - vflowy
          vy2 = vy * vy
          do i3 = n3a, n3b
            vx = vx0 + (0.5 + i3) * dvx - vflowx
            N_local = N_local + uxybar(idx)
            KE_local = KE_local + (vx * vx + vy2) * uxybar(idx)
            idx = idx + 1
          end do
        end do
      endif
      call MPI_ALLREDUCE(KE_local, KE, 1, MPI_DOUBLE_PRECISION,
     *                   MPI_SUM, comm, ierr)
      call MPI_ALLREDUCE(N_local, N, 1, MPI_DOUBLE_PRECISION,
     *                   MPI_SUM, comm, ierr)
      T = 0.5 * KE / N
      vth = sqrt(T)
      return
      end
