c
c Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
c See the top-level LICENSE file for details.
c Produced at the Lawrence Livermore National Laboratory
c
c SPDX-License-Identifier: MIT
c
c Fortran functions called by CollisionOperator.
c
      subroutine computeDistFuncAvg(
     &     nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b,
     &     n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b,
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
     &     vflowx, vflowy,
     &     velocities,
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
      real vflowx, vflowy
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real u(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real vth
      integer comm
c
c.. declarations of local variables
      real vx, vy
      real uxybar, T, N, N_local, KE, KE_local
      integer NxNy, i1, i2, i3, i4, ierr
c
      N_local = 0.0
      KE_local = 0.0
      NxNy = (n1b-n1a+1)*(n2b-n2a+1)
      do i4 = n4a, n4b
        do i3 = n3a, n3b
          vx = velocities(i3, i4, 0) - vflowx
          vy = velocities(i3, i4, 1) - vflowy
          uxybar = 0.0
          do i2 = n2a, n2b
            do i1 = n1a, n1b
              uxybar = uxybar + u(i1, i2, i3, i4)
            end do
          end do
          uxybar = uxybar / NxNy
          N_local = N_local + uxybar
          KE_local = KE_local + (vx * vx + vy * vy) * uxybar
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
     &     vel_space_size,
     &     vflowx, vflowy,
     &     velocities,
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
      real vflowx, vflowy
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real uxybar(1:vel_space_size)
      real vth
      integer is_config_space_head_node, comm
c
c.. declarations of local variables
      real vx, vy
      real T, N, N_local, KE, KE_local
      integer idx, i1, i2, i3, i4, ierr
c
      N_local = 0.0
      KE_local = 0.0
      if (is_config_space_head_node .eq. 1) then
        idx = 1
        do i4 = n4a, n4b
          do i3 = n3a, n3b
            vx = velocities(i3, i4, 0) - vflowx
            vy = velocities(i3, i4, 1) - vflowy
            N_local = N_local + uxybar(idx)
            KE_local = KE_local + (vx * vx + vy * vy) * uxybar(idx)
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
