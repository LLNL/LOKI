c
c Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
c See the top-level LICENSE file for details.
c Produced at the Lawrence Livermore National Laboratory
c
c SPDX-License-Identifier: MIT
c
c Fortran functions called by Poisson.
c
      subroutine neutralizeCharge4D( 
     *     md1a,md1b,md2a,md2b,
     *     n1a,n1b,n2a,n2b,
     *     rhs,
     *     comm )
c
c.. function to neutralize charge densities
      implicit none
      include 'mpif.h'
c
      integer md1a,md1b,md2a,md2b
      integer n1a,n1b,n2a,n2b
      integer comm
      real rhs( md1a:md1b,md2a:md2b )
c
      integer i1,i2
      integer count_local,count,ierr
      double precision sum_local,sum
c
c.. Currently, "parallel" solves of Poisson's equation are achieved by
c.. communicating each Poisson processor's part of the charge density
c.. to the other Poisson processors so that each may do the full serial
c.. solve.  Hence the charge neutralization is always done on the
c.. serialized charge density and this routine does not need to reduce
c.. the integrated charge or count.  When a true parallel Poisson solve
c.. is in place the commented out code to acheieve the reductions will
c.. need to be activated.
c      count_local = 0
c      sum_local = 0.0
      count = 0
      sum = 0.0
c
      do i2 = n2a,n2b
      do i1 = n1a,n1b
c        sum_local = sum_local + rhs(i1,i2)
c        count_local = count_local + 1
        sum = sum + rhs(i1,i2)
        count = count + 1
      end do
      end do
c
c      call MPI_ALLREDUCE( sum_local, sum, 1, MPI_DOUBLE_PRECISION,
c     *                    MPI_SUM, comm, ierr )
c      call MPI_ALLREDUCE( count_local, count, 1, MPI_INTEGER,
c     *                    MPI_SUM, comm, ierr )
      sum = sum / count

      do i2 = n2a,n2b
      do i1 = n1a,n1b
        rhs(i1,i2) = rhs(i1,i2) - sum
      end do
      end do
c
      return
      end
c
c+++++++++++
c
      subroutine computeEFieldFromPotential( 
     *     nd1a,nd1b,nd2a,nd2b,
     *     n1a,n1b,n2a,n2b,
     *     solution_order,
     *     em_vars_dim,
     *     dx,emVars,phi)
c
c.. function to neutralize charge densities
      implicit none
c
      integer nd1a,nd1b,nd2a,nd2b
      integer n1a,n1b,n2a,n2b
      integer solution_order, em_vars_dim
      real dx( * )
      real emVars(  nd1a:nd1b,nd2a:nd2b,1:em_vars_dim )
      real phi( nd1a:nd1b,nd2a:nd2b )
c
      integer i1,i2
c
      if( solution_order .eq. 4 ) then
        ! 4th order
        do i2 = n2a,n2b
        do i1 = n1a,n1b
          emVars(i1,i2,1) = (phi(i1-2,i2)-8.0*phi(i1-1,i2)+
     *         8.0*phi(i1+1,i2)-phi(i1+2,i2))/(12.0*dx(1))
          emVars(i1,i2,2) = (phi(i1,i2-2)-8.0*phi(i1,i2-1)+
     *         8.0*phi(i1,i2+1)-phi(i1,i2+2))/(12.0*dx(2))

        end do
        end do
      else
        ! 6th order
        do i2 = n2a,n2b
        do i1 = n1a,n1b
          emVars(i1,i2,1) = (
     *         -1.0  *phi(i1-3,i2)
     *         +9.0  *phi(i1-2,i2)
     *         -45.0 *phi(i1-1,i2)
     *         +45.0 *phi(i1+1,i2)
     *         -9.0  *phi(i1+2,i2)
     *         +1.0  *phi(i1+3,i2))/(60.0*dx(1))
          emVars(i1,i2,2) = (
     *         -1.0  *phi(i1,i2-3)
     *         +9.0  *phi(i1,i2-2)
     *         -45.0 *phi(i1,i2-1)
     *         +45.0 *phi(i1,i2+1)
     *         -9.0  *phi(i1,i2+2)
     *         +1.0  *phi(i1,i2+3))/(60.0*dx(2))
        end do
        end do
      end if
c
      return
      end
c
c+++++++++++
c
