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
      count_local = 0
      sum_local = 0.0
c
      do i2 = n2a,n2b
      do i1 = n1a,n1b
        sum_local = sum_local + rhs(i1,i2)
        count_local = count_local + 1
      end do
      end do
c
      call MPI_ALLREDUCE( sum_local, sum, 1, MPI_DOUBLE_PRECISION,
     *                    MPI_SUM, comm, ierr )
      call MPI_ALLREDUCE( count_local, count, 1, MPI_INTEGER,
     *                    MPI_SUM, comm, ierr )
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
     *     dx,
     *     Ex,Ey,phi,solution_order)
c
c.. function to neutralize charge densities
      implicit none
c
      integer nd1a,nd1b,nd2a,nd2b
      integer n1a,n1b,n2a,n2b
      integer solution_order
      real dx( * )
      real Ex(  nd1a:nd1b,nd2a:nd2b )
      real Ey(  nd1a:nd1b,nd2a:nd2b )
      real phi( nd1a:nd1b,nd2a:nd2b )
c
      integer i1,i2
c
      if( solution_order .eq. 4 ) then
        ! 4th order
        do i2 = n2a,n2b
        do i1 = n1a,n1b
          Ex(i1,i2) = (phi(i1-2,i2)-8.0*phi(i1-1,i2)+
     *         8.0*phi(i1+1,i2)-phi(i1+2,i2))/(12.0*dx(1))
          Ey(i1,i2) = (phi(i1,i2-2)-8.0*phi(i1,i2-1)+
     *         8.0*phi(i1,i2+1)-phi(i1,i2+2))/(12.0*dx(2))

        end do
        end do
      else
        ! 6th order
        do i2 = n2a,n2b
        do i1 = n1a,n1b
          Ex(i1,i2) = (
     *         -1.0  *phi(i1-3,i2)
     *         +9.0  *phi(i1-2,i2)
     *         -45.0 *phi(i1-1,i2)
     *         +45.0 *phi(i1+1,i2)
     *         -9.0  *phi(i1+2,i2)
     *         +1.0  *phi(i1+3,i2))/(60.0*dx(1))
          Ey(i1,i2) = (
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
