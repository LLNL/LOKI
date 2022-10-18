c
c Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
c See the top-level LICENSE file for details.
c Produced at the Lawrence Livermore National Laboratory
c
c SPDX-License-Identifier: MIT
c
c Fortran function to convert cell centered average values to face centered
c average values at 4th order.
c
      subroutine centAvg_to_faceAvg4D( nd1a,nd1b,nd2a,nd2b,
     *                                 nd3a,nd3b,nd4a,nd4b,
     *                                 nf1a,nf1b,nf2a,nf2b,
     *                                 nf3a,nf3b,nf4a,nf4b,
     *                                 dir,cell,face )
c
c.. function to convert cell centered average values to face centered average
c     values at 4th order. We assume a uniform grid. 
      implicit none
c
c.. declarations of incoming variables
      integer nd1a,nd1b,nd2a,nd2b
      integer nd3a,nd3b,nd4a,nd4b
      integer nf1a,nf1b,nf2a,nf2b
      integer nf3a,nf3b,nf4a,nf4b
      integer dir
      real cell( nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b )
      real face( nf1a:nf1b,nf2a:nf2b,nf3a:nf3b,nf4a:nf4b )
c
c.. declarations of local variables
      integer f1a,f1b
      integer i1,i2,i3,i4
c
      f1a = nf1a+2
      f1b = nf1b-2
c
      if( dir.eq.1 ) then
        do i4 = nf4a,nf4b    
        do i3 = nf3a,nf3b    
        do i2 = nf2a,nf2b    
        do i1 = f1a,f1b  
          face(i1,i2,i3,i4) = (7.0*(cell(i1-1,i2,i3,i4)+
     *                              cell(i1,i2,i3,i4))-
     *       (cell(i1-2,i2,i3,i4)+cell(i1+1,i2,i3,i4)))/12.0
        end do
        end do
        end do
        end do
      else if( dir.eq.2 ) then
        do i4 = nf4a,nf4b    
        do i3 = nf3a,nf3b    
        do i2 = nf2a,nf2b    
        do i1 = f1a,f1b  
c        do i3 = nf3a,nf3b    
c        do i2 = nf2a,nf2b    
c        do i1 = f1a,f1b  
c        do i4 = nf4a,nf4b    
          face(i1,i2,i3,i4) = (7.0*(cell(i4,i1-1,i2,i3)+
     *                              cell(i4,i1,i2,i3))-
     *       (cell(i4,i1-2,i2,i3)+cell(i4,i1+1,i2,i3)))/12.0
        end do
        end do
        end do
        end do
      else if( dir.eq.3 ) then
        do i4 = nf4a,nf4b    
        do i3 = nf3a,nf3b    
        do i2 = nf2a,nf2b    
        do i1 = f1a,f1b  
c        do i2 = nf2a,nf2b    
c        do i1 = f1a,f1b  
c        do i4 = nf4a,nf4b    
c        do i3 = nf3a,nf3b    
          face(i1,i2,i3,i4) = (7.0*(cell(i3,i4,i1-1,i2)+
     *                              cell(i3,i4,i1,i2))-
     *       (cell(i3,i4,i1-2,i2)+cell(i3,i4,i1+1,i2)))/12.0
        end do
        end do
        end do
        end do
      else
        do i4 = nf4a,nf4b    
        do i3 = nf3a,nf3b    
        do i2 = nf2a,nf2b    
        do i1 = f1a,f1b  
c        do i4 = nf4a,nf4b    
c        do i3 = nf3a,nf3b    
c        do i2 = nf2a,nf2b    
c        do i1 = f1a,f1b  
          face(i1,i2,i3,i4) = (7.0*(cell(i2,i3,i4,i1-1)+
     *                              cell(i2,i3,i4,i1))-
     *       (cell(i2,i3,i4,i1-2)+cell(i2,i3,i4,i1+1)))/12.0
        end do
        end do
        end do
        end do
      end if
c
      return
      end
