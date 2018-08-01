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
