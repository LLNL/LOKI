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
c Fortran functions called by KineticSpecies.
c
      subroutine xpby4d( 
     &     x,
     &     nx1lo,nx1hi,nx2lo,nx2hi,nx3lo,nx3hi,nx4lo,nx4hi,
     &     y,
     &     ny1lo,ny1hi,ny2lo,ny2hi,ny3lo,ny3hi,ny4lo,ny4hi,
     &     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     &     b )
c
      implicit none
c
      integer nx1lo,nx1hi,nx2lo,nx2hi,nx3lo,nx3hi,nx4lo,nx4hi
      integer ny1lo,ny1hi,ny2lo,ny2hi,ny3lo,ny3hi,ny4lo,ny4hi
      integer n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b
      real x( nx1lo:nx1hi,nx2lo:nx2hi,nx3lo:nx3hi,nx4lo:nx4hi )
      real y( ny1lo:ny1hi,ny2lo:ny2hi,ny3lo:ny3hi,ny4lo:ny4hi )
      real b
c
      integer i1,i2,i3,i4
c
      do i4 = n4a,n4b
       do i3 = n3a,n3b
        do i2 = n2a,n2b
         do i1 = n1a,n1b
          x(i1,i2,i3,i4) = x(i1,i2,i3,i4) + b * y(i1,i2,i3,i4)
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
      subroutine setphasespacevel4D(
     &     vel3,vel4,
     &     nv1a,nv1b,nv2a,nv2b,nv3a,nv3b,nv4a,nv4b,
     &     accel,
     &     na1a,na1b,na2a,na2b,na3a,na3b,na4a,na4b,
     &     axmax,aymax )
c
c.. declarations of incoming variables
c
      implicit none
      integer nv1a,nv1b,nv2a,nv2b
      integer nv3a,nv3b,nv4a,nv4b
      integer na1a,na1b,na2a,na2b
      integer na3a,na3b,na4a,na4b
      real accel( na1a:na1b,na2a:na2b,0:1 )
      real vel3( nv3a:nv3b+1,nv4a:nv4b,nv1a:nv1b,nv2a:nv2b )
      real vel4( nv4a:nv4b+1,nv1a:nv1b,nv2a:nv2b,nv3a:nv3b )
      real axmax,aymax
c
c.. declarations of local variables
c
      integer i1,i2,i3,i4

      axmax = 0.0
      do i2 = na2a,na2b
      do i1 = na1a,na1b
      do i4 = nv4a,nv4b
      do i3 = nv3a,nv3b+1
        vel3(i3,i4,i1,i2) = accel(i1,i2,0)
        axmax = max(axmax,abs(accel(i1,i2,0)))
      end do
      end do
      end do
      end do

      aymax = 0.0
      do i3 = nv3a,nv3b
      do i2 = na2a,na2b
      do i1 = na1a,na1b
      do i4 = nv4a,nv4b+1
        vel4(i4,i1,i2,i3) = accel(i1,i2,1)
        aymax = max(aymax,abs(accel(i1,i2,1)))
      end do
      end do
      end do
      end do
      
      return
      end
c
c **************
c
      subroutine setphasespacevelmaxwell4D(
     &     vel3, vel4,
     &     nv1a, nv1b, nv2a, nv2b, nv3a, nv3b, nv4a, nv4b,
     &     xlo, xhi, dx,
     &     charge_per_mass,
     &     em_vars,
     &     vz,
     &     axmax, aymax)
c
c.. function to compute lorentz force E + v x B and propagate
c   velocities to faces
      implicit none
c
c.. declarations of incoming variables
      integer nv1a, nv1b, nv2a, nv2b, nv3a, nv3b, nv4a, nv4b
      real xlo(1:4), xhi(1:4), dx(1:4)
      real charge_per_mass
      real vel3(nv3a:nv3b+1, nv4a:nv4b, nv1a:nv1b, nv2a:nv2b)
      real vel4(nv4a:nv4b+1, nv1a:nv1b, nv2a:nv2b, nv3a:nv3b)
      real em_vars(nv1a:nv1b, nv2a:nv2b, 1:6)
      real vz(nv1a:nv1b, nv2a:nv2b, 1:1)
      real axmax, aymax
c
c.. declaration of local variables
      integer i1, i2, i3, i4
      real vx, vy, vx0, vy0, dvx, dvy
      real accel
c
c.. x component
      vx0 = xlo(3)
      vy0 = xlo(4)
      dvx = dx(3)
      dvy = dx(4)
      axmax = 0.0
      do i3 = nv3a, nv3b+1
      do i4 = nv4a, nv4b
        vy = vy0+(0.5+i4)*dvy
        do i1 = nv1a, nv1b
        do i2 = nv2a, nv2b
          accel = charge_per_mass*(em_vars(i1, i2, 1) +
     *       vy * em_vars(i1, i2, 6) -
     *       vz(i1, i2, 1)*em_vars(i1, i2, 5))
          vel3(i3, i4, i1, i2) = accel
          axmax = max(axmax, abs(accel))
        end do
        end do
      end do
      end do
c
c.. y component
      aymax = 0.0
      do i4 = nv4a, nv4b+1
      do i1 = nv1a, nv1b
      do i2 = nv2a, nv2b
      do i3 = nv3a, nv3b
        vx = vx0+(0.5+i3)*dvx
        accel = charge_per_mass*(em_vars(i1, i2, 2) +
     *     vz(i1, i2, 1) * em_vars(i1, i2, 4) -
     *     vx * em_vars(i1, i2, 6))
        vel4(i4, i1, i2, i3) = accel
        aymax = max(aymax, abs(accel))
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
      subroutine artVisFlux4D(
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     nf1a,nf1b,nf2a,nf2b,nf3a,nf3b,nf4a,nf4b,
     *     u,flux,dx,dir )
c
      implicit none
c
c.. declarations of incoming variables
      integer nd1a,nd1b,nd2a,nd2b
      integer nd3a,nd3b,nd4a,nd4b
      integer nf1a,nf1b,nf2a,nf2b
      integer nf3a,nf3b,nf4a,nf4b
      integer dir
      real u( nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b )
      real flux( nf1a:nf1b,nf2a:nf2b,nf3a:nf3b,nf4a:nf4b )
      real dx
c
c.. declarations of local variables
c
c.. declarations of local variables
      integer f1a,f1b
      integer i1,i2,i3,i4
      real mu
c
      mu = 1.0e-1
c
      f1a = nf1a+2
      f1b = nf1b-2
c
      mu = mu*dx
c
      if( dir.eq.1 ) then
        do i4 = nf4a,nf4b    
        do i3 = nf3a,nf3b    
        do i2 = nf2a,nf2b    
        do i1 = f1a,f1b 
          flux(i1,i2,i3,i4) = flux(i1,i2,i3,i4)+mu*(
     *       -u(i1-2,i2,i3,i4)+
     *       3.0*u(i1-1,i2,i3,i4)-
     *       3.0*u(i1,i2,i3,i4)+
     *       u(i1+1,i2,i3,i4))
        end do
        end do
        end do
        end do
      else if( dir.eq.2 ) then
        do i4 = nf4a,nf4b    
        do i3 = nf3a,nf3b    
        do i2 = nf2a,nf2b    
        do i1 = f1a,f1b 
          flux(i1,i2,i3,i4) = flux(i1,i2,i3,i4)+mu*(
     *       -u(i4,i1-2,i2,i3)+
     *       3.0*u(i4,i1-1,i2,i3)-
     *       3.0*u(i4,i1,i2,i3)+
     *       u(i4,i1+1,i2,i3))
        end do
        end do
        end do
        end do
      else if( dir.eq.3 ) then
        do i4 = nf4a,nf4b    
        do i3 = nf3a,nf3b    
        do i2 = nf2a,nf2b    
        do i1 = f1a,f1b 
          flux(i1,i2,i3,i4) = flux(i1,i2,i3,i4)+mu*(
     *       -u(i3,i4,i1-2,i2)+
     *       3.0*u(i3,i4,i1-1,i2)-
     *       3.0*u(i3,i4,i1,i2)+
     *       u(i3,i4,i1+1,i2))
        end do
        end do
        end do
        end do
      else if( dir.eq.4 ) then
        do i4 = nf4a,nf4b    
        do i3 = nf3a,nf3b    
        do i2 = nf2a,nf2b    
        do i1 = f1a,f1b 
          flux(i1,i2,i3,i4) = flux(i1,i2,i3,i4)+mu*(
     *       -u(i2,i3,i4,i1-2)+
     *       3.0*u(i2,i3,i4,i1-1)-
     *       3.0*u(i2,i3,i4,i1)+
     *       u(i2,i3,i4,i1+1))
        end do
        end do
        end do
        end do
c
      end if !dir
c
      return
      end
c
c ++++++++++++++
c
      subroutine SKLimit4D(
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     nf1a,nf1b,nf2a,nf2b,nf3a,nf3b,nf4a,nf4b,
     *     dir,cell,face,
     *     vel,temp )
c
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
      real vel(  nf1a:nf1b,nf2a:nf2b,nf3a:nf3b,nf4a:nf4b )
      real temp( nf1a:nf1b,nf2a:nf2b,nf3a:nf3b,nf4a:nf4b )
c
c.. declarations of local variables
      integer f1a,f1b
      integer i1,i2,i3,i4
      real skC,dc,dl,dr,sg,dLim,val(2)
c
      skC = 1.25
c
      f1a = nf1a+3
      f1b = nf1b-3
c
c.. first part of the limiter here ... this constrains the approximations at cell faces
      if( dir.eq.1 ) then
        do i4 = nf4a,nf4b    
        do i3 = nf3a,nf3b    
        do i2 = nf2a,nf2b    
        do i1 = f1a,f1b  
          if( face(i1,i2,i3,i4).lt.
     *         min(cell(i1-1,i2,i3,i4),cell(i1,i2,i3,i4)).or.
     *        face(i1,i2,i3,i4).gt.
     *         max(cell(i1-1,i2,i3,i4),cell(i1,i2,i3,i4)) ) then
            dc = cell(i1-1,i2,i3,i4)-
     *           2.0*face(i1,i2,i3,i4)+cell(i1,i2,i3,i4)
            dl = cell(i1-2,i2,i3,i4)-
     *           2.0*cell(i1-1,i2,i3,i4)+cell(i1,i2,i3,i4)
            dr = cell(i1-1,i2,i3,i4)-
     *           2.0*cell(i1,i2,i3,i4)+cell(i1+1,i2,i3,i4)

            if( dc.gt.0.0 .and. dl.gt.0.0 .and. dr.gt.0.0 ) then
              sg = 1.0
            else if( dc.lt.0.0 .and. dl.lt.0.0 .and. dr.lt.0.0 ) then
              sg = -1.0
            else
              sg = 0.0
            end if
            dLim = sg*(min(skC*min(abs(dl),abs(dr)),abs(dc)))
            face(i1,i2,i3,i4) = 
     *           0.5*(cell(i1-1,i2,i3,i4)+cell(i1,i2,i3,i4)-dLim)
          end if
        end do
        end do
        end do
        end do
c
      else if( dir.eq.2 ) then
        do i4 = nf4a,nf4b    
        do i3 = nf3a,nf3b 
        do i2 = nf2a,nf2b    
        do i1 = f1a,f1b 
          if( face(i1,i2,i3,i4).lt.
     *         min(cell(i2,i3,i4,i1-1),cell(i2,i3,i4,i1)).or.
     *        face(i1,i2,i3,i4).gt.
     *         max(cell(i2,i3,i4,i1-1),cell(i2,i3,i4,i1)) ) then
            dc = cell(i2,i3,i4,i1-1)-
     *           2.0*face(i1,i2,i3,i4)+cell(i2,i3,i4,i1)
            dl = cell(i2,i3,i4,i1-2)-
     *           2.0*cell(i2,i3,i4,i1-1)+cell(i2,i3,i4,i1)
            dr = cell(i2,i3,i4,i1-1)-
     *           2.0*cell(i2,i3,i4,i1)+cell(i2,i3,i4,i1+1)

            if( dc.gt.0.0 .and. dl.gt.0.0 .and. dr.gt.0.0 ) then
              sg = 1.0
            else if( dc.lt.0.0 .and. dl.lt.0.0 .and. dr.lt.0.0 ) then
              sg = -1.0
            else
              sg = 0.0
            end if
            dLim = sg*(min(skC*min(abs(dl),abs(dr)),abs(dc)))
            face(i1,i2,i3,i4) = 
     *           0.5*(cell(i2,i3,i4,i1-1)+cell(i2,i3,i4,i1)-dLim)
          end if
        end do
        end do
        end do
        end do
      else if( dir.eq.3 ) then
        do i4 = nf4a,nf4b    
        do i3 = nf3a,nf3b 
        do i2 = nf2a,nf2b    
        do i1 = f1a,f1b 
          if( face(i1,i2,i3,i4).lt.
     *         min(cell(i3,i4,i1-1,i2),cell(i3,i4,i1,i2)).or.
     *        face(i1,i2,i3,i4).gt.
     *         max(cell(i3,i4,i1-1,i2),cell(i3,i4,i1,i2)) ) then
            dc = cell(i3,i4,i1-1,i2)-
     *           2.0*face(i1,i2,i3,i4)+cell(i3,i4,i1,i2)
            dl = cell(i3,i4,i1-2,i2)-
     *           2.0*cell(i3,i4,i1-1,i2)+cell(i3,i4,i1,i2)
            dr = cell(i3,i4,i1-1,i2)-
     *           2.0*cell(i3,i4,i1,i2)+cell(i3,i4,i1+1,i2)

            if( dc.gt.0.0 .and. dl.gt.0.0 .and. dr.gt.0.0 ) then
              sg = 1.0
            else if( dc.lt.0.0 .and. dl.lt.0.0 .and. dr.lt.0.0 ) then
              sg = -1.0
            else
              sg = 0.0
            end if
            dLim = sg*(min(skC*min(abs(dl),abs(dr)),abs(dc)))
            face(i1,i2,i3,i4) = 
     *           0.5*(cell(i3,i4,i1-1,i2)+cell(i3,i4,i1,i2)-dLim)
          end if
        end do
        end do
        end do
        end do
      else if( dir.eq.4 ) then
        do i4 = nf4a,nf4b    
        do i3 = nf3a,nf3b 
        do i2 = nf2a,nf2b    
        do i1 = f1a,f1b 
          if( face(i1,i2,i3,i4).lt.
     *         min(cell(i4,i1-1,i2,i3),cell(i4,i1,i2,i3)).or.
     *        face(i1,i2,i3,i4).gt.
     *         max(cell(i4,i1-1,i2,i3),cell(i4,i1,i2,i3)) ) then
            dc = cell(i4,i1-1,i2,i3)-
     *           2.0*face(i1,i2,i3,i4)+cell(i4,i1,i2,i3)
            dl = cell(i4,i1-2,i2,i3)-
     *           2.0*cell(i4,i1-1,i2,i3)+cell(i4,i1,i2,i3)
            dr = cell(i4,i1-1,i2,i3)-
     *           2.0*cell(i4,i1,i2,i3)+cell(i4,i1+1,i2,i3)

            if( dc.gt.0.0 .and. dl.gt.0.0 .and. dr.gt.0.0 ) then
              sg = 1.0
            else if( dc.lt.0.0 .and. dl.lt.0.0 .and. dr.lt.0.0 ) then
              sg = -1.0
            else
              sg = 0.0
            end if
            dLim = sg*(min(skC*min(abs(dl),abs(dr)),abs(dc)))
            face(i1,i2,i3,i4) = 
     *           0.5*(cell(i4,i1-1,i2,i3)+cell(i4,i1,i2,i3)-dLim)
          end if
        end do
        end do
        end do
        end do
      end if
c
c.. second part of limiter here ... this does the parabolic reconstruction and upwind determination
      if( dir.eq.1 ) then
        do i4 = nf4a,nf4b    
        do i3 = nf3a,nf3b 
        do i2 = nf2a,nf2b
        do i1 = f1a,f1b
          if( vel(i1,i2,i3,i4).gt.0.0 ) then
            ! we are really looking at cell i1-1,i2,i3,i4
            call ppmFit4D( cell(i1-3,i2,i3,i4),cell(i1-2,i2,i3,i4),
     *           cell(i1-1,i2,i3,i4),cell(i1,i2,i3,i4),
     *           cell(i1+1,i2,i3,i4),face(i1-1,i2,i3,i4),
     *           face(i1,i2,i3,i4),val )
            temp(i1,i2,i3,i4) = val(1)
          else
            ! we are really looking at cell i1,i2
            call ppmFit4D( cell(i1-2,i2,i3,i4),cell(i1-1,i2,i3,i4),
     *           cell(i1,i2,i3,i4),cell(i1+1,i2,i3,i4),
     *           cell(i1+2,i2,i3,i4),face(i1,i2,i3,i4),
     *           face(i1+1,i2,i3,i4),val )
            temp(i1,i2,i3,i4) = val(2)
          end if
        end do
        end do
        end do
        end do
c
      else if( dir.eq.2 ) then
        do i4 = nf4a,nf4b    
        do i3 = nf3a,nf3b 
        do i2 = nf2a,nf2b
        do i1 = f1a,f1b
          if( vel(i1,i2,i3,i4).gt.0.0 ) then
            ! we are really looking at cell i1-1,i2,i3,i4
            call ppmFit4D( cell(i2,i3,i4,i1-3),cell(i2,i3,i4,i1-2),
     *           cell(i2,i3,i4,i1-1),cell(i2,i3,i4,i1),
     *           cell(i2,i3,i4,i1+1),face(i1-1,i2,i3,i4),
     *           face(i1,i2,i3,i4),val )
            temp(i1,i2,i3,i4) = val(1)
          else
            ! we are really looking at cell i1,i2
            call ppmFit4D( cell(i2,i3,i4,i1-2),cell(i2,i3,i4,i1-1),
     *           cell(i2,i3,i4,i1),cell(i2,i3,i4,i1+1),
     *           cell(i2,i3,i4,i1+2),face(i1,i2,i3,i4),
     *           face(i1+1,i2,i3,i4),val )
            temp(i1,i2,i3,i4) = val(2)
          end if
        end do
        end do
        end do
        end do
c
      else if( dir.eq.3 ) then
        do i4 = nf4a,nf4b    
        do i3 = nf3a,nf3b 
        do i2 = nf2a,nf2b
        do i1 = f1a,f1b
          if( vel(i1,i2,i3,i4).gt.0.0 ) then
            ! we are really looking at cell i1-1,i2,i3,i4
            call ppmFit4D( cell(i3,i4,i1-3,i2),cell(i3,i4,i1-2,i2),
     *           cell(i3,i4,i1-1,i2),cell(i3,i4,i1,i2),
     *           cell(i3,i4,i1+1,i2),face(i1-1,i2,i3,i4),
     *           face(i1,i2,i3,i4),val )
            temp(i1,i2,i3,i4) = val(1)
          else
            ! we are really looking at cell i1,i2
            call ppmFit4D( cell(i3,i4,i1-2,i2),cell(i3,i4,i1-1,i2),
     *           cell(i3,i4,i1,i2),cell(i3,i4,i1+1,i2),
     *           cell(i3,i4,i1+2,i2),face(i1,i2,i3,i4),
     *           face(i1+1,i2,i3,i4),val )
            temp(i1,i2,i3,i4) = val(2)
          end if
        end do
        end do
        end do
        end do
      else if( dir.eq.4 ) then
        do i4 = nf4a,nf4b    
        do i3 = nf3a,nf3b 
        do i2 = nf2a,nf2b
        do i1 = f1a,f1b
          if( vel(i1,i2,i3,i4).gt.0.0 ) then
            ! we are really looking at cell i1-1,i2,i3,i4
            call ppmFit4D( cell(i4,i1-3,i2,i3),cell(i4,i1-2,i2,i3),
     *           cell(i4,i1-1,i2,i3),cell(i4,i1,i2,i3),
     *           cell(i4,i1+1,i2,i3),face(i1-1,i2,i3,i4),
     *           face(i1,i2,i3,i4),val )
            temp(i1,i2,i3,i4) = val(1)
          else
            ! we are really looking at cell i1,i2
            call ppmFit4D( cell(i4,i1-2,i2,i3),cell(i4,i1-1,i2,i3),
     *           cell(i4,i1,i2,i3),cell(i4,i1+1,i2,i3),
     *           cell(i4,i1+2,i2,i3),face(i1,i2,i3,i4),
     *           face(i1+1,i2,i3,i4),val )
            temp(i1,i2,i3,i4) = val(2)
          end if
        end do
        end do
        end do
        end do
c
      end if ! dir
c
      do i4 = nf4a,nf4b    
      do i3 = nf3a,nf3b 
      do i2 = nf2a,nf2b
      do i1 = f1a,f1b
        face(i1,i2,i3,i4) = temp(i1,i2,i3,i4)
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
      subroutine ppmFit4D(
     *     um2,um1,u0,up1,up2,
     *     fl,fr,val )
c
      implicit none
c
c.. declarations of incoming variables
      real um2,um1,u0,up1,up2,fl,fr,val(2)
c
c.. declarations of local variables
      real skC,dc,dl,dr,sg,dLim,dj
c
      skC = 1.25
      val(1) = fr
      val(2) = fl
c
      if( (fr-u0)*(u0-fl).le.0.0 .and.
     *    (um1-u0)*(u0-up1).le.0.0 ) then
        dj = 4.0*(fl-2.0*u0+fr)
        dc = um1-2.0*u0+up1
        dl = um2-2.0*um1+u0
        dr = u0-2.0*up1+up2

        if( dj.gt.0.0 .and. dc.gt.0.0 .and. 
     *      dl.gt.0.0 .and. dr.gt.0.0 ) then
          sg = 1.0
        else if( dj.lt.0.0 .and. dc.lt.0.0 .and. 
     *           dl.lt.0.0 .and. dr.lt.0.0 ) then
          sg = -1.0
        else
          sg = 0.0
        end if
        dLim = sg*(min(skC*min(min(abs(dl),abs(dr)),abs(dc)),
     *     abs(dj)))
        if( abs(dj).gt.1.e-10 ) then
          val(1) = u0+(fr-u0)*dLim/dj
          val(2) = u0+(fl-u0)*dLim/dj
        else
          val(1) = u0
          val(2) = u0
        end if
      else if( abs(fr-u0).gt.2.0*abs(fl-u0) ) then
        if( 4.0*abs(fr-2.0*u0+fl).lt.abs(up2-2.0*u0+um2)/4.0 ) then
          val(1) = fr
        else
          val(1) = u0-2.0*(fl-u0)
        end if
      else if( abs(fl-u0).gt.2.0*abs(fr-u0) ) then
        if( 4.0*abs(fr-2.0*u0+fl).lt.abs(up2-2.0*u0+um2)/4.0 ) then
          val(2) = fl
        else
          val(2) = u0-2.0*(fr-u0)
        end if
      end if
c
      return
      end
c
c ++++++++++++++
c
      subroutine WENO43Avg4D(
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     nf1a,nf1b,nf2a,nf2b,nf3a,nf3b,nf4a,nf4b,
     *     dir,cell,vel,face )
c
c.. function to call WENO averaging code.
c  
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
      real vel(  nf1a:nf1b,nf2a:nf2b,nf3a:nf3b,nf4a:nf4b )
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
          call WENO43Fit4D( cell(i1-2,i2,i3,i4),cell(i1-1,i2,i3,i4),
     *                      cell(i1,i2,i3,i4),  cell(i1+1,i2,i3,i4),
     *                      face(i1,i2,i3,i4),  vel(i1,i2,i3,i4) )
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
          call WENO43Fit4D( cell(i4,i1-2,i2,i3),cell(i4,i1-1,i2,i3),
     *                      cell(i4,i1,i2,i3),  cell(i4,i1+1,i2,i3),
     *                      face(i1,i2,i3,i4),  vel(i1,i2,i3,i4) )
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
          call WENO43Fit4D( cell(i3,i4,i1-2,i2),cell(i3,i4,i1-1,i2),
     *                      cell(i3,i4,i1,i2),  cell(i3,i4,i1+1,i2),
     *                      face(i1,i2,i3,i4),  vel(i1,i2,i3,i4) )
        end do
        end do
        end do
        end do
      else
        do i4 = nf4a,nf4b    
        do i3 = nf3a,nf3b    
        do i2 = nf2a,nf2b    
        do i1 = f1a,f1b
c        do i1 = f1a,f1b
c        do i4 = nf4a,nf4b    
c        do i3 = nf3a,nf3b    
c        do i2 = nf2a,nf2b    
          call WENO43Fit4D( cell(i2,i3,i4,i1-2),cell(i2,i3,i4,i1-1),
     *                      cell(i2,i3,i4,i1),  cell(i2,i3,i4,i1+1),
     *                      face(i1,i2,i3,i4),  vel(i1,i2,i3,i4) )
        end do
        end do
        end do
        end do
      end if
c
      return
      end
c
c ++++++++++++++
c
      subroutine WENO43Fit4D(
     *     um2,um1,u0,up1,face,vel )
c
c.. function to get WENO averaged face average values from cell average values
c     This routine should be 4th order for smooth and 3rd order for non-smooth data. We assume a uniform grid. 
      implicit none
c
c.. declarations of incoming variables
      real um2,um1,u0,up1,face,vel
c
c.. declarations of local variables
      real eps
      real fl,fr,bl,br,al,ar,wl,wr
      real c1l,c2l,c1r,c2r
      real wmax,wmin
      real tmp
c
      eps = 1.e-6
c
      ! get left and right 3rd order approximations
c      fl = (-um2+5.0*um1+2.0*u0)/6.0
c      fr = (2.0*um1+5.0*u0-up1)/6.0
      tmp = 1.0/6.0
      fl = tmp*(-um2+5.0*um1+2.0*u0)
      fr = tmp*(2.0*um1+5.0*u0-up1)

      ! get smoothness indicators
      c1l = u0-2.0*um1+um2
      c2l = u0-um2
      c1r = up1-2.0*u0+um1
      c2r = up1-um1
c      bl  = 4.0*(c1l**2)/3.0+0.5*c1l*c2l+0.25*c2l**2
c      br  = 4.0*(c1r**2)/3.0-0.5*c1r*c2r+0.25*c2r**2
      bl  = 8.0*tmp*(c1l**2)+0.5*c1l*c2l+0.25*c2l**2
      br  = 8.0*tmp*(c1r**2)-0.5*c1r*c2r+0.25*c2r**2

      ! get weights
      al = 1.0/((eps+bl)**2)
      ar = 1.0/((eps+br)**2)
c      wl = al/(al+ar)
c      wr = ar/(al+ar)
      tmp = 1.0/(al+ar)
      wl = tmp*al
      wr = tmp*ar

      ! perform mapping of the weights (mapped weno as in Henrick JCP 2005)
      al = wl*(0.75+wl*(wl-1.5))
      ar = wr*(0.75+wr*(wr-1.5))
c      wl = al/(al+ar)
c      wr = ar/(al+ar)
      tmp = 1.0/(al+ar)
      wl = tmp*al
      wr = tmp*ar

      wmax = max(wl,wr)
      wmin = min(wl,wr)
      if( vel.gt.0.0 ) then
        wl = wmax
        wr = wmin
      else
        wl = wmin
        wr = wmax
      end if

      face = (wl*fl+wr*fr)
c
      return
      end

ccccccc
ccccccc
c
c ++++++++++++++
c
      subroutine WENO65Avg4D(
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     nf1a,nf1b,nf2a,nf2b,nf3a,nf3b,nf4a,nf4b,
     *     dir,cell,vel,face )
c
c.. function to call 6/5 BWENO averaging code.
c  
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
      real vel(  nf1a:nf1b,nf2a:nf2b,nf3a:nf3b,nf4a:nf4b )
c
c.. declarations of local variables
      integer f1a,f1b
      integer i1,i2,i3,i4
c
      f1a = nf1a+3
      f1b = nf1b-3
c
      if( dir.eq.1 ) then
        do i4 = nf4a,nf4b    
        do i3 = nf3a,nf3b    
        do i2 = nf2a,nf2b    
        do i1 = f1a,f1b
          call WENO65Fit4D( 
     *         cell(i1-3,i2,i3,i4),
     *         cell(i1-2,i2,i3,i4),
     *         cell(i1-1,i2,i3,i4),
     *         cell(i1,i2,i3,i4),  
     *         cell(i1+1,i2,i3,i4),
     *         cell(i1+2,i2,i3,i4),
     *         face(i1,i2,i3,i4),  
     *         vel(i1,i2,i3,i4) )
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
          call WENO65Fit4D( 
     *         cell(i4,i1-3,i2,i3),
     *         cell(i4,i1-2,i2,i3),
     *         cell(i4,i1-1,i2,i3),
     *         cell(i4,i1,i2,i3),  
     *         cell(i4,i1+1,i2,i3),
     *         cell(i4,i1+2,i2,i3),
     *         face(i1,i2,i3,i4),  
     *         vel(i1,i2,i3,i4) )
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
          call WENO65Fit4D( 
     *         cell(i3,i4,i1-3,i2),
     *         cell(i3,i4,i1-2,i2),
     *         cell(i3,i4,i1-1,i2),
     *         cell(i3,i4,i1,i2),  
     *         cell(i3,i4,i1+1,i2),
     *         cell(i3,i4,i1+2,i2),
     *         face(i1,i2,i3,i4),  
     *         vel(i1,i2,i3,i4) )
        end do
        end do
        end do
        end do
      else
        do i4 = nf4a,nf4b    
        do i3 = nf3a,nf3b    
        do i2 = nf2a,nf2b    
        do i1 = f1a,f1b
c        do i1 = f1a,f1b
c        do i4 = nf4a,nf4b    
c        do i3 = nf3a,nf3b    
c        do i2 = nf2a,nf2b    
          call WENO65Fit4D( 
     *         cell(i2,i3,i4,i1-3),
     *         cell(i2,i3,i4,i1-2),
     *         cell(i2,i3,i4,i1-1),
     *         cell(i2,i3,i4,i1),  
     *         cell(i2,i3,i4,i1+1),
     *         cell(i2,i3,i4,i1+2),
     *         face(i1,i2,i3,i4),  
     *         vel(i1,i2,i3,i4) )
        end do
        end do
        end do
        end do
      end if
c
      return
      end
c
c ++++++++++++++
c
      subroutine WENO65Fit4D(
     *     um3,um2,um1,u0,up1,up2,face,vel )
c
c.. function to get 6/5 WENO averaged face values
c     This routine should be 6th order for smooth and 5rd order for non-smooth data. We assume a uniform grid. 
      implicit none
c
c.. declarations of incoming variables
      real um3,um2,um1,u0,up1,up2,face,vel
c
c.. declarations of local variables
      real eps
      real fl,fr,bl,br,al,ar,wl,wr
      real wmax,wmin
c
      eps = 1.e-10
c
      ! get left and right 5th order approximations
      fl = ( 2.0*um3-13.0*um2+47.0*um1+27.0*u0 -3.0*up1)/60.0
      fr = (-3.0*um2+27.0*um1+47.0*u0 -13.0*up1+2.0*up2)/60.0

      ! get smoothness indicators
      bl = 0.5489E4 / 0.105E3 * um1 ** 2 + (-0.2242428E7 * u0 - 0.1887
     *108E7 * um2 + 0.410226E6 * um3 + 0.557646E6 * up1) * um1 / 0.30240
     *E5 + 0.75329E5 / 0.3780E4 * um2 ** 2 + (0.1259696E7 * u0 - 0.27531
     *8E6 * um3 - 0.302534E6 * up1) * um2 / 0.30240E5 + 0.33727E5 / 0.30
     *240E5 * um3 ** 2 + (-0.264314E6 * u0 + 0.61952E5 * up1) * um3 / 0.
     *30240E5 + 0.106409E6 / 0.3780E4 * u0 ** 2 - 0.227749E6 / 0.15120E5
     * * u0 * up1 + 0.69217E5 / 0.30240E5 * up1 ** 2

      br = 0.106409E6 / 0.3780E4 * um1 ** 2 + (-0.2242428E7 * u0 - 0.4
     *55498E6 * um2 + 0.1259696E7 * up1 - 0.264314E6 * up2) * um1 / 0.30
     *240E5 + 0.69217E5 / 0.30240E5 * um2 ** 2 + (0.557646E6 * u0 - 0.30
     *2534E6 * up1 + 0.61952E5 * up2) * um2 / 0.30240E5 + 0.75329E5 / 0.
     *3780E4 * up1 ** 2 + (-0.1887108E7 * u0 - 0.275318E6 * up2) * up1 /
     * 0.30240E5 + 0.5489E4 / 0.105E3 * u0 ** 2 + 0.68371E5 / 0.5040E4 *
     * u0 * up2 + 0.33727E5 / 0.30240E5 * up2 ** 2

      ! get weights
      al = 1.0/((eps+bl)**2)
      ar = 1.0/((eps+br)**2)
      wl = al/(al+ar)
      wr = ar/(al+ar)

      ! perform mapping of the weights (mapped weno as in Henrick JCP 2005)
c      al = 16.0*(wl-0.5)^5+0.5
c      ar = 16.0*(wr-0.5)^5+0.5
      al = wl*(0.75+wl*(wl-1.5))
      ar = wr*(0.75+wr*(wr-1.5))
      wl = al/(al+ar)
      wr = ar/(al+ar)

      wmax = max(wl,wr)
      wmin = min(wl,wr)
      if( vel.gt.0.0 ) then
        wl = wmax
        wr = wmin
      else
        wl = wmin
        wr = wmax
      end if

      face = (wl*fl+wr*fr)
c
      return
      end
ccccccc
ccccccc
c
c ++++++++++++++
c
      subroutine accumfluxdiv4D(
     *     rhs,
     *     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     *     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     *     fluxx1,fluxx2,fluxx3,fluxx4,
     *     xlo, xhi, deltax)
c
c.. add updates ...
      implicit none
c
c.. declarations of incoming variables 
      integer nd1a,nd1b,nd2a,nd2b
      integer nd3a,nd3b,nd4a,nd4b
      integer n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b
c
c      real rhs ( n1a:n1b,n2a:n2b,n3a:n3b,n4a:n4b )
      real rhs ( nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b )
      real fluxx1( nd1a:nd1b+1,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b )
      real fluxx2( nd2a:nd2b+1,nd3a:nd3b,nd4a:nd4b,nd1a:nd1b )
      real fluxx3( nd3a:nd3b+1,nd4a:nd4b,nd1a:nd1b,nd2a:nd2b )
      real fluxx4( nd4a:nd4b+1,nd1a:nd1b,nd2a:nd2b,nd3a:nd3b )
c
      real xlo(1:4),xhi(1:4),deltax(1:4)
c.. declarations of local variables
      integer i1,i2,i3,i4
      real dx, dy, dvx, dvy
      real temp
c
      dx = deltax(1)
      dy = deltax(2)
      dvx = deltax(3)
      dvy = deltax(4)
      do i4 = n4a,n4b
      do i3 = n3a,n3b
      do i2 = n2a,n2b
      do i1 = n1a,n1b
        temp = -(fluxx1(i1+1,i2,i3,i4)-fluxx1(i1,i2,i3,i4))/dx
     *         -(fluxx2(i2+1,i3,i4,i1)-fluxx2(i2,i3,i4,i1))/dy
     *         -(fluxx3(i3+1,i4,i1,i2)-fluxx3(i3,i4,i1,i2))/dvx
     *         -(fluxx4(i4+1,i1,i2,i3)-fluxx4(i4,i1,i2,i3))/dvx
        rhs(i1,i2,i3,i4) = temp
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
      subroutine setAccelerationBCs4D( 
     &     u,
     &     ng1a,ng1b,ng2a,ng2b,ng3a,ng3b,ng4a,ng4b,
     &     nl1a,nl1b,nl2a,nl2b,nl3a,nl3b,nl4a,nl4b,
     &     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     &     solution_order,
     &     vel3,vel4,
     &     xlo, xhi, dx,
     &     ic_phase_param )
c
c.. boundary conditions 
c
c.. declarations of incoming variables
      implicit none
      integer ng1a,ng1b,ng2a,ng2b
      integer ng3a,ng3b,ng4a,ng4b
      integer nl1a,nl1b,nl2a,nl2b
      integer nl3a,nl3b,nl4a,nl4b
      integer n1a,n1b,n2a,n2b
      integer n3a,n3b,n4a,n4b
      integer solution_order
      real u(nl1a:nl1b,nl2a:nl2b,nl3a:nl3b,nl4a:nl4b)
      real vel3( nl3a:nl3b+1,nl4a:nl4b,nl1a:nl1b,nl2a:nl2b )
      real vel4( nl4a:nl4b+1,nl1a:nl1b,nl2a:nl2b,nl3a:nl3b )
      real xlo(1:4), xhi(1:4), dx(1:4)
      real ic_phase_param(*)
c
c.. declarations of local variables
      integer nghosts
      integer i1,i2,i3,i4
      real alpha,beta,vx0,vy0,vflowinitx,vflowinity
      real vxlo,vylo,dvx,dvy
      real fgt1,fgt2,fgt3
      real fgb1,fgb2,fgb3
      real vx,vy
      real pi
c
      real fx
      real fy
c
      real one,four
c
      one = 1.0
      four = 4.0
c
      pi  = four*atan(one)
c
      alpha      = ic_phase_param(1)
      beta       = ic_phase_param(2)
      vx0        = ic_phase_param(3)
      vy0        = ic_phase_param(4)
      vflowinitx = ic_phase_param(5)
      vflowinity = ic_phase_param(6)
c
      vxlo = xlo(3)
      vylo = xlo(4)
      dvx = dx(3)
      dvy = dx(4)
      if (solution_order .eq. 4) then
        nghosts = 2
      else
        nghosts = 3
      end if
c
c.. set i3 boundaries
      if ((ng3b-nghosts .eq. n3b) .or. (ng3a+nghosts .eq. n3a)) then
      do i4 = nl4a,nl4b
cccc
cccc  sample distribution on top
        i3 = n3b
        vx = vxlo+(0.5+(i3+1))*dvx
        vy = vylo+(0.5+i4)*dvy
        fx = exp(-0.5*((alpha*(vx-vx0-vflowinitx))**2))
        fy = exp(-0.5*((beta*(vy-vy0-vflowinity))**2))
        fgt1 = fx*fy
        fgt1 = alpha*beta/(2.0*pi)*fgt1
c
        vx = vxlo+(0.5+(i3+2))*dvx
        fx = exp(-0.5*((alpha*(vx-vx0-vflowinitx))**2))
        fgt2 = fx*fy
        fgt2 = alpha*beta/(2.0*pi)*fgt2
c
        if (solution_order .eq. 6) then
          vx = vxlo+(0.5+(i3+3))*dvx
          fx = exp(-0.5*((alpha*(vx-vx0-vflowinitx))**2))
          fgt3 = fx*fy
          fgt3 = alpha*beta/(2.0*pi)*fgt3
        end if
cccc
cccc  sample distribution on bottom
        i3 = n3a
        vx = vxlo+(0.5+(i3-1))*dvx
        fx = exp(-0.5*((alpha*(vx-vx0-vflowinitx))**2))
        fgb1 = fx*fy
        fgb1 = alpha*beta/(2.0*pi)*fgb1
c
        vx = vxlo+(0.5+(i3-2))*dvx
        fx = exp(-0.5*((alpha*(vx-vx0-vflowinitx))**2))
        fgb2 = fx*fy
        fgb2 = alpha*beta/(2.0*pi)*fgb2
c
        if (solution_order .eq. 6) then
          vx = vxlo+(0.5+(i3-3))*dvx
          fx = exp(-0.5*((alpha*(vx-vx0-vflowinitx))**2))
          fgb3 = fx*fy
          fgb3 = alpha*beta/(2.0*pi)*fgb3
        end if
cccc
c
        if( ng3b-nghosts.eq.n3b ) then
        do i2 = nl2a,nl2b
        do i1 = nl1a,nl1b
          ! top
          if( vel3(n3b+1,i4,i1,i2) .ge. 0.0 ) then
            u(i1,i2,n3b+1,i4) = 3.0*u(i1,i2,n3b,i4)-
     *           3.0*u(i1,i2,n3b-1,i4)+u(i1,i2,n3b-2,i4)
            u(i1,i2,n3b+2,i4) = 3.0*u(i1,i2,n3b+1,i4)-
     *           3.0*u(i1,i2,n3b,i4)+u(i1,i2,n3b-1,i4)
            if (solution_order .eq. 6) then
              u(i1,i2,n3b+3,i4) = 3.0*u(i1,i2,n3b+2,i4)-
     *             3.0*u(i1,i2,n3b+1,i4)+u(i1,i2,n3b,i4)
            end if
          else
            u(i1,i2,n3b+1,i4) = fgt1
            u(i1,i2,n3b+2,i4) = fgt2
            if (solution_order .eq. 6) then
              u(i1,i2,n3b+3,i4) = fgt3
            end if
          end if
        end do
        end do
        end if

        if( ng3a+nghosts.eq.n3a ) then
        do i2 = nl2a,nl2b
        do i1 = nl1a,nl1b

          ! bottom
          if( vel3(n3a,i4,i1,i2) .gt. 0.0 ) then
            u(i1,i2,n3a-1,i4) = fgb1
            u(i1,i2,n3a-2,i4) = fgb2
            if (solution_order .eq. 6) then
              u(i1,i2,n3a-3,i4) = fgb3
            end if
          else
            u(i1,i2,n3a-1,i4) = 3.0*u(i1,i2,n3a,i4)-
     *           3.0*u(i1,i2,n3a+1,i4)+u(i1,i2,n3a+2,i4)
            u(i1,i2,n3a-2,i4) = 3.0*u(i1,i2,n3a-1,i4)-
     *           3.0*u(i1,i2,n3a,i4)+u(i1,i2,n3a+1,i4)
            if (solution_order .eq. 6) then
              u(i1,i2,n3a-3,i4) = 3.0*u(i1,i2,n3a-2,i4)-
     *             3.0*u(i1,i2,n3a-1,i4)+u(i1,i2,n3a,i4)
            end if
          end if
        end do
        end do
        end if

      end do
      end if
c
c.. set i4 boundaries
      if ((ng4b-nghosts.eq.n4b) .or. (ng4a+nghosts .eq. n4a)) then
      do i3 = nl3a,nl3b
cccc
cccc  sample distribution on top
        i4 = n4b
        vx = vxlo+(0.5+i3)*dvx
        vy = vylo+(0.5+(i4+1))*dvy
        fx = exp(-0.5*((alpha*(vx-vx0-vflowinitx))**2))
        fy = exp(-0.5*((beta*(vy-vy0-vflowinity))**2))
        fgt1 = fx*fy
        fgt1 = alpha*beta/(2.0*pi)*fgt1
c
        vy = vylo+(0.5+(i4+2))*dvy
        fy = exp(-0.5*((beta*(vy-vy0-vflowinity))**2))
        fgt2 = fx*fy
        fgt2 = alpha*beta/(2.0*pi)*fgt2
c
        if (solution_order .eq. 6) then
          vy = vylo+(0.5+(i4+3))*dvy
          fy = exp(-0.5*((beta*(vy-vy0-vflowinity))**2))
          fgt3 = fx*fy
          fgt3 = alpha*beta/(2.0*pi)*fgt3
        end if
cccc
cccc  sample distribution on bottom
        i4 = n4a
        vy = vylo+(0.5+(i4-1))*dvy
        fy = exp(-0.5*((beta*(vy-vy0-vflowinity))**2))
        fgb1 = fx*fy
        fgb1 = alpha*beta/(2.0*pi)*fgb1
c
        vy = vylo+(0.5+(i4-2))*dvy
        fy = exp(-0.5*((beta*(vy-vy0-vflowinity))**2))
        fgb2 = fx*fy
        fgb2 = alpha*beta/(2.0*pi)*fgb2
c
        if (solution_order .eq. 6) then
          vy = vylo+(0.5+(i4-3))*dvy
          fy = exp(-0.5*((beta*(vy-vy0-vflowinity))**2))
          fgb3 = fx*fy
          fgb3 = alpha*beta/(2.0*pi)*fgb3
        end if
c
        if( ng4b-nghosts.eq.n4b ) then
        do i2 = nl2a,nl2b
        do i1 = nl1a,nl1b
          ! top
          if( vel4(n4b+1,i1,i2,i3) .ge. 0.0 ) then
            u(i1,i2,i3,n4b+1) = 3.0*u(i1,i2,i3,n4b)-
     *           3.0*u(i1,i2,i3,n4b-1)+u(i1,i2,i3,n4b-2)
            u(i1,i2,i3,n4b+2) = 3.0*u(i1,i2,i3,n4b+1)-
     *           3.0*u(i1,i2,i3,n4b)+u(i1,i2,i3,n4b-1)
            if (solution_order .eq. 6) then
              u(i1,i2,i3,n4b+3) = 3.0*u(i1,i2,i3,n4b+2)-
     *             3.0*u(i1,i2,i3,n4b+1)+u(i1,i2,i3,n4b)
            end if
          else
            u(i1,i2,i3,n4b+1) = fgt1
            u(i1,i2,i3,n4b+2) = fgt2
            if (solution_order .eq. 6) then
              u(i1,i2,i3,n4b+3) = fgt3
            end if
          end if
        end do
        end do
        end if

        if( ng4a+nghosts.eq.n4a ) then
        do i2 = nl2a,nl2b
        do i1 = nl1a,nl1b
          ! bottom
          if( vel4(n4a,i1,i2,i3) .gt. 0.0 ) then
            u(i1,i2,i3,n4a-1) = fgb1
            u(i1,i2,i3,n4a-2) = fgb2
            if (solution_order .eq. 6) then
              u(i1,i2,i3,n4a-3) = fgb3
            end if
          else
            u(i1,i2,i3,n4a-1) = 3.0*u(i1,i2,i3,n4a)-
     *           3.0*u(i1,i2,i3,n4a+1)+u(i1,i2,i3,n4a+2)
            u(i1,i2,i3,n4a-2) = 3.0*u(i1,i2,i3,n4a-1)-
     *           3.0*u(i1,i2,i3,n4a)+u(i1,i2,i3,n4a+1)
            if (solution_order .eq. 6) then
              u(i1,i2,i3,n4a-3) = 3.0*u(i1,i2,i3,n4a-2)-
     *             3.0*u(i1,i2,i3,n4a-1)+u(i1,i2,i3,n4a)
            end if
          end if
        end do
        end do
        end if

      end do
      end if
c
      return
      end
c
c +++++++++++++
c
      subroutine setAdvectionBCs4D( 
     &     u,
     &     ng1a,ng1b,ng2a,ng2b,ng3a,ng3b,ng4a,ng4b,
     &     nl1a,nl1b,nl2a,nl2b,nl3a,nl3b,nl4a,nl4b,
     &     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     &     solution_order,
     &     vel1,vel2,
     &     xlo, xhi, dx,
     &     ic_phase_param,
     &     xPeriodic, yPeriodic )
c
c.. boundary conditions 
c
c.. declarations of incoming variables
      implicit none
      integer ng1a,ng1b,ng2a,ng2b
      integer ng3a,ng3b,ng4a,ng4b
      integer nl1a,nl1b,nl2a,nl2b
      integer nl3a,nl3b,nl4a,nl4b
      integer n1a,n1b,n2a,n2b
      integer n3a,n3b,n4a,n4b
      integer solution_order
      integer xPeriodic,yPeriodic
      real u(nl1a:nl1b,nl2a:nl2b,nl3a:nl3b,nl4a:nl4b)
      real vel1( nl1a:nl1b+1,nl2a:nl2b,nl3a:nl3b,nl4a:nl4b )
      real vel2( nl2a:nl2b+1,nl3a:nl3b,nl4a:nl4b,nl1a:nl1b )
      real xlo(1:4), xhi(1:4), dx(1:4)
      real ic_phase_param(*)
c
c.. declarations of local variables
      integer nghosts
      integer i1,i2,i3,i4
      real alpha,beta,vx0,vy0,vflowinitx,vflowinity
      real vxlo,vylo,dvx,dvy
      real fg1
      real vx,vy
      real pi
c
      real fx
      real fy
c
      real one,four
c
      one = 1.0
      four = 4.0
c
      pi  = four*atan(one)
c
      alpha      = ic_phase_param(1)
      beta       = ic_phase_param(2)
      vx0        = ic_phase_param(3)
      vy0        = ic_phase_param(4)
      vflowinitx = ic_phase_param(5)
      vflowinity = ic_phase_param(6)
c
      vxlo = xlo(3)
      vylo = xlo(4)
      dvx = dx(3)
      dvy = dx(4)
      if (solution_order .eq. 4) then
        nghosts = 2
      else
        nghosts = 3
      end if
c.. set i1 boundaries (if not periodic)
      if ((xPeriodic .ne. 1) .and.
     *    ((ng1b-nghosts .eq. n1b) .or. (ng1a+nghosts .eq. n1a))) then
        do i4 = nl4a,nl4b
          vy = vylo+(0.5+i4)*dvy
        do i3 = nl3a,nl3b
          ! sample the background distribution
          vx = vxlo+(0.5+i3)*dvx

          fx = exp(-0.5*((alpha*(vx-vx0-vflowinitx))**2))
          fy = exp(-0.5*((beta*(vy-vy0-vflowinity))**2))
          fg1 = fx*fy
          fg1 = alpha*beta/(2.0*pi)*fg1

          if( ng1b-nghosts.eq.n1b ) then
          do i2 = nl2a,nl2b
            ! top
            if( vel1(n1b+1,i2,i3,i4) .ge. 0.0 ) then
              u(n1b+1,i2,i3,i4) = 3.0*u(n1b,i2,i3,i4)-
     *             3.0*u(n1b-1,i2,i3,i4)+u(n1b-2,i2,i3,i4)
              u(n1b+2,i2,i3,i4) = 3.0*u(n1b+1,i2,i3,i4)-
     *             3.0*u(n1b,i2,i3,i4)+u(n1b-1,i2,i3,i4)
              if (solution_order .eq. 6) then
                u(n1b+3,i2,i3,i4) = 3.0*u(n1b+2,i2,i3,i4)-
     *               3.0*u(n1b+1,i2,i3,i4)+u(n1b,i2,i3,i4)
              end if
            else
              ! sample the background distribution
              u(n1b+1,i2,i3,i4) = fg1
              u(n1b+2,i2,i3,i4) = fg1
              if (solution_order .eq. 6) then
                u(n1b+3,i2,i3,i4) = fg1
              end if
            end if
          end do
          end if

          if( ng1a+nghosts.eq.n1a ) then
          do i2 = nl2a,nl2b
            ! bottom
            if( vel1(n1a,i2,i3,i4) .gt. 0.0 ) then
              ! sample the background distribution
              u(n1a-1,i2,i3,i4) = fg1
              u(n1a-2,i2,i3,i4) = fg1
              if (solution_order .eq. 6) then
                u(n1a-3,i2,i3,i4) = fg1
              end if
            else
              u(n1a-1,i2,i3,i4) = 3.0*u(n1a,i2,i3,i4)-
     *             3.0*u(n1a+1,i2,i3,i4)+u(n1a+2,i2,i3,i4)
              u(n1a-2,i2,i3,i4) = 3.0*u(n1a-1,i2,i3,i4)-
     *             3.0*u(n1a,i2,i3,i4)+u(n1a+1,i2,i3,i4)
              if (solution_order .eq. 6) then
                u(n1a-3,i2,i3,i4) = 3.0*u(n1a-2,i2,i3,i4)-
     *               3.0*u(n1a-1,i2,i3,i4)+u(n1a,i2,i3,i4)
              end if
            end if
              
          end do
          end if
        end do
        end do
      end if
c
c.. set i2 boundaries (if not periodic)
      if ((yPeriodic .ne. 1) .and.
     *    (( ng2b-nghosts.eq.n2b ) .or. (ng2a+nghosts .eq. n2a))) then
        do i4 = nl4a,nl4b
          vy = vylo+(0.5+i4)*dvy
        do i3 = nl3a,nl3b
          vx = vxlo+(0.5+i3)*dvx

          fx = exp(-0.5*((alpha*(vx-vx0-vflowinitx))**2))
          fy = exp(-0.5*((beta*(vy-vy0-vflowinity))**2))
          fg1 = fx*fy
          fg1 = alpha*beta/(2.0*pi)*fg1

          if( ng2b-nghosts.eq.n2b ) then
          do i1 = nl1a,nl1b
            ! top
            if( vel2(n2b+1,i3,i4,i1) .ge. 0.0 ) then
              u(i1,n2b+1,i3,i4) = 3.0*u(i1,n2b,i3,i4)-
     *             3.0*u(i1,n2b-1,i3,i4)+u(i1,n2b-2,i3,i4)
              u(i1,n2b+2,i3,i4) = 3.0*u(i1,n2b+1,i3,i4)-
     *             3.0*u(i1,n2b,i3,i4)+u(i1,n2b-1,i3,i4)
              if (solution_order .eq. 6) then
                u(i1,n2b+3,i3,i4) = 3.0*u(i1,n2b+2,i3,i4)-
     *               3.0*u(i1,n2b+1,i3,i4)+u(i1,n2b,i3,i4)
              end if
            else
              ! sample the background distribution
              u(i1,n2b+1,i3,i4) = fg1
              u(i1,n2b+2,i3,i4) = fg1
              if (solution_order .eq. 6) then
                u(i1,n2b+3,i3,i4) = fg1
              end if
            end if
          end do
          end if

          if( ng2a+nghosts.eq.n2a ) then
          do i1 = nl1a,nl1b
            ! bottom
            if( vel2(n2a,i3,i4,i1) .gt. 0.0 ) then
              ! sample the background distribution
              u(i1,n2a-1,i3,i4) = fg1
              u(i1,n2a-2,i3,i4) = fg1
              if (solution_order .eq. 6) then
                u(i1,n2a-3,i3,i4) = fg1
              end if
            else
              u(i1,n2a-1,i3,i4) = 3.0*u(i1,n2a,i3,i4)-
     *             3.0*u(i1,n2a+1,i3,i4)+u(i1,n2a+2,i3,i4)
              u(i1,n2a-2,i3,i4) = 3.0*u(i1,n2a-1,i3,i4)-
     *             3.0*u(i1,n2a,i3,i4)+u(i1,n2a+1,i3,i4)
              if (solution_order .eq. 6) then
                u(i1,n2a-3,i3,i4) = 3.0*u(i1,n2a-2,i3,i4)-
     *               3.0*u(i1,n2a-1,i3,i4)+u(i1,n2a,i3,i4)
              end if
            end if
          end do
          end if
        end do
        end do
      end if
c
      return
      end
c
c +++++++++++++
c
      subroutine setAccelerationBCs4DJB(
     &     u,
     &     ng1a,ng1b,ng2a,ng2b,ng3a,ng3b,ng4a,ng4b,
     &     nl1a,nl1b,nl2a,nl2b,nl3a,nl3b,nl4a,nl4b,
     &     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     &     solution_order, nGhost,
     &     vel3,vel4,
     &     xlo, xhi, dx,
     &     ic_phase_param )
c
c.. boundary conditions 
c
c.. declarations of incoming variables
      implicit none
      integer ng1a,ng1b,ng2a,ng2b
      integer ng3a,ng3b,ng4a,ng4b
      integer nl1a,nl1b,nl2a,nl2b
      integer nl3a,nl3b,nl4a,nl4b
      integer n1a,n1b,n2a,n2b
      integer n3a,n3b,n4a,n4b
      integer solution_order, nGhost
      real u(nl1a:nl1b,nl2a:nl2b,nl3a:nl3b,nl4a:nl4b)
      real vel3( nl3a:nl3b+1,nl4a:nl4b,nl1a:nl1b,nl2a:nl2b )
      real vel4( nl4a:nl4b+1,nl1a:nl1b,nl2a:nl2b,nl3a:nl3b )
      real xlo(1:4), xhi(1:4), dx(1:4)
      real ic_phase_param(*)
c
c.. declarations of local variables
      integer i1,i2,i3,i4
      integer ig,extrapEq,iSten
      real eCoeffs(1:6,1:6)
      real alpha,beta,vx0,vy0,vflowinitx,vflowinity
      real vxlo,vylo,vxhi,vyhi
      real dvx,dvy
      real fsamp
      real vx,vy
      real pi
c
      real fx
      real fy
c
      real one,four
c
      one = 1.0
      four = 4.0
c
      pi  = four*atan(one)
c
c fill in matrix of extrapolation coefficients (formula,coefficient)
      eCoeffs(1,1) = 1.0

      eCoeffs(2,1) =  2.0
      eCoeffs(2,2) = -1.0

      eCoeffs(3,1) =  3.0
      eCoeffs(3,2) = -3.0
      eCoeffs(3,3) =  1.0

      eCoeffs(4,1) =  4.0
      eCoeffs(4,2) = -6.0
      eCoeffs(4,3) =  4.0
      eCoeffs(4,4) = -1.0

      eCoeffs(5,1) =  5.0
      eCoeffs(5,2) = -10.0
      eCoeffs(5,3) =  10.0
      eCoeffs(5,4) = -5.0
      eCoeffs(5,5) =  1.0

      eCoeffs(6,1) =  6.0
      eCoeffs(6,2) = -15.0
      eCoeffs(6,3) =  20.0
      eCoeffs(6,4) = -15.0
      eCoeffs(6,5) =  6.0
      eCoeffs(6,6) = -1.0
c
      alpha      = ic_phase_param(1)
      beta       = ic_phase_param(2)
      vx0        = ic_phase_param(3)
      vy0        = ic_phase_param(4)
      vflowinitx = ic_phase_param(5)
      vflowinity = ic_phase_param(6)
c
      vxlo = xlo(3)
      vylo = xlo(4)
      vxhi = xhi(3)
      vyhi = xhi(4)
      dvx = dx(3)
      dvy = dx(4)
c
c.. set lower i3 boundary
      ! check if I am the first interior pt (using a tolerance sacaled by dx)
      i3 = n3a
      vx = vxlo+(0.5+i3-0.5)*dvx
c      write(6,*) vx,vx-vxlo,dvx
      if( abs(vx-vxlo) .lt. 0.5*dvx ) then
c        write(6,*) 'vxlo'
        do i4 = nl4a,nl4b
        do i2 = nl2a,nl2b
        do i1 = nl1a,nl1b
          ! check if outflow or inflow
          if( vel3(i3,i4,i1,i2) .gt. 0.0 ) then
            ! inflow so sample distribution fcn
            ! loop over all ghosts
            do ig = 1,nGhost
              vx = vxlo+(0.5+i3-ig)*dvx
              vy = vylo+(0.5+i4)*dvy
              fx = exp(-0.5*((alpha*(vx-vx0-vflowinitx))**2))
              fy = exp(-0.5*((beta*(vy-vy0-vflowinity))**2))
              fsamp = fx*fy
              fsamp = alpha*beta/(2.0*pi)*fsamp
              u(i1,i2,i3-ig,i4) = fsamp
            end do
          else
            ! outflow so extrapolate
            ! determine the extrapolation condition so that we do
            ! not accidentally reach into ghosts or parallel ghosts
            extrapEq = min( n3b-n3a+1,solution_order )
            !write(6,*) extrapEq
            do ig = 1,nGhost
              u(i1,i2,i3-ig,i4) = 0.0
              do iSten = 1,extrapEq
                u(i1,i2,i3-ig,i4) = u(i1,i2,i3-ig,i4)
     *                +eCoeffs(extrapEq,iSten)*
     *               u(i1,i2,i3-ig+iSten,i4)
              end do
            end do
          end if
        end do
        end do
        end do
      end if
c
c.. set upper i3 boundary
      ! check if I am the last interior pt (using a tolerance sacaled by dx)
      i3 = n3b
      vx = vxlo+(0.5+i3+0.5)*dvx
      if( abs(vx-vxhi) .lt. 0.5*dvx ) then
c        write(6,*) 'vxhi'
        do i4 = nl4a,nl4b
        do i2 = nl2a,nl2b
        do i1 = nl1a,nl1b
          ! check if outflow or inflow
          if( vel3(i3+1,i4,i1,i2) .lt. 0.0 ) then
            ! inflow so sample distribution fcn
            ! loop over all ghosts
            do ig = 1,nGhost
              vx = vxlo+(0.5+i3+ig)*dvx
              vy = vylo+(0.5+i4)*dvy
              fx = exp(-0.5*((alpha*(vx-vx0-vflowinitx))**2))
              fy = exp(-0.5*((beta*(vy-vy0-vflowinity))**2))
              fsamp = fx*fy
              fsamp = alpha*beta/(2.0*pi)*fsamp
              u(i1,i2,i3+ig,i4) = fsamp
            end do
          else
            ! outflow so extrapolate
            ! determine the extrapolation condition so that we do
            ! not accidentally reach into ghosts or parallel ghosts
            extrapEq = min( n3b-n3a+1,solution_order )
            do ig = 1,nGhost
              u(i1,i2,i3+ig,i4) = 0.0
              do iSten = 1,extrapEq
                u(i1,i2,i3+ig,i4) = u(i1,i2,i3+ig,i4)
     *               +eCoeffs(extrapEq,iSten)*
     *               u(i1,i2,i3+ig-iSten,i4)
              end do
            end do
          end if
        end do
        end do
        end do
      end if
c
c.. set lower i4 boundary
      ! check if I am the first interior pt (using a tolerance sacaled by dx)
      i4 = n4a
      vy = vylo+(0.5+i4-0.5)*dvy
      if( abs(vy-vylo) .lt. 0.5*dvy ) then
c        write(6,*) 'vylo'
        do i3 = nl3a,nl3b
        do i2 = nl2a,nl2b
        do i1 = nl1a,nl1b
          ! check if outflow or inflow
          if( vel4(i4,i1,i2,i3) .gt. 0.0 ) then
            ! inflow so sample distribution fcn
            ! loop over all ghosts
            do ig = 1,nGhost
              vx = vxlo+(0.5+i3)*dvx
              vy = vylo+(0.5+i4-ig)*dvy
              fx = exp(-0.5*((alpha*(vx-vx0-vflowinitx))**2))
              fy = exp(-0.5*((beta*(vy-vy0-vflowinity))**2))
              fsamp = fx*fy
              fsamp = alpha*beta/(2.0*pi)*fsamp
              u(i1,i2,i3,i4-ig) = fsamp
            end do
          else
            ! outflow so extrapolate
            ! determine the extrapolation condition so that we do
            ! not accidentally reach into ghosts or parallel ghosts
            extrapEq = min( n4b-n4a+1,solution_order )
            !write(6,*) extrapEq
            do ig = 1,nGhost
              u(i1,i2,i3,i4-ig) = 0.0
              do iSten = 1,extrapEq
                u(i1,i2,i3,i4-ig) = u(i1,i2,i3,i4-ig)
     *               +eCoeffs(extrapEq,iSten)*
     *               u(i1,i2,i3,i4-ig+iSten)
              end do
            end do
          end if
        end do
        end do
        end do
      end if
c
c.. set upper i4 boundary
      ! check if I am the last interior pt (using a tolerance sacaled by dx)
      i4 = n4b
      vy = vylo+(0.5+i4+0.5)*dvy
      if( abs(vy-vyhi) .lt. 0.5*dvy ) then
c        write(6,*) 'vyhi'
        do i3 = nl3a,nl3b
        do i2 = nl2a,nl2b
        do i1 = nl1a,nl1b
          ! check if outflow or inflow
          if( vel4(i4+1,i1,i2,i3) .lt. 0.0 ) then
            ! inflow so sample distribution fcn
            ! loop over all ghosts
            do ig = 1,nGhost
              vx = vxlo+(0.5+i3)*dvx
              vy = vylo+(0.5+i4+ig)*dvy
              fx = exp(-0.5*((alpha*(vx-vx0-vflowinitx))**2))
              fy = exp(-0.5*((beta*(vy-vy0-vflowinity))**2))
              fsamp = fx*fy
              fsamp = alpha*beta/(2.0*pi)*fsamp
              u(i1,i2,i3,i4+ig) = fsamp
            end do
          else
            ! outflow so extrapolate
            ! determine the extrapolation condition so that we do
            ! not accidentally reach into ghosts or parallel ghosts
            extrapEq = min( n4b-n4a+1,solution_order )
            do ig = 1,nGhost
              u(i1,i2,i3,i4+ig) = 0.0
              do iSten = 1,extrapEq
                u(i1,i2,i3,i4+ig) = u(i1,i2,i3,i4+ig)
     *               +eCoeffs(extrapEq,iSten)*
     *               u(i1,i2,i3,i4+ig-iSten)
              end do
            end do
          end if
        end do
        end do
        end do
      end if
c
      return
      end
c
c +++++++++++++
c
      subroutine setAdvectionBCs4DJB(
     &     u,
     &     ng1a,ng1b,ng2a,ng2b,ng3a,ng3b,ng4a,ng4b,
     &     nl1a,nl1b,nl2a,nl2b,nl3a,nl3b,nl4a,nl4b,
     &     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     &     solution_order, nGhost,
     &     vel1,vel2,
     &     xloV, xhiV, dxV,
     &     ic_phase_param,
     &     xPeriodic, yPeriodic )
c
c.. boundary conditions 
c
c.. declarations of incoming variables
      implicit none
      integer ng1a,ng1b,ng2a,ng2b
      integer ng3a,ng3b,ng4a,ng4b
      integer nl1a,nl1b,nl2a,nl2b
      integer nl3a,nl3b,nl4a,nl4b
      integer n1a,n1b,n2a,n2b
      integer n3a,n3b,n4a,n4b
      integer solution_order, nGhost
      integer xPeriodic,yPeriodic
      real u(nl1a:nl1b,nl2a:nl2b,nl3a:nl3b,nl4a:nl4b)
      real vel1( nl1a:nl1b+1,nl2a:nl2b,nl3a:nl3b,nl4a:nl4b )
      real vel2( nl2a:nl2b+1,nl3a:nl3b,nl4a:nl4b,nl1a:nl1b )
      real xloV(1:4), xhiV(1:4), dxV(1:4)
      real ic_phase_param(*)
c
c.. declrations of local variables
      integer i1,i2,i3,i4
      integer ig,extrapEq,iSten
      real eCoeffs(1:6,1:6)
      real alpha,beta,vx0,vy0,vflowinitx,vflowinity
      real xlo,xhi,ylo,yhi
      real dx,dy
      real vxlo,vylo,dvx,dvy
      real fsamp
      real x,y,vx,vy
      real pi
c
      real fx
      real fy
c
      real one,four
c
      one = 1.0
      four = 4.0
c
      pi  = four*atan(one)
c
c fill in matrix of extrapolation coefficients (formula,coefficient)
      eCoeffs(1,1) = 1.0

      eCoeffs(2,1) =  2.0
      eCoeffs(2,2) = -1.0

      eCoeffs(3,1) =  3.0
      eCoeffs(3,2) = -3.0
      eCoeffs(3,3) =  1.0

      eCoeffs(4,1) =  4.0
      eCoeffs(4,2) = -6.0
      eCoeffs(4,3) =  4.0
      eCoeffs(4,4) = -1.0

      eCoeffs(5,1) =  5.0
      eCoeffs(5,2) = -10.0
      eCoeffs(5,3) =  10.0
      eCoeffs(5,4) = -5.0
      eCoeffs(5,5) =  1.0

      eCoeffs(6,1) =  6.0
      eCoeffs(6,2) = -15.0
      eCoeffs(6,3) =  20.0
      eCoeffs(6,4) = -15.0
      eCoeffs(6,5) =  6.0
      eCoeffs(6,6) = -1.0
c
      alpha      = ic_phase_param(1)
      beta       = ic_phase_param(2)
      vx0        = ic_phase_param(3)
      vy0        = ic_phase_param(4)
      vflowinitx = ic_phase_param(5)
      vflowinity = ic_phase_param(6)

      xlo  = xloV(1)
      ylo  = xloV(2)
      vxlo = xloV(3)
      vylo = xloV(4)
      
      xhi  = xhiV(1)
      yhi  = xhiV(2)
      
      dx  = dxV(1)
      dy  = dxV(2)
      dvx = dxV(3)
      dvy = dxV(4)
c
c.. set lower i1 boundary
      i1 = n1a
      x  = xlo+(0.5+i1+0.5)*dx
      if( abs(x-xlo) .lt. 0.5*dx .and. xPeriodic .ne. 1 ) then
        do i4 = nl4a,nl4b
        do i3 = nl3a,nl3b
        do i2 = nl2a,nl2b
          ! check if outflow or inflow
          if( vel1(i1,i2,i3,i4) .lt. 0.0 ) then
            ! inflow so sample distribution fcn
            ! loop over all ghosts
            do ig = 1,nGhost
              vx = vxlo+(0.5+i3)*dvx
              vy = vylo+(0.5+i4)*dvy
              fx = exp(-0.5*((alpha*(vx-vx0-vflowinitx))**2))
              fy = exp(-0.5*((beta*(vy-vy0-vflowinity))**2))
              fsamp = fx*fy
              fsamp = alpha*beta/(2.0*pi)*fsamp
              u(i1-ig,i2,i3,i4) = fsamp
            end do
          else
            ! outflow so extrapolate
            ! determine the extrapolation condition so that we do
            ! not accidentally reach into ghosts or parallel ghosts
            extrapEq = min( n1b-n1a+1,solution_order )
            do ig = 1,nGhost
              u(i1-ig,i2,i3,i4) = 0.0
              do iSten = 1,extrapEq
                u(i1-ig,i2,i3,i4) = u(i1-ig,i2,i3,i4)
     *               +eCoeffs(extrapEq,iSten)*
     *               u(i1-ig+iSten,i2,i3,i4)
              end do
            end do
          end if
        end do
        end do
        end do
      end if
c     
c.. set upper i1 boundary
      i1 = n1b
      x  = xlo+(0.5+i1+0.5)*dx
      if( abs(x-xhi) .lt. 0.5*dx .and. xPeriodic .ne. 1 ) then
        do i4 = nl4a,nl4b
        do i3 = nl3a,nl3b
        do i2 = nl2a,nl2b
          ! check if outflow or inflow
          if( vel1(i1+1,i2,i3,i4) .lt. 0.0 ) then
            ! inflow so sample distribution fcn
            ! loop over all ghosts
            do ig = 1,nGhost
              vx = vxlo+(0.5+i3)*dvx
              vy = vylo+(0.5+i4)*dvy
              fx = exp(-0.5*((alpha*(vx-vx0-vflowinitx))**2))
              fy = exp(-0.5*((beta*(vy-vy0-vflowinity))**2))
              fsamp = fx*fy
              fsamp = alpha*beta/(2.0*pi)*fsamp
              u(i1+ig,i2,i3,i4) = fsamp
            end do
          else
            ! outflow so extrapolate
            ! determine the extrapolation condition so that we do
            ! not accidentally reach into ghosts or parallel ghosts
            extrapEq = min( n1b-n1a+1,solution_order )
            do ig = 1,nGhost
              u(i1+ig,i2,i3,i4) = 0.0
              do iSten = 1,extrapEq
                u(i1+ig,i2,i3,i4) = u(i1+ig,i2,i3,i4)
     *               +eCoeffs(extrapEq,iSten)*
     *               u(i1+ig-iSten,i2,i3,i4)
              end do
            end do
          end if
        end do
        end do
        end do
      end if
c
c.. set lower i2 boundary
      i2 = n2a
      y  = ylo+(0.5+i2+0.5)*dy
      if( abs(y-ylo) .lt. 0.5*dy .and. yPeriodic .ne. 1 ) then
        do i4 = nl4a,nl4b
        do i3 = nl3a,nl3b
        do i1 = nl1a,nl1b
          ! check if outflow or inflow
          if( vel2(i2,i3,i4,i1) .lt. 0.0 ) then
            ! inflow so sample distribution fcn
            ! loop over all ghosts
            do ig = 1,nGhost
              vx = vxlo+(0.5+i3)*dvx
              vy = vylo+(0.5+i4)*dvy
              fx = exp(-0.5*((alpha*(vx-vx0-vflowinitx))**2))
              fy = exp(-0.5*((beta*(vy-vy0-vflowinity))**2))
              fsamp = fx*fy
              fsamp = alpha*beta/(2.0*pi)*fsamp
              u(i1,i2-ig,i3,i4) = fsamp
            end do
          else
            ! outflow so extrapolate
            ! determine the extrapolation condition so that we do
            ! not accidentally reach into ghosts or parallel ghosts
            extrapEq = min( n2b-n2a+1,solution_order )
            do ig = 1,nGhost
              u(i1,i2-ig,i3,i4) = 0.0
              do iSten = 1,extrapEq
                u(i1,i2-ig,i3,i4) = u(i1,i2-ig,i3,i4)
     *               +eCoeffs(extrapEq,iSten)*
     *               u(i1,i2-ig+iSten,i3,i4)
              end do
            end do
          end if
        end do
        end do
        end do
      end if
c     
c.. set upper i2 boundary
      i2 = n2b
      y  = ylo+(0.5+i2+0.5)*dy
      if( abs(y-yhi) .lt. 0.5*dy .and. yPeriodic .ne. 1 ) then
        do i4 = nl4a,nl4b
        do i3 = nl3a,nl3b
        do i1 = nl1a,nl1b
          ! check if outflow or inflow
          if( vel2(i2+1,i3,i4,i1) .lt. 0.0 ) then
            ! inflow so sample distribution fcn
            ! loop over all ghosts
            do ig = 1,nGhost
              vx = vxlo+(0.5+i3)*dvx
              vy = vylo+(0.5+i4)*dvy
              fx = exp(-0.5*((alpha*(vx-vx0-vflowinitx))**2))
              fy = exp(-0.5*((beta*(vy-vy0-vflowinity))**2))
              fsamp = fx*fy
              fsamp = alpha*beta/(2.0*pi)*fsamp
              u(i1,i2+ig,i3,i4) = fsamp
            end do
          else
            ! outflow so extrapolate
            ! determine the extrapolation condition so that we do
            ! not accidentally reach into ghosts or parallel ghosts
            extrapEq = min( n2b-n2a+1,solution_order )
            do ig = 1,nGhost
              u(i1,i2+ig,i3,i4) = 0.0
              do iSten = 1,extrapEq
                u(i1,i2+ig,i3,i4) = u(i1,i2+ig,i3,i4)
     *               +eCoeffs(extrapEq,iSten)*
     *               u(i1,i2+ig-iSten,i3,i4)
              end do
            end do
          end if
        end do
        end do
        end do
      end if
c
      return
      end
c
c +++++++++++++
c
      subroutine zero_Ghost4D(
     & u,
     & n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     & nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b )
c
      implicit none
      integer nd1a,nd1b,nd2a,nd2b
      integer nd3a,nd3b,nd4a,nd4b
      integer n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b
      real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      integer i1,i2,i3,i4
c
c .. i1 left
      do i4=nd4a,nd4b
      do i3=nd3a,nd3b
      do i2=nd2a,nd2b
      do i1=nd1a,n1a-1
        u(i1,i2,i3,i4) = 0.0
      end do
      end do
      end do
      end do
c.. i1 right
      do i4=nd4a,nd4b
      do i3=nd3a,nd3b
      do i2=nd2a,nd2b
      do i1=n1b+1,nd1b
        u(i1,i2,i3,i4) = 0.0
      end do
      end do
      end do
      end do
c
c .. i2 left
      do i4=nd4a,nd4b
      do i3=nd3a,nd3b
      do i2=nd2a,n2a-1
      do i1=nd1a,nd1b
        u(i1,i2,i3,i4) = 0.0
      end do
      end do
      end do
      end do
c.. i2 right
      do i4=nd4a,nd4b
      do i3=nd3a,nd3b
      do i2=n2b+1,nd2b
      do i1=nd1a,nd1b
        u(i1,i2,i3,i4) = 0.0
      end do
      end do
      end do
      end do
c
c .. i3 left
      do i4=nd4a,nd4b
      do i3=nd3a,n3a-1
      do i2=nd2a,nd2b
      do i1=nd1a,nd1b
        u(i1,i2,i3,i4) = 0.0
      end do
      end do
      end do
      end do
c.. i3 right
      do i4=nd4a,nd4b
      do i3=n3b+1,nd3b
      do i2=nd2a,nd2b
      do i1=nd1a,nd1b
        u(i1,i2,i3,i4) = 0.0
      end do
      end do
      end do
      end do
c
c .. i4 left
      do i4=nd4a,n4a-1
      do i3=nd3a,nd3b
      do i2=nd2a,nd2b
      do i1=nd1a,nd1a
        u(i1,i2,i3,i4) = 0.0
      end do
      end do
      end do
      end do
c.. i4 right
      do i4=n4b+1,nd4b
      do i3=nd3a,nd3b
      do i2=nd2a,nd2b
      do i1=nd1a,nd1b
        u(i1,i2,i3,i4) = 0.0
      end do
      end do
      end do
      end do

      return
      end
c
c +++++++++++++
c
      subroutine computeadvectionfluxes4D( 
     &     flux1,flux2,
     &     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     &     vel1,vel2,
     &     face1,face2,
     &     u,
     &     dx,
     &     solution_order)
c
c.. evaluate just the advection fluxes 
      implicit none
c
c.. declarations of incoming variables 
      integer nd1a,nd1b,nd2a,nd2b
      integer nd3a,nd3b,nd4a,nd4b
      integer solution_order
      real u( nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b )

      real vel1( nd1a:nd1b+1,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b )
      real vel2( nd2a:nd2b+1,nd3a:nd3b,nd4a:nd4b,nd1a:nd1b )

      real face1( nd1a:nd1b+1,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b )
      real face2( nd2a:nd2b+1,nd3a:nd3b,nd4a:nd4b,nd1a:nd1b )

      real flux1( nd1a:nd1b+1,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b )
      real flux2( nd2a:nd2b+1,nd3a:nd3b,nd4a:nd4b,nd1a:nd1b )

      real dx(1:4)
c
c.. compute face averages 
      if( .true. ) then
        if( solution_order .eq. 4 ) then
          ! WENO, as in 4/3 that I cooked up
          call WENO43Avg4D( nd1a,nd1b,nd2a,nd2b,
     *                      nd3a,nd3b,nd4a,nd4b,
     *                      nd1a,nd1b+1,nd2a,nd2b,
     *                      nd3a,nd3b,nd4a,nd4b,
     *                      1,u,vel1,face1 )
          call WENO43Avg4D( nd1a,nd1b,nd2a,nd2b,
     *                      nd3a,nd3b,nd4a,nd4b,
     *                      nd2a,nd2b+1,nd3a,nd3b,
     *                      nd4a,nd4b,nd1a,nd1b,
     *                      2,u,vel2,face2 )
        else
          ! BWENO, as in 6/5 that I cooked up
          call WENO65Avg4D( nd1a,nd1b,nd2a,nd2b,
     *                      nd3a,nd3b,nd4a,nd4b,
     *                      nd1a,nd1b+1,nd2a,nd2b,
     *                      nd3a,nd3b,nd4a,nd4b,
     *                      1,u,vel1,face1 )
          call WENO65Avg4D( nd1a,nd1b,nd2a,nd2b,
     *                      nd3a,nd3b,nd4a,nd4b,
     *                      nd2a,nd2b+1,nd3a,nd3b,
     *                      nd4a,nd4b,nd1a,nd1b,
     *                      2,u,vel2,face2 )
        end if
      else
        ! unlimited 4th order fluxes
        call centAvg_to_faceAvg4D( nd1a,nd1b,nd2a,nd2b,
     *                             nd3a,nd3b,nd4a,nd4b,
     *                             nd1a,nd1b+1,nd2a,nd2b,
     *                             nd3a,nd3b,nd4a,nd4b,
     *                             1,u,face1 )
        call centAvg_to_faceAvg4D( nd1a,nd1b,nd2a,nd2b,
     *                             nd3a,nd3b,nd4a,nd4b,
     *                             nd2a,nd2b+1,nd3a,nd3b,
     *                             nd4a,nd4b,nd1a,nd1b,
     *                             2,u,face2 )
      end if
c
c.. ppm limiter
      if( .false. ) then
        call SKLimit4D( nd1a,nd1b,nd2a,nd2b,
     *                  nd3a,nd3b,nd4a,nd4b,
     *                  nd1a,nd1b+1,nd2a,nd2b,
     *                  nd3a,nd3b,nd4a,nd4b,
     *                  1,u,face1,vel1,flux1 )
        call SKLimit4D( nd1a,nd1b,nd2a,nd2b,
     *                  nd3a,nd3b,nd4a,nd4b,
     *                  nd2a,nd2b+1,nd3a,nd3b,
     *                  nd4a,nd4b,nd1a,nd1b,
     *                  2,u,face2,vel2,flux2 )

      end if
c
c.. compute fluxes
      call computeFlux4D( nd1a,nd1b+1,nd2a,nd2b,
     *                    nd3a,nd3b,nd4a,nd4b,
     *                    face1,vel1,flux1 )
      call computeFlux4D( nd2a,nd2b+1,nd3a,nd3b,
     *                    nd4a,nd4b,nd1a,nd1b,
     *                    face2,vel2,flux2 )
c
      if( .false. ) then
        call artVisFlux4D( nd1a,nd1b,nd2a,nd2b,
     *                     nd3a,nd3b,nd4a,nd4b,
     *                     nd1a,nd1b+1,nd2a,nd2b,
     *                     nd3a,nd3b,nd4a,nd4b,
     *                     u,flux1,dx(1),1 )
        call artVisFlux4D( nd1a,nd1b,nd2a,nd2b,
     *                     nd3a,nd3b,nd4a,nd4b,
     *                     nd2a,nd2b+1,nd3a,nd3b,
     *                     nd4a,nd4b,nd1a,nd1b,
     *                     u,flux2,dx(2),2 )
      end if
c
      return
      end
c
c +++++++++++++
c
      subroutine computeadvectionderivatives4D( 
     &     rhs,f,
     &     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     &     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     &     vel1,vel2,
     &     deltax,
     &     solution_order)
c
c.. evaluate just the advection derivatives
      implicit none
c
c.. declarations of incoming variables 
      integer nd1a,nd1b,nd2a,nd2b
      integer nd3a,nd3b,nd4a,nd4b
      integer n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b
      integer solution_order
      real f(    nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b )
      real rhs(  nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b )

      real vel1( nd1a:nd1b+1,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b )
      real vel2( nd2a:nd2b+1,nd3a:nd3b,nd4a:nd4b,nd1a:nd1b )

      real deltax(1:4)
c
c.. declarations of local variables
      integer i1,i2,i3,i4
      real uRightx,uLeftx
      real uRighty,uLefty
      real vx,vy
      real dx,dy
c
c.. compute x and y derivatives
      dx = deltax(1)
      dy = deltax(2)
      do i4 = n4a,n4b
        vy = vel2(n2a,n3a,i4,n1a)

        do i3 = n3a,n3b
          vx = vel1(n1a,n2a,i3,n4a)

          if( solution_order .eq. 4 ) then
            do i2 = n2a,n2b
              i1 = n1a
              call WENO43Fit4D( f(i1-2,i2,i3,i4),f(i1-1,i2,i3,i4),
     *             f(i1,i2,i3,i4),f(i1+1,i2,i3,i4),
     *             uLeftx,vx )

              do i1 = n1a,n1b
                call WENO43Fit4D( f(i1-1,i2,i3,i4),f(i1,i2,i3,i4),
     *                            f(i1+1,i2,i3,i4),f(i1+2,i2,i3,i4),
     *                            uRightx,vx )

                rhs(i1,i2,i3,i4) = 
c     *                -vx*(uRightx-uLeftx)/dx
     *                -(vx*uRightx-vx*uLeftx)/dx
                uLeftx = uRightx

              end do
            end do
            do i1 = n1a,n1b
              i2 = n2a
              call WENO43Fit4D( f(i1,i2-2,i3,i4),f(i1,i2-1,i3,i4),
     *                         f(i1,i2,i3,i4),f(i1,i2+1,i3,i4),
     *                         uLefty,vy )

              do i2 = n2a,n2b

                call WENO43Fit4D( f(i1,i2-1,i3,i4),f(i1,i2,i3,i4),
     *                           f(i1,i2+1,i3,i4),f(i1,i2+2,i3,i4),
     *                           uRighty,vy )
               
                rhs(i1,i2,i3,i4) = rhs(i1,i2,i3,i4)
c     *               -vy*(uRighty-uLefty)/dy
     *               -(vy*uRighty-vy*uLefty)/dy
                uLefty = uRighty

              end do
            end do
          else ! 6th order
            do i2 = n2a,n2b
              i1 = n1a
              call WENO65Fit4D( 
     *             f(i1-3,i2,i3,i4),
     *             f(i1-2,i2,i3,i4),
     *             f(i1-1,i2,i3,i4),
     *             f(i1,i2,i3,i4),
     *             f(i1+1,i2,i3,i4),
     *             f(i1+2,i2,i3,i4),
     *             uLeftx,vx )

              do i1 = n1a,n1b
                call WENO65Fit4D( 
     *               f(i1-2,i2,i3,i4),
     *               f(i1-1,i2,i3,i4),
     *               f(i1,i2,i3,i4),
     *               f(i1+1,i2,i3,i4),
     *               f(i1+2,i2,i3,i4),
     *               f(i1+3,i2,i3,i4),
     *               uRightx,vx )

                rhs(i1,i2,i3,i4) = 
c     *                -vx*(uRightx-uLeftx)/dx
     *               -(vx*uRightx-vx*uLeftx)/dx
                uLeftx = uRightx

              end do
            end do
            do i1 = n1a,n1b
              i2 = n2a
              call WENO65Fit4D( 
     *             f(i1,i2-3,i3,i4),
     *             f(i1,i2-2,i3,i4),
     *             f(i1,i2-1,i3,i4),
     *             f(i1,i2,i3,i4),
     *             f(i1,i2+1,i3,i4),
     *             f(i1,i2+2,i3,i4),
     *             uLefty,vy )

              do i2 = n2a,n2b
                call WENO65Fit4D( 
     *               f(i1,i2-2,i3,i4),
     *               f(i1,i2-1,i3,i4),
     *               f(i1,i2,i3,i4),
     *               f(i1,i2+1,i3,i4),
     *               f(i1,i2+2,i3,i4),
     *               f(i1,i2+3,i3,i4),
     *               uRighty,vy )
               
                rhs(i1,i2,i3,i4) = rhs(i1,i2,i3,i4)
c     *               -vy*(uRighty-uLefty)/dy
     *               -(vy*uRighty-vy*uLefty)/dy
                uLefty = uRighty

              end do
            end do
          end if
         
        end do
      end do
c
      return
      end
c
c ++++++++++++++
c
      subroutine computeaccelerationderivatives4D( 
     &     rhs,f,
     &     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     &     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     &     vel3,vel4,
     &     dx,
     &     solution_order)
c
c.. evaluate just the acceleration derivatives
      implicit none
c
c.. declarations of incoming variables 
      integer nd1a,nd1b,nd2a,nd2b
      integer nd3a,nd3b,nd4a,nd4b
      integer n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b
      integer solution_order
      real f( nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b )
      real rhs ( nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b )

      real vel3( nd3a:nd3b+1,nd4a:nd4b,nd1a:nd1b,nd2a:nd2b )
      real vel4( nd4a:nd4b+1,nd1a:nd1b,nd2a:nd2b,nd3a:nd3b )

      real dx(1:4)
c
c.. declarations of local variables
      integer i1,i2,i3,i4
      real uRightvx,uLeftvx
      real uRightvy,uLeftvy
      real ax,ay
      real dvx,dvy
c
c.. compute vx and vy derivatives
      dvx = dx(3)
      dvy = dx(4)
      do i2 = n2a,n2b
c        ay = vel4(n4a,n1a,i2,n3a)

        do i1 = n1a,n1b
cc          ax = vel3(n3a,n4a,i1,n2a)
c          ax = vel3(n3a,n4a,i1,i2)
c          ay = vel4(n4a,i1,i2,n3a)

          if( solution_order .eq. 4 ) then
            do i4 = n4a,n4b
              i3 = n3a
              ax = vel3(i3,i4,i1,i2)
              call WENO43Fit4D( f(i1,i2,i3-2,i4),f(i1,i2,i3-1,i4),
     *                          f(i1,i2,i3,i4),f(i1,i2,i3+1,i4),
     *                          uLeftvx,ax )

              do i3 = n3a,n3b
                ax = vel3(i3,i4,i1,i2)
                call WENO43Fit4D( f(i1,i2,i3-1,i4),f(i1,i2,i3,i4),
     *                            f(i1,i2,i3+1,i4),f(i1,i2,i3+2,i4),
     *                            uRightvx,ax )

                rhs(i1,i2,i3,i4) = rhs(i1,i2,i3,i4)
c     *               -ax*(uRightvx-uLeftvx)/dvx
     *               -(ax*uRightvx-ax*uLeftvx)/dvx
                uLeftvx = uRightvx
              end do
            end do

            do i3 = n3a,n3b
              i4 = n4a
              ay = vel4(i4,i1,i2,i3)
              call WENO43Fit4D( f(i1,i2,i3,i4-2),f(i1,i2,i3,i4-1),
     *                          f(i1,i2,i3,i4),f(i1,i2,i3,i4+1),
     *                          uLeftvy,ay )

              do i4 = n4a,n4b
                ay = vel4(i4,i1,i2,i3)
                call WENO43Fit4D( f(i1,i2,i3,i4-1),f(i1,i2,i3,i4),
     *                            f(i1,i2,i3,i4+1),f(i1,i2,i3,i4+2),
     *                            uRightvy,ay )
           
                rhs(i1,i2,i3,i4) = rhs(i1,i2,i3,i4)
c     *               -ay*(uRightvy-uLeftvy)/dvy
     *               -(ay*uRightvy-ay*uLeftvy)/dvy
                uLeftvy = uRightvy
              end do
            end do
          else ! 6th order 
            do i4 = n4a,n4b
              i3 = n3a
              ax = vel3(i3,i4,i1,i2)
              call WENO65Fit4D( 
     *             f(i1,i2,i3-3,i4),
     *             f(i1,i2,i3-2,i4),
     *             f(i1,i2,i3-1,i4),
     *             f(i1,i2,i3,i4),
     *             f(i1,i2,i3+1,i4),
     *             f(i1,i2,i3+2,i4),
     *             uLeftvx,ax )

              do i3 = n3a,n3b
                ax = vel3(i3,i4,i1,i2)
                call WENO65Fit4D( 
     *               f(i1,i2,i3-2,i4),
     *               f(i1,i2,i3-1,i4),
     *               f(i1,i2,i3,i4),
     *               f(i1,i2,i3+1,i4),
     *               f(i1,i2,i3+2,i4),
     *               f(i1,i2,i3+3,i4),
     *               uRightvx,ax )

                rhs(i1,i2,i3,i4) = rhs(i1,i2,i3,i4)
c     *               -ax*(uRightvx-uLeftvx)/dvx
     *               -(ax*uRightvx-ax*uLeftvx)/dvx
                uLeftvx = uRightvx
              end do
            end do

            do i3 = n3a,n3b
              i4 = n4a
              ay = vel4(i4,i1,i2,i3)
              call WENO65Fit4D( 
     *             f(i1,i2,i3,i4-3),
     *             f(i1,i2,i3,i4-2),
     *             f(i1,i2,i3,i4-1),
     *             f(i1,i2,i3,i4),
     *             f(i1,i2,i3,i4+1),
     *             f(i1,i2,i3,i4+2),
     *             uLeftvy,ay )

              do i4 = n4a,n4b
                ay = vel4(i4,i1,i2,i3)
                call WENO65Fit4D( 
     *               f(i1,i2,i3,i4-2),
     *               f(i1,i2,i3,i4-1),
     *               f(i1,i2,i3,i4),
     *               f(i1,i2,i3,i4+1),
     *               f(i1,i2,i3,i4+2),
     *               f(i1,i2,i3,i4+3),
     *               uRightvy,ay )

                rhs(i1,i2,i3,i4) = rhs(i1,i2,i3,i4)
c     *               -ay*(uRightvy-uLeftvy)/dvy
     *               -(ay*uRightvy-ay*uLeftvy)/dvy
                uLeftvy = uRightvy
              end do
            end do
          end if

        end do
      end do
c     
      return
      end
c
c ++++++++++++++
c
      subroutine computeaccelerationfluxes4D( 
     &     flux3,flux4,
     &     nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,
     &     vel3,vel4,
     &     face3,face4,
     &     u,
     &     dx,
     &     solution_order)
c
c.. evaluate just the acceleration fluxes
      implicit none
c
c.. declarations of incoming variables 
      integer nd1a,nd1b,nd2a,nd2b
      integer nd3a,nd3b,nd4a,nd4b
      integer solution_order
      real u( nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b )

      real vel3( nd3a:nd3b+1,nd4a:nd4b,nd1a:nd1b,nd2a:nd2b )
      real vel4( nd4a:nd4b+1,nd1a:nd1b,nd2a:nd2b,nd3a:nd3b )

      real face3( nd3a:nd3b+1,nd4a:nd4b,nd1a:nd1b,nd2a:nd2b )
      real face4( nd4a:nd4b+1,nd1a:nd1b,nd2a:nd2b,nd3a:nd3b )

      real flux3( nd3a:nd3b+1,nd4a:nd4b,nd1a:nd1b,nd2a:nd2b )
      real flux4( nd4a:nd4b+1,nd1a:nd1b,nd2a:nd2b,nd3a:nd3b )

      real dx(1:4)
c
c.. compute face averages 
      if( .true. ) then
        if( solution_order .eq. 4 ) then
          ! WENO, as in 4/3 that I cooked up
          call WENO43Avg4D( nd1a,nd1b,nd2a,nd2b,
     *                      nd3a,nd3b,nd4a,nd4b,
     *                      nd3a,nd3b+1,nd4a,nd4b,
     *                      nd1a,nd1b,nd2a,nd2b,
     *                      3,u,vel3,face3 )
          call WENO43Avg4D( nd1a,nd1b,nd2a,nd2b,
     *                      nd3a,nd3b,nd4a,nd4b,
     *                      nd4a,nd4b+1,nd1a,nd1b,
     *                      nd2a,nd2b,nd3a,nd3b,
     *                      4,u,vel4,face4 )
        else
          ! BWENO, as in 6/5 that I cooked up
          call WENO65Avg4D( nd1a,nd1b,nd2a,nd2b,
     *                      nd3a,nd3b,nd4a,nd4b,
     *                      nd3a,nd3b+1,nd4a,nd4b,
     *                      nd1a,nd1b,nd2a,nd2b,
     *                      3,u,vel3,face3 )
          call WENO65Avg4D( nd1a,nd1b,nd2a,nd2b,
     *                      nd3a,nd3b,nd4a,nd4b,
     *                      nd4a,nd4b+1,nd1a,nd1b,
     *                      nd2a,nd2b,nd3a,nd3b,
     *                      4,u,vel4,face4 )
        end if
      else
        ! unlimited 4th order fluxes
        call centAvg_to_faceAvg4D( nd1a,nd1b,nd2a,nd2b,
     *                             nd3a,nd3b,nd4a,nd4b,
     *                             nd3a,nd3b+1,nd4a,nd4b,
     *                             nd1a,nd1b,nd2a,nd2b,
     *                             3,u,face3 )
        call centAvg_to_faceAvg4D( nd1a,nd1b,nd2a,nd2b,
     *                             nd3a,nd3b,nd4a,nd4b,
     *                             nd4a,nd4b+1,nd1a,nd1b,
     *                             nd2a,nd2b,nd3a,nd3b,
     *                             4,u,face4 )
      end if
c
c.. ppm limiter
      if( .false. ) then
        call SKLimit4D( nd1a,nd1b,nd2a,nd2b,
     *                  nd3a,nd3b,nd4a,nd4b,
     *                  nd3a,nd3b+1,nd4a,nd4b,
     *                  nd1a,nd1b,nd2a,nd2b,
     *                  3,u,face3,vel3,flux3 )
        call SKLimit4D( nd1a,nd1b,nd2a,nd2b,
     *                  nd3a,nd3b,nd4a,nd4b,
     *                  nd4a,nd4b+1,nd1a,nd1b,
     *                  nd2a,nd2b,nd3a,nd3b,
     *                  4,u,face4,vel4,flux4 )
      end if
c
c.. compute fluxes
      call computeFlux4D( nd3a,nd3b+1,nd4a,nd4b,
     *                    nd1a,nd1b,nd2a,nd2b,
     *                    face3,vel3,flux3 )
      call computeFlux4D( nd4a,nd4b+1,nd1a,nd1b,
     *                    nd2a,nd2b,nd3a,nd3b,
     *                    face4,vel4,flux4 )
c
      if( .false. ) then
        call artVisFlux4D( nd1a,nd1b,nd2a,nd2b,
     *                     nd3a,nd3b,nd4a,nd4b,
     *                     nd3a,nd3b+1,nd4a,nd4b,
     *                     nd1a,nd1b,nd2a,nd2b,
     *                     u,flux3,dx(3),3 )
        call artVisFlux4D( nd1a,nd1b,nd2a,nd2b,
     *                     nd3a,nd3b,nd4a,nd4b,
     *                     nd4a,nd4b+1,nd1a,nd1b,
     *                     nd2a,nd2b,nd3a,nd3b,
     *                     u,flux4,dx(4),4 )
      end if
c
      return
      end
c
c ++++++++++++++
c
      subroutine computeFlux4D(
     *     nf1a,nf1b,nf2a,nf2b,nf3a,nf3b,nf4a,nf4b,
     *     face,vel,flux )
c
      implicit none
c
c.. declarations of incoming variables
      integer nf1a,nf1b,nf2a,nf2b
      integer nf3a,nf3b,nf4a,nf4b
      real face( nf1a:nf1b,nf2a:nf2b,nf3a:nf3b,nf4a:nf4b )
      real vel( nf1a:nf1b,nf2a:nf2b,nf3a:nf3b,nf4a:nf4b )
      real flux( nf1a:nf1b,nf2a:nf2b,nf3a:nf3b,nf4a:nf4b )
c
c.. declarations of local variables
      integer f1a,f1b,f2a,f2b,f3a,f3b,f4a,f4b
      integer i1,i2,i3,i4,k
c
      f1a = nf1a+2
      f1b = nf1b-2
      f2a = nf2a+2
      f2b = nf2b-2
      f3a = nf3a+2
      f3b = nf3b-2
      f4a = nf4a+2
      f4b = nf4b-2
c
      do i4 = f4a,f4b    
      do i3 = f3a,f3b    
      do i2 = f2a,f2b    
      do i1 = f1a,f1b
        flux(i1,i2,i3,i4) = vel(i1,i2,i3,i4)*face(i1,i2,i3,i4)
      end do
      end do
      end do
      end do
      
      return
      end
c
c **************
c
      subroutine computecurrents(
     &     nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b,
     &     n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b,
     &     xlo, xhi, dx,
     &     u,
     &     vz,
     &     Jx, Jy, Jz)
c
c.. function to obtain local currents
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
      real xlo(1:4), xhi(1:4), dx(1:4)
      real u(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real vz(nd1a:nd1b, nd2a:nd2b, 1:1)
      real Jx(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real Jy(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real Jz(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      real vx, vy, vx0, vy0, dvx, dvy
c
      vx0 = xlo(3)
      vy0 = xlo(4)
      dvx = dx(3)
      dvy = dx(4)
      do i4 = n4a, n4b
        vy = vy0+(0.5+i4)*dvy
      do i3 = n3a, n3b
        vx = vx0+(0.5+i3)*dvx
      do i2 = n2a, n2b
      do i1 = n1a, n1b

        Jx(i1, i2, i3, i4) = u(i1, i2, i3, i4)*vx

        Jy(i1, i2, i3, i4) = u(i1, i2, i3, i4)*vy

        Jz(i1, i2, i3, i4) = u(i1, i2, i3, i4)*vz(i1, i2, 1)
      end do
      end do
      end do
      end do
c
      return
      end
c
c **************
c
      subroutine computeke4d(
     &     nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b,
     &     n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b,
     &     xlo, xhi, dx,
     &     u,
     &     ke)
c
c.. function to obtain local currents
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
      real xlo(1:4), xhi(1:4), dx(1:4)
      real u(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real ke(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      real vx, vy, vx0, vy0, dvx, dvy, vy2, v2
c
      vx0 = xlo(3)
      vy0 = xlo(4)
      dvx = dx(3)
      dvy = dx(4)
      do i4 = n4a, n4b
        vy = vy0+(0.5+i4)*dvy
        vy2 = vy*vy
      do i3 = n3a, n3b
        vx = vx0+(0.5+i3)*dvx
        v2 = vx*vx+vy2
      do i2 = n2a, n2b
      do i1 = n1a, n1b

        ke(i1, i2, i3, i4) = 0.5*u(i1, i2, i3, i4)*v2
      end do
      end do
      end do
      end do
c
      return
      end
c
c **************
c
      subroutine computeke2d(
     &     nd1a, nd1b, nd2a, nd2b,
     &     n1a, n1b, n2a, n2b,
     &     xlo, xhi, dx,
     &     ke_2d,
     &     ke)
c
c.. function to obtain local currents
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b
      integer n1a, n1b, n2a, n2b
      real xlo(1:4), xhi(1:4), dx(1:4)
      real ke_2d(nd1a:nd1b, nd2a:nd2b)
      real ke
c
c.. declarations of local variables
      integer i1, i2
c
      ke = 0.0
      do i2 = n2a, n2b
      do i1 = n1a, n1b

        ke = ke + ke_2d(i1, i2)
      end do
      end do
      ke = ke*dx(1)*dx(2)
c
      return
      end
c
c **************
c
      subroutine computekeflux4d(
     &     nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b,
     &     n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b,
     &     ng1a, ng1b, ng2a, ng2b, ng3a, ng3b, ng4a, ng4b,
     &     xlo, xhi, dx,
     &     face_flux1,
     &     face_flux2,
     &     face_vel1,
     &     face_vel2,
     &     ke_flux,
     &     dir)
c
c.. function to compute energy flux at physical boundary
      implicit none
c
c.. declaration of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
      integer ng1a, ng1b, ng2a, ng2b, ng3a, ng3b, ng4a, ng4b
      real xlo(1:4), xhi(1:4), dx(1:4)
      real face_flux1(nd1a:nd1b+1, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real face_flux2(nd2a:nd2b+1, nd3a:nd3b, nd4a:nd4b, nd1a:nd1b)
      real face_vel1(nd1a:nd1b+1, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real face_vel2(nd2a:nd2b+1, nd3a:nd3b, nd4a:nd4b, nd1a:nd1b)
      real ke_flux(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      integer dir
c
c.. declaration of local variables
      integer i1, i2, i3, i4
      real vx, vy, fluxx, fluxy, vx0, vy0, dvx, dvy, vx2, vy2
c
      vx0 = xlo(3)
      vy0 = xlo(4)
      dvx = dx(3)
      dvy = dx(4)
      if (dir .eq. 0) then
        do i4 = n4a, n4b
          vy = vy0+(i4+0.5)*dvy
          vy2 = vy*vy
        do i3 = n3a, n3b
        do i2 = n2a, n2b
        do i1 = n1a, n1b
          if (i1 .eq. ng1b) then
            vx = face_vel1(i1+1, i2, i3, i4)
            fluxx = face_flux1(i1+1, i2, i3, i4)
          else
            vx = face_vel1(i1, i2, i3, i4)
            fluxx = face_flux1(i1, i2, i3, i4)
          endif
          ke_flux(i1, i2, i3, i4) = 0.5*fluxx*(vx*vx+vy2)
        end do
        end do
        end do
        end do
      else if (dir .eq. 1) then
        do i4 = n4a, n4b
        do i3 = n3a, n3b
          vx = vx0+(i3+0.5)*dvx
          vx2 = vx*vx
        do i2 = n2a, n2b
        do i1 = n1a, n1b
          if (i2 .eq. ng2b) then
            vy = face_vel2(i2+1, i3, i4, i1)
            fluxy = face_flux2(i2+1, i3, i4, i1)
          else
            vy = face_vel2(i2, i3, i4, i1)
            fluxy = face_flux2(i2, i3, i4, i1)
          endif
          ke_flux(i1, i2, i3, i4) = 0.5*fluxy*(vx2+vy*vy)
        end do
        end do
        end do
        end do
      endif
c
      return
      end
c
c **************
c
      subroutine computekeflux2d(
     &     nd1a, nd1b, nd2a, nd2b,
     &     n1a, n1b, n2a, n2b,
     &     ng1a, ng1b, ng2a, ng2b,
     &     xlo, xhi, dx,
     &     ke_flux_2d,
     &     ke_flux,
     &     dir, side)
c
c.. function to obtain local currents
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b
      integer n1a, n1b, n2a, n2b
      integer ng1a, ng1b, ng2a, ng2b
      real xlo(1:4), xhi(1:4), dx(1:4)
      real ke_flux_2d(nd1a:nd1b, nd2a:nd2b)
      real ke_flux
      integer dir, side
c
c.. declarations of local variables
      integer i1, i2
      integer i1lo, i1hi, i2lo, i2hi
      real ddir
      integer dosum
c
      dosum = 1
      if (dir .eq. 0) then
        i2lo = n2a
        i2hi = n2b
        ddir = dx(2)
        if (side .eq. 0) then
          if (n1a .eq. ng1a) then
            i1lo = n1a
            i1hi = n1a
          else
            dosum = 0
          endif
        else
          if (n1b .eq. ng1b) then
            i1lo = n1b
            i1hi = n1b
          else
            dosum = 0
          endif
        endif
      else
        i1lo = n1a
        i1hi = n1b
        ddir = dx(1)
        if (side .eq. 0) then
          if (n2a .eq. ng2a) then
            i2lo = n2a
            i2hi = n2a
          else
            dosum = 0
          endif
        else
          if (n2b .eq. ng2b) then
            i2lo = n2b
            i2hi = n2b
          else
            dosum = 0
          endif
        endif
      endif
      ke_flux = 0.0
      if (dosum .eq. 1) then
        do i2 = i2lo, i2hi
        do i1 = i1lo, i1hi

          ke_flux = ke_flux + ke_flux_2d(i1, i2)
        end do
        end do
        ke_flux = ke_flux*ddir
      endif
c
      return
      end
c
c **************
c
      subroutine computekevelspaceflux(
     &     nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b,
     &     n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b,
     &     ng1a, ng1b, ng2a, ng2b, ng3a, ng3b, ng4a, ng4b,
     &     xlo, xhi, dx,
     &     face_flux3,
     &     face_flux4,
     &     ke_flux,
     &     mass,
     &     side,
     &     dir)
c
c.. function to compute energy flux at velocity boundary
      implicit none
c
c.. declaration of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
      integer ng1a, ng1b, ng2a, ng2b, ng3a, ng3b, ng4a, ng4b
      real xlo(1:4), xhi(1:4), dx(1:4)
      real face_flux3(nd3a:nd3b+1, nd4a:nd4b, nd1a:nd1b, nd2a:nd2b)
      real face_flux4(nd4a:nd4b+1, nd1a:nd1b, nd2a:nd2b, nd3a:nd3b)
      real ke_flux(nd1a:nd1b, nd2a:nd2b)
      real mass
      integer side, dir
c
c.. declaration of local variables
      integer i1, i2, i3, i4
      real vx, vy, ddir, vx0, vy0, dvx, dvy, vx2, vy2, v2
c
      vx0 = xlo(3)
      vy0 = xlo(4)
      dvx = dx(3)
      dvy = dx(4)
      if (dir .eq. 2) then
        ddir = dx(4)
        if (side .eq. 0) then
          i3 = ng3a
        else
          i3 = ng3b+1
        endif
        vx = vx0+i3*dvx
        vx2 = vx*vx
        do i4 = n4a, n4b
          vy = vy0+(i4+0.5)*dvy
          v2 = vx2+vy*vy
        do i2 = n2a, n2b
        do i1 = n1a, n1b
          ke_flux(i1, i2) = ke_flux(i1, i2)+
     *      0.5*mass*face_flux3(i3, i4, i1, i2)*v2*ddir
        end do
        end do
        end do
      else if (dir .eq. 3) then
        ddir = dx(3)
        if (side .eq. 0) then
          i4 = ng4a
        else
          i4 = ng4b+1
        endif
        vy = vy0+i4*dvy
        vy2 = vy*vy
        do i3 = n3a, n3b
          vx = vx0+(i3+0.5)*dvx
          v2 = vx*vx+vy2
        do i2 = n2a, n2b
        do i1 = n1a, n1b
          ke_flux(i1, i2) = ke_flux(i1, i2)+
     *      0.5*mass*face_flux4(i4, i1, i2, i3)*v2*ddir
        end do
        end do
        end do
      endif
c
      return
      end
c
c **************
c
      subroutine appendkrook(
     &     nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b,
     &     n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b,
     &     xlo, xhi, deltax,
     &     xlo_krook, xhi_krook,
     &     krookPower, krookCoeff,
     &     u,
     &     rhs)
c
c.. function to append krook layer damping to distribution function
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
      real xlo(1:4), xhi(1:4), deltax(1:4)
      real xlo_krook(1:2), xhi_krook(1:2)
      real krookPower, krookCoeff
      real u(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real rhs(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
c
c.. declarations of local variables
      real g1, g2, g3
      real w1, w2, w3
      real x1, x2, x3, x4
      real vxt1, vxt2, vxt3
      real vyt1, vyt2, vyt3
      real fx1, fx2, fx3
      real fy1, fy2, fy3
      real integral1, integral2, integral3
      real f0, pi, xi, nu, eta, nux, nuy
      real xa, xb, ya, yb
      integer i1, i2, i3, i4
      real x0, y0, vx0, vy0, dx, dy, dvx, dvy, xmax, ymax
      real x0_krook, y0_krook, xmax_krook, ymax_krook
c
      real one, four
c
      one = 1.0
      four = 4.0
c
      pi = four*atan(one)
c
      ! define Gaussian quadrature points
      g1 = -0.5*sqrt(3.0/5.0)
      g2 = 0.0
      g3 = 0.5*sqrt(3.0/5.0)
c
      ! define Gaussian quadrature weights
      w1 = 5.0/18.0
      w2 = 8.0/18.0
      w3 = 5.0/18.0
c
      x0 = xlo(1)
      y0 = xlo(2)
      vx0 = xlo(3)
      vy0 = xlo(4)
      dx = deltax(1)
      dy = deltax(2)
      dvx = deltax(3)
      dvy = deltax(4)
      xmax = xhi(1)
      ymax = xhi(2)
      x0_krook = xlo_krook(1)
      y0_krook = xlo_krook(2)
      xmax_krook = xhi_krook(1)
      ymax_krook = xhi_krook(2)
      do i4 = n4a, n4b
        x4 = vy0+(0.5+i4)*dvy
        do i3 = n3a, n3b
          x3 = vx0+(0.5+i3)*dvx

          vxt1 = x3+g1*dvx
          vxt2 = x3+g2*dvx
          vxt3 = x3+g3*dvx

          vyt1 = x4+g1*dvy
          vyt2 = x4+g2*dvy
          vyt3 = x4+g3*dvy

          fx1 = exp(-0.5*(vxt1**2))
          fx2 = exp(-0.5*(vxt2**2))
          fx3 = exp(-0.5*(vxt3**2))

          fy1 = exp(-0.5*(vyt1**2))
          fy2 = exp(-0.5*(vyt2**2))
          fy3 = exp(-0.5*(vyt3**2))

          integral1 = w1*(fx1*fy1)+w2*(fx2*fy1)+w3*(fx3*fy1)
          integral2 = w1*(fx1*fy2)+w2*(fx2*fy2)+w3*(fx3*fy2)
          integral3 = w1*(fx1*fy3)+w2*(fx2*fy3)+w3*(fx3*fy3)
          f0 = w1*integral1+w2*integral2+w3*integral3
          f0 = f0/(2.0*pi)

          do i2 = n2a, n2b
            x2 = y0+(0.5+i2)*dy

            if (x2 .lt. y0_krook) then
              ya = y0_krook
              yb = y0
              eta = (x2-ya)/(yb-ya)
              nuy = eta**krookPower
            else if (x2 .gt. ymax_krook) then
              ya = ymax_krook
              yb = ymax
              eta = (x2-ya)/(yb-ya)
              nuy = eta**krookPower
            else
              nuy = 0.0
            end if

            do i1 = n1a, n1b
              x1 = x0+(0.5+i1)*dx

              if (x1 .lt. x0_krook) then
                xa = x0_krook
                xb = x0
                xi = (x1-xa)/(xb-xa)
                nux = xi**krookPower
              else if (x1 .gt. xmax_krook) then
                xa = xmax_krook
                xb = xmax
                xi = (x1-xa)/(xb-xa)
                nux = xi**krookPower
              else
                nux = 0.0
              end if

              nu = krookCoeff*((1.0-nuy)*nux+nuy)

              rhs(i1, i2, i3, i4) =
     *          rhs(i1, i2, i3, i4)-nu*(u(i1, i2, i3, i4)-f0)

            end do
          end do
        end do
      end do
c
      return
      end
