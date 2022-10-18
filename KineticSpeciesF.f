c
c Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
c See the top-level LICENSE file for details.
c Produced at the Lawrence Livermore National Laboratory
c
c SPDX-License-Identifier: MIT
c
c Fortran functions called by KineticSpecies.
c
      subroutine xpby4d( 
     &     x,
     &     y,
     &     b,
     &     nd1lo,nd1hi,nd2lo,nd2hi,nd3lo,nd3hi,nd4lo,nd4hi,
     &     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b )
c
      implicit none
c
      integer nd1lo,nd1hi,nd2lo,nd2hi,nd3lo,nd3hi,nd4lo,nd4hi
      integer n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b
      real x( nd1lo:nd1hi,nd2lo:nd2hi,nd3lo:nd3hi,nd4lo:nd4hi )
      real y( nd1lo:nd1hi,nd2lo:nd2hi,nd3lo:nd3hi,nd4lo:nd4hi )
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
     &     ni1a,ni1b,ni2a,ni2b,ni3a,ni3b,ni4a,ni4b,
     &     vxface_velocities,
     &     vyface_velocities,
     &     normalization, bz_const,
     &     accel,
     &     na1a,na1b,na2a,na2b,
     &     axmax,aymax )
c
c.. declarations of incoming variables
c
      implicit none
      integer nv1a,nv1b,nv2a,nv2b
      integer nv3a,nv3b,nv4a,nv4b
      integer ni1a,ni1b,ni2a,ni2b
      integer ni3a,ni3b,ni4a,ni4b
      integer na1a,na1b,na2a,na2b
      real vxface_velocities( nv3a:nv3b+1,nv4a:nv4b,0:1 )
      real vyface_velocities( nv3a:nv3b,nv4a:nv4b+1,0:1 )
      real normalization, bz_const !IEO
      real accel( na1a:na1b,na2a:na2b,0:1 )
      real vel3( nv3a:nv3b+1,nv4a:nv4b,nv1a:nv1b,nv2a:nv2b )
      real vel4( nv4a:nv4b+1,nv1a:nv1b,nv2a:nv2b,nv3a:nv3b )
      real axmax,aymax
c
c.. declarations of local variables
c
      integer i1,i2,i3,i4
      real vx, vy !IEO

      axmax = 0.0
      do i4 = nv4a,nv4b
      do i3 = nv3a,nv3b+1
        vy  = vxface_velocities(i3, i4, 1)
      do i2 = na2a,na2b
      do i1 = na1a,na1b
        vel3(i3,i4,i1,i2) = accel(i1,i2,0) + 
     *    normalization * vy * bz_const
        if ((i1 .ge. ni1a .and. i1 .le. ni1b) .and.
     *      (i2 .ge. ni2a .and. i2 .le. ni2b) .and.
     *      (i3 .ge. ni3a .and. i3 .le. ni3b+1) .and.
     *      (i4 .ge. ni4a .and. i4 .le. ni4b)) then
          axmax = max(axmax,abs(vel3(i3,i4,i1,i2)))
        end if
      end do
      end do
      end do
      end do


      aymax = 0.0
      do i3 = nv3a,nv3b
      do i2 = na2a,na2b
      do i1 = na1a,na1b
      do i4 = nv4a,nv4b+1
        vx  = vyface_velocities(i3, i4, 0)
        vel4(i4,i1,i2,i3) = accel(i1,i2,1) - 
     *    normalization * vx * bz_const
        if ((i1 .ge. ni1a .and. i1 .le. ni1b) .and.
     *      (i2 .ge. ni2a .and. i2 .le. ni2b) .and.
     *      (i3 .ge. ni3a .and. i3 .le. ni3b) .and.
     *      (i4 .ge. ni4a .and. i4 .le. ni4b+1)) then
          aymax = max(aymax,abs(vel4(i4,i1,i2,i3)))
        end if
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
     &     ni1a, ni1b, ni2a, ni2b, ni3a, ni3b, ni4a, ni4b,
     &     vxface_velocities,
     &     vyface_velocities,
     &     normalization, bz_const,
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
      integer ni1a, ni1b, ni2a, ni2b, ni3a, ni3b, ni4a, ni4b
      real vxface_velocities(nv3a:nv3b+1, nv4a:nv4b, 0:1)
      real vyface_velocities(nv3a:nv3b, nv4a:nv4b+1, 0:1)
      real normalization, bz_const !IEO
      real vel3(nv3a:nv3b+1, nv4a:nv4b, nv1a:nv1b, nv2a:nv2b)
      real vel4(nv4a:nv4b+1, nv1a:nv1b, nv2a:nv2b, nv3a:nv3b)
      real em_vars(nv1a:nv1b, nv2a:nv2b, 1:6)
      real vz(nv1a:nv1b, nv2a:nv2b)
      real axmax, aymax
c
c.. declaration of local variables
      integer i1, i2, i3, i4
      real vx, vy
      real accel
c
c.. x component
      axmax = 0.0
      do i3 = nv3a, nv3b+1
      do i4 = nv4a, nv4b
        vy  = vxface_velocities(i3, i4, 1)
        do i1 = nv1a, nv1b
        do i2 = nv2a, nv2b
          accel = normalization*(em_vars(i1, i2, 1) +
     *       vy * em_vars(i1, i2, 6) +
     *       vy * bz_const - !IEO
     *       vz(i1, i2)*em_vars(i1, i2, 5))
          vel3(i3, i4, i1, i2) = accel
          if ((i1 .ge. ni1a .and. i1 .le. ni1b) .and.
     *        (i2 .ge. ni2a .and. i2 .le. ni2b) .and.
     *        (i3 .ge. ni3a .and. i3 .le. ni3b+1) .and.
     *        (i4 .ge. ni4a .and. i4 .le. ni4b)) then
            axmax = max(axmax, abs(accel))
          end if
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
        vx  = vyface_velocities(i3, i4, 0)
        accel = normalization*(em_vars(i1, i2, 2) +
     *     vz(i1, i2) * em_vars(i1, i2, 4) -
     *     vx * em_vars(i1, i2, 6) -
     *     vx * bz_const) !IEO
        vel4(i4, i1, i2, i3) = accel
        if ((i1 .ge. ni1a .and. i1 .le. ni1b) .and.
     *      (i2 .ge. ni2a .and. i2 .le. ni2b) .and.
     *      (i3 .ge. ni3a .and. i3 .le. ni3b) .and.
     *      (i4 .ge. ni4a .and. i4 .le. ni4b+1)) then
          aymax = max(aymax, abs(accel))
        end if
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
      eps = 1.e-10
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
     *     deltax)
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
      real deltax(1:4)
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
     &     vel3,vel4,ic )
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
      integer*8 ic
c
c.. declarations of local variables
      real initialconditionatpoint
      integer nghosts
      integer i1,i2,i3,i4,ig,loc
c
      if (solution_order .eq. 4) then
        nghosts = 2
      else
        nghosts = 3
      end if
c
c.. set i3 boundaries
      if ((ng3b-nghosts .eq. n3b) .or. (ng3a+nghosts .eq. n3a)) then
      do i4 = nl4a,nl4b
c
        if( ng3b-nghosts.eq.n3b ) then
        do i2 = nl2a,nl2b
        do i1 = nl1a,nl1b
          ! top
          if( vel3(n3b+1,i4,i1,i2) .ge. 0.0 ) then
            do ig = 1,nghosts
              u(i1,i2,n3b+ig,i4) = 3.0*u(i1,i2,n3b+ig-1,i4)-
     *          3.0*u(i1,i2,n3b+ig-2,i4)+u(i1,i2,n3b+ig-3,i4)
            end do
          else
            do ig = 1,nghosts
              loc = n3b+ig
              u(i1,i2,loc,i4) = initialconditionatpoint(ic,i1,i2,loc,i4)
            end do
          end if
        end do
        end do
        end if

        if( ng3a+nghosts.eq.n3a ) then
        do i2 = nl2a,nl2b
        do i1 = nl1a,nl1b
          ! bottom
          if( vel3(n3a,i4,i1,i2) .gt. 0.0 ) then
            do ig = 1,nghosts
              loc = n3a-ig
              u(i1,i2,loc,i4) = initialconditionatpoint(ic,i1,i2,loc,i4)
            end do
          else
            do ig = 1,nghosts
              u(i1,i2,n3a-ig,i4) = 3.0*u(i1,i2,n3a-ig+1,i4)-
     *          3.0*u(i1,i2,n3a-ig+2,i4)+u(i1,i2,n3a-ig+3,i4)
            end do
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
c
        if( ng4b-nghosts.eq.n4b ) then
        do i2 = nl2a,nl2b
        do i1 = nl1a,nl1b
          ! top
          if( vel4(n4b+1,i1,i2,i3) .ge. 0.0 ) then
            do ig = 1,nghosts
              u(i1,i2,i3,n4b+ig) = 3.0*u(i1,i2,i3,n4b+ig-1)-
     *          3.0*u(i1,i2,i3,n4b+ig-2)+u(i1,i2,i3,n4b+ig-3)
            end do
          else
            do ig = 1,nghosts
              loc = n4b+ig
              u(i1,i2,i3,loc) = initialconditionatpoint(ic,i1,i2,i3,loc)
            end do
          end if
        end do
        end do
        end if

        if( ng4a+nghosts.eq.n4a ) then
        do i2 = nl2a,nl2b
        do i1 = nl1a,nl1b
          ! bottom
          if( vel4(n4a,i1,i2,i3) .gt. 0.0 ) then
            do ig = 1,nghosts
              loc = n4a-ig
              u(i1,i2,i3,loc) = initialconditionatpoint(ic,i1,i2,i3,loc)
            end do
          else
            do ig = 1,nghosts
              u(i1,i2,i3,n4a-ig) = 3.0*u(i1,i2,i3,n4a-ig+1)-
     *          3.0*u(i1,i2,i3,n4a-ig+2)+u(i1,i2,i3,n4a-ig+3)
            end do
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
     &     vel1, vel2,
     &     xPeriodic, yPeriodic, ic )
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
      integer*8 ic
c
c.. declarations of local variables
      real initialconditionatpoint
      integer nghosts
      integer i1,i2,i3,i4,ig,loc
c
      if (solution_order .eq. 4) then
        nghosts = 2
      else
        nghosts = 3
      end if

c.. set i1 boundaries (if not periodic)
      if ((xPeriodic .ne. 1) .and.
     *    ((ng1b-nghosts .eq. n1b) .or. (ng1a+nghosts .eq. n1a))) then
        do i4 = nl4a,nl4b
        do i3 = nl3a,nl3b
          if( ng1b-nghosts.eq.n1b ) then
          do i2 = nl2a,nl2b
            ! top
            if( vel1(n1b+1,i2,i3,i4) .ge. 0.0 ) then
              do ig = 1,nghosts
                u(n1b+ig,i2,i3,i4) = 3.0*u(n1b+ig-1,i2,i3,i4)-
     *            3.0*u(n1b+ig-2,i2,i3,i4)+u(n1b+ig-3,i2,i3,i4)
              end do
            else
              ! sample the background distribution
              do ig = 1,nghosts
                loc = n1b+ig
                u(loc,i2,i3,i4) =
     *            initialconditionatpoint(ic,loc,i2,i3,i4)
              end do
            end if
          end do
          end if

          if( ng1a+nghosts.eq.n1a ) then
          do i2 = nl2a,nl2b
            ! bottom
            if( vel1(n1a,i2,i3,i4) .gt. 0.0 ) then
              ! sample the background distribution
              do ig = 1,nghosts
                loc = n1a-ig
                u(loc,i2,i3,i4) =
     *            initialconditionatpoint(ic,loc,i2,i3,i4)
              end do
            else
              do ig = 1,nghosts
                u(n1a-ig,i2,i3,i4) = 3.0*u(n1a-ig+1,i2,i3,i4)-
     *            3.0*u(n1a-ig+2,i2,i3,i4)+u(n1a-ig+3,i2,i3,i4)
              end do
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
        do i3 = nl3a,nl3b
          if( ng2b-nghosts.eq.n2b ) then
          do i1 = nl1a,nl1b
            ! top
            if( vel2(n2b+1,i3,i4,i1) .ge. 0.0 ) then
              do ig = 1,nghosts
                u(i1,n2b+ig,i3,i4) = 3.0*u(i1,n2b+ig-1,i3,i4)-
     *            3.0*u(i1,n2b+ig-2,i3,i4)+u(i1,n2b+ig-3,i3,i4)
              end do
            else
              ! sample the background distribution
              do ig = 1,nghosts
                loc = n2b+ig
                u(i1,loc,i3,i4) =
     *            initialconditionatpoint(ic,i1,loc,i3,i4)
              end do
            end if
          end do
          end if

          if( ng2a+nghosts.eq.n2a ) then
          do i1 = nl1a,nl1b
            ! bottom
            if( vel2(n2a,i3,i4,i1) .gt. 0.0 ) then
              ! sample the background distribution
              do ig = 1,nghosts
                loc = n2a-ig
                u(i1,loc,i3,i4) =
     *            initialconditionatpoint(ic,i1,loc,i3,i4)
              end do
            else
              do ig = 1,nghosts
                u(i1,n2a-ig,i3,i4) = 3.0*u(i1,n2a-ig+1,i3,i4)-
     *            3.0*u(i1,n2a-ig+2,i3,i4)+u(i1,n2a-ig+3,i3,i4)
              end do
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
     &     xlo, xhi, deltax,
     &     ic)
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
      real xlo(1:4), xhi(1:4), deltax(1:4)
      integer*8 ic
c
c.. declarations of local variables
      real initialconditionatpoint
      integer i1,i2,i3,i4,loc
      integer ig,extrapEq,iSten
      real eCoeffs(1:6,1:6)
      real vxlo,vylo,vxhi,vyhi
      real dvx,dvy
      real vx,vy
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
      vxlo = xlo(3)
      vylo = xlo(4)
      vxhi = xhi(3)
      vyhi = xhi(4)
      dvx = deltax(3)
      dvy = deltax(4)
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
              loc = i3-ig
              u(i1,i2,loc,i4) = initialconditionatpoint(ic,i1,i2,loc,i4)
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
              loc = i3+ig
              u(i1,i2,loc,i4) = initialconditionatpoint(ic,i1,i2,loc,i4)
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
              loc = i4-ig
              u(i1,i2,i3,loc) = initialconditionatpoint(ic,i1,i2,i3,loc)
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
              loc = i4+ig
              u(i1,i2,i3,loc) = initialconditionatpoint(ic,i1,i2,i3,loc)
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
     &     nl1a,nl1b,nl2a,nl2b,nl3a,nl3b,nl4a,nl4b,
     &     n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,
     &     solution_order, nGhost,
     &     vel1, vel2,
     &     xloV, xhiV, dxV,
     &     xPeriodic, yPeriodic,
     &     ic)
c
c.. boundary conditions 
c
c.. declarations of incoming variables
      implicit none
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
      integer*8 ic
c
c.. declartions of local variables
      real initialconditionatpoint
      integer i1,i2,i3,i4,loc
      integer ig,extrapEq,iSten
      real eCoeffs(1:6,1:6)
      real xlo,xhi,ylo,yhi
      real dx,dy
      real x,y
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
      xlo  = xloV(1)
      ylo  = xloV(2)
      
      xhi  = xhiV(1)
      yhi  = xhiV(2)
      
      dx  = dxV(1)
      dy  = dxV(2)
c
c.. set lower i1 boundary
      i1 = n1a
      x  = xlo+(0.5+i1-0.5)*dx
      if( abs(x-xlo) .lt. 0.5*dx .and. xPeriodic .ne. 1 ) then
        do i4 = nl4a,nl4b
        do i3 = nl3a,nl3b
        do i2 = nl2a,nl2b
          ! check if outflow or inflow
          if( vel1(i1,i2,i3,i4) .gt. 0.0 ) then
            ! inflow so sample distribution fcn
            ! loop over all ghosts
            do ig = 1,nGhost
              loc = i1-ig
              u(loc,i2,i3,i4) = initialconditionatpoint(ic,loc,i2,i3,i4)
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
              loc = i1+ig
              u(loc,i2,i3,i4) = initialconditionatpoint(ic,loc,i2,i3,i4)
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
      y  = ylo+(0.5+i2-0.5)*dy
      if( abs(y-ylo) .lt. 0.5*dy .and. yPeriodic .ne. 1 ) then
        do i4 = nl4a,nl4b
        do i3 = nl3a,nl3b
        do i1 = nl1a,nl1b
          ! check if outflow or inflow
          if( vel2(i2,i3,i4,i1) .gt. 0.0 ) then
            ! inflow so sample distribution fcn
            ! loop over all ghosts
            do ig = 1,nGhost
              loc = i2-ig
              u(i1,loc,i3,i4) = initialconditionatpoint(ic,i1,loc,i3,i4)
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
              loc = i2+ig
              u(i1,loc,i3,i4) = initialconditionatpoint(ic,i1,loc,i3,i4)
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
     &     vel1, vel2,
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
        do i3 = n3a,n3b
          vx = vel1(n1a,n2a,i3,i4)
          vy = vel2(n2a,i3,i4,n1a)

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

              ax = vel3(n3a,i4,i1,i2) !IEO
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

              ay = vel4(n4a,i1,i2,i3) !IEO
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

              ax = vel3(n3a,i4,i1,i2) !IEO
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

              ay = vel4(n4a,i1,i2,i3) !IEO
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
     &     velocities,
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
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real u(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real vz(nd1a:nd1b, nd2a:nd2b)
      real Jx(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real Jy(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real Jz(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      real vx, vy
c
      do i4 = n4a, n4b
      do i3 = n3a, n3b
        vx = velocities(i3, i4, 0)
        vy = velocities(i3, i4, 1)
      do i2 = n2a, n2b
      do i1 = n1a, n1b

        Jx(i1, i2, i3, i4) = u(i1, i2, i3, i4)*vx

        Jy(i1, i2, i3, i4) = u(i1, i2, i3, i4)*vy

        Jz(i1, i2, i3, i4) = u(i1, i2, i3, i4)*vz(i1, i2)
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
      subroutine computeke(
     &     nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b,
     &     n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b,
     &     xlo, xhi, deltax,
     &     u,
     &     mass,
     &     velocities,
     &     ke, ke_x, ke_y, px, py)
c
c.. function to obtain local currents
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
      real xlo(1:4), xhi(1:4), deltax(1:4)
      real u(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real mass, ke, ke_x, ke_y, px, py
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      real vx, vy, dx, dy, dvx, dvy, vx2, vy2, this_u
c
      dx  = deltax(1)
      dy  = deltax(2)
      dvx = deltax(3)
      dvy = deltax(4)
      do i4 = n4a, n4b
      do i3 = n3a, n3b
        vx = velocities(i3, i4, 0)
        vx2 = vx*vx
        vy = velocities(i3, i4, 1)
        vy2 = vy*vy
      do i2 = n2a, n2b
      do i1 = n1a, n1b

        this_u = u(i1, i2, i3, i4)
        ke_x = ke_x + 0.5*this_u*vx2
        ke_y = ke_y + 0.5*this_u*vy2
        px = px + this_u*vx
        py = py + this_u*vy
      end do
      end do
      end do
      end do
      ke_x = ke_x*mass*dx*dy*dvx*dvy
      ke_y = ke_y*mass*dx*dy*dvx*dvy
      px = px*mass*dx*dy*dvx*dvy
      py = py*mass*dx*dy*dvx*dvy
      ke = ke_x+ke_y
c
      return
      end
c
c **************
c
      subroutine computekemaxwell(
     &     nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b,
     &     n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b,
     &     xlo, xhi, deltax,
     &     u,
     &     mass,
     &     velocities, vz_in,
     &     ke, ke_x, ke_y)
c
c.. function to obtain local currents
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
      real xlo(1:4), xhi(1:4), deltax(1:4)
      real u(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real vz_in(nd1a:nd1b, nd2a:nd2b)
      real mass, ke, ke_x, ke_y
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      real vx, vy, vz, dx, dy, dvx, dvy, vx2, vy2, vz2,ke_z, this_u
c
      ke_z = 0.0
      dx  = deltax(1)
      dy  = deltax(2)
      dvx = deltax(3)
      dvy = deltax(4)
      do i4 = n4a, n4b
      do i3 = n3a, n3b
        vx = velocities(i3, i4, 0)
        vx2 = vx*vx
        vy = velocities(i3, i4, 1)
        vy2 = vy*vy
      do i2 = n2a, n2b
      do i1 = n1a, n1b

        vz = vz_in(i1, i2)
        vz2 = vz*vz
        this_u = u(i1, i2, i3, i4)
        ke_x = ke_x + 0.5*this_u*vx2
        ke_y = ke_y + 0.5*this_u*vy2
        ke_z = ke_z + 0.5*this_u*vz2
      end do
      end do
      end do
      end do
      ke_x = ke_x*mass*dx*dy*dvx*dvy
      ke_y = ke_y*mass*dx*dy*dvx*dvy
      ke_z = ke_z*mass*dx*dy*dvx*dvy
      ke = ke_x+ke_y+ke_z
c
      return
      end
c
c **************
c
      subroutine computekeedot(
     &     nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b,
     &     n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b,
     &     xlo, xhi, dx,
     &     u,
     &     charge,
     &     velocities,
     &     ext_efield,
     &     ke_e_dot)
c
c.. function to obtain local currents
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
      real xlo(1:4), xhi(1:4), dx(1:4)
      real u(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real ext_efield(nd1a:nd1b, nd2a:nd2b, 0:1)
      real charge, ke_e_dot
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      real vx
c
      do i4 = n4a, n4b
      do i3 = n3a, n3b
        vx = velocities(i3, i4, 0)
      do i2 = n2a, n2b
      do i1 = n1a, n1b
        ke_e_dot = ke_e_dot + ext_efield(i1, i2, 0)*vx*u(i1, i2, i3, i4)
      end do
      end do
      end do
      end do
      ke_e_dot = ke_e_dot*charge*dx(1)*dx(2)*dx(3)*dx(4)
c
      return
      end
c
c **************
c
      subroutine computemom4d(
     &     nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b,
     &     n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b,
     &     velocities,
     &     u, momx, momy, ke, ent)
c
c.. function to obtain local currents
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real u(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real momx(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real momy(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real ke(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real ent(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      real vx, vy, rlogu, utemp, v2
c
      do i4 = n4a, n4b
      do i3 = n3a, n3b
        vx = velocities(i3, i4, 0)
        vy = velocities(i3, i4, 1)
        v2 = vx*vx+vy*vy
      do i2 = n2a, n2b
         do i1 = n1a, n1b
            utemp = u(i1, i2, i3, i4)
            momx(i1, i2, i3, i4) = utemp*vx
            momy(i1, i2, i3, i4) = utemp*vy
            ke(i1, i2, i3, i4) = 0.5*u(i1, i2, i3, i4)*v2
            rlogu = log(abs(utemp) + 1.0e-15)
            ent(i1, i2, i3, i4) = ent(i1, i2, i3, i4) - utemp * rlogu
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
      subroutine integrate2d(
     &     nd1a, nd1b, nd2a, nd2b,
     &     n1a, n1b, n2a, n2b,
     &     dx,
     &     integrand_2d,
     &     integral)
c
c.. function to obtain local currents
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b
      integer n1a, n1b, n2a, n2b
      real dx(1:4)
      real integrand_2d(nd1a:nd1b, nd2a:nd2b)
      real integral
c
c.. declarations of local variables
      integer i1, i2
c
      integral = 0.0
      do i2 = n2a, n2b
      do i1 = n1a, n1b

        integral = integral + integrand_2d(i1, i2)
      end do
      end do
      integral = integral*dx(1)*dx(2)
c
      return
      end
c
c **************
c
      subroutine computemom2d(
     &     nd1a, nd1b, nd2a, nd2b,
     &     n1a, n1b, n2a, n2b,
     &     dx,
     &     momx_2d, momy_2d, ke_2d, ent_2d,
     &     momx, momy, ke, ent)
C
C.. Function to obtain local currents
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b
      integer n1a, n1b, n2a, n2b
      real dx(1:4)
      real momx_2d(nd1a:nd1b, nd2a:nd2b)
      real momx
      real momy_2d(nd1a:nd1b, nd2a:nd2b)
      real momy
      real ke_2d(nd1a:nd1b, nd2a:nd2b)
      real ke
      real ent_2d(nd1a:nd1b, nd2a:nd2b)
      real ent
c
c.. declarations of local variables
      integer i1, i2
c
      momx = 0.0
      momy = 0.0
      ke = 0.0
      ent = 0.0
      do i2 = n2a, n2b
      do i1 = n1a, n1b
        momx = momx + momx_2d(i1, i2)
        momy = momy + momy_2d(i1, i2)
        ke = ke + ke_2d(i1, i2)
        ent = ent + ent_2d(i1, i2)
      end do
      end do
      momx = momx*dx(1)*dx(2)
      momy = momy*dx(1)*dx(2)
      ke = ke*dx(1)*dx(2)
      ent = ent*dx(1)*dx(2)
c
      return
      end
c
c **************
c
      subroutine computekeflux(
     &     nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b,
     &     n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b,
     &     ng1a, ng1b, ng2a, ng2b, ng3a, ng3b, ng4a, ng4b,
     &     dx,
     &     face_flux1,
     &     face_flux2,
     &     face_flux3,
     &     face_flux4,
     &     velocities,
     &     vxface_velocities,
     &     vyface_velocities,
     &     dir,
     &     side,
     &     mass,
     &     ke_flux)
c
c.. function to compute energy flux at physical boundary
      implicit none
c
c.. declaration of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
      integer ng1a, ng1b, ng2a, ng2b, ng3a, ng3b, ng4a, ng4b
      real dx(1:4)
      real face_flux1(nd1a:nd1b+1, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real face_flux2(nd2a:nd2b+1, nd3a:nd3b, nd4a:nd4b, nd1a:nd1b)
      real face_flux3(nd3a:nd3b+1, nd4a:nd4b, nd1a:nd1b, nd2a:nd2b)
      real face_flux4(nd4a:nd4b+1, nd1a:nd1b, nd2a:nd2b, nd3a:nd3b)
      real velocities(nd3a:nd3b, nd4a:nd4b, 0:1)
      real vxface_velocities(nd3a:nd3b+1, nd4a:nd4b, 0:1)
      real vyface_velocities(nd3a:nd3b, nd4a:nd4b+1, 0:1)
      integer dir, side
      real mass, ke_flux
c
c.. declaration of local variables
      integer i1, i2, i3, i4
      real ddir
      integer faceidx, dosum
      real vx, vy, fluxx, fluxy, fluxvx, fluxvy
      real v2
c
      dosum = 1
      if (dir .eq. 0) then
        ddir = dx(2)*dx(3)*dx(4)
        if (side .eq. 0) then
          if (n1a .eq. ng1a) then
            faceidx = n1a
          else
            dosum = 0
          endif
        else
          if (n1b .eq. ng1b) then
            faceidx = n1b+1
          else
            dosum = 0
          endif
        endif
        if (dosum .eq. 1) then
          do i4 = n4a, n4b
          do i3 = n3a, n3b
            vx = velocities(i3, i4, 0)
            vy = velocities(i3, i4, 1)
            v2 = vx*vx+vy*vy
          do i2 = n2a, n2b
            fluxx = face_flux1(faceidx, i2, i3, i4)
            ke_flux = ke_flux + 0.5*fluxx*v2
          end do
          end do
          end do
        endif
      else if (dir .eq. 1) then
        ddir = dx(1)*dx(3)*dx(4)
        if (side .eq. 0) then
          if (n2a .eq. ng2a) then
            faceidx = n2a
          else
            dosum = 0
          endif
        else
          if (n2b .eq. ng2b) then
            faceidx = n2b+1
          else
            dosum = 0
          endif
        endif
        if (dosum .eq. 1) then
          do i4 = n4a, n4b
          do i3 = n3a, n3b
            vx = velocities(i3, i4, 0)
            vy = velocities(i3, i4, 1)
            v2 = vx*vx+vy*vy
          do i1 = n1a, n1b
            fluxy = face_flux2(faceidx, i3, i4, i1)
            ke_flux = ke_flux + 0.5*fluxy*v2
          end do
          end do
          end do
        endif
      else if (dir .eq. 2) then
        ddir = dx(1)*dx(2)*dx(4)
        if (side .eq. 0) then
          if (n3a .eq. ng3a) then
            faceidx = n3a
          else
            dosum = 0
          endif
        else
          if (n3b .eq. ng3b) then
            faceidx = n3b + 1
          else
            dosum = 0
          endif
        endif
        if (dosum .eq. 1) then
          do i4 = n4a, n4b
            vx = vxface_velocities(faceidx, i4, 0)
            vy = vxface_velocities(faceidx, i4, 1)
            v2 = vx*vx+vy*vy
          do i2 = n2a, n2b
          do i1 = n1a, n1b
            fluxvx = face_flux3(faceidx, i4, i1, i2)
            ke_flux = ke_flux + 0.5*fluxvx*v2
          end do
          end do
          end do
        endif
      else if (dir .eq. 3) then
        ddir = dx(1)*dx(2)*dx(3)
        if (side .eq. 0) then
          if (n4a .eq. ng4a) then
            faceidx = n4a
          else
            dosum = 0
          endif
        else
          if (n4b .eq. ng4b) then
            faceidx = n4b + 1
          else
            dosum = 0
          endif
        endif
        if (dosum .eq. 1) then
          do i3 = n3a, n3b
            vx = vyface_velocities(i3, faceidx, 0)
            vy = vyface_velocities(i3, faceidx, 1)
            v2 = vx*vx+vy*vy
          do i2 = n2a, n2b
          do i1 = n1a, n1b
            fluxvy = face_flux4(faceidx, i1, i2, i3)
            ke_flux = ke_flux + 0.5*fluxvy*v2
          end do
          end do
          end do
        endif
      endif
      ke_flux = ke_flux*mass*ddir
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
     &     dx,
     &     face_flux3,
     &     face_flux4,
     &     ke_flux,
     &     mass,
     &     vxface_velocities,
     &     vyface_velocities,
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
      real dx(1:4)
      real face_flux3(nd3a:nd3b+1, nd4a:nd4b, nd1a:nd1b, nd2a:nd2b)
      real face_flux4(nd4a:nd4b+1, nd1a:nd1b, nd2a:nd2b, nd3a:nd3b)
      real ke_flux(nd1a:nd1b, nd2a:nd2b)
      real mass
      real vxface_velocities(nd3a:nd3b+1, nd4a:nd4b, 0:1)
      real vyface_velocities(nd3a:nd3b, nd4a:nd4b+1, 0:1)
      integer side, dir
c
c.. declaration of local variables
      integer i1, i2, i3, i4, dosum
      real vx, vy, ddir, v2
c
      dosum = 1
      if (dir .eq. 2) then
        ddir = dx(4)
        if (side .eq. 0) then
          if (n3a .eq. ng3a) then
            i3 = ng3a
          else
            dosum = 0
          endif
        else
          if (n3b .eq. ng3b) then
            i3 = ng3b+1
          else
            dosum = 0
          endif
        endif
        if (dosum .eq. 1) then
          do i4 = n4a, n4b
            vx = vxface_velocities(i3, i4, 0)
            vy = vxface_velocities(i3, i4, 1)
            v2 = vx*vx+vy*vy
          do i2 = n2a, n2b
          do i1 = n1a, n1b
            ke_flux(i1, i2) = ke_flux(i1, i2)+
     *        0.5*mass*face_flux3(i3, i4, i1, i2)*v2*ddir
          end do
          end do
          end do
        endif
      else if (dir .eq. 3) then
        ddir = dx(3)
        if (side .eq. 0) then
          if (n4a .eq. ng4a) then
            i4 = ng4a
          else
            dosum = 0
          endif
        else
          if (n4b .eq. ng4b) then
            i4 = ng4b+1
          else
            dosum = 0
          endif
        endif
        if (dosum .eq. 1) then
          do i3 = n3a, n3b
            vx = vyface_velocities(i3, i4, 0)
            vy = vyface_velocities(i3, i4, 1)
            v2 = vx*vx+vy*vy
          do i2 = n2a, n2b
          do i1 = n1a, n1b
            ke_flux(i1, i2) = ke_flux(i1, i2)+
     *        0.5*mass*face_flux4(i4, i1, i2, i3)*v2*ddir
          end do
          end do
          end do
        endif
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
     &     dt,
     &     ic,
     &     nu, u, rhs)
c
c.. function to append krook layer damping to distribution function
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b, nd3a, nd3b, nd4a, nd4b
      integer n1a, n1b, n2a, n2b, n3a, n3b, n4a, n4b
      real dt
      integer*8 ic
      real nu(nd1a:nd1b, nd2a:nd2b)
      real u(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real rhs(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
c
c.. declarations of local variables
      real f0, initialconditionatpoint
      integer i1, i2, i3, i4
c
      do i4 = n4a, n4b
        do i3 = n3a, n3b
          do i2 = n2a, n2b
            do i1 = n1a, n1b
              if (nu(i1, i2) .ne. 0.0) then
                f0 = initialconditionatpoint(ic, i1, i2, i3, i4)
                rhs(i1, i2, i3, i4) = rhs(i1, i2, i3, i4)-
     *            nu(i1, i2)/dt*(u(i1, i2, i3, i4)-f0)
              end if

            end do
          end do
        end do
      end do
c
      return
      end
