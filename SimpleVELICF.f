c
c Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
c See the top-level LICENSE file for details.
c Produced at the Lawrence Livermore National Laboratory
c
c SPDX-License-Identifier: MIT
c
c Fortran functions called by SimpleVELIC.
c
      subroutine setsimplevelic(vz,
     *                          nd1a, nd1b, nd2a, nd2b,
     *                          xlo, xhi, dx,
     *                          num_waves,
     *                          amp, kx, ky, phi)
c
c.. function to set solution
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b
      real vz(nd1a:nd1b, nd2a:nd2b)
      real xlo(1:2), xhi(1:2), dx(1:2)
      integer num_waves
      real amp(1:num_waves), kx(1:num_waves), ky(1:num_waves)
      real phi(1:num_waves)
c
c.. declarations of local variables
      integer i1, i2, w
      real x1, x2
c
      do w = 1,num_waves
        do i2 = nd2a,nd2b
          x2 = xlo(2)+(i2+0.5)*dx(2)
          do i1 = nd1a,nd1b
            x1 = xlo(1)+(i1+0.5)*dx(1)
            vz(i1, i2) = vz(i1, i2) +
     *        amp(w)*cos(kx(w)*x1+ky(w)*x2+phi(w))
          end do
        end do
      end do
c
      return
      end
