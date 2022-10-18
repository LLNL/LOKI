c
c Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
c See the top-level LICENSE file for details.
c Produced at the Lawrence Livermore National Laboratory
c
c SPDX-License-Identifier: MIT
c
c Fortran functions called by SimpleEMIC.
c
      subroutine setsimpleemic(nd1a, nd1b, nd2a, nd2b,
     *                         u,
     *                         xlo, xhi, dx,
     *                         num_waves,
     *                         xamp, yamp, zamp, kx, ky, phi,
     *                         iparams)
c
c.. function to set solution
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b
      real u(nd1a:nd1b, nd2a:nd2b, 1:6)
      real xlo(1:2), xhi(1:2), dx(1:2)
      integer num_waves
      real xamp(1:num_waves), yamp(1:num_waves), zamp(1:num_waves)
      real kx(1:num_waves), ky(1:num_waves), phi(1:num_waves)
      integer iparams(*)
c
c.. declarations of local variables
      integer i1, i2, w
      real x1, x2, env
      integer field

      field    = iparams(1)
      do w = 1,num_waves
        do i2 = nd2a,nd2b
          x2 = xlo(2)+(i2+0.5)*dx(2)
          do i1 = nd1a,nd1b
            x1 = xlo(1)+(i1+0.5)*dx(1)
            env = cos(kx(w)*x1+ky(w)*x2+phi(w))
            u(i1, i2, field)   = u(i1, i2, field) + xamp(w)*env
            u(i1, i2, field+1) = u(i1, i2, field+1) + yamp(w)*env
            u(i1, i2, field+2) = u(i1, i2, field+2) + zamp(w)*env
          end do
        end do
      end do
c
      return
      end
