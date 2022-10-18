c
c Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
c See the top-level LICENSE file for details.
c Produced at the Lawrence Livermore National Laboratory
c
c SPDX-License-Identifier: MIT
c
c Computes specified relativistic momentum.
c
      subroutine getMomentum(
     &     vx, vy, mass, c, px, py)
c
c.. function to compute velocity
      implicit none
c
c.. declarations of incoming variables
      real vx, vy, mass, c, px, py
c
c.. declarations of local variables
      real vdivc
c
      vdivc = (vx*vx+vy*vy)/(c*c)
      px = mass*vx / sqrt(1.0 - vdivc)
      py = mass*vy / sqrt(1.0 - vdivc)
      return
      end
      
