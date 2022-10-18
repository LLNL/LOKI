c
c Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
c See the top-level LICENSE file for details.
c Produced at the Lawrence Livermore National Laboratory
c
c SPDX-License-Identifier: MIT
c
c Computes specified relativistic velocity.
c
      subroutine getVelocity(
     &     px, py, mass, c, vx, vy)
c
c.. function to compute velocity
      implicit none
c
c.. declarations of incoming variables
      real px, py, mass, c, vx, vy
c
c.. declarations of local variables
      real pdivc
c
      pdivc = (px*px+py*py)/(c*c)
      vx = px / sqrt(mass*mass + pdivc)
      vy = py / sqrt(mass*mass + pdivc)
      return
      end
      
