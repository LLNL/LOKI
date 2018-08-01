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
c Fortran functions called by External2DIC.
c
      subroutine setIC_external_2D(nd1a, nd1b, nd2a, nd2b,
     *                             nd3a, nd3b, nd4a, nd4b,
     *                             ng1a, ng1b, ng2a, ng2b,
     *                             ng3a, ng3b, ng4a, ng4b,
     *                             u,
     *                             u_spatial_init,
     *                             ic_phase_param,
     *                             xlo, xhi, deltax,
     *                             n_cells1,
     *                             n_cells2,
     *                             n_cells3,
     *                             n_cells4)
c
c.. function to set solution
      implicit none
c
c.. declarations of incoming variables
      integer nd1a, nd1b, nd2a, nd2b
      integer nd3a, nd3b, nd4a, nd4b
      integer ng1a, ng1b, ng2a, ng2b
      integer ng3a, ng3b, ng4a, ng4b
      integer n_cells1, n_cells2, n_cells3, n_cells4
      real u(nd1a:nd1b, nd2a:nd2b, nd3a:nd3b, nd4a:nd4b)
      real u_spatial_init(ng1a:ng1b, ng2a:ng2b)
      real ic_phase_param(*)
      real xlo(1:4), xhi(1:4), deltax(1:4)
c
c.. declarations of local variables
      integer i1, i2, i3, i4
      real xmin, ymin, vxmin, vymin, dx, dy, dvx, dvy
      real x1, x2, x3, x4
      real alpha, beta, vx0, vy0, vflowx, vflowy, kx, ky, flow_vel_phi
c
      real pi, one, four, vx, vy
c
      one = 1.0
      four = 4.0
c
      pi = four*atan(one)
c
      xmin   = xlo(1)
      ymin   = xlo(2)
      vxmin  = xlo(3)
      vymin  = xlo(4)
      dx     = deltax(1)
      dy     = deltax(2)
      dvx    = deltax(3)
      dvy    = deltax(4)
c
      alpha        = ic_phase_param(1)
      beta         = ic_phase_param(2)
      vx0          = ic_phase_param(3)
      vy0          = ic_phase_param(4)
      vflowx       = ic_phase_param(5)
      vflowy       = ic_phase_param(6)
      kx           = ic_phase_param(7)
      ky           = ic_phase_param(8)
      flow_vel_phi = ic_phase_param(9)
c
      do i4 = nd4a, nd4b
        x4 = vymin+(i4+0.5)*dvy
        do i3 = nd3a, nd3b
          x3 = vxmin+(i3+0.5)*dvx
          do i2 = nd2a, nd2b
            x2 = ymin+(i2+0.5)*dy
            do i1 = nd1a, nd1b
              x1 = xmin+(i1+0.5)*dx

              vx = vx0*cos(kx*x1+ky*x2+flow_vel_phi)+vflowx
              vy = vy0*cos(kx*x1+ky*x2+flow_vel_phi)+vflowy
              u(i1, i2, i3, i4) =
     *           alpha*beta/(2.0*pi)*exp(-0.5*((alpha*(x3-vx))**2+
     *           (beta*(x4-vy))**2))*u_spatial_init(i1, i2)
              
            end do
          end do
        end do
      end do
c
      return
      end
