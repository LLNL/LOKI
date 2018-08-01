/*************************************************************************
 *
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * Written by Jeffrey Banks banksj3@rpi.edu (Rensselaer Polytechnic Institute,
 * Amos Eaton 301, 110 8th St., Troy, NY 12180); Jeffrey Hittinger
 * hittinger1@llnl.gov, William Arrighi arrighi2@llnl.gov, Richard Berger
 * berger5@llnl.gov, Thomas Chapman chapman29@llnl.gov (LLNL, P.O Box 808,
 * Livermore, CA 94551); Stephan Brunner stephan.brunner@epfl.ch (Ecole
 * Polytechnique Federale de Lausanne, EPFL SB SPC-TH, PPB 312, Station 13,
 * CH-1015 Lausanne, Switzerland).
 * CODE-744849
 *
 * All rights reserved.
 *
 * This file is part of Loki.  For details, see.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 ************************************************************************/
#include "CollisionOperator.H"
#include "CollisionOperatorF.H"

namespace Loki {

CollisionOperator::CollisionOperator()
   : m_vthermal_method(INPUT_VTHERMAL),
     m_dist_func_avg(0)
{
}

CollisionOperator::~CollisionOperator()
{
   if (m_dist_func_avg != 0) {
      delete [] m_dist_func_avg;
   }
}

real
CollisionOperator::computeVthLocal(
   const RealArray& a_u,
   const tbox::Box& a_local_box,
   const tbox::Box& a_interior_box,
   const ProblemDomain& a_domain,
   real a_vflowx,
   real a_vflowy)
{
   // Compute the mean spatial distribution locally and use it in a global
   // calculation of vthermal.
   real vth;
   FORT_COMPUTE_VTH_LOCAL(BOX4D_TO_FORT(a_local_box),
      BOX4D_TO_FORT(a_interior_box),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      a_vflowx,
      a_vflowy,
      m_comm,
      *a_u.getDataPointer(),
      vth);
   return vth;
}

real
CollisionOperator::computeVthGlobal(
   const RealArray& a_u,
   const tbox::Box& a_local_box,
   const tbox::Box& a_interior_box,
   const ProblemDomain& a_domain,
   real a_vflowx,
   real a_vflowy)
{
   // Compute the mean spatial distribution globally.
   if (m_dist_func_avg == 0) {
      m_vel_space_size =
         (a_interior_box.upper(2)-a_interior_box.lower(2)+1)*
         (a_interior_box.upper(3)-a_interior_box.lower(3)+1);
      m_dist_func_avg = new real [m_vel_space_size];

      m_config_space_size =
         a_domain.box().numberCells(0)*a_domain.box().numberCells(1);

      int color = a_interior_box.lower(3)*a_domain.box().numberCells(2)+
                  a_interior_box.lower(2);

      int comm_id;
      MPI_Comm_rank(m_comm, &comm_id);
      const int status = MPI_Comm_split(m_comm,
         color,
         comm_id,
         &m_config_space_comm);
      if (status != MPI_SUCCESS) {
         OV_ABORT("Configuration space splitting of MPI communicator failed");
      }

      int config_space_comm_id;
      MPI_Comm_rank(m_config_space_comm, &config_space_comm_id);
      m_is_config_space_head_node = config_space_comm_id == 0 ? 1 : 0;
   }
   FORT_COMPUTE_DIST_FUNC_AVG(BOX4D_TO_FORT(a_local_box),
      BOX4D_TO_FORT(a_interior_box),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      m_vel_space_size,
      m_config_space_size,
      *a_u.getDataPointer(),
      *m_dist_func_avg);
   MPI_Allreduce(MPI_IN_PLACE,
      m_dist_func_avg,
      m_vel_space_size,
      MPI_DOUBLE,
      MPI_SUM,
      m_config_space_comm);

   // Use the mean spatial distribution in a global calculation of vthermal.
   real vth;
   FORT_COMPUTE_VTH_GLOBAL(BOX4D_TO_FORT(a_local_box),
      BOX4D_TO_FORT(a_interior_box),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      m_vel_space_size,
      a_vflowx,
      a_vflowy,
      *m_dist_func_avg,
      m_is_config_space_head_node,
      m_comm,
      vth);
   return vth;
}

} // end namespace Loki
