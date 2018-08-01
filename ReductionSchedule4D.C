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
#include "ReductionSchedule4D.H"

#include "BoxOps.H"
#include "KernelOp.H"
#include "TimerManager.H"
#include "Loki_Utilities.H"

namespace Loki {

const int ReductionSchedule4D::s_MAX_DIM = 4;

ReductionSchedule4D::ReductionSchedule4D(
   const RealArray& a_local_src_array,
   const tbox::Box& a_local_interior_box,
   const tbox::Box& a_domain_box,
   const Array<bool>& a_collapse_dir,
   const tbox::Pointer<ReductionOp>& a_reduction,
   const Range& a_processor_range,
   const MPI_Comm& a_comm)
   : m_local_src_array(a_local_src_array),
     m_collapse_dir(a_collapse_dir),
     m_is_in_proc_range(false),
     m_reduction_comm(a_comm),
     m_reduction(a_reduction)
{
   const int my_id(std::max(0, Communication_Manager::My_Process_Number));
   m_is_in_proc_range = isInRange(my_id, a_processor_range);

   // Construct sub-communicators for each group of processors dealing with
   // the same piece of configuration space.
   constructCommunicator(a_local_interior_box, a_domain_box, a_comm);
}


ReductionSchedule4D::~ReductionSchedule4D()
{
}


void
ReductionSchedule4D::execute(
   RealArray& a_local_dst_array)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("reduction to phase");
   // do sum reduction on phase space nodes
   if (m_is_in_proc_range) {
      if (a_local_dst_array.numberOfDimensions() !=
          m_local_src_array.numberOfDimensions()) {
         OV_ABORT("Local dest array and local src array dims do not match");
      }

      a_local_dst_array = 0;
      m_reduction->applyGlobal(a_local_dst_array,
         m_local_src_array,
         m_reduction_comm,
         ReductionOp::REDUCE_TO_ALL_NODES);
   }
   timers->stopTimer("reduction to phase");
}


//////// PRIVATE METHODS //////////////////////////////////////////////////


void
ReductionSchedule4D::constructCommunicator(
   const tbox::Box& a_local_interior_box,
   const tbox::Box& a_domain_box,
   const MPI_Comm& a_comm)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("reduction to phase");

   // Find the groups of processors that deal with the same piece of
   // configuration space and make a sub-communicator for each group.  The
   // color variable is essentially an offset to the first point of
   // configuration space dealt with by this processor.  So color is the same
   // for all processors dealing with the same piece of configuration space.
   int color(-1);
   if (m_is_in_proc_range) {
      int ibeg(-1);
      for (int d(0); d < a_domain_box.getDim(); ++d) {
         if (!m_collapse_dir[d]) {
            ibeg = d;
         }
      }
      color = a_local_interior_box.lower(ibeg);
      for (int d(ibeg-1); d >= 0; --d) {
         if (!m_collapse_dir[d]) {
            color *= a_domain_box.numberCells(d);
            color += a_local_interior_box.lower(d);
         }
      }

      int comm_id;
      MPI_Comm_rank(a_comm, &comm_id);
      const int status(MPI_Comm_split(a_comm,
         color,
         comm_id,
         &m_reduction_comm));
      if (status != MPI_SUCCESS) {
         OV_ABORT("Configuration space splitting of MPI communicator failed");
      }
   }
   timers->stopTimer("reduction to phase");
}

} // end namespace Loki
