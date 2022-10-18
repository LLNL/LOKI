/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "ReductionSchedule4D.H"

#include "TimerManager.H"
#include "Loki_Utilities.H"

namespace Loki {

const int ReductionSchedule4D::s_MAX_DIM = 4;

ReductionSchedule4D::ReductionSchedule4D(
   const ParallelArray& a_local_src_array,
   const ParallelArray::Box& a_local_interior_box,
   const tbox::Box& a_domain_box,
   const deque<bool>& a_collapse_dir,
   double a_dv,
   int a_proc_lo,
   int a_proc_hi,
   const MPI_Comm& a_comm)
   : m_local_src_array(a_local_src_array),
     m_dv(a_dv),
     m_is_in_proc_range(false),
     m_reduction_comm(a_comm)
{
   m_is_in_proc_range =
      isInRange(Loki_Utilities::s_my_id, a_proc_lo, a_proc_hi);

   // Construct sub-communicators for each group of processors dealing with
   // the same piece of configuration space.
   constructCommunicator(a_local_interior_box,
      a_domain_box,
      a_collapse_dir,
      a_comm);
}


ReductionSchedule4D::~ReductionSchedule4D()
{
}


void
ReductionSchedule4D::execute(
   ParallelArray& a_local_dst_array)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("reduction to phase");
   // do sum reduction on phase space nodes
   if (m_is_in_proc_range) {
      if (a_local_dst_array.dim() != m_local_src_array.dim()) {
         LOKI_ABORT("Local dest array and local src array dims do not match");
      }

      a_local_dst_array = 0;
      Loki_Utilities::getSums(&(*m_local_src_array.getData()),
         &(*a_local_dst_array.getData()),
         m_local_src_array.dataBox().size(),
         REDUCE_TO_ALL_NODES,
         m_reduction_comm);
      a_local_dst_array *= m_dv;
   }
   timers->stopTimer("reduction to phase");
}


//////// PRIVATE METHODS //////////////////////////////////////////////////


void
ReductionSchedule4D::constructCommunicator(
   const ParallelArray::Box& a_local_interior_box,
   const tbox::Box& a_domain_box,
   const deque<bool>& a_collapse_dir,
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
         if (!a_collapse_dir[d]) {
            ibeg = d;
         }
      }
      color = a_local_interior_box.lower(ibeg);
      for (int d(ibeg-1); d >= 0; --d) {
         if (!a_collapse_dir[d]) {
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
         LOKI_ABORT("Configuration space splitting of MPI communicator failed");
      }

      m_commSize = 0;
      MPI_Comm_size(m_reduction_comm, &m_commSize);
   }
   timers->stopTimer("reduction to phase");
}

} // end namespace Loki
