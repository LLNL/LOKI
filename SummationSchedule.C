/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "SummationSchedule.H"

#include "TimerManager.H"
#include "Loki_Utilities.H"

namespace Loki {

const int SummationSchedule::s_TAG = 4380;

SummationSchedule::SummationSchedule(
   const ParallelArray& a_local_src_array,
   const ParallelArray& a_dist_func,
   const tbox::Box& a_domain_box,
   const MPI_Comm& a_comm)
   : m_local_src_array(a_local_src_array),
     m_dist_func(a_dist_func),
     m_local_dst_array(a_local_src_array),
     m_is_in_proc_range(false),
     m_summation_comm(a_comm),
     m_summation_comm_id(-1),
     m_need_communication_pattern(true)
{
   m_is_in_proc_range = isInRange(Loki_Utilities::s_my_id,
      m_dist_func.procLo(),
      m_dist_func.procHi());

   // Construct sub-communicators for each group of processors dealing with
   // the same piece of configuration space.
   constructCommunicator(a_local_src_array.interiorBox(), a_domain_box, a_comm);
}


SummationSchedule::~SummationSchedule()
{
}


void
SummationSchedule::execute(
   ParallelArray& a_global_dst_array)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("summation");

   // Sum data on phase space nodes onto the head nodes.
   m_local_dst_array = 0;
   if (m_is_in_proc_range) {
      if (a_global_dst_array.dim() != m_local_src_array.dim()) {
         LOKI_ABORT("global dst array dim and local src array dim do not match");
      }

      Loki_Utilities::getSums(&(*m_local_src_array.getData()),
         &(*m_local_dst_array.getData()),
         m_local_src_array.dataBox().size(),
         REDUCE_TO_HEAD_NODE,
         m_summation_comm);
   }

   // If this processor owns part of the global array or is a head node then
   // it must be involved in the data communication.  For each head node,
   // transfer data in m_local_dst_array to the corresponding part of
   // a_global_dst_array.
   if (m_is_in_proc_range && isHeadNode()) {
      if (m_need_communication_pattern) {
         computeSendTargets(a_global_dst_array);
         m_need_communication_pattern = false;
      }
      sendHeadProcDataToDestProcs(a_global_dst_array);
   }
   else if (a_global_dst_array.procLo() <= Loki_Utilities::s_my_id &&
            Loki_Utilities::s_my_id <= a_global_dst_array.procHi()){
      if (m_need_communication_pattern) {
         computeRecvTargets(a_global_dst_array);
         m_need_communication_pattern = false;
      }
      destProcsRecvHeadProcData(a_global_dst_array);
   }

   timers->stopTimer("summation");
}


//////// PRIVATE METHODS //////////////////////////////////////////////////


void
SummationSchedule::constructCommunicator(
   const ParallelArray::Box& a_local_interior_box,
   const tbox::Box& a_domain_box,
   const MPI_Comm& a_comm)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("summation");

   // Find the groups of processors that deal with the same piece of
   // configuration space and make a sub-communicator for each group.  The
   // color variable is essentially an offset to the first point of
   // configuration space dealt with by this processor.  So color is the same
   // for all processors dealing with the same piece of configuration space.
   if (m_is_in_proc_range) {
      int color = a_local_interior_box.lower(0) +
         a_domain_box.numberCells(0)*a_local_interior_box.lower(1);

      int comm_id;
      MPI_Comm_rank(a_comm, &comm_id);
      const int status(MPI_Comm_split(a_comm,
         color,
         comm_id,
         &m_summation_comm));
      MPI_Comm_rank(m_summation_comm, &m_summation_comm_id);
      if (status != MPI_SUCCESS) {
         LOKI_ABORT("Configuration space splitting of MPI communicator failed");
      }
   }
   timers->stopTimer("summation");
}


void
SummationSchedule::computeSendTargets(
   ParallelArray& a_global_dst_array)
{
   // Determine which destination targets a head node overlaps.
   const ParallelArray::Box& src_local_box = m_local_src_array.localBox();
   for (int i = a_global_dst_array.procLo();
        i <= a_global_dst_array.procHi();
        ++i) {
      const ParallelArray::Box dst_local_box = a_global_dst_array.localBox(i);
      if (!dst_local_box.intersect(src_local_box).empty()) {
         m_send_targets.push_back(i);
      }
   }
}


void
SummationSchedule::computeRecvTargets(
   ParallelArray& a_global_dst_array)
{
   // Determine which heads nodes overlap a destination target and how much
   // data will be sent by each head node.
   const vector<int>& src_dim_partitions = m_dist_func.getDimPartitions();
   vector<int> src_idx_rank(m_dist_func.dim(), 0);
   const ParallelArray::Box& dest_proc_local_box =
      a_global_dst_array.localBox();
   for (int i = 0; i < src_dim_partitions[0]; ++i) {
      src_idx_rank[0] = i;
      for (int j = 0; j < src_dim_partitions[1]; ++j) {
         src_idx_rank[1] = j;
         int this_head_proc = m_dist_func.getGlobalRank(src_idx_rank);
         const ParallelArray::Box this_head_proc_local_box_unreduced =
            m_dist_func.localBox(this_head_proc);
         int dst_dist_dim = a_global_dst_array.distDim();
         int dst_dim = a_global_dst_array.dim();
         ParallelArray::Box this_head_proc_local_box(dst_dim);
         for (int dim = 0; dim < dst_dist_dim; ++dim) {
            this_head_proc_local_box.lower(dim) =
               this_head_proc_local_box_unreduced.lower(dim);
            this_head_proc_local_box.upper(dim) =
               this_head_proc_local_box_unreduced.upper(dim);
         }
         for (int dim = dst_dist_dim; dim < dst_dim; ++dim) {
            this_head_proc_local_box.lower(dim) =
               a_global_dst_array.localBox().lower(dim);
            this_head_proc_local_box.upper(dim) =
               a_global_dst_array.localBox().upper(dim);
         }
         if (!dest_proc_local_box.intersect(this_head_proc_local_box).empty()) {
            m_recv_targets.push_back(this_head_proc);
            m_recv_local_boxes.push_back(this_head_proc_local_box);
            const ParallelArray::Box this_head_proc_recv_box_unreduced =
               m_dist_func.dataBox(this_head_proc);
            ParallelArray::Box this_head_proc_recv_box(dst_dim);
            for (int dim = 0; dim < dst_dist_dim; ++dim) {
               this_head_proc_recv_box.lower(dim) =
                  this_head_proc_recv_box_unreduced.lower(dim);
               this_head_proc_recv_box.upper(dim) =
                  this_head_proc_recv_box_unreduced.upper(dim);
            }
            for (int dim = dst_dist_dim; dim < dst_dim; ++dim) {
               this_head_proc_recv_box.lower(dim) =
                  a_global_dst_array.localBox().lower(dim);
               this_head_proc_recv_box.upper(dim) =
                  a_global_dst_array.localBox().upper(dim);
            }
            m_recv_boxes.push_back(this_head_proc_recv_box);
         }
      }
   }
}


void
SummationSchedule::sendHeadProcDataToDestProcs(
   ParallelArray& a_global_dst_array)
{
   // Post sends of m_local_dst_array's data to each destination target with
   // which a head node overlaps.
   int num_sends = static_cast<int>(m_send_targets.size());
   vector<MPI_Request> sendReqs(num_sends);
   double* send_buffer = m_local_dst_array.getData();
   int buffer_size = m_local_dst_array.dataBox().size();
   for (int i = 0; i < num_sends; ++i) {
      MPI_Isend(send_buffer,
         buffer_size,
         MPI_DOUBLE,
         m_send_targets[i],
         s_TAG,
         MPI_COMM_WORLD,
         &sendReqs[i]);
   }

   // Verify that all the sends have occurred.
   vector<MPI_Status> sendStatus(num_sends);
   MPI_Waitall(num_sends, &sendReqs[0], &sendStatus[0]);
}


void
SummationSchedule::destProcsRecvHeadProcData(
   ParallelArray& a_global_dst_array)
{
   // Issue the receives for each head node that sends data to this destination
   // target.
   int num_recvs = static_cast<int>(m_recv_targets.size());
   vector<double*> receive_buffers(num_recvs);
   vector<MPI_Request> recv_reqs(num_recvs);
   for (int i = 0; i < num_recvs; ++i) {
      int buffer_size = m_recv_boxes[i].size();
      receive_buffers[i] = new double [buffer_size];
      MPI_Irecv(receive_buffers[i],
         buffer_size,
         MPI_DOUBLE,
         m_recv_targets[i],
         s_TAG,
         MPI_COMM_WORLD,
         &recv_reqs[i]);
   }

   // Process receives.
   const ParallelArray::Box& dest_proc_local_box =
      a_global_dst_array.localBox();
   for (int i = 0; i < num_recvs; ++i) {
      // Find a completed receive.
      int recv_idx;
      MPI_Status stat;
      MPI_Waitany(num_recvs, &recv_reqs[0], &recv_idx, &stat);

      // Get the receive buffer with which this receive is associated.
      double* this_recv_buffer = receive_buffers[recv_idx];
      const ParallelArray::Box& this_recv_box = m_recv_boxes[recv_idx];
      const ParallelArray::Box& this_local_box = m_recv_local_boxes[recv_idx];

      // Fill a_global_dst_array with the received data.
      int i0lo = max(dest_proc_local_box.lower(0), this_local_box.lower(0));
      int i0hi = min(dest_proc_local_box.upper(0), this_local_box.upper(0));
      int i1lo = max(dest_proc_local_box.lower(1), this_local_box.lower(1));
      int i1hi = min(dest_proc_local_box.upper(1), this_local_box.upper(1));
      int n0 = this_recv_box.numberOfCells(0);
      int i0buff = i0lo-this_recv_box.lower(0);
      for (int i0 = i0lo; i0 <= i0hi; ++i0) {
         int i1buff = i1lo-this_recv_box.lower(1);
         for (int i1 = i1lo; i1 <= i1hi; ++i1) {
            a_global_dst_array(i0, i1) = this_recv_buffer[i1buff*n0 + i0buff];
            ++i1buff;
         }
         ++i0buff;
      }
      delete [] this_recv_buffer;
   }
}

} // end namespace Loki
