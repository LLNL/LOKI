/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "ContractionSchedule.H"
#include "Loki_Utilities.H"
#include "TimerManager.H"

namespace Loki {

const int ContractionSchedule::s_TAG = 4383;

ContractionSchedule::ContractionSchedule(
   const ParallelArray& a_global_array,
   const ParallelArray& a_local_src_array,
   const tbox::Box& a_domain_box,
   const MPI_Comm& a_comm)
   : m_global_array(a_global_array),
     m_local_src_array(a_local_src_array),
     m_is_in_proc_range(false),
     m_2D_comm_id(-1),
     m_need_communication_pattern(true)
{
   if (m_local_src_array.dim() != 2 && m_local_src_array.dim() != 3) {
      LOKI_ABORT("Attemping to use ContractionSchedule on arrays that are not 2D or 3D");
   }

   m_is_in_proc_range = m_global_array.procLo() <= Loki_Utilities::s_my_id &&
      Loki_Utilities::s_my_id <= m_global_array.procHi();

   // Construct communicator for the 4D processors that deal with the same
   // piece of configuration space.  We don't need the communicators per se but
   // we need to know which processor on each communicator is the "head node",
   // the one processor with rank 0 on each communicator.  These are the
   // processors that send data.
   constructCommunicator(m_global_array.interiorBox(), a_domain_box, a_comm);
}


ContractionSchedule::~ContractionSchedule()
{
}


void
ContractionSchedule::execute(
   ParallelArray& a_global_dst_array)
{
   if (m_local_src_array.dim() != a_global_dst_array.dim()) {
      LOKI_ABORT("Source and destination arrays have different dimensions.");
   }

   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("contraction");

   // If this processor is a configuration space head processor then it sends
   // its part of m_local_src_array to the configuration space processors.  The
   // configuration space processors receive this data and put it in
   // a_global_dst_array.
   int config_proc_lo = a_global_dst_array.procLo();
   int config_proc_hi = a_global_dst_array.procHi();
   if (m_is_in_proc_range && isHeadNode()) {
      if (m_need_communication_pattern) {
         computeSendTargets(a_global_dst_array);
         m_need_communication_pattern = false;
      }
      sendPhaseSpaceDataToConfigSpaceProcs(config_proc_lo, config_proc_hi);
   }
   else if (config_proc_lo <= Loki_Utilities::s_my_id &&
            Loki_Utilities::s_my_id <= config_proc_hi) {
      if (m_need_communication_pattern) {
         computeRecvTargets(a_global_dst_array);
         m_need_communication_pattern = false;
      }
      configSpaceRecvPhaseSpaceData(a_global_dst_array);
   }

   timers->stopTimer("contraction");
}


//////// PRIVATE METHODS //////////////////////////////////////////////////


void
ContractionSchedule::constructCommunicator(
   const ParallelArray::Box& a_global_interior_box,
   const tbox::Box& a_domain_box,
   const MPI_Comm& a_comm)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("contraction");

   // Find the groups of processors that deal with the same piece of
   // configuration space.  Make a sub-communicator for each group.  The process
   // with rank 0 on each sub-communicator is the "head node" and the one that
   // communicates that piece of the configuration space data from the local to
   // the global array.
   int color(-1);
   if (m_is_in_proc_range) {
      color = a_global_interior_box.lower(1)*a_domain_box.numberCells(0) +
              a_global_interior_box.lower(0);
      int comm_id;
      MPI_Comm_rank(a_comm, &comm_id);
      MPI_Comm config_space_comm;
      const int status(MPI_Comm_split(a_comm,
         color,
         comm_id,
         &config_space_comm));
      if (status != MPI_SUCCESS) {
         LOKI_ABORT("Configuration space splitting of MPI communicator failed");
      }
      MPI_Comm_rank(config_space_comm, &m_2D_comm_id);
   }
   timers->stopTimer("contraction");
}


void
ContractionSchedule::computeSendTargets(
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
ContractionSchedule::computeRecvTargets(
   ParallelArray& a_global_dst_array)
{
   // Determine which heads nodes overlap a destination target and how much
   // data will be sent by each head node. 
   const vector<int>& src_dim_partitions = m_global_array.getDimPartitions();
   vector<int> src_idx_rank(m_global_array.dim(), 0);
   const ParallelArray::Box& dest_proc_local_box =
      a_global_dst_array.localBox();
   for (int i = 0; i < src_dim_partitions[0]; ++i) {
      src_idx_rank[0] = i;
      for (int j = 0; j < src_dim_partitions[1]; ++j) {
         src_idx_rank[1] = j;
         int this_head_proc = m_global_array.getGlobalRank(src_idx_rank);
         const ParallelArray::Box this_head_proc_local_box_unreduced =
            m_global_array.localBox(this_head_proc);
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
               m_global_array.dataBox(this_head_proc);
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
ContractionSchedule::sendPhaseSpaceDataToConfigSpaceProcs(
   int a_config_proc_lo,
   int a_config_proc_hi)
{
   // Post sends of m_local_src_array's data to each destination target with
   // which a head node overlaps.
   int num_sends = static_cast<int>(m_send_targets.size());
   vector<MPI_Request> send_reqs(num_sends);
   const double* send_buffer = m_local_src_array.getData();
   int buffer_size = m_local_src_array.dataBox().size();
   for (int i = 0; i < num_sends; ++i) {
      MPI_Isend(send_buffer,
         buffer_size,
         MPI_DOUBLE,
         m_send_targets[i],
         s_TAG,
         MPI_COMM_WORLD,
         &send_reqs[i]);
   }

   // Verify that all the sends have occurred.
   vector<MPI_Status> send_status(num_sends);
   MPI_Waitall(num_sends, &send_reqs[0], &send_status[0]);
}


void
ContractionSchedule::configSpaceRecvPhaseSpaceData(
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
   
      // Get the receive bufffer with which the receive is associated.
      double* this_recv_buffer = receive_buffers[recv_idx];
      ParallelArray::Box& this_recv_box = m_recv_boxes[recv_idx];
      const ParallelArray::Box& this_local_box = m_recv_local_boxes[recv_idx];
   
      // Fill a_global_dst_array with the received data.
      if (this_recv_box.dim() == 2) {
         int i0lo = max(dest_proc_local_box.lower(0), this_local_box.lower(0));
         int i0hi = min(dest_proc_local_box.upper(0), this_local_box.upper(0));
         int i1lo = max(dest_proc_local_box.lower(1), this_local_box.lower(1));
         int i1hi = min(dest_proc_local_box.upper(1), this_local_box.upper(1));
         int n0 = this_recv_box.numberOfCells(0);
         int i0buff = i0lo-this_recv_box.lower(0);
         for (int i0 = i0lo; i0 <= i0hi; ++i0) {
            int i1buff = i1lo-this_recv_box.lower(1);
            for (int i1 = i1lo; i1 <= i1hi; ++i1) {
               a_global_dst_array(i0, i1) =
                  this_recv_buffer[i1buff*n0 + i0buff];
               ++i1buff;
            }
            ++i0buff;
         }
      }
      else {
         int i0lo = max(dest_proc_local_box.lower(0), this_local_box.lower(0));
         int i0hi = min(dest_proc_local_box.upper(0), this_local_box.upper(0));
         int i1lo = max(dest_proc_local_box.lower(1), this_local_box.lower(1));
         int i1hi = min(dest_proc_local_box.upper(1), this_local_box.upper(1));
         int i2lo = max(dest_proc_local_box.lower(2), this_local_box.lower(2));
         int i2hi = min(dest_proc_local_box.upper(2), this_local_box.upper(2));
         int n0 = this_recv_box.numberOfCells(0);
         int n1 = this_recv_box.numberOfCells(1);
         int i0buff = i0lo-this_recv_box.lower(0);
         for (int i0 = i0lo; i0 <= i0hi; ++i0) {
            int i1buff = i1lo-this_recv_box.lower(1);
            for (int i1 = i1lo; i1 <= i1hi; ++i1) {
               int i2buff = i2lo-this_recv_box.lower(2);
               for (int i2 = i2lo; i2 <= i2hi; ++i2) {
                  a_global_dst_array(i0, i1, i2) =
                     this_recv_buffer[i2buff*n0*n1 + i1buff*n0 + i0buff];
                  ++i2buff;
               }
               ++i1buff;
            }
            ++i0buff;
         }
      }

      // Delete the receive buffer.
      delete [] this_recv_buffer;
   }
}

} // end namespace Loki
