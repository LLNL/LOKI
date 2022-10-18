/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "ReductionSchedule.H"

#include "TimerManager.H"
#include "Loki_Utilities.H"

namespace Loki {

const int ReductionSchedule::s_MAX_DIM = 4;
const int ReductionSchedule::s_TAG = 4381;

ReductionSchedule::ReductionSchedule(
   const ParallelArray& a_src_array,
   const ProblemDomain& a_domain,
   const deque<bool>& a_collapse_dir,
   int a_nghosts,
   double a_dv,
   double a_weight,
   MPI_Comm& a_comm)
   : m_src_array(a_src_array),
     m_src_interior_box(a_src_array.interiorBox()),
     m_domain_box(a_domain.box()),
     m_collapse_dir(a_collapse_dir),
     m_dv(a_dv),
     m_weight(a_weight),
     m_dst_dim(0),
     m_is_in_proc_range(false),
     m_reduction_comm(a_comm),
     m_reduction_comm_id(-1),
     m_need_communication_pattern(true)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("reduction");

   if (m_src_array.dim() > s_MAX_DIM) {
      LOKI_ABORT("Source array exceeds maximum dimension.");
   }
   m_is_in_proc_range = isInRange(Loki_Utilities::s_my_id,
      m_src_array.procLo(),
      m_src_array.procHi());

   // This is a misnomer as it only computes the dimension of dst array.
   computeSrcAndDstDimensions();

   // Construct sub-communicators for each group of processors dealing with
   // the same piece of configuration space.
   constructCommunicator(a_comm);

   // Create the local array defined on the reduced base space of a_src_array.
   allocateReducedSrcArray(a_nghosts);

   timers->stopTimer("reduction");
}


ReductionSchedule::~ReductionSchedule()
{
}


void
ReductionSchedule::execute(
   ParallelArray& a_dst_array)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("reduction");

   // do sum reduction on phase space nodes
   if (m_is_in_proc_range) {
      // Perform the reduction of dimension.
      reduceDimension();

      // Now sum the reduced results onto the head nodes.
      sumToHeadProcs();

      // Now multiply the summed reduced results on the head processors by m_dv
      // and m_weight.
      if (isHeadNode()) {
         m_reduced_src_array *= m_dv;
         m_reduced_src_array *= m_weight;
      }
   }

   // If this processor owns part of the global array or is a head node then
   // it must be involved in the data communication.  For each head node,
   // transfer data in m_reduced_src_array to the corresponding part of
   // a_dst_array.
   if (m_is_in_proc_range && isHeadNode()) {
      if (m_need_communication_pattern) {
         computeSendTargets(a_dst_array);
         m_need_communication_pattern = false;
      }
      sendHeadProcDataToDestProcs(a_dst_array);
   }
   else if (a_dst_array.procLo() <= Loki_Utilities::s_my_id &&
            Loki_Utilities::s_my_id <= a_dst_array.procHi()) {
      if (m_need_communication_pattern) {
         computeRecvTargets(a_dst_array);
         m_need_communication_pattern = false;
      }
      destProcsRecvHeadProcData(a_dst_array);
   }

   timers->stopTimer("reduction");
}


//////// PRIVATE METHODS //////////////////////////////////////////////////


void
ReductionSchedule::computeSrcAndDstDimensions()
{
   // The destination dimension is the source dimension minus the number of
   // dimensions being collapsed (integrate over).
   int src_dim = m_src_interior_box.dim();
   if (src_dim != m_collapse_dir.size()) {
      LOKI_ABORT("Incongruent dimensions in ReductionSchedule construction");
   }
   for (int dir(0); dir < src_dim; ++dir) {
      if (!m_collapse_dir[dir]) {
         ++m_dst_dim;
         m_kept_dims.push_back(dir);
      }
   }
   if (m_dst_dim == src_dim) {
      LOKI_ABORT("No reduction in ReductionSchedule construction");
   }
}


void
ReductionSchedule::constructCommunicator(
   MPI_Comm& a_comm)
{
   // Find the groups of processors that deal with the same piece of
   // configuration space and make a sub-communicator for each group.  The
   // color variable is essentially an offset to the first point of
   // configuration space dealt with by this processor.  So color is the same
   // for all processors dealing with the same piece of configuration space.
   int color(-1);
   if (m_is_in_proc_range) {
      int ibeg(-1);
      for (int d(0); d < m_domain_box.getDim(); ++d) {
         if (!m_collapse_dir[d]) {
            ibeg = d;
         }
      }
      color = m_src_interior_box.lower(ibeg);
      for (int d(ibeg-1); d >= 0; --d) {
         if (!m_collapse_dir[d]) {
            color *= m_domain_box.numberCells(d);
            color += m_src_interior_box.lower(d);
         }
      }

      int comm_id;
      MPI_Comm_rank(a_comm, &comm_id);
      const int status(MPI_Comm_split(a_comm,
         color,
         comm_id,
         &m_reduction_comm));
      MPI_Comm_rank(m_reduction_comm, &m_reduction_comm_id);
      if (status != MPI_SUCCESS) {
         LOKI_ABORT("Configuration space splitting of MPI communicator failed");
      }
   }
}


void
ReductionSchedule::allocateReducedSrcArray(
   int a_nghosts)
{
   // The reduced source array is a local array defined on the collapsed base
   // space of m_src_array.
   ParallelArray::Box dst_box(m_dst_dim);
   vector<int> num_global_cells(m_dst_dim);
   int dst_dim = 0;
   for (int d = 0; d < m_src_array.dim(); ++d) {
      if (!m_collapse_dir[d]) {
         dst_box.lower(dst_dim) = m_src_interior_box.lower(d);
         dst_box.upper(dst_dim) = m_src_interior_box.upper(d);
         num_global_cells[dst_dim] = m_domain_box.numberCells(d);
         ++dst_dim;
      }
   }

   m_reduced_src_array.partition(dst_box,
      m_dst_dim,
      a_nghosts,
      num_global_cells);
}


void
ReductionSchedule::reduceDimension()
{
   switch (m_src_array.dim()) {
      case 1:
         switch (m_dst_dim) {
            case 0:
               sum_reduce_1d_to_0d();
               break;
         }
         break;
      case 2:
         switch (m_dst_dim) {
            case 0:
               sum_reduce_2d_to_0d();
               break;
            case 1:
               sum_reduce_2d_to_1d();
              break;
         }
         break;
      case 3:
         switch (m_dst_dim) {
            case 0:
               sum_reduce_3d_to_0d();
               break;
            case 1:
               sum_reduce_3d_to_1d();
               break;
            case 2:
               sum_reduce_3d_to_2d();
               break;
         }
         break;
      case 4:
         switch (m_dst_dim) {
            case 0:
               sum_reduce_4d_to_0d();
               break;
            case 1:
               sum_reduce_4d_to_1d();
               break;
            case 2:
               sum_reduce_4d_to_2d();
               break;
            case 3:
               sum_reduce_4d_to_3d();
               break;
         }
      break;
   }
}


void
ReductionSchedule::sum_reduce_1d_to_0d()
{
   m_reduced_src_array = 0;
   int i0lo = m_src_interior_box.lower(0);
   int i0hi = m_src_interior_box.upper(0);
   tbox::IntVector i(tbox::Dimension(1));
   for (i[0] = i0lo; i[0] <= i0hi; ++i[0]) {
      m_reduced_src_array += m_src_array(i[0]);
   }
}


void
ReductionSchedule::sum_reduce_2d_to_0d()
{
   m_reduced_src_array = 0;
   int i0lo = m_src_interior_box.lower(0);
   int i0hi = m_src_interior_box.upper(0);
   int i1lo = m_src_interior_box.lower(1);
   int i1hi = m_src_interior_box.upper(1);
   tbox::IntVector i(tbox::Dimension(2));
   for (i[1] = i1lo; i[1] <= i1hi; ++i[1]) {
      for (i[0] = i0lo; i[0] <= i0hi; ++i[0]) {
         m_reduced_src_array += m_src_array(i[0], i[1]);
      }
   }
}


void
ReductionSchedule::sum_reduce_2d_to_1d()
{
   m_reduced_src_array = 0;
   int i0lo = m_src_interior_box.lower(0);
   int i0hi = m_src_interior_box.upper(0);
   int i1lo = m_src_interior_box.lower(1);
   int i1hi = m_src_interior_box.upper(1);
   tbox::IntVector i(tbox::Dimension(2));
   for (i[1] = i1lo; i[1] <= i1hi; ++i[1]) {
      for (i[0] = i0lo; i[0] <= i0hi; ++i[0]) {
         m_reduced_src_array(i[m_kept_dims[0]]) += m_src_array(i[0], i[1]);
      }
   }
}


void
ReductionSchedule::sum_reduce_3d_to_0d()
{
   m_reduced_src_array = 0;
   int i0lo = m_src_interior_box.lower(0);
   int i0hi = m_src_interior_box.upper(0);
   int i1lo = m_src_interior_box.lower(1);
   int i1hi = m_src_interior_box.upper(1);
   int i2lo = m_src_interior_box.lower(2);
   int i2hi = m_src_interior_box.upper(2);
   tbox::IntVector i(tbox::Dimension(3));
   for (i[2] = i2lo; i[2] <= i2hi; ++i[2]) {
      for (i[1] = i1lo; i[1] <= i1hi; ++i[1]) {
         for (i[0] = i0lo; i[0] <= i0hi; ++i[0]) {
            m_reduced_src_array += m_src_array(i[0], i[1], i[2]);
         }
      }
   }
}


void
ReductionSchedule::sum_reduce_3d_to_1d()
{
   m_reduced_src_array = 0;
   int i0lo = m_src_interior_box.lower(0);
   int i0hi = m_src_interior_box.upper(0);
   int i1lo = m_src_interior_box.lower(1);
   int i1hi = m_src_interior_box.upper(1);
   int i2lo = m_src_interior_box.lower(2);
   int i2hi = m_src_interior_box.upper(2);
   tbox::IntVector i(tbox::Dimension(3));
   for (i[2] = i2lo; i[2] <= i2hi; ++i[2]) {
      for (i[1] = i1lo; i[1] <= i1hi; ++i[1]) {
         for (i[0] = i0lo; i[0] <= i0hi; ++i[0]) {
            m_reduced_src_array(i[m_kept_dims[0]]) +=
               m_src_array(i[0], i[1], i[2]);
         }
      }
   }
}


void
ReductionSchedule::sum_reduce_3d_to_2d()
{
   m_reduced_src_array = 0;
   int i0lo = m_src_interior_box.lower(0);
   int i0hi = m_src_interior_box.upper(0);
   int i1lo = m_src_interior_box.lower(1);
   int i1hi = m_src_interior_box.upper(1);
   int i2lo = m_src_interior_box.lower(2);
   int i2hi = m_src_interior_box.upper(2);
   tbox::IntVector i(tbox::Dimension(3));
   for (i[2] = i2lo; i[2] <= i2hi; ++i[2]) {
      for (i[1] = i1lo; i[1] <= i1hi; ++i[1]) {
         for (i[0] = i0lo; i[0] <= i0hi; ++i[0]) {
            m_reduced_src_array(i[m_kept_dims[0]], i[m_kept_dims[1]]) +=
               m_src_array(i[0], i[1], i[2]);
         }
      }
   }
}


void
ReductionSchedule::sum_reduce_4d_to_0d()
{
   m_reduced_src_array = 0;
   int i0lo = m_src_interior_box.lower(0);
   int i0hi = m_src_interior_box.upper(0);
   int i1lo = m_src_interior_box.lower(1);
   int i1hi = m_src_interior_box.upper(1);
   int i2lo = m_src_interior_box.lower(2);
   int i2hi = m_src_interior_box.upper(2);
   int i3lo = m_src_interior_box.lower(3);
   int i3hi = m_src_interior_box.upper(3);
   tbox::IntVector i(tbox::Dimension(4));
   for (i[3] = i3lo; i[3] <= i3hi; ++i[3]) {
      for (i[2] = i2lo; i[2] <= i2hi; ++i[2]) {
         for (i[1] = i1lo; i[1] <= i1hi; ++i[1]) {
            for (i[0] = i0lo; i[0] <= i0hi; ++i[0]) {
               m_reduced_src_array += m_src_array(i[0], i[1], i[2], i[3]);
            }
         }
      }
   }
}


void
ReductionSchedule::sum_reduce_4d_to_1d()
{
   m_reduced_src_array = 0;
   int i0lo = m_src_interior_box.lower(0);
   int i0hi = m_src_interior_box.upper(0);
   int i1lo = m_src_interior_box.lower(1);
   int i1hi = m_src_interior_box.upper(1);
   int i2lo = m_src_interior_box.lower(2);
   int i2hi = m_src_interior_box.upper(2);
   int i3lo = m_src_interior_box.lower(3);
   int i3hi = m_src_interior_box.upper(3);
   tbox::IntVector i(tbox::Dimension(4));
   for (i[3] = i3lo; i[3] <= i3hi; ++i[3]) {
      for (i[2] = i2lo; i[2] <= i2hi; ++i[2]) {
         for (i[1] = i1lo; i[1] <= i1hi; ++i[1]) {
            for (i[0] = i0lo; i[0] <= i0hi; ++i[0]) {
               m_reduced_src_array(i[m_kept_dims[0]]) +=
                  m_src_array(i[0], i[1], i[2], i[3]);
            }
         }
      }
   }
}


void
ReductionSchedule::sum_reduce_4d_to_2d()
{
   m_reduced_src_array = 0;
   int i0lo = m_src_interior_box.lower(0);
   int i0hi = m_src_interior_box.upper(0);
   int i1lo = m_src_interior_box.lower(1);
   int i1hi = m_src_interior_box.upper(1);
   int i2lo = m_src_interior_box.lower(2);
   int i2hi = m_src_interior_box.upper(2);
   int i3lo = m_src_interior_box.lower(3);
   int i3hi = m_src_interior_box.upper(3);
   tbox::IntVector i(tbox::Dimension(4));
   for (i[3] = i3lo; i[3] <= i3hi; ++i[3]) {
      for (i[2] = i2lo; i[2] <= i2hi; ++i[2]) {
         for (i[1] = i1lo; i[1] <= i1hi; ++i[1]) {
            for (i[0] = i0lo; i[0] <= i0hi; ++i[0]) {
               m_reduced_src_array(i[m_kept_dims[0]], i[m_kept_dims[1]]) +=
                  m_src_array(i[0], i[1], i[2], i[3]);
            }
         }
      }
   }
}


void
ReductionSchedule::sum_reduce_4d_to_3d()
{
   m_reduced_src_array = 0;
   int i0lo = m_src_interior_box.lower(0);
   int i0hi = m_src_interior_box.upper(0);
   int i1lo = m_src_interior_box.lower(1);
   int i1hi = m_src_interior_box.upper(1);
   int i2lo = m_src_interior_box.lower(2);
   int i2hi = m_src_interior_box.upper(2);
   int i3lo = m_src_interior_box.lower(3);
   int i3hi = m_src_interior_box.upper(3);
   tbox::IntVector i(tbox::Dimension(4));
   for (i[3] = i3lo; i[3] <= i3hi; ++i[3]) {
      for (i[2] = i2lo; i[2] <= i2hi; ++i[2]) {
         for (i[1] = i1lo; i[1] <= i1hi; ++i[1]) {
            for (i[0] = i0lo; i[0] <= i0hi; ++i[0]) {
               m_reduced_src_array(i[m_kept_dims[0]], i[m_kept_dims[1]], i[m_kept_dims[2]]) +=
                  m_src_array(i[0], i[1], i[2], i[3]);
            }
         }
      }
   }
}


void
ReductionSchedule::sumToHeadProcs()
{
   if (isHeadNode()) {
      MPI_Reduce(MPI_IN_PLACE,
         m_reduced_src_array.getData(),
         m_reduced_src_array.dataBox().size(),
         MPI_DOUBLE,
         MPI_SUM,
         0,
         m_reduction_comm);
   }
   else {
      MPI_Reduce(m_reduced_src_array.getData(),
         m_reduced_src_array.getData(),
         m_reduced_src_array.dataBox().size(),
         MPI_DOUBLE,
         MPI_SUM,
         0,
         m_reduction_comm);
   }
}


void
ReductionSchedule::computeSendTargets(
   ParallelArray& a_dst_array)
{
   // Determine which destination targets a head node overlaps.
   const ParallelArray::Box& src_local_box = m_reduced_src_array.localBox();
   for (int i = a_dst_array.procLo(); i <= a_dst_array.procHi(); ++i) {
      const ParallelArray::Box dst_local_box = a_dst_array.localBox(i);
      if (!dst_local_box.intersect(src_local_box).empty()) {
         m_send_targets.push_back(i);
      }
   }
}


void
ReductionSchedule::computeRecvTargets(
   ParallelArray& a_dst_array)
{
   // Determine which heads nodes overlap a destination target and how much
   // data will be sent by each head node.
   const vector<int>& src_dim_partitions = m_src_array.getDimPartitions();
   vector<int> ihi(s_MAX_DIM, 1);
   for (int i = 0; i < m_dst_dim; ++i) {
      ihi[m_kept_dims[i]] = src_dim_partitions[m_kept_dims[i]];
   }
   vector<int> src_idx_rank(s_MAX_DIM, 0);
   const ParallelArray::Box& dest_proc_local_box = a_dst_array.localBox();
   for (int i = 0; i < ihi[0]; ++i) {
      src_idx_rank[0] = i;
      for (int j = 0; j < ihi[1]; ++j) {
         src_idx_rank[1] = j;
         for (int k = 0; k < ihi[2]; ++k) {
            src_idx_rank[2] = k;
            for (int l = 0; l < ihi[3]; ++l) {
               src_idx_rank[3] = l;
               int this_head_proc = m_src_array.getGlobalRank(src_idx_rank);
               const ParallelArray::Box this_head_proc_local_box_unreduced =
                  m_src_array.localBox(this_head_proc);
               ParallelArray::Box this_head_proc_local_box(m_dst_dim);
               for (int dim = 0; dim < m_dst_dim; ++dim) {
                  this_head_proc_local_box.lower(dim) =
                     this_head_proc_local_box_unreduced.lower(m_kept_dims[dim]);
                  this_head_proc_local_box.upper(dim) =
                     this_head_proc_local_box_unreduced.upper(m_kept_dims[dim]);
               }
               if (!dest_proc_local_box.intersect(this_head_proc_local_box).empty()) {
                  m_recv_targets.push_back(this_head_proc);
                  m_recv_local_boxes.push_back(this_head_proc_local_box);
                  const ParallelArray::Box this_head_proc_recv_box_unreduced =
                     m_src_array.dataBox(this_head_proc);
                  ParallelArray::Box this_head_proc_recv_box(m_dst_dim);
                  for (int dim = 0; dim < m_dst_dim; ++dim) {
                     this_head_proc_recv_box.lower(dim) =
                        this_head_proc_recv_box_unreduced.lower(m_kept_dims[dim]);
                     this_head_proc_recv_box.upper(dim) =
                        this_head_proc_recv_box_unreduced.upper(m_kept_dims[dim]);
                  }
                  m_recv_boxes.push_back(this_head_proc_recv_box);
               }
            }
         }
      }
   }
}


void
ReductionSchedule::sendHeadProcDataToDestProcs(
   ParallelArray& a_dst_array)
{
   // Post sends of m_reduced_src_array's data to each destination target with
   // which a head node overlaps.
   int num_sends = static_cast<int>(m_send_targets.size());
   vector<MPI_Request> sendReqs(num_sends);
   double* send_buffer = m_reduced_src_array.getData();
   int buffer_size = m_reduced_src_array.dataBox().size();
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
ReductionSchedule::destProcsRecvHeadProcData(
   ParallelArray& a_dst_array)
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
   const ParallelArray::Box& dest_proc_local_box = a_dst_array.localBox();
   for (int i = 0; i < num_recvs; ++i) {
      // Find a completed receive.
      int recv_idx;
      MPI_Status stat;
      MPI_Waitany(num_recvs, &recv_reqs[0], &recv_idx, &stat);

      // Get the receive buffer with which this receive is associated.
      double* this_recv_buffer = receive_buffers[recv_idx];
      const ParallelArray::Box& this_recv_box = m_recv_boxes[recv_idx];
      const ParallelArray::Box& this_local_box = m_recv_local_boxes[recv_idx];

      // Fill a_dst_array with the received data.
      if (m_dst_dim == 0) {
         a_dst_array += this_recv_buffer[0];
      }
      else if (m_dst_dim == 1) {
         int i0lo = max(dest_proc_local_box.lower(0), this_local_box.lower(0));
         int i0hi = min(dest_proc_local_box.upper(0), this_local_box.upper(0));
         int i0buff = i0lo-this_recv_box.lower(0);
         for (int i0 = i0lo; i0 <= i0hi; ++i0) {
            a_dst_array(i0) = this_recv_buffer[i0buff];
            ++i0buff;
         }
      }
      else if (m_dst_dim == 2) {
         int i0lo = max(dest_proc_local_box.lower(0), this_local_box.lower(0));
         int i0hi = min(dest_proc_local_box.upper(0), this_local_box.upper(0));
         int i1lo = max(dest_proc_local_box.lower(1), this_local_box.lower(1));
         int i1hi = min(dest_proc_local_box.upper(1), this_local_box.upper(1));
         int n0 = this_recv_box.numberOfCells(0);
         int i0buff = i0lo-this_recv_box.lower(0);
         for (int i0 = i0lo; i0 <= i0hi; ++i0) {
            int i1buff = i1lo-this_recv_box.lower(1);
            for (int i1 = i1lo; i1 <= i1hi; ++i1) {
               a_dst_array(i0, i1) = this_recv_buffer[i1buff*n0 + i0buff];
               ++i1buff;
            }
            ++i0buff;
         }
      }
      else if (m_dst_dim == 3) {
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
                  a_dst_array(i0, i1, i2) =
                     this_recv_buffer[((i2buff*n1)+i1buff)*n0+i0buff];
                  ++i2buff;
               }
               ++i1buff;
            }
            ++i0buff;
         }
      }
      delete [] this_recv_buffer;
   }
}

} // end namespace Loki
