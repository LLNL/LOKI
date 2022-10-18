/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "ExpansionSchedule.H"
#include "Directions.H"
#include "Loki_Utilities.H"
#include "TimerManager.H"

namespace Loki {

const int ExpansionSchedule::s_TAG = 4382;

ExpansionSchedule::ExpansionSchedule(
   const ParallelArray& a_phase_space_dist,
   const ProblemDomain& a_domain,
   MPI_Comm a_comm)
   : m_phase_space_dist(a_phase_space_dist),
     m_need_communication_pattern(true)
{
   if (m_phase_space_dist.procLo() <= Loki_Utilities::s_my_id &&
       Loki_Utilities::s_my_id <= m_phase_space_dist.procHi()) {
      constructExpansionCommunicator(a_domain, a_comm);
   }
}


ExpansionSchedule::~ExpansionSchedule()
{
}


void
ExpansionSchedule::execute(
   const ParallelArray& a_global_config_array,
   ParallelArray& a_local_phase_array)
{
   if (a_global_config_array.dim() != a_local_phase_array.dim()) {
      LOKI_ABORT("Source and destination arrays have different dimensions.");
   }
   if (a_global_config_array.dim() != 2 && a_global_config_array.dim() != 3) {
      LOKI_ABORT("Attemping to use ExpansionSchedule on arrays that are not 2D or 3D");
   }

   // If this processor owns part of the global configuration space array or is
   // a phase space processor it must be involved in the data communication.
   // Transfer the data in a_global_config_array to the corresponding part of
   // a_local_phase_array.
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("phys to phase");

   // The processors owning the global configuration space array send data to
   // the head processors owning the local phase space arrays.
   if (a_global_config_array.procLo() <= Loki_Utilities::s_my_id &&
       Loki_Utilities::s_my_id <= a_global_config_array.procHi()) {
      if (m_need_communication_pattern) {
         computeHeadProcTargets(a_global_config_array);
         m_need_communication_pattern = false;
      }
      sendConfigSpaceDataToHeadProcTargets(a_global_config_array);
   }
   else if (m_phase_space_dist.procLo() <= Loki_Utilities::s_my_id &&
            Loki_Utilities::s_my_id <= m_phase_space_dist.procHi()) {
      if (m_expansion_comm_id == 0) {
         if (m_need_communication_pattern) {
            computeConfigRecvTargets(a_global_config_array,
               a_local_phase_array);
            m_need_communication_pattern = false;
         }
         headProcsRecvConfigSpaceData(a_global_config_array,
            a_local_phase_array);
      }
   }

   // The head processors owning the local phase space arrays send data to the
   // non-head processors owning the local phase space arrays.
   if (m_phase_space_dist.procLo() <= Loki_Utilities::s_my_id &&
       Loki_Utilities::s_my_id <= m_phase_space_dist.procHi()) {
      MPI_Bcast(a_local_phase_array.getData(),
         a_local_phase_array.dataBox().size(),
         MPI_DOUBLE,
         0,
         m_expansion_comm);
   }

   timers->stopTimer("phys to phase");
}


void
ExpansionSchedule::computeHeadProcTargets(
   const ParallelArray& a_global_config_array)
{
   int dist_dim = a_global_config_array.distDim();
   int dim = a_global_config_array.dim();
   const ParallelArray::Box config_space_box = a_global_config_array.localBox();
   const vector<int>& dim_partitions = m_phase_space_dist.getDimPartitions();
   int num_vel_space_partitions = dim_partitions[2]*dim_partitions[3];
   int this_phase_space_proc = m_phase_space_dist.procLo();
   for (int i = 0; i < dim_partitions[0]*dim_partitions[1]; ++i) {
      const ParallelArray::Box phase_space_box =
         m_phase_space_dist.localBox(this_phase_space_proc);
      ParallelArray::Box phase_space_box_reduced(dim);
      for (int d = 0; d < dist_dim; ++d) {
         phase_space_box_reduced.lower(d) = phase_space_box.lower(d);
         phase_space_box_reduced.upper(d) = phase_space_box.upper(d);
      }
      for (int d = dist_dim; d < dim; ++d) {
         phase_space_box_reduced.lower(d) = config_space_box.lower(d);
         phase_space_box_reduced.upper(d) = config_space_box.upper(d);
      }
      const ParallelArray::Box send_box = config_space_box.intersect(
         phase_space_box_reduced);
      if (!send_box.empty()) {
         m_send_boxes.push_back(send_box);
         m_send_targets.push_back(this_phase_space_proc);
      }
      this_phase_space_proc += num_vel_space_partitions;
   }
}


void
ExpansionSchedule::computeConfigRecvTargets(
   const ParallelArray& a_global_config_array,
   const ParallelArray& a_local_phase_array)
{
   for (int i = a_global_config_array.procLo();
        i <= a_global_config_array.procHi(); ++ i) {
      ParallelArray::Box recv_box = a_local_phase_array.localBox().intersect(
         a_global_config_array.localBox(i));
      if (!recv_box.empty()) {
         m_recv_targets.push_back(i);
         m_recv_boxes.push_back(recv_box);
      }
   }
}


void
ExpansionSchedule::sendConfigSpaceDataToHeadProcTargets(
   const ParallelArray& a_global_config_array)
{
   // Post sends.
   int num_sends = static_cast<int>(m_send_targets.size());
   vector<MPI_Request> sendReqs(num_sends);
   vector<double*> send_buffers(num_sends);
   int which_buffer = 0;
   int which_send = 0;
   for (int i = 0; i < num_sends; ++i) {
      const ParallelArray::Box& this_send_box = m_send_boxes[i];
      int this_buffer_size = this_send_box.size();
      double* this_send_buffer = new double [this_buffer_size];
      if (a_global_config_array.dim() == 2) {
         int buf_idx = 0;
         for (int y = this_send_box.lower(1);
              y <= this_send_box.upper(1); ++y) {
            for (int x = this_send_box.lower(0);
                 x <= this_send_box.upper(0); ++x) {
               this_send_buffer[buf_idx++] = a_global_config_array(x, y);
            }
         }
      }
      else {
         int buf_idx = 0;
         for (int comp = this_send_box.lower(2);
              comp <= this_send_box.upper(2); ++comp) {
            for (int y = this_send_box.lower(1);
                 y <= this_send_box.upper(1); ++y) {
               for (int x = this_send_box.lower(0);
                    x <= this_send_box.upper(0); ++x) {
                  this_send_buffer[buf_idx++] =
                     a_global_config_array(x, y, comp);
               }
            }
         }
      }
      send_buffers[i] = this_send_buffer;
      MPI_Isend(this_send_buffer,
         this_buffer_size,
         MPI_DOUBLE,
         m_send_targets[i],
         s_TAG,
         MPI_COMM_WORLD,
         &sendReqs[which_send++]);
   }

   // Verify that all the sends have occurred.
   for (int dest = 0; dest < num_sends; ++dest) {
      int send_idx;
      MPI_Status stat;
      MPI_Waitany(num_sends, &sendReqs[0], &send_idx, &stat);
      delete [] send_buffers[send_idx];
   }
}


void
ExpansionSchedule::headProcsRecvConfigSpaceData(
   const ParallelArray& a_global_config_array,
   ParallelArray& a_local_phase_array)
{
   // Post receives.
   int num_recvs = static_cast<int>(m_recv_targets.size());
   vector<MPI_Request> recvReqs(num_recvs);
   vector<double*> recv_buffers(num_recvs);
   for (int src = 0; src < num_recvs; ++src) {
      int buffer_size = m_recv_boxes[src].size();
      double* this_recv_buffer = new double [buffer_size];
      recv_buffers[src] = this_recv_buffer;
      MPI_Irecv(this_recv_buffer,
         buffer_size,
         MPI_DOUBLE,
         m_recv_targets[src],
         s_TAG,
         MPI_COMM_WORLD,
         &recvReqs[src]);
   }

   // Process receives as they complete.
   for (int i = 0; i < num_recvs; ++i) {
      // Find a completed receive.
      int recv_idx;
      MPI_Status stat;
      MPI_Waitany(num_recvs, &recvReqs[0], &recv_idx, &stat);

      // Get the receive bufffer with which the receive is associated.
      double* this_recv_buffer = recv_buffers[recv_idx];
      const ParallelArray::Box& this_recv_box = m_recv_boxes[recv_idx];

      // Fill a_local_phase_array with the received data.
      if (this_recv_box.dim() == 2) {
         int buff_idx = 0;
         for (int j = this_recv_box.lower(1);
              j <= this_recv_box.upper(1); ++j) {
            for (int i = this_recv_box.lower(0);
                 i <= this_recv_box.upper(0); ++i) {
               a_local_phase_array(i, j) = this_recv_buffer[buff_idx++];
            }
         }
      }
      else {
         int buff_idx = 0;
         for (int k = this_recv_box.lower(2);
              k <= this_recv_box.upper(2); ++k) {
            for (int j = this_recv_box.lower(1);
                 j <= this_recv_box.upper(1); ++j) {
               for (int i = this_recv_box.lower(0);
                    i <= this_recv_box.upper(0); ++i) {
                  a_local_phase_array(i, j, k) = this_recv_buffer[buff_idx++];
               }
            }
         }
      }

      // Delete the receive buffer.
      delete [] this_recv_buffer;
   }
}


void
ExpansionSchedule::constructExpansionCommunicator(
   const ProblemDomain& a_domain,
   MPI_Comm a_comm)
{
   const ParallelArray::Box& interior_box = m_phase_space_dist.interiorBox();
   int config_space_id =
      interior_box.lower(X2)*a_domain.box().numberCells(X1)+
      interior_box.lower(X1);
   int comm_id;
   MPI_Comm_rank(a_comm, &comm_id);
   const int status = MPI_Comm_split(a_comm,
      config_space_id,
      comm_id,
      &m_expansion_comm);
   if (status != MPI_SUCCESS) {
      LOKI_ABORT("Configuration space splitting of MPI communicator failed");
   }
   MPI_Comm_rank(m_expansion_comm, &m_expansion_comm_id);
}

} // end namespace Loki
