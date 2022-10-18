/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _CONTRACTION_SCHEDULE_H_
#define _CONTRACTION_SCHEDULE_H_

#include "ParallelArray.H"
#include "tbox/Box.H"

#include <mpi.h>

namespace Loki {

/**
 * Performs the communication and other operations necessary to communicate a
 * user supplied local 2D array from the 4D processors to the 2D processor(s).
 * Pretty much the inverse of ExpansionSchedule except that only the head 4D
 * processors are involved.
 */
class ContractionSchedule
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_global_array A global 4D array which describes the
    *                           decomposition of the 4D processors.
    * @param[in] a_local_src_array The local 2D array to be contracted.
    * @param[in] a_domain_box Span of the problem domain.
    * @param[in] a_comm Communicator corresponding to processor_range
    */
   ContractionSchedule(
      const ParallelArray& a_global_array,
      const ParallelArray& a_local_src_array,
      const tbox::Box& a_domain_box,
      const MPI_Comm& a_comm);

   /**
    * @brief Destructor.
    */
   ~ContractionSchedule();

   /**
    * @brief Perform the contraction.
    *
    * @param[in] a_global_dst_array The global array for the contraction.
    */
   void
   execute(
      ParallelArray& a_global_dst_array);

private:
   // Unimplemented default constructor.
   ContractionSchedule();

   // Unimplemented copy constructor.
   ContractionSchedule(
      const ContractionSchedule& other);

   // Unimplemented assignment operator.
   ContractionSchedule&
   operator = (
      const ContractionSchedule& rhs);

   // Construct communicators for each set of processors that deal with the
   // same piece of 2D configuration space.
   void
   constructCommunicator(
      const ParallelArray::Box& a_global_interior_box,
      const tbox::Box& a_domain_box,
      const MPI_Comm& a_comm);

   // Precomputes the destination processors to send data to.
   void
   computeSendTargets(
      ParallelArray& a_global_dst_array);

   // Precomputes the head nodes to receive data from.
   void
   computeRecvTargets(
      ParallelArray& a_global_dst_array);

   // Sends each head node's data to each configuraion space processor.
   void
   sendPhaseSpaceDataToConfigSpaceProcs(
      int a_config_proc_lo,
      int a_config_proc_hi);

   // Each configuration space processor receives data from each head node.
   void
   configSpaceRecvPhaseSpaceData(
      ParallelArray& a_global_dst_array);

   // Returns true if this is a head node.
   bool
   isHeadNode()
   {
      return m_2D_comm_id == 0;
   }

   // The global 4D array which describes the decomposition of the 4D
   // processors.
   const ParallelArray& m_global_array;

   // Local source of contraction.
   const ParallelArray& m_local_src_array;

   bool m_is_in_proc_range;

   // Processor id of groups of 4D processors all dealing with same piece of 2D
   // configuration space.
   int m_2D_comm_id;

   // If true, the communication pattern must be precompute.
   bool m_need_communication_pattern;

   // The processors that head nodes send data to.
   vector<int> m_send_targets;

   // The head nodes that the reduction destination receives data from.
   vector<int> m_recv_targets;

   // The box describing the local part of the data sent from each head node.
   vector<ParallelArray::Box> m_recv_local_boxes;

   // The box describing the data sent from each head node.
   vector<ParallelArray::Box> m_recv_boxes;

   // The global ranks of the sending processes.
   vector<int> m_src_heads;

   // Tag for send of phase space 2D data to configuration space processors.
   static const int s_TAG;
};

} // end namespace Loki

#endif
