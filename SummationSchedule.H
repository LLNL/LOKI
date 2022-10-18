/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _SUMMATION_SCHEDULE_H_
#define _SUMMATION_SCHEDULE_H_

#include "ParallelArray.H"
#include "Directions.H"
#include "tbox/Box.H"
#include <mpi.h>

namespace Loki {

/**
 * Performs the communication and other operations necessary to accomplish a
 * summation of a user supplied local array into a global array.  This differs
 * from a ReductionSchedule in that the arrays are already reduced.  All we
 * want to do is to sum onto the head nodes the data across all the 4D
 * processors dealing with a common piece of configuration space.  Then
 * communicate that data to a global 2D destination.
 */
class SummationSchedule
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_local_src_array The local array to be summed.
    * @param[in] a_dist_func A global 4D array which describes the decomposition
    *                        of the 4D processors.
    * @param[in] a_domain_box Global problem domain box.
    * @param[in] a_comm Communicator corresponding to processor_range.
    */
   SummationSchedule(
      const ParallelArray& a_local_src_array,
      const ParallelArray& a_dist_func,
      const tbox::Box& a_domain_box,
      const MPI_Comm& a_comm);

   /**
    * @brief Destructor.
    */
   ~SummationSchedule();

   /**
    * @brief Perform the summation.
    *
    * @param[in] a_global_dst_array The global array for the summation.
    */
   void
   execute(
      ParallelArray& a_global_dst_array);

private:
   // Unimplemented default constructor.
   SummationSchedule();

   // Unimplemented copy constructor.
   SummationSchedule(
      const SummationSchedule& other);

   // Unimplemented assignment operator.
   SummationSchedule&
   operator = (
      const SummationSchedule& rhs);

   // Construct communicator over which summation takes place.
   void
   constructCommunicator(
      const ParallelArray::Box& a_local_interior_box,
      const tbox::Box& a_domain_box,
      const MPI_Comm& a_comm);

   // Returns true if a_proc_id is in the range of processors
   // a_processor_range.
   static
   bool
   isInRange(
      int a_proc_id,
      int a_proc_lo,
      int a_proc_hi)
   {
      return ((a_proc_id >= a_proc_lo) && (a_proc_id <= a_proc_hi));
   }

   // Returns true if this is a head node.
   bool
   isHeadNode()
   {
      return m_summation_comm_id == 0;
   }

   // Precomputes the destination processors to send data to.
   void
   computeSendTargets(
      ParallelArray& a_global_dst_array);

   // Precomputes the head nodes to receive data from.
   void
   computeRecvTargets(
      ParallelArray& a_global_dst_array);

   // Send the summed data from the head processors to the destination
   // processors.
   void
   sendHeadProcDataToDestProcs(
      ParallelArray& a_global_dst_array);

   // The destination processors receive the summed data from the head
   // processors.
   void
   destProcsRecvHeadProcData(
      ParallelArray& a_global_dst_array);

   // Local source of summation.
   const ParallelArray& m_local_src_array;

   // Distribution function for description of decomposition of species
   // summation is for.
   const ParallelArray& m_dist_func;

   // Local result of summation.
   ParallelArray m_local_dst_array;

   // True if this processor is involved in the reduction.
   bool m_is_in_proc_range;

   // The communicator containing all processors dealing with the same piece of
   // configuration space.
   MPI_Comm m_summation_comm;

   // The processor ID on m_reduction_comm.
   int m_summation_comm_id;

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

   // Tag for send of source array to destination array.
   static const int s_TAG;

   // Indicates if the global reduction should only be sent to the head or all
   // nodes on the communicator.
   enum ReductionTargets {
     REDUCE_TO_ALL_NODES = -1,
     REDUCE_TO_HEAD_NODE
   };
};

} // end namespace Loki

#endif
