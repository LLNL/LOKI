/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _REDUCTION_SCHEDULE_H_
#define _REDUCTION_SCHEDULE_H_

#include "ParallelArray.H"
#include "Directions.H"
#include "ProblemDomain.H"
#include "tbox/Dimension.H"
#include <deque>

namespace Loki {

/**
 * Performs the communication and other operations necessary to accomplish a
 * user defined reduction of dimension on a user supplied local array.
 */
class ReductionSchedule
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_src_array The array to be reduced.
    * @param[in] a_domain Global problem domain.
    * @param[in] a_collapse_dir True in each dimension to be reduced.
    * @param[in] a_nghosts Number of ghosts in a_src_array.
    * @param[in] a_dv Integration volume differential.
    * @param[in] a_weight Weight to apply to reduction.
    * @param[in] a_comm Communicator corresponding to processor_range
    */
   ReductionSchedule(
      const ParallelArray& a_src_array,
      const ProblemDomain& a_domain,
      const deque<bool>& a_collapse_dir,
      int a_nghosts,
      double a_dv,
      double a_weight,
      MPI_Comm& a_comm);

   /**
    * @brief Destructor.
    */
   ~ReductionSchedule();

   /**
    * @brief Perform the reduction.
    *
    * @param[in] a_dst_array The result of the reduction.
    */
   void
   execute(
      ParallelArray& a_dst_array);

private:
   // Unimplemented default constructor.
   ReductionSchedule();

   // Unimplemented copy constructor.
   ReductionSchedule(
      const ReductionSchedule& other);

   // Unimplemented assignment operator.
   ReductionSchedule&
   operator = (
      const ReductionSchedule& rhs);

   // Computes dimension of source and destination of reduction.
   void
   computeSrcAndDstDimensions();

   // Allocates the result of the dimension reduction of m_src_array.
   void
   allocateReducedSrcArray(
      int a_nghosts);

   // Construct communicator over which reduction takes place.
   void
   constructCommunicator(
      MPI_Comm& a_comm);

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
      return m_reduction_comm_id == 0;
   }

   // Looks at src and dst dimensions and calls appropriate sum_reduce* method
   // to reduce dimension.
   void
   reduceDimension();

   // Perform diimension reduction from 1 to 0 dimensions.
   void
   sum_reduce_1d_to_0d();

   // Perform diimension reduction from 2 to 0 dimensions.
   void
   sum_reduce_2d_to_0d();

   // Perform diimension reduction from 2 to 1 dimensions.
   void
   sum_reduce_2d_to_1d();

   // Perform diimension reduction from 3 to 0 dimensions.
   void
   sum_reduce_3d_to_0d();

   // Perform diimension reduction from 3 to 1 dimensions.
   void
   sum_reduce_3d_to_1d();

   // Perform diimension reduction from 3 to 2 dimensions.
   void
   sum_reduce_3d_to_2d();

   // Perform diimension reduction from 4 to 0 dimensions.
   void
   sum_reduce_4d_to_0d();

   // Perform diimension reduction from 4 to 1 dimensions.
   void
   sum_reduce_4d_to_1d();

   // Perform diimension reduction from 4 to 2 dimensions.
   void
   sum_reduce_4d_to_2d();

   // Perform diimension reduction from 4 to 3 dimensions.
   void
   sum_reduce_4d_to_3d();

   // Sum the reduced array to the head processors for the reduced dimension.
   void
   sumToHeadProcs();

   // Precomputes the destination processors to send data to.
   void
   computeSendTargets(
      ParallelArray& a_dst_array);

   // Precomputes the head nodes to receive data from.
   void
   computeRecvTargets(
      ParallelArray& a_dst_array);

   // Send the summed, reduced data from the head processors to the destination
   // processors.
   void
   sendHeadProcDataToDestProcs(
      ParallelArray& a_dst_array);

   // The destination processors receive the summed, reduced data from the head
   // processors.
   void
   destProcsRecvHeadProcData(
      ParallelArray& a_dst_array);

   // Local source of reduction.
   const ParallelArray& m_src_array;

    // m_src_array after dimension reduction.
   ParallelArray m_reduced_src_array;

   // Range over which reduction occurs.
   ParallelArray::Box m_src_interior_box;

   // Box describing global problem domain.
   tbox::Box m_domain_box;

   // True in each direction to be reduced.
   deque<bool> m_collapse_dir;

   // Integration volume differential.
   double m_dv;

   // Weight to apply to reduction.
   double m_weight;

   // The source array dimension minus the number of dimensions being collapsed.
   int m_dst_dim;

   // The dimensions that are not collapsed.
   vector<int> m_kept_dims;

   // True if this processor is involved in the reduction.
   bool m_is_in_proc_range;

   // The communicator containing all processors dealing with the same piece
   // of configuration space.
   MPI_Comm m_reduction_comm;

   // The processor ID on m_reduction_comm.
   int m_reduction_comm_id;

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

   // Maximum dimension reduction may be applied to.
   static const int s_MAX_DIM;

   // Tag for send of reduced source array to destination array.
   static const int s_TAG;
};

} // end namespace Loki

#endif
