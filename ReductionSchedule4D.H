/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _REDUCTION_SCHEDULE4D_H_
#define _REDUCTION_SCHEDULE4D_H_

#include "Directions.H"
#include "tbox/Box.H"
#include "tbox/Pointer.H"
#include "tbox/Dimension.H"
#include "ParallelArray.H"

#include <mpi.h>

namespace Loki {

/**
 * Performs the communication and other operations necessary to accomplish a
 * user defined reduction of dimension on a user supplied local array.  This
 * differs from ReductionSchedule in that the result remains on the 4D
 * processors and is not communicated to the 2D processor(s).
 */
class ReductionSchedule4D
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_local_src_array The local array to be reduced.
    * @param[in] a_local_interior_box Range of the reduction.
    * @param[in] a_domain_box Global problem domain box.
    * @param[in] a_collapse_dir The directions the reduction takes place on.
    * @param[in] a_dv Cell volume over reduced dimensions.
    * @param[in] a_proc_lo Lowest rank processor involved in reduction.
    * @param[in] a_proc_hi Highest rank processor involved in reduction.
    * @param[in] a_comm Communicator corresponding to range of processors.
    */
   ReductionSchedule4D(
      const ParallelArray& a_local_src_array,
      const ParallelArray::Box& a_local_interior_box,
      const tbox::Box& a_domain_box,
      const deque<bool>& a_collapse_dir,
      const double a_dv,
      int a_proc_lo,
      int a_proc_hi,
      const MPI_Comm& a_comm);

   /**
    * @brief Destructor.
    */
   ~ReductionSchedule4D();

   /**
    * @brief Perform the reduction.
    *
    * @param[in] a_local_dst_array The local array for the reduction.
    */
   void
   execute(
      ParallelArray& a_local_dst_array);

   int
   GetCommSize()
   {
      return m_commSize;
   }
  
private:
   // Unimplemented default constructor.
   ReductionSchedule4D();

   // Unimplemented copy constructor.
   ReductionSchedule4D(
      const ReductionSchedule4D& other);

   // Unimplemented assignment operator.
   ReductionSchedule4D&
   operator = (
      const ReductionSchedule4D& rhs);

   // Construct communicator over which reduction takes place.
   void
   constructCommunicator(
      const ParallelArray::Box& a_local_interior_box,
      const tbox::Box& a_domain_box,
      const deque<bool>& a_collapse_dir,
      const MPI_Comm& a_comm);

   // Returns true if a_proc_id is in the range of supplied processors.
   static
   bool
   isInRange(
      int a_proc_id,
      int a_proc_lo,
      int a_proc_hi)
   {
      return ((a_proc_id >= a_proc_lo) && (a_proc_id <= a_proc_hi));
   }

   // Local source of reduction.
   const ParallelArray& m_local_src_array;

   // Cell volume over reduced dimensions.
   double m_dv;

   bool m_is_in_proc_range;

   MPI_Comm m_reduction_comm;
   int m_commSize;

   // Maximum dimension reduction may be applied to.
   const static int s_MAX_DIM;

   // Indicates if the global reduction should only be sent to the head or all
   // nodes on the communicator.
   enum ReductionTargets {
     REDUCE_TO_ALL_NODES = -1,
     REDUCE_TO_HEAD_NODE
   };
};

} // end namespace Loki

#endif
