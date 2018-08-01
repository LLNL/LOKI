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
#ifndef _SUMMATION_SCHEDULE_H_
#define _SUMMATION_SCHEDULE_H_

#include "Overture.h"
#include "Directions.H"
#include "ProblemDomain.H"
#include "tbox/Pointer.H"
#include "tbox/Dimension.H"

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
    * @param[in] a_local_interior_box Range of the summation.
    * @param[in] a_domain_box Global problem domain box.
    * @param[in] a_processor_range Processors involved in summation.
    * @param[in] a_comm Communicator corresponding to processor_range
    */
   SummationSchedule(
      const RealArray& a_local_src_array,
      const tbox::Box& a_local_interior_box,
      const tbox::Box& a_domain_box,
      const Range& a_processor_range,
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
      realArray& a_global_dst_array);

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

   // Compute range of local destination.
   void
   computeLocalDstIndexArray(
      const tbox::Box& a_local_interior_box);

   // Compute range of global destination.
   void
   computeGlobalDstIndexArray(
      const tbox::Box& a_domain_box);

   // Constructs set of processors that hold local reduced data.
   void
   constructSetOfProcessors(
      const Range& a_processor_range);

   // Construct communicator over which summation takes place.
   void
   constructCommunicator(
      const tbox::Box& a_local_interior_box,
      const tbox::Box& a_domain_box,
      const MPI_Comm& a_comm);

   // Returns true if a_proc_id is in the range of processors
   // a_processor_range.
   static
   bool
   isInRange(
      int a_proc_id,
      const Range& a_processor_range)
   {
      return ((a_proc_id >= a_processor_range.getBase()) &&
              (a_proc_id <= a_processor_range.getBound()));
   }

   // Returns true if this is a head node.
   bool
   isHeadNode()
   {
      return m_summation_comm_id == 0;
   }

   // Returns true if this processor sends or receives data.
   bool
   thisProcessorSendsOrReceivesData(
      const bool& a_this_processor_will_receive_data)
   {
      return (m_this_processor_has_data || a_this_processor_will_receive_data);
   }

   // Local source of summation.
   const RealArray& m_local_src_array;

   // Local result of summation.
   RealArray m_local_dst_array;

   // The range of data in m_local_dst_array to be communicated to the global
   // destination array.
   Array<Index> m_local_dst_index;

   // The range of data in the global destination array that receives data from
   // m_local_dst_array.
   Array<Index> m_global_dst_index;

   // The processor IDs based on m_comm of the head nodes.
   IntegerArray m_proc_set;

   // True if this processor is involved in the reduction.
   bool m_is_in_proc_range;

   // True if this processor has data.
   bool m_this_processor_has_data;

   // Processor ID on the communicator in the constructor.
   int m_comm_id;

   // The communicator containing all processors dealing with the same piece of
   // configuration space.
   MPI_Comm m_summation_comm;

   // The processor ID on m_reduction_comm.
   int m_summation_comm_id;
};

} // end namespace Loki

#endif
