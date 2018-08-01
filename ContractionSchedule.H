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
#ifndef _CONTRACTION_SCHEDULE_H_
#define _CONTRACTION_SCHEDULE_H_

#include "Overture.h"
#include "Directions.H"
#include "ProblemDomain.H"
#include "tbox/Pointer.H"
#include "tbox/Dimension.H"

namespace Loki {

/**
 * Performs the communication and other operations necessary to communicate a
 * user supplied local array from the 4D processors to the 2D processor(s).
 * Pretty much the inverse of ExpansionSchedule except that here there is no
 * need to communicate among the group of 4D processors covering the same
 * piece of 2D space.
 */
class ContractionSchedule
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_local_src_array The local array to be contracted.
    * @param[in] a_local_interior_box Range of the contraction.
    * @param[in] a_global_box The Range of data to be sent to the global array.
    * @param[in] a_domain Global problem domain.
    * @param[in] a_processor_range Processors involved in contraction.
    * @param[in] a_comm Communicator corresponding to processor_range
    */
   ContractionSchedule(
      const RealArray& a_local_src_array,
      const tbox::Box& a_local_interior_box,
      const tbox::Box& a_global_box,
      const ProblemDomain& a_domain,
      const Range& a_processor_range,
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
      realArray& a_global_dst_array);

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

   // Compute range of local source.
   void
   computeLocalSrcIndexArray();

   // Compute range of global destination.
   void
   computeGlobalDstIndexArray();

   // Construct communicators for each set of processors that deal with the
   // same piece of 2D configuration space.
   void
   constructCommunicator();

   // Constructs set of processors that deal with the same piece of 2D
   // configuration space.
   void
   constructSetOfProcessors(
      const Range& a_processor_range);

   // Returns true if a_proc_id is in the range of processors a_processor_range.
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
      return m_2D_comm_id == 0;
   }

   // Returns true if this processor sends or receives data.
   bool
   thisProcessorSendsOrReceivesData(
      const bool& a_this_processor_will_receive_data)
   {
      return (m_this_processor_has_data || a_this_processor_will_receive_data);
   }

   // Local source of contraction.
   const RealArray& m_local_src_array;

   // Range over which contraction occurs.
   tbox::Box m_local_interior_box;

   // Box describing the Range of data transferred to the global array.
   tbox::Box m_global_box;

   // Box describing global problem domain.
   tbox::Box m_domain_box;

   Array<Index> m_local_src_index;

   Array<Index> m_global_dst_index;

   IntegerArray m_proc_set;

   bool m_is_in_proc_range;

  // True if this processor has data 
   bool m_this_processor_has_data;

   MPI_Comm m_comm;

   int m_comm_id;

   // Processor id of groups of 4D processors all dealing with same piece of 2D
   // configuration space.
   int m_2D_comm_id;
};

} // end namespace Loki

#endif
