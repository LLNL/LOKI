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
#ifndef _REDUCTION_SCHEDULE4D_H_
#define _REDUCTION_SCHEDULE4D_H_

#include "Overture.h"
#include "Directions.H"
#include "ProblemDomain.H"
#include "tbox/Pointer.H"
#include "tbox/Dimension.H"
#include "ReductionOp.H"

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
    * @param[in] a_reduction Specific reduction to be performed.
    * @param[in] a_processor_range Processors involved in reduction.
    * @param[in] a_comm Communicator corresponding to processor_range
    */
   ReductionSchedule4D(
      const RealArray& a_local_src_array,
      const tbox::Box& a_local_interior_box,
      const tbox::Box& a_domain_box,
      const Array<bool>& a_collapse_dir,
      const tbox::Pointer<ReductionOp>& a_reduction,
      const Range& a_processor_range,
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
      RealArray& a_local_dst_array);

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

   // Local source of reduction.
   const RealArray& m_local_src_array;

   // True in each direction to be reduced.
   Array<bool> m_collapse_dir;

   bool m_is_in_proc_range;

   MPI_Comm m_reduction_comm;

   // User supplied reduction operation.
   tbox::Pointer<ReductionOp> m_reduction;

   // Maximum dimension reduction may be applied to.
   const static int s_MAX_DIM;
};

} // end namespace Loki

#endif
