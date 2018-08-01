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
#include "ReductionSchedule.H"

#include "BoxOps.H"
#include "KernelOp.H"
#include "TimerManager.H"
#include "Loki_Utilities.H"

namespace Loki {

const int ReductionSchedule::s_MAX_DIM = 4;

ReductionSchedule::ReductionSchedule(
   const RealArray& a_local_src_array,
   const tbox::Box& a_local_interior_box,
   const ProblemDomain& a_domain,
   const tbox::Pointer<Reduction>& a_reduction,
   const Range& a_processor_range,
   const MPI_Comm& a_comm)
   : m_local_src_array(a_local_src_array),
     m_local_interior_box(a_local_interior_box),
     m_domain_box(a_domain.box()),
     m_collapse_dir(a_reduction->collapseDirection()),
     m_dst_dim(0),
     m_is_in_proc_range(false),
     m_this_processor_has_data(false),
     m_comm(a_comm),
     m_comm_id(-1),
     m_reduction_comm(a_comm),
     m_reduction_comm_id(-1),
     m_reduction(a_reduction)
{
   const int my_id(std::max(0, Communication_Manager::My_Process_Number));
   m_is_in_proc_range = isInRange(my_id, a_processor_range);

   // This is a misnomer as it only computes the dimension of dst array.
   computeSrcAndDstDimensions();

   // Compute the range in each dimension of the global destination array.
   computeLocalDstIndexArray();

   // Compute the range in each dimension of the global destination array.
   computeGlobalDstIndexArray();

   // Construct sub-communicators for each group of processors dealing with
   // the same piece of configuration space.
   constructCommunicator();

   // Construct the list of processors that hold data to send to the 2D
   // processor.
   constructSetOfProcessors(a_processor_range);

   // Determine if this processor holds data to send to the 2D processor.
   const int myid(Communication_Manager::My_Process_Number);
   for (int p(0); p < m_proc_set.getLength(0); ++p) {
      if (m_proc_set(p) == myid) {
         m_this_processor_has_data = true;
         break;
      }
   }

   allocateLocalDestinationArray(m_local_dst_array);
}


ReductionSchedule::~ReductionSchedule()
{
}


void
ReductionSchedule::execute(
   realArray& a_global_dst_array,
   const KernelOp& a_kernel_op)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("reduction");
   // do sum reduction on phase space nodes
   if (m_is_in_proc_range) {
      m_reduction->apply(m_local_dst_array,
         m_local_src_array,
         m_local_interior_box,
         a_kernel_op,
         m_reduction_comm);
   }
   else {
      m_local_dst_array.redim(0);
   }
   // now we take the local_dst_array on the head nodes and form a global dst
   // array across these head nodes

   // The processors that own the global array will receive data.  See if this
   // processor is one of them.
   const int myid(Communication_Manager::My_Process_Number);
   const intSerialArray& copyToProcessorSet(
      a_global_dst_array.getPartition().getProcessorSet());
   bool this_processor_will_receive_data(false);
   for (int p(0); p < copyToProcessorSet.getLength(0); ++p) {
      if (copyToProcessorSet(p) == myid) {
         this_processor_will_receive_data = true;
         break;
      }
   }

   // If this processor owns part of the global array or is a head node then
   // it must be involved in the data communication.  For each processor in
   // m_proc_set, the head nodes, transfer data in m_local_src_array defined by
   // m_local_src_index to the location in a_global_dst_array defined by
   // m_global_dst_index.
   if (thisProcessorSendsOrReceivesData(this_processor_will_receive_data)) {
      CopyArray::copyArray(m_local_dst_array,
         m_local_dst_index.dataPtr(),
         m_proc_set,
         a_global_dst_array,
         m_global_dst_index.dataPtr());
   }

   timers->stopTimer("reduction");
}


//////// PRIVATE METHODS //////////////////////////////////////////////////


void
ReductionSchedule::computeSrcAndDstDimensions()
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("reduction");

   // The destination dimension is the source dimension minus the number of
   // dimensions being collapsed (integrate over).
   const tbox::Dimension& src_dim(m_local_interior_box.getDim());
   if (src_dim != m_collapse_dir.length()) {
      OV_ABORT("Incongruent dimensions in ReductionSchedule construction");
   }
   for (int dir(0); dir < src_dim; ++dir) {
      if (!m_collapse_dir[dir]) {
         ++m_dst_dim;
      }
   }
   if (m_dst_dim == src_dim) {
      OV_ABORT("No reduction in ReductionSchedule construction");
   }
   timers->stopTimer("reduction");
}


void
ReductionSchedule::computeLocalDstIndexArray()
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("reduction");

   // Similar to computeSecAndDstDimensions.  The Ranges of each dimension of
   // the local destination are the Ranges of the local interior box.  These
   // are the source dimensions that are not collapsed.
   const tbox::Dimension& src_dim(m_local_interior_box.getDim());
   m_local_dst_index.resize(src_dim);
   if (m_is_in_proc_range) {
      for (int dir(0); dir < src_dim; ++dir) {
         if (m_collapse_dir[dir]) {
            m_local_dst_index[dir] = Range(0, 0);
         }
         else {
            m_local_dst_index[dir] =
               BoxOps::range(m_local_interior_box, dir);
         }
      }
   }
   else {
      for (int dir(0); dir < src_dim; ++dir) {
         m_local_dst_index[dir] = Range(0, -1);
      }
   }
   timers->stopTimer("reduction");
}


void
ReductionSchedule::computeGlobalDstIndexArray()
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("reduction");

   // The Ranges of each dimension of the global destination array are the
   // dimensions of the domain box.  These are the source dimensions that are
   // not collapsed.
   m_global_dst_index.resize(s_MAX_DIM);
   {
      for (int dir(0); dir < s_MAX_DIM; ++dir) {
         if (!m_collapse_dir[dir]) {
            m_global_dst_index[dir] = BoxOps::range(m_domain_box, dir);
         }
         else {
            m_global_dst_index[dir] = Range(0, 0);
         }
      }
   }
   timers->stopTimer("reduction");
}


void
ReductionSchedule::constructSetOfProcessors(
   const Range& a_processor_range)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("reduction");

   // Find all the processors that are head nodes.
   int n(a_processor_range.length());
   std::vector<int> is_head_node_lcl(n, 0);
   if (m_comm_id >= 0) {
      is_head_node_lcl[m_comm_id] = isHeadNode() ? 1 : 0;
   }
   std::vector<int> is_head_node(n, 0);
   MPI_Allreduce(&(is_head_node_lcl[0]),
      &(is_head_node[0]),
      n,
      MPI_INT,
      MPI_SUM,
      MPI_COMM_WORLD);

   // Now get the processor ids based on m_comm of the head nodes.
   int count(0);
   for (int np(0); np < a_processor_range.length(); ++np) {
      count += is_head_node[np];
   }
   m_proc_set.resize(count);
   int next(0);
   for (int np(0); np < a_processor_range.length(); ++np) {
      if (is_head_node[np] == 1) {
         m_proc_set(next) = a_processor_range.getBase() + np;
         ++next;
      }
   }
   timers->stopTimer("reduction");
}


void
ReductionSchedule::allocateLocalDestinationArray(
   RealArray& a_local_dst_array)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("reduction");

   // The local destination array has non-empty Ranges in the uncollapsed
   // dimensions.  These Ranges have been determined already and are in
   // m_local_dst_index.
   const tbox::Box dst_box(BoxOps::indexArrayToBox(m_local_dst_index));
   std::vector<Range> r(s_MAX_DIM);
   for (int d(0); d < s_MAX_DIM; ++d) {
      if (m_collapse_dir[d]) {
         r[d] = Range(0, 0);
      }
      else {
         r[d] = BoxOps::range(dst_box, d);
      }
   }
   if (dst_box.getDim() == tbox::Dimension(4)) {
      a_local_dst_array.redim(r[0], r[1], r[2], r[3]);
   }
   else {
      OV_ABORT("Not implemented for D>4!");
   }
   a_local_dst_array = 0;
   timers->stopTimer("reduction");
}


void
ReductionSchedule::constructCommunicator()
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("reduction");

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
      color = m_local_interior_box.lower(ibeg);
      for (int d(ibeg-1); d >= 0; --d) {
         if (!m_collapse_dir[d]) {
            color *= m_domain_box.numberCells(d);
            color += m_local_interior_box.lower(d);
         }
      }

      MPI_Comm_rank(m_comm, &m_comm_id);
      const int status(MPI_Comm_split(m_comm,
         color,
         m_comm_id,
         &m_reduction_comm));
      MPI_Comm_rank(m_reduction_comm, &m_reduction_comm_id);
      if (status != MPI_SUCCESS) {
         OV_ABORT("Configuration space splitting of MPI communicator failed");
      }
   }
   timers->stopTimer("reduction");
}

} // end namespace Loki
