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
#include "SummationSchedule.H"

#include "BoxOps.H"
#include "KernelOp.H"
#include "TimerManager.H"
#include "Loki_Utilities.H"
#include "ReductionOp.H"

namespace Loki {

SummationSchedule::SummationSchedule(
   const RealArray& a_local_src_array,
   const tbox::Box& a_local_interior_box,
   const tbox::Box& a_domain_box,
   const Range& a_processor_range,
   const MPI_Comm& a_comm)
   : m_local_src_array(a_local_src_array),
     m_local_dst_array(a_local_src_array),
     m_is_in_proc_range(false),
     m_this_processor_has_data(false),
     m_comm_id(-1),
     m_summation_comm(a_comm),
     m_summation_comm_id(-1)
{
   const int my_id(std::max(0, Communication_Manager::My_Process_Number));
   m_is_in_proc_range = isInRange(my_id, a_processor_range);

   // Compute the range in each dimension of the global destination array.
   computeLocalDstIndexArray(a_local_interior_box);

   // Compute the range in each dimension of the global destination array.
   computeGlobalDstIndexArray(a_domain_box);

   // Construct sub-communicators for each group of processors dealing with
   // the same piece of configuration space.
   constructCommunicator(a_local_interior_box, a_domain_box, a_comm);

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
}


SummationSchedule::~SummationSchedule()
{
}


void
SummationSchedule::execute(
   realArray& a_global_dst_array)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("summation");

   // Sum data on phase space nodes onto the head nodes.
   m_local_dst_array = 0;
   if (m_is_in_proc_range) {
      if (a_global_dst_array.numberOfDimensions() !=
          m_local_src_array.numberOfDimensions()) {
         OV_ABORT("global dst array dim and local src array dim do not match");
      }

      ParallelUtility::getSums(&(*m_local_src_array.getDataPointer()),
         &(*m_local_dst_array.getDataPointer()),
         m_local_src_array.elementCount(),
         ReductionOp::REDUCE_TO_HEAD_NODE,
         m_summation_comm);
   }

   // Now copy data from the phase space head nodes to the field defined on
   // the configuration space processor(s).

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

   timers->stopTimer("summation");
}


//////// PRIVATE METHODS //////////////////////////////////////////////////


void
SummationSchedule::computeLocalDstIndexArray(
   const tbox::Box& a_local_interior_box)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("summation");

   // The Ranges of each dimension of the local destination are the Ranges of
   // the local interior box.
   const tbox::Dimension& src_dim(a_local_interior_box.getDim());
   m_local_dst_index.resize(src_dim);
   if (m_is_in_proc_range) {
      for (int dir(0); dir < src_dim; ++dir) {
         m_local_dst_index[dir] = BoxOps::range(a_local_interior_box, dir);
      }
   }
   else {
      for (int dir(0); dir < src_dim; ++dir) {
         m_local_dst_index[dir] = Range(0, -1);
      }
   }
   timers->stopTimer("summation");
}


void
SummationSchedule::computeGlobalDstIndexArray(
   const tbox::Box& a_domain_box)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("summation");

   // The Ranges of each dimension of the global destination array are the
   // dimensions of the domain box.
   const tbox::Dimension& src_dim(a_domain_box.getDim());
   m_global_dst_index.resize(src_dim);
   for (int dir(0); dir < src_dim; ++dir) {
      m_global_dst_index[dir] = BoxOps::range(a_domain_box, dir);
   }
   timers->stopTimer("summation");
}


void
SummationSchedule::constructCommunicator(
   const tbox::Box& a_local_interior_box,
   const tbox::Box& a_domain_box,
   const MPI_Comm& a_comm)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("summation");

   // Find the groups of processors that deal with the same piece of
   // configuration space and make a sub-communicator for each group.  The
   // color variable is essentially an offset to the first point of
   // configuration space dealt with by this processor.  So color is the same
   // for all processors dealing with the same piece of configuration space.
   if (m_is_in_proc_range) {
      int color = a_local_interior_box.lower(0) +
         a_domain_box.numberCells(0)*a_local_interior_box.lower(1);

      MPI_Comm_rank(a_comm, &m_comm_id);
      const int status(MPI_Comm_split(a_comm,
         color,
         m_comm_id,
         &m_summation_comm));
      MPI_Comm_rank(m_summation_comm, &m_summation_comm_id);
      if (status != MPI_SUCCESS) {
         OV_ABORT("Configuration space splitting of MPI communicator failed");
      }
   }
   timers->stopTimer("summation");
}


void
SummationSchedule::constructSetOfProcessors(
   const Range& a_processor_range)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("summation");

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

   // Now get the processor IDs based on the communicator passed into
   // constructCommunicator of the head nodes.
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
   timers->stopTimer("summation");
}

} // end namespace Loki
