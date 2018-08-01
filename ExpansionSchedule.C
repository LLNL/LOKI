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
#include "ExpansionSchedule.H"
#include "Directions.H"

namespace Loki {

ExpansionSchedule::ExpansionSchedule(
   const tbox::Box& a_box,
   int a_number_of_species,
   int a_config_space_id,
   const MPI_Comm& a_comm) :
     m_is_phase_space_processor(false),
     m_phase_copy_to(Communication_Manager::numberOfProcessors()),
     m_expansion_comm(a_comm)
{
   // Get the phase space extents of each processor.,
   int my_id = Communication_Manager::My_Process_Number;
   const int n_procs(Communication_Manager::numberOfProcessors());
   int buf[2*CDIM+2];
   for (int i = 0; i < n_procs; ++i) {
      if (my_id == i) {
         for (int j = 0; j < (CDIM+1); ++j) {
            buf[2*j] = a_box.lower(j);
            buf[2*j+1] = a_box.upper(j);
         }
      }
      MPI_Bcast(buf, 2*CDIM+2, MPI_INT, i, MPI_COMM_WORLD);
      m_phase_copy_to[i].processor = i;
      m_phase_copy_to[i].setBounds(buf[0],
         buf[1],
         buf[2],
         buf[3],
         buf[4],
         buf[5]);
   }

   // Determine if this processor deals with phase space.
   if (!m_phase_copy_to[my_id].isEmpty()) {
      m_is_phase_space_processor = true;
   }

   // Only allow the phase space head processors to have extents defining the
   // range of data copied to them as the first step of the communication is
   // from configuration space processors to the phase space head processors.
   // If each species is distributed on only 1 processor then this step is not
   // necessary as there is only 1 phase space processor and by definition it is
   // the head processor.  Note that there is an implicit assumption that there
   // is only 1 Maxwell or Poisson processor.
   if (n_procs > 1+a_number_of_species) {
      for (int i = 0; i < n_procs;) {
         if (m_phase_copy_to[i].isEmpty()) {
            ++i;
         }
         else {
            const IndexBox& head_index_box = m_phase_copy_to[i];
            for (int j = i+1; j < n_procs; ++j) {
               IndexBox& other_index_box = m_phase_copy_to[j];
               bool is_foot_proc = true;
               for (int k = 0; k < CDIM; ++k) {
                  if (head_index_box.base(k) != other_index_box.base(k) ||
                      head_index_box.bound(k) != other_index_box.bound(k)) {
                     is_foot_proc = false;
                     break;
                  }
               }
               if (is_foot_proc) {
                  other_index_box.setBounds(0, -1, 0, -1, 0, -1);
                  if (j == n_procs - 1) {
                     i = n_procs;
                  }
               }
               else {
                  if (j == n_procs - 1) {
                     i = n_procs;
                  }
                  else {
                     i = j;
                  }
                  break;
               }
            }
         }
      }
   }

   // Now create a sub-communicator for each set of phase space processors that
   // deal with the same piece of configuration space.
   if (m_is_phase_space_processor) {
      int comm_id;
      MPI_Comm_rank(a_comm, &comm_id);
      const int status = MPI_Comm_split(a_comm,
         a_config_space_id,
         comm_id,
         &m_expansion_comm);
      if (status != MPI_SUCCESS) {
         OV_ABORT("Configuration space splitting of MPI communicator failed");
      }
   }
}


ExpansionSchedule::~ExpansionSchedule()
{
}


void
ExpansionSchedule::execute(
   const realArray& a_global_config_array,
   RealArray& a_local_phase_array)
{
   // Communicate a_global_config_array which is defined on configuration space
   // processors to a_local_phase_array on the phase space head processors.
   // Then each head processor broadcasts a_local_phase_array to the other
   // phase space processors that deal with that same piece of configuration
   // space.

   // The source range is defined by a_global_config_array.
   Index src_index[PDIM];
   for (int d(0); d < CDIM+1; ++d) {
      int na = min(a_global_config_array.getBase(d), 0);
      int nb = max(a_global_config_array.getBound(d), 0);
      src_index[d] = Range(na, nb);
   }
   src_index[3] = Range(0, 0);

   // Do the communication
   CopyArray::copyArray(a_global_config_array,
      src_index,
      m_phase_copy_to.dataPtr(),
      a_local_phase_array);

   // Each phase space head node communicates its a_local_phase_array to
   // a_local_phase_array from one of its related phase space processors.
   if (isPhaseSpaceProcessor()) {
      broadCast(a_local_phase_array, 0, m_expansion_comm);
   }
}

} // end namespace Loki
