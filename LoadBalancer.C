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
#include "LoadBalancer.H"

namespace Loki
{

LoadBalancer::LoadBalancer()
   :  m_num_procs(Communication_Manager::numberOfProcessors())
{
}


LoadBalancer::~LoadBalancer()
{
}


MPI_Comm
LoadBalancer::createCommunicator(
   const Range& a_range)
{
   // Create a communicator consisting of the processors with ranks in a_range.
   std::vector<int> ranks(a_range.length());
   for (int i(0); i < static_cast<int>(ranks.size()); ++i) {
      ranks[i] = i + a_range.getBase();
   }

   MPI_Group world_group(MPI_GROUP_NULL);
   MPI_Comm_group(MPI_COMM_WORLD, &world_group);

   int result(MPI_SUCCESS);
   MPI_Group group(MPI_GROUP_NULL);
   result = MPI_Group_incl(world_group,
      static_cast<int>(ranks.size()),
      &(ranks[0]),
      &group);
   if (result != MPI_SUCCESS) {
      OV_ABORT("Creating MPI communicator group failed");
   }

   MPI_Comm comm(MPI_COMM_NULL);
   result = MPI_Comm_create(MPI_COMM_WORLD, group, &comm);
   if (result != MPI_SUCCESS) {
      OV_ABORT("Creating MPI communicator failed");
   }

   return comm;
}


void
LoadBalancer::balance(
   const std::vector<Load*>& a_loads)
{
   // Let's first make a list of the loads.  The load balancing process it a
   // 2 step process so as we balance things in the first step we want to
   // remove them from the continer of loade.  Hence, We want a list.
   std::list<Load*> load_list;
   for (int i(0); i < static_cast<int>(a_loads.size()); ++i) {
      load_list.push_back(a_loads[i]);
   }

   // First partition those loads with fixed numbers of processors.  This is
   // essentially the Poisson or Maxwell load.
   int next_proc(0);
   std::list<Load*>::iterator llit(load_list.begin());
   while (llit != load_list.end()) {
      if ((*(llit))->fixedNumberOfProcessors()) {
         const int n((*(llit))->numberOfProcessors());
         if (n + next_proc < m_num_procs) {
            Range range(next_proc, next_proc + n - 1);
            MPI_Comm comm(createCommunicator(range));
            (*(llit))->createPartition(range, comm);
            next_proc += n;
            llit = load_list.erase(llit);
         }
         else {
            OV_ABORT("Too few processors to meet request!");
         }
      }
      else {
         ++llit;
      }
   }

   // Now load balance everything else in a proportional way.

   // Count up the number of loads remaining and sum up their total work.
   const int n_procs_left(m_num_procs - next_proc);
   float net_work(0.0);
   int count(0);
   for (llit = load_list.begin(); llit != load_list.end(); ++llit) {
      net_work += (*llit)->netCost();
      ++count;
   }

   if (count > n_procs_left) {
      OV_ABORT("Too few processors to meet request!");
   }

   // Figure out the relative work of each load and give each its share of the
   // processors that are left.
   for (llit = load_list.begin(); llit != load_list.end(); ++llit) {
      const float relative_work =
         ((*llit)->netCost() * static_cast<float>(n_procs_left)) / net_work;
      const int n_proportional(static_cast<int>(floor(relative_work + 0.5)));
      const int n_left(m_num_procs - next_proc);
      const int n = std::min(n_proportional, n_left);
      if (n > 0) {
         Range range(next_proc, next_proc + n - 1);
         MPI_Comm comm(createCommunicator(range));
         (*llit)->createPartition(range, comm);
         next_proc += n;
      }
      else {
         OV_ABORT("Too few processors to meet request!");
      }
   }
}

} // end namespace Loki
