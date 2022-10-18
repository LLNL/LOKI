/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "LoadBalancer.H"
#include "Loki_Defines.H"
#include "Loki_Utilities.H"

namespace Loki
{


void
LoadBalancer::balance(
   const vector<Load*>& a_loads)
{
   // Let's first make a list of the loads.  The load balancing process it a
   // 2 step process so as we balance things in the first step we want to
   // remove them from the continer of loade.  Hence, We want a list.
   list<Load*> load_list;
   for (int i(0); i < static_cast<int>(a_loads.size()); ++i) {
      load_list.push_back(a_loads[i]);
   }

   // First partition those loads with fixed numbers of processors.  This is
   // essentially the Poisson or Maxwell load.
   int next_proc(0);
   list<Load*>::iterator llit(load_list.begin());
   while (llit != load_list.end()) {
      if ((*(llit))->fixedNumberOfProcessors()) {
         const int n((*(llit))->numberOfProcessors());
         if (n + next_proc < Loki_Utilities::s_num_procs) {
            int proc_lo = next_proc;
            int proc_hi = next_proc + n -1;
            MPI_Comm comm(createCommunicator(proc_lo, proc_hi));
            (*(llit))->createPartition(proc_lo, proc_hi, comm);
            next_proc += n;
            llit = load_list.erase(llit);
         }
         else {
            LOKI_ABORT("Too few processors to meet request!");
         }
      }
      else {
         ++llit;
      }
   }

   // Now load balance everything else in a proportional way.

   // Count up the number of loads remaining and sum up their total work.
   const int n_procs_left(Loki_Utilities::s_num_procs - next_proc);
   float net_work(0.0);
   int count(0);
   for (llit = load_list.begin(); llit != load_list.end(); ++llit) {
      net_work += (*llit)->netCost();
      ++count;
   }

   if (count > n_procs_left) {
      LOKI_ABORT("Too few processors to meet request!");
   }

   // Figure out the relative work of each load and give each its share of the
   // processors that are left.
   for (llit = load_list.begin(); llit != load_list.end(); ++llit) {
      const float relative_work =
         ((*llit)->netCost() * static_cast<float>(n_procs_left)) / net_work;
      const int n_proportional(static_cast<int>(floor(relative_work + 0.5)));
      const int n_left(Loki_Utilities::s_num_procs - next_proc);
      const int n = min(n_proportional, n_left);
      if (n > 0) {
         int proc_lo = next_proc;
         int proc_hi = next_proc + n -1;
         MPI_Comm comm(createCommunicator(proc_lo, proc_hi));
         (*llit)->createPartition(proc_lo, proc_hi, comm);
         next_proc += n;
      }
      else {
         LOKI_ABORT("Too few processors to meet request!");
      }
   }
}


MPI_Comm
LoadBalancer::createCommunicator(
   int a_proc_lo,
   int a_proc_hi)
{
   // Create a communicator consisting of the processors with ranks between
   // a_proc_lo and a_proc_hi, inclusive.
   vector<int> ranks(a_proc_hi-a_proc_lo+1);
   for (int i(0); i < static_cast<int>(ranks.size()); ++i) {
      ranks[i] = i + a_proc_lo;
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
      LOKI_ABORT("Creating MPI communicator group failed");
   }

   MPI_Comm comm(MPI_COMM_NULL);
   result = MPI_Comm_create(MPI_COMM_WORLD, group, &comm);
   if (result != MPI_SUCCESS) {
      LOKI_ABORT("Creating MPI communicator failed");
   }

   return comm;
}

} // end namespace Loki
