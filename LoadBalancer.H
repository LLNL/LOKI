/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _LOAD_BALANCER_H_
#define _LOAD_BALANCER_H_

#include <list>
#include <vector>
#include <mpi.h>

#include "Load.H"

using namespace std;

namespace Loki
{

/**
 * A class that distributes work of a collection of abstract Loads across
 * processors.
 */
class LoadBalancer 
{
public:
   /**
    * @brief Balance a collection of loads across processors.
    *
    * @param[in] a_loads What needs to be load balanced.
    */
   static
   void
   balance(
      const vector<Load*>& a_loads);

private:
   // Unimplemented default constructor.
   LoadBalancer();

   // Unimplemented destructor.
   ~LoadBalancer();

   // Unimplemented copy constructor.
   LoadBalancer(
      const LoadBalancer& other);

   // Unimplemented assignment operator.
   LoadBalancer&
   operator = (
      const LoadBalancer& rhs);

   // Create communicator for a load on a range of processors.
   static
   MPI_Comm
   createCommunicator(
      int a_proc_lo,
      int a_proc_hi);
};

} // end namespace Loki

#endif
