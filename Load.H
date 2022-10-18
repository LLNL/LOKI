/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _LOAD_H_
#define _LOAD_H_

#include <mpi.h>

namespace Loki {

/**
 * Abstract base class defining the interface to a computational load.
 */
class Load
{
public:
   /**
    * @brief Constructor.
    */
   Load();

   /**
    * @brief Destructor.
    */
   virtual
   ~Load();

   /**
    * @brief Returns the computational cost of this load.
    */
   virtual
   float
   netCost() const = 0;

   /**
    * @brief Returns the number of processors over which this load is
    *        distributed.
    */
   virtual
   int
   numberOfProcessors() const = 0 ;

   /**
    * @brief Returns true if this load is distributed over a user specified
    *        fixed number of processors.
    */
   virtual
   bool
   fixedNumberOfProcessors() const = 0;

   /**
    * @brief Partitions this load over the necessary number of processors.
    *
    * @param[in] a_proc_lo Lowest rank processor distributed over.
    * @param[in] a_proc_hi Highest rank processor distributed over.
    * @param[in] a_comm Corresponding communicator.
    */
   virtual
   void
   createPartition(
      int a_proc_lo,
      int a_proc_hi,
      const MPI_Comm& a_comm) = 0;

   /**
    * @brief Returns true if this load is distributed to processor with
    *        supplied id.
    *
    * @param[in] a_proc_id ID of processor of interest.
    */
   virtual
   bool
   isInRange(
      int a_proc_id) const = 0;

   /**
    * @brief Print information about processors this load is distributed over.
    */
   virtual
   void
   printDecomposition() const = 0;

private:
   // Unimplemented copy constructor.
   Load(
      const Load& other);

   // Unimplemented assignment operator.
   Load& 
   operator = (
      const Load& rhs);
};

} // end namespace Loki

#endif
