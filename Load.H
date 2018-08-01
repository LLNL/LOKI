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
#ifndef _LOAD_H_
#define _LOAD_H_

#include "Overture.h"

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
    * distributed.
    */
   virtual
   int
   numberOfProcessors() const = 0 ;

   /**
    * @brief Returns true if this load is distributed over a user specified
    * fixed number of processors.
    */
   virtual
   bool
   fixedNumberOfProcessors() const = 0;

   /**
    * @brief Partitions this load over the necessary number of processors.
    *
    * @param[in] a_range Processor range distributed over.
    * @param[in] a_comm Corresponding communicator.
    */
   virtual
   void
   createPartition(
      const Range& a_range,
      const MPI_Comm& a_comm) = 0;

   /**
    * @brief Returns true if this load is distributed to processor with
    * supplied id.
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
