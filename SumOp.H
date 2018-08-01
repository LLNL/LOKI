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
#ifndef _SUM_OP_H_
#define _SUM_OP_H_

#include "Overture.h"
#include "ParallelUtility.h"
#include "Directions.H"
#include "tbox/Dimension.H"
#include "ReductionOp.H"

#define MAX_DIMENSION 4

namespace Loki {

/**
 * Implementation of ReductionOp interface that sums over specified directions.
 */
class SumOp : public ReductionOp
{
public:
   /**
    * @brief Constructor.
    */
   SumOp();

   /**
    * @brief Destructor.
    */
   virtual
   ~SumOp();

   /**
    * @brief Sum a_source over the directions specified in a_collpse_dir and
    * place result in a_result.
    *
    * @param[out] a_result Result of summation.
    * @param[in] a_source Source of summation.
    * @param[in] a_box Box defining range of summation.
    * @param[in] a_collpse_dir Specifies directions to sum over.
    */
   virtual
   void
   applyLocal(
      RealArray& a_result,
      const RealArray& a_source,
      const tbox::Box& a_box,
      const Array<bool>& a_collpse_dir) const;

   /**
    * @brief Do sum reduction of a_source using communicator a_comm placing
    * result in a_result.
    *
    * @param[out] a_result Final result of sum reduction.
    * @param[in] a_source Intermediate result of applyLocal.
    * @param[in] a_comm Communicator on which a_source is further reduced.
    * @param[in] a_reduction_target Reduce to all or head node of a_comm.
    */
   virtual
   void
   applyGlobal(
      RealArray& a_result,
      const RealArray& a_source,
      MPI_Comm& a_comm,
      ReductionTargets a_reduction_target) const;

private:

   // Unimplemented copy constructor.
   SumOp(
      const SumOp& other);

   // Unimplemented assignment operator.
   SumOp&
   operator = (
      const SumOp& rhs);
};

} // end namespace Loki

#endif
