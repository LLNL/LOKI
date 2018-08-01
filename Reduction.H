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
#ifndef _REDUCTION_H_
#define _REDUCTION_H_

#include "Overture.h"
#include "Array.H"
#include "tbox/Dimension.H"
#include "tbox/Box.H"
#include "tbox/Pointer.H"
#include "BoxOps.H"
#include "KernelOp.H"
#include "ReductionOp.H"

namespace Loki {

/**
 * Implements a dimension reduction of a RealArray over a user supplied set of
 * directions using a user supplied operator.
 */
class Reduction
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_dim The dimension before reduction.
    * @param[in] a_collapse_direction The directions the reduction takes place
    *                                 on.
    * @param[in] a_reduction_operator The type of reduction to apply.
    */
   Reduction(
      const tbox::Dimension& a_dim,
      const Array<bool>& a_collapse_direction,
      const tbox::Pointer<ReductionOp>& a_reduction_operator);

   /**
    * @brief Destructor.
    */
   ~Reduction();

   /**
    * @brief Apply the reduction to a_in placing results in a_out.
    *
    * @param[out] a_out Result of reduction.
    * @param[in] a_in Array to reduce.
    * @param[in] a_box Space reduction takes place on.
    * @param[in] a_kernel_operator Pre and post processing of kernel. 
    * @param[in] a_comm Communicator for global reduction.
    */
   void
   apply(
      RealArray& a_out,
      const RealArray& a_in,
      const tbox::Box& a_box,
      const KernelOp& a_kernel_operator,
      MPI_Comm& a_comm)
   {
      if (m_dim != a_box.getDim()) {
         OV_ABORT("a_box dimension != Dimension before reduction");
      }
      if (m_dim != a_in.numberOfDimensions()) {
         OV_ABORT("a_in dimension != Dimension before reduction");
      }

      RealArray temp(a_in);
      temp = 0;
      a_kernel_operator.applyLocal(temp, a_in);

      RealArray tempLocal(a_out);
      tempLocal = 0;
      m_reduction_op->applyLocal(tempLocal, temp, a_box, m_collpse_dir);
      m_reduction_op->applyGlobal(a_out,
         tempLocal,
         a_comm,
         ReductionOp::REDUCE_TO_HEAD_NODE);

      a_kernel_operator.applyGlobal(a_out);
   }

   /**
    * @brief Returns true for each direction to be reduced.
    */
   const Array<bool>&
   collapseDirection() const
   {
      return m_collpse_dir;
   }

private:
   // Unimplemented default constructor.
   Reduction();

   // Unimplemented copy constructor.
   Reduction(
      const Reduction& other);

   // Unimplemented assignment operator.
   Reduction&
   operator = (
      const Reduction& rhs);

   // Dimension over which reduction occurs.
   tbox::Dimension m_dim;

   // User supplied reduction operator.
   tbox::Pointer<ReductionOp> m_reduction_op;

   // True for each direction to be reduced.
   Array<bool> m_collpse_dir;
};

} // end namespace Loki

#endif
