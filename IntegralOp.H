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
#ifndef _INTEGRAL_OP_H_
#define _INTEGRAL_OP_H_

#include "Overture.h"
#include "SumOp.H"

namespace Loki {

/**
 * A specific type of SumOp that forms integral f(v) dv where f is constant
 * over each dv element and each dv element is the same.  This reduces to
 * sum(i) (fi) * dv.  Therefore SumOp's local operation which performs the
 * local sum is inherited, and the global operation is reimplemented to first
 * perform SumOp's reduction of the local sums followed by multiplication by
 * dv.
 */
class IntegralOp : public SumOp
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_measure Integral over independent variable used by
    * applyGlobal.
    */
   explicit IntegralOp(
      double a_measure);

   /**
    * @brief Destructor.
    */
   virtual
   ~IntegralOp();

   /**
    * @brief Call SumOp's method and multiply result by m_measure.
    *
    * @param[out] a_result Final result of reduction.
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
   // Unimplemented default constructor.
   IntegralOp();

   // Unimplemented copy constructor.
   IntegralOp(
      const IntegralOp& other);

   // Unimplemented assignment operator.
   IntegralOp&
   operator = (
      const IntegralOp& rhs);

   // Integral over independent variable used by applyGlobal.
   double m_measure;
};

} // end namespace Loki

#endif
