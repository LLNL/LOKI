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
#ifndef _KERNEL_OP_H_
#define _KERNEL_OP_H_

#include "Overture.h"
#include "tbox/Box.H"

namespace Loki {

/**
 * KernelOp is an abstract base class defining the interface to an operator
 * that transforms a function.
 */
class KernelOp
{
public:
   /**
    * @brief Constructor.
    */
   KernelOp();

   /**
    * @brief Destructor.
    */
   virtual
   ~KernelOp();

   /**
    * @brief Perform the local part of the transformation.
    *
    * @param[out] a_result The local result.
    * @param[in] a_source The local source.
    */
   virtual
   void
   applyLocal(
      RealArray& a_result,
      const RealArray& a_source) const = 0;

   /**
    * @brief Perform the global part of the transformation.
    *
    * @param[out] a_result The global result.
    */
   virtual
   void
   applyGlobal(
      RealArray& a_result) const = 0;

private:
   // Unimplemented copy constructor.
   KernelOp(
      const KernelOp& rhs);

   // Unimplemented assignment operator.
   KernelOp&
   operator = (
      const KernelOp& rhs);
};


/**
 * DensityKernal implements the interface of KernelOp for the transformation of
 * a distribution function of entities to a distribution function of some
 * underlying property like charge or mass.
 */
class DensityKernel : public KernelOp
{
public:
   /**
    * @brief Constructor.
    */
   explicit DensityKernel(
      double a_mass);

   // Copy constructor.
   DensityKernel(
      const DensityKernel& other);

   /**
    * @brief Destructor.
    */
   ~DensityKernel();

   /**
    * @brief Perform the local part of the transformation.
    *
    * For DensityKernel this just assigns source to result.
    *
    * @param[out] a_result The local result.
    * @param[in] a_source The local source.
    */
   virtual
   void
   applyLocal(
      RealArray& a_result,
      const RealArray& a_source) const;

   /**
    * @brief Perform the global part of the transformation.
    *
    * For DensityKernel this just multiplies a_result by m_mass.
    *
    * @param[out] a_result The global result.
    */
   virtual
   void
   applyGlobal(
      RealArray& a_result) const;

private:
   // Unimplemented default constructor.
   DensityKernel();

   // Unimplemented assignment operator.
   DensityKernel&
   operator = (
      const DensityKernel& rhs);

   // "Mass" to multiply by in applyGlobal.
   real m_mass;

};

} // end namespace Loki

#endif
