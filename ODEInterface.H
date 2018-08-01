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
#ifndef _ODE_INTEGRATOR_H_
#define _ODE_INTEGRATOR_H_

#include "Overture.h"

namespace Loki {

/**
 * Abstract base class that defines the interface to an ODE operator.  The
 * template parameter is the entity which must be time integrated.
 */
template <class TSoln>
class ODEInterface
{
public:
   /**
    * @brief Default constructor.
    */
   ODEInterface()
   {
   }

   /**
    * @brief Destructor.
    */
   virtual
   ~ODEInterface()
   {
   }

   /**
    * @brief Defines the TSoln data holder for data without ghost cells.
    *
    * @param[out] a_rhs Data holder for rhs (without ghost cells).
    * @param[in] a_prototype Prototype data holder.
    */
   virtual
   void
   defineRHSData(
      TSoln& a_rhs,
      const TSoln& a_prototype) = 0;

   /**
    * @brief Defines the TSoln data holder for data with ghost cells.
    *
    * @param[out] a_soln Data holder for solution with ghost cells.
    * @param[in] a_prototype Prototype data holder.
    */
   virtual
   void
   defineSolnData(
      TSoln& a_soln,
      const TSoln& a_prototype) = 0;

   /**
    * @brief Copies all data from srcSoln to dstSoln.
    *
    * @param[out] a_dstSoln Destination data holder.
    * @param[in] a_srcSoln Source data holder.
    */
   virtual
   void
   copySolnData(
      TSoln& a_dstSoln,
      const TSoln& a_srcSoln) = 0;

   // FIXME: soln should be const for an ODE but we are solving a PDE which
   // requires that its boundary conditions be set.  Instead of making a
   // copy of the soln inside evalRHS we are making soln non-const.
   /**
    * @brief Given the current time and solution, evaluates the RHS of the ODE.
    *
    * @param[out] a_rhs Data holder for rhs (without ghost cells).
    * @param[in] a_soln Current solution with which to evaluate RHS.
    * @param[in] a_time Current time.
    * @param[in] a_dt Current simulation time step.
    * @param[in] a_stage RK stage number.
    */
   virtual
   void
   evalRHS(
      TSoln& a_rhs,
      TSoln& a_soln,
      real   a_time,
      real   a_dt,
      int    a_stage) = 0;

   /**
    * @brief Checks for conformity between soln and protoSoln.
    *
    * @param[in] a_soln Data holder.
    * @param[in] a_protoSoln Prototype soln data holder.
    */
   virtual
   bool
   validSolnData(
      const TSoln& a_soln,
      const TSoln& a_protoSoln) = 0;

   /**
    * @brief Checks for conformity between soln and protoSoln ignoring ghost
    * cells.
    *
    * @param[in] a_rhs Data holder.
    * @param[in] a_protoSoln Prototype soln data holder.
    */
   virtual
   bool
   validRHSData(
      const TSoln& a_rhs,
      const TSoln& a_protoSoln) = 0;

   /**
    * @brief Sets the TSoln data to zero.
    *
    * @param[out] a_soln Data holder.
    */
   virtual
   void
   zeroSolnData(
      TSoln& a_soln) = 0;

   /**
    * @brief Adds (scale * increment) to the TSoln data holder.
    *
    * @param[out] a_soln Solution data to which is added.
    * @param[in] a_increment Amount to add.
    * @param[in] a_scale Multiplier by which to scale increment.
    * @param[in] a_final_rk If true, this has been called as part of the final
    *                       RK update and tracking particle positions need to
    *                       be limited.
    */
   virtual
   void
   addSolnData(
      TSoln&       a_soln,
      const TSoln& a_increment,
      real         a_scale = 1.0,
      bool         a_final_rk = false) = 0;

   /**
    * @brief Callback to allow for additional operations after stage advance.
    *
    * @param[in,out] a_soln Data holder.
    * @param[in] a_stage Which solution stage.
    */
   virtual
   void
   postStageAdvance(
      TSoln& a_soln,
      int a_stage) = 0;

private:
   // Unimplemented copy constructor.
   ODEInterface(
      const ODEInterface<TSoln>& other);

   // Unimplemented assignment operator.
   ODEInterface&
   operator = (
      const ODEInterface<TSoln>& rhs);

};

} // end namespace Loki

#endif
