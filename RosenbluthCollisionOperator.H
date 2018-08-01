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
#ifndef _ROSENBLUTH_COLLISION_OPERATOR_H_
#define _ROSENBLUTH_COLLISION_OPERATOR_H_

#include "Overture.h"
#include "OvertureTypes.h"
#include "ParmParse.H"

#include "CollisionOperator.H"
#include "tbox/Pointer.H"
#include "ReductionSchedule4D.H"

#include <string>

namespace Loki {

/**
 * Implementation of CollisionOperator interface for the Rosenbluth collision
 * operator.
 */
class RosenbluthCollisionOperator : public CollisionOperator
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_pp Where the input parameters are to be found.
    * @param[in] a_solution_order The order of the solution.
    */
   RosenbluthCollisionOperator(
      ParmParse& a_pp,
      int a_solution_order);

   /**
    * @brief Destructor.
    */
   virtual
   ~RosenbluthCollisionOperator();

   //@{

   //! @name Implementation of CollisionOperator interface.

   /**
    * @brief Append the collision operator to the Vlasov rhs given the current
    *        distribution function.
    *
    * @param[out] a_rhs The distribution function of the rhs of the Vlasov
    *                   equation.
    * @param[in] a_u The distribution function at the current RK stage.
    * @param[in] a_local_box Box describing the bounds of a_rhs and a_u.
    * @param[in] a_interior_box Box describing the interior of a_rhs and a_u.
    * @param[in] a_domain Information about the problem domain.
    * @param[in] a_vflowx Current flow velocity in x.
    * @param[in] a_vflowy Current flow velocity in y.
    */
   virtual
   void
   evaluate(
      RealArray& a_rhs,
      const RealArray& a_u,
      const tbox::Box& a_local_box,
      const tbox::Box& a_interior_box,
      const ProblemDomain& a_domain,
      real a_vflowx,
      real a_vflowy);

   /**
    * @brief Compute the real part of lambda for the time step limit.
    *
    * @param[in] a_dv The smaller velocity resolution.
    * @return The real part of lambda for the time step limit.
    */
   virtual
   real
   computeRealLam(
      real a_dv) const;

   /**
    * @brief Print the operator's parameters.
    */
   virtual
   void
   printParameters() const;

   //@}

   /**
    * @brief Returns true if name describes this class' type.
    *
    * @param[in] a_name Name of a type of collision operator.
    *
    * @return true if name describes this class' type.
    */
   static
   bool
   isType(
      const aString& a_name);

private:
   // Unimplemented default constructor.
   RosenbluthCollisionOperator();

   // Unimplemented copy constructor.
   RosenbluthCollisionOperator(
      const RosenbluthCollisionOperator& other);

   // Unimplemented assignment operator.
   RosenbluthCollisionOperator&
   operator = (
      const RosenbluthCollisionOperator& rhs);

   // Compute constant vx, vy, Khat, and Hhat terms used by operator.
   void
   computeConstants(
      const tbox::Box& a_interior_box,
      const ProblemDomain& a_domain);

   // Construct the various momenta and kinetic energy array data members used
   // in the back reaction calculation.
   void
   constructBackReactionArrays(
      const tbox::Box& a_local_box);

   // Construct the momenta and kinetic energy reduction schedules used in the
   // back reaction calculation.
   void
   constructBackReactionReductionSchedules(
      const tbox::Box& a_interior_box,
      const ProblemDomain& a_domain);

   // Computes momenta and kinetic energy and their 2D reductions used in the
   // back reaction calculation.
   void
   computeBackReactionMoments(
      const tbox::Box& a_local_box,
      const tbox::Box& a_interior_box,
      const ProblemDomain& a_domain);

   // Parse the input parameters.
   void
   parseParameters(
      ParmParse& a_pp);

   // The name of this type of collision operator.
   static const aString s_CLASS_NAME;

   // Indices of specific parameters in m_dparameters.
   enum {
      VCEIL,
      VFLOOR,
      VTHERMAL,
      VTHERMAL_DT,
      NU,
      ALPHA,
      MASSR,
      VFLOWX,
      VFLOWY,
      NUM_DPARAMS
   };

   // Indices of specific parameters in m_iparameters.
   enum {
      SOLN_ORDER,
      NUM_GHOSTS,
      NUM_IPARAMS
   };


   // Flag indicating if the back reaction should be computed.
   bool m_back_reaction;

   // Flag indicating that the moments need to be computed.
   bool m_compute_moments;

   // Flag indicating that the constants need to be computed.
   bool m_compute_constants;

   // Persistent 4D (unreduced) x momentum.
   RealArray m_cMomx;

   // Persistent 4D (unreduced) y momentum.
   RealArray m_cMomy;

   // Persistent 4D (unreduced) kinetic energy.
   RealArray m_cKE;

   // Persistent 4D maxwellian distribution.
   RealArray m_fM;

   // Intermediate 2D (locally reduced) x momentum.
   RealArray m_rMomx;

   // Intermediate 2D (locally reduced) y momentum.
   RealArray m_rMomy;

   // Intermediate 2D (locally reduced) kinetic energy.
   RealArray m_rKE;

   // Persistent 2D (globally reduced) x momentum denominator term.
   RealArray m_IMomxD;

   // Persistent 2D (globally reduced) y momentum denominator term.
   RealArray m_IMomyD;

   // Persistent 2D (globally reduced) kinetic energy denominator term.
   RealArray m_IKED;

   // Persistent 2D (globally reduced) x momentum numerator term.
   RealArray m_IMomxN;

   // Persistent 2D (globally reduced) y momentum numerator term.
   RealArray m_IMomyN;

   // Persistent 2D (globally reduced) kinetic energy numerator term.
   RealArray m_IKEN;

   // Persistent Rosenbluth potential K/v**2.
   RealArray m_khat;

   // Persistent Rosenbluth potential H/v**2.
   RealArray m_hhat;

   // Persistent x velocities for Rosenbluth collision calculation.
   RealArray m_vx;

   // Persistent y velocities for Rosenbluth collision calculation.
   RealArray m_vy;

   // x momentum reduction schedule
   tbox::Pointer<ReductionSchedule4D> m_x_mom_reduction;

   // y momentum reduction schedule
   tbox::Pointer<ReductionSchedule4D> m_y_mom_reduction;

   // kinetic energy reduction schedule
   tbox::Pointer<ReductionSchedule4D> m_ke_reduction;
};

} // end namespace Loki

#endif
