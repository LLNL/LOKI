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
#ifndef _PITCH_ANGLE_COLLISION_OPERATOR_H_
#define _PITCH_ANGLE_COLLISION_OPERATOR_H_

#include "Overture.h"
#include "OvertureTypes.h"
#include "ParmParse.H"

#include "CollisionOperator.H"

#include <string>

namespace Loki {

/**
 * Implementation of CollisionOperator interface for a pitch angle collision.
 */
class PitchAngleCollisionOperator : public CollisionOperator
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_pp Where the input parameters are to be found.
    * @param[in] a_solution_order The order of the solution.
    */
   PitchAngleCollisionOperator(
      ParmParse& a_pp,
      int a_solution_order);

   /**
    * @brief Destructor.
    */
   virtual
   ~PitchAngleCollisionOperator();

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
   PitchAngleCollisionOperator();

   // Unimplemented copy constructor.
   PitchAngleCollisionOperator(
      const PitchAngleCollisionOperator& other);

   // Unimplemented assignment operator.
   PitchAngleCollisionOperator&
   operator = (
      const PitchAngleCollisionOperator& rhs);

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
      VFLOWX,
      VFLOWY,
      NUM_DPARAMS
   };

   // Indices of specific parameters in m_iparameters.
   enum {
      ICONS,
      SOLN_ORDER,
      NUM_IPARAMS
   };
};

} // end namespace Loki

#endif
