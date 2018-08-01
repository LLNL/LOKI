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
#ifndef _COLLISION_OPERATOR_H_
#define _COLLISION_OPERATOR_H_

#include "Overture.h"
#include "OvertureTypes.h"

#include "tbox/Box.H"
#include "ProblemDomain.H"

namespace Loki {

/**
 * CollisionOperator is an abstract base class for all CollisionOperator
 * implementations.
 */
class CollisionOperator
{
public:
   /**
    * @brief Default Constructor.
    */
   CollisionOperator();

   /**
    * @brief Destructor.
    */
   virtual
   ~CollisionOperator();

   /**
    * @brief Set the range of processors and their communicator for the
    *        KineticSpecies associated with this operator.
    *
    * @param[in] a_processor_range The range of processors for the
    *                              KineticSpecies associated with this
    *                              operator.
    * @param[in] a_comm The communicator for the KineticSpecies associated
    *                   with this operator.
    */
   void
   setCommunicationInfo(
      const Range& a_processor_range,
      const MPI_Comm& a_comm)
   {
      m_processor_range = a_processor_range;
      m_comm = a_comm;
   }

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
      real a_vflowy) = 0;

   /**
    * @brief Compute the real part of lambda for the time step limit.
    *
    * @param[in] a_dv The smaller velocity resolution.
    * @return The real part of lambda for the time step limit.
    */
   virtual
   real
   computeRealLam(
      real a_dv) const = 0;

   /**
    * @brief Print the operator's parameters.
    */
   virtual
   void
   printParameters() const = 0;

protected:
   /**
    * @brief Compute vthermal using an algorithm that is spatially local.  The
    * mean spatial distribution is computed locally and then used in a global
    * computation of vthermal.
    *
    * @param[in] a_u The distribution function at the current RK stage.
    * @param[in] a_local_box Box describing the bounds of a_rhs and a_u.
    * @param[in] a_interior_box Box describing the interior of a_rhs and a_u.
    * @param[in] a_domain Information about the problem domain.
    * @param[in] a_vflowx Current flow velocity in x.
    * @param[in] a_vflowy Current flow velocity in y.
    *
    * @return The computed vthermal.
    */
   real
   computeVthLocal(
      const RealArray& a_u,
      const tbox::Box& a_local_box,
      const tbox::Box& a_interior_box,
      const ProblemDomain& a_domain,
      real a_vflowx,
      real a_vflowy);

   /**
    * @brief Compute vthermal using an algorithm that is spatially global.  The
    * mean spatial distribution is computed globally and then used in a global
    * computation of vthermal.
    *
    * @param[in] a_u The distribution function at the current RK stage.
    * @param[in] a_local_box Box describing the bounds of a_rhs and a_u.
    * @param[in] a_interior_box Box describing the interior of a_rhs and a_u.
    * @param[in] a_domain Information about the problem domain.
    * @param[in] a_vflowx Current flow velocity in x.
    * @param[in] a_vflowy Current flow velocity in y.
    *
    * @return The computed vthermal.
    */
   real
   computeVthGlobal(
      const RealArray& a_u,
      const tbox::Box& a_local_box,
      const tbox::Box& a_interior_box,
      const ProblemDomain& a_domain,
      real a_vflowx,
      real a_vflowy);

   // The range of processors for the KineticSpecies associated with this
   // operator.
   Range m_processor_range;

   // The communicator for the KineticSpecies associated with this operator.
   MPI_Comm m_comm;

   // The double valued parameters necessary for the evaluation of this
   // collision operator.
   RealArray m_dparameters;

   // The integer valued parameters necessary for the evaluation of this
   // collision operator.
   IntegerArray m_iparameters;

   // Options for how to determine vthermal.
   enum {
      INPUT_VTHERMAL,
      LOCAL_VTHERMAL,
      GLOBAL_VTHERMAL
   };

   // How vthermal will be determined.
   // INPUT_VTHERMAL is the default to match historical usage.
   // LOCAL_VTHERMAL means to use a local computation of vthermal.
   // GLOBAL_VTHERMAL means to use a global computation of vthermal.
  int m_vthermal_method;

   // Total number of points in local velocity space.
   int m_vel_space_size;

   // Total number of points in global configuration space.
   int m_config_space_size;

   // Average of distribution function over configuration space.
   real* m_dist_func_avg;

   // Communicator for all KineticSpecies sharing the same part of
   // configuration space.
   MPI_Comm m_config_space_comm;

   // If 1, then this processor is a head node WRT configuration space.
   int m_is_config_space_head_node;

private:
   // Unimplemented copy constructor.
   CollisionOperator(
      const CollisionOperator& other);

   // Unimplemented assignment operator.
   CollisionOperator&
   operator = (
      const CollisionOperator& rhs);
};

} // end namespace Loki

#endif
