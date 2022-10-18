/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _COLLISION_OPERATOR_H_
#define _COLLISION_OPERATOR_H_

#include "ProblemDomain.H"
#include <vector>

namespace Loki {

class KineticSpecies;

/**
 * CollisionOperator is an abstract base class for all CollisionOperator
 * implementations.
 */
class CollisionOperator
{
public:
   /**
    * @brief The number of types of collision diagnostics--x and y momentum, KE,
    *        and entropy.
    */
   static const int s_DIAGNOSTIC_WORK_SIZE;

   /**
    * @brief Constructor.
    * @param[in] a_domain Information about the problem domain.
    * @param[in] a_num_dparams The number of double valued parameters,
    * @param[in] a_num_iparams The number of integer valued parameters,
    */
   CollisionOperator(
      const ProblemDomain& a_domain,
      int a_num_dparams,
      int a_num_iparams);

   /**
    * @brief Destructor.
    */
   virtual
   ~CollisionOperator();

   /**
    * @brief Initialize operator given information about supplies species.
    *
    * @param[in] a_species The species this operator is assiciated with.
    */
   virtual
   void
   initialize(
      KineticSpecies* a_species);

   /**
    * @brief Append the collision operator to the Vlasov rhs given the current
    *        distribution function.
    *
    * @param[out] a_rhs The KineticSpecies of the rhs of the Vlasov equation.
    * @param[in] a_u The KineticSpecies at the current RK stage.
    * @param[in] a_dt Current simulation time step.
    * @param[in] a_last_rk_stage True if called from the last RK stage.
    */
   virtual
   void
   evaluate(
      KineticSpecies& a_rhs,
      const KineticSpecies& a_u,
      double a_dt,
      bool a_last_rk_stage) = 0;
  
   /**
    * @brief Compute the real part of lambda for the time step limit.
    *
    * @param[in] a_u The KineticSpecies at the current RK stage.
    * @return The real part of lambda for the time step limit.
    */
   virtual
   double
   computeRealLam(
      KineticSpecies& a_u) = 0;

   /**
    * @brief Print the operator's parameters.
    */
   virtual
   void
   printParameters() const = 0;

   /**
    * @brief Copy collision diagnostic fields into supplied arrays of fields.
    *
    * @param[out] a_d Collision diagnostic fields to output.
    */
   virtual
   void
   copyDiagnosticFields(
      vector<ParallelArray>& a_d) const = 0;
  
protected:
   // Unimplemented default Constructor.
   CollisionOperator();

   // The lowest rank processor for the KineticSpecies associated with this
   // operator.
   int m_proc_lo;

   // The highest rank processor for the KineticSpecies associated with this
   // operator.
   int m_proc_hi;

   // The communictor for the KineticSpecies associated with this operator.
   MPI_Comm m_comm;

   // Box describing the interior plus ghosts/boundary layer of the
   // KineticSpecies this operator is associated with on this processor.
   ParallelArray::Box m_data_box;

   // Box describing the extent of the interior of the KineticSpecies this
   // operator is associated with on this processor.
   ParallelArray::Box m_interior_box;

   // The global problem domain.
   ProblemDomain m_domain;

   // For each "velocity cell", the local velocities at the "velocity cell"
   // center.
   ParallelArray m_velocities;

   // The double valued parameters necessary for the evaluation of this
   // collision operator.
   vector<double> m_dparameters;

   // The integer valued parameters necessary for the evaluation of this
   // collision operator.
   vector<int> m_iparameters;

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
