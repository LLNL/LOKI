/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _PITCH_ANGLE_COLLISION_OPERATOR_H_
#define _PITCH_ANGLE_COLLISION_OPERATOR_H_

#include "LokiInputParser.H"
#include "CollisionOperator.H"
#include "tbox/Pointer.H"
#include "ReductionSchedule4D.H"

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
    * @param[in] a_species The species this operator is associated with.
    */
   PitchAngleCollisionOperator(
      LokiInputParser& a_pp,
      const KineticSpecies* a_species);

   /**
    * @brief Destructor.
    */
   virtual
   ~PitchAngleCollisionOperator();

   //@{

   //! @name Implementation of CollisionOperator interface.

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
      bool a_last_rk_stage);

   /**
    * @brief Compute the real part of lambda for the time step limit.
    *
    * @param[in] a_u The KineticSpecies at the current RK stage.
    * @return The real part of lambda for the time step limit.
    */
   virtual
   double
   computeRealLam(
      KineticSpecies& a_u);

   /**
    * @brief Print the operator's parameters.
    */
   virtual
   void
   printParameters() const;

   /**
    * @brief Copy collision diagnostic fields into supplied arrays of fields.
    *
    * @param[out] a_d Collision diagnostic fields to output.
    */
   virtual
   void
   copyDiagnosticFields(
      vector<ParallelArray>& a_d) const;

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
      const string& a_name);

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

   // Construct fundamental arrays used by operator.
   void
   constructArrays();

   // Construct the momenta and kinetic energy reduction schedules.
   void
   constructReductionSchedules();

   // Parse the input parameters.
   void
   parseParameters(
      LokiInputParser& a_pp,
      const KineticSpecies* a_species);

   // The name of this type of collision operator.
   static const string s_CLASS_NAME;

   // Indices of specific parameters in m_dparameters.
   enum {
      VFLOOR,
      VTHERMAL_DT,
      NU,
      NUM_DPARAMS
   };

   // Indices of specific parameters in m_iparameters.
   enum {
      ICONS,
      SOLN_ORDER,
      DO_RELATIVITY,
      NUM_IPARAMS
   };

   // Coordinates of the lower bounds of the range of the collision in velocity
   // space.
   vector<double> m_range_lo;

   // Coordinates of the upper bounds of the range of the collision in velocity
   // space.
   vector<double> m_range_hi;

   // Locally reduced 2D density.
   ParallelArray m_rN;

   // Locally reduced 2D x momentum.
   ParallelArray m_rGammax;

   // Locally reduced 2D y momentum.
   ParallelArray m_rGammay;

   // Locally reduced 2D kinetic energy.
   ParallelArray m_rsKEC;

   // Globally reduced 2D density.
   ParallelArray m_IN;

   // Globally reduced 2D x momentum.
   ParallelArray m_IGammax;

   // Globally reduced 2D y momentum.
   ParallelArray m_IGammay;

   // Globally reduced 2D kinetic energy.
   ParallelArray m_IsKEC;

   // Globally reduced 2D x flow velocity.
   ParallelArray m_IVx;

   // Globally reduced 2D y flow velocity.
   ParallelArray m_IVy;

   // Globally reduced 2D v_thermal.
   ParallelArray m_IVth;

   // Density reduction schedule
   tbox::Pointer<ReductionSchedule4D> m_N_reduction;

   // x momentum reduction schedule
   tbox::Pointer<ReductionSchedule4D> m_Gammax_reduction;

   // y momentum reduction schedule
   tbox::Pointer<ReductionSchedule4D> m_Gammay_reduction;

   // Kinetic energy reduction schedule
   tbox::Pointer<ReductionSchedule4D> m_skec_reduction;
};

} // end namespace Loki

#endif
