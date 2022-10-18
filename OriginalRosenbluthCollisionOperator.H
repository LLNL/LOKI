/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _ORIGINAL_ROSENBLUTH_COLLISION_OPERATOR_H_
#define _ORIGINAL_ROSENBLUTH_COLLISION_OPERATOR_H_

#include "LokiInputParser.H"
#include "CollisionOperator.H"
#include "tbox/Pointer.H"
#include "ReductionSchedule4D.H"

#include <string>

namespace Loki {

/**
 * Implementation of CollisionOperator interface for the Rosenbluth collision
 * operator.
 */
class OriginalRosenbluthCollisionOperator : public CollisionOperator
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_pp Where the input parameters are to be found.
    * @param[in] a_species The species this operator is associated with.
    */
   OriginalRosenbluthCollisionOperator(
      LokiInputParser& a_pp,
      const KineticSpecies* a_species);

   /**
    * @brief Destructor.
    */
   virtual
   ~OriginalRosenbluthCollisionOperator();

   //@{

   /**
    * @brief Initialize operator given information about supplies species.
    *
    * @param[in] a_species The species this operator is assiciated with.
    */
   virtual
   void
   initialize(
      KineticSpecies* a_species);

   //! @name Implementation of CollisionOperator interface.

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
   OriginalRosenbluthCollisionOperator();

   // Unimplemented copy constructor.
   OriginalRosenbluthCollisionOperator(
      const OriginalRosenbluthCollisionOperator& other);

   // Unimplemented assignment operator.
   OriginalRosenbluthCollisionOperator&
   operator = (
      const OriginalRosenbluthCollisionOperator& rhs);

  /**
    * @brief Compute vthermal using an algorithm that is spatially local.  The
    *        mean spatial distribution is computed locally and then used in a
    *        global computation of vthermal.
    *
    * @param[in] a_u The distribution function at the current RK stage.
    * @param[in] a_vflowx Current flow velocity in x.
    * @param[in] a_vflowy Current flow velocity in y.
    *
    * @return The computed vthermal.
    */
   double
   computeVthLocal(
      const ParallelArray& a_u,
      double a_vflowx,
      double a_vflowy);

   /**
    * @brief Compute vthermal using an algorithm that is spatially global.  The
    *        mean spatial distribution is computed globally and then used in a
    *        global computation of vthermal.
    *
    * @param[in] a_u The distribution function at the current RK stage.
    * @param[in] a_vflowx Current flow velocity in x.
    * @param[in] a_vflowy Current flow velocity in y.
    *
    * @return The computed vthermal.
    */
   double
   computeVthGlobal(
      const ParallelArray& a_u,
      double a_vflowx,
      double a_vflowy);

   // Construct arrays always needed by operator.
   void
   constructArrays();

   // Compute constant diffusion tensor term used by operator.
   void
   computeConstants();

   // Construct the various momenta and kinetic energy array data members used
   // in the back reaction calculation.
   void
   constructBackReactionArrays();

   // Construct the momenta and kinetic energy reduction schedules used in the
   // back reaction calculation.
   void
   constructBackReactionReductionSchedules();

   // Computes momenta and kinetic energy and their 2D reductions used in the
   // back reaction calculation.
   void
   computeBackReactionMoments();

   // Parse the input parameters.
   void
   parseParameters(
      LokiInputParser& a_pp,
      const KineticSpecies* a_species);

   // The name of this type of collision operator.
   static const string s_CLASS_NAME;

   // Indices of specific parameters in m_dparameters.
   enum {
      VTHERMAL,
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
      KERNEL_TYPE,
      DO_RELATIVITY,
      NUM_IPARAMS
   };

   // Options for how to determine vthermal.
   enum {
      INPUT_VTHERMAL,
      LOCAL_VTHERMAL,
      GLOBAL_VTHERMAL
   };

   // Total number of points in local velocity space.
   int m_vel_space_size;

   // Total number of points in global configuration space.
   int m_config_space_size;

   // Average of distribution function over configuration space.
   double* m_dist_func_avg;

   // Communicator for all KineticSpecies sharing the same part of
   // configuration space.
   MPI_Comm m_config_space_comm;

   // If 1, then this processor is a head node WRT configuration space.
   int m_is_config_space_head_node;

   // How vthermal will be determined.
   // INPUT_VTHERMAL is the default to match historical usage.
   // LOCAL_VTHERMAL means to use a local computation of vthermal.
   // GLOBAL_VTHERMAL means to use a global computation of vthermal.
   int m_vthermal_method;

   // Coordinates of the lower bounds of the range of the collision in velocity
   // space.
   vector<double> m_range_lo;

   // Coordinates of the upper bounds of the range of the collision in velocity
   // space.
   vector<double> m_range_hi;

   // Flag indicating if the back reaction should be computed.
   bool m_back_reaction;

   // Persistent 4D (unreduced) x momentum.
   ParallelArray m_cMomx;

   // Persistent 4D (unreduced) y momentum.
   ParallelArray m_cMomy;

   // Persistent 4D (unreduced) kinetic energy.
   ParallelArray m_cKE;

   // Persistent 4D maxwellian distribution.
   ParallelArray m_fM;

   // Intermediate 2D (locally reduced) x momentum.
   ParallelArray m_rMomx;

   // Intermediate 2D (locally reduced) y momentum.
   ParallelArray m_rMomy;

   // Intermediate 2D (locally reduced) kinetic energy.
   ParallelArray m_rKE;

   // Persistent 2D (globally reduced) x momentum denominator term.
   ParallelArray m_IMomxD;

   // Persistent 2D (globally reduced) y momentum denominator term.
   ParallelArray m_IMomyD;

   // Persistent 2D (globally reduced) kinetic energy denominator term.
   ParallelArray m_IKED;

   // Persistent 2D (globally reduced) x momentum numerator term.
   ParallelArray m_IMomxN;

   // Persistent 2D (globally reduced) y momentum numerator term.
   ParallelArray m_IMomyN;

   // Persistent 2D (globally reduced) kinetic energy numerator term.
   ParallelArray m_IKEN;

   // Persistent diffusion tensor.
   ParallelArray m_d;

   // x momentum reduction schedule
   tbox::Pointer<ReductionSchedule4D> m_x_mom_reduction;

   // y momentum reduction schedule
   tbox::Pointer<ReductionSchedule4D> m_y_mom_reduction;

   // kinetic energy reduction schedule
   tbox::Pointer<ReductionSchedule4D> m_ke_reduction;
};

} // end namespace Loki

#endif
