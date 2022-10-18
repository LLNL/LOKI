/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _ROSENBLUTH_COLLISION_OPERATOR_H_
#define _ROSENBLUTH_COLLISION_OPERATOR_H_

#include "LokiInputParser.H"
#include "CollisionOperator.H"
#include "Loki_Utilities.H"
#include "tbox/Pointer.H"
#include "ReductionSchedule4D.H"

#include <map>
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
    * @param[in] a_species The species this operator is associated with.
    */
   RosenbluthCollisionOperator(
      LokiInputParser& a_pp,
      const KineticSpecies* a_species);

   /**
    * @brief Destructor.
    */
   virtual
   ~RosenbluthCollisionOperator();

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
   /**
    * Class to encapsulate the interspecies communication of the Rosenbluth
    * collision operator.  The interspecies communication is a 3 step process
    * and relies on a few conventions.  First, the "head process" for each
    * species is just the lowest rank processor for that species.  Secondly,
    * each species is partitioned in all 4 dimensions but in general there are
    * multiple processors associated with a given partition of configuration
    * space.  Third, interspecies communication consists of data reduced over
    * velocity space.  This means that the data for a given partition of
    * configuration space is identical on all processors associated with this
    * partition.  For a given configuration space partition we elect to use the
    * data from the processor that contains (0,0) in velocity space.
    *
    * The 3 steps of the interspecies communication are as follows:
    * 1)  The head processor of each species collects the communicated data for
    *     each configuration space partition in order to form the communicated
    *     data for the entire configuration space.
    * 2)  The head processor of each species communicates the data it has just
    *     formed to each other species' head processor.  Now each species' head
    *     processor has the data for all the other species.
    * 3)  Looping over each species interspecies collision, the data for the
    *     other species in this interaction is communicated from the head
    *     processor to all the other processors.
    *
    * The nature of the data communicated for the forward and back reactions are
    * different so there are differences in the details of the 3 steps outlined
    * above in these 2 cases.  For the forward reaction, each species computes
    * and communicates to all the other species 4 fields describing itself.  For
    * the back reaction, each species computes 3 back reaction numerators for
    * each other species with which it interacts--3*(number of species-1) fields
    * in all.  Unlike the 4 forward reaction fields which are all needed by each
    * other species, the back reaction terms are specific to each other species.
    * Hence, for the back reaction only the 3 fields specific to each other
    * species are communicated to that other species.
    */
   class RosenbluthCommunicator
   {
   public:
      /**
       * @brief Constructor
       *
       * @param[in] a_species KineticSpecies the operator is associated with.
       * @param[in] a_maxNumFields Maximum number of fields communicated.
       * @param[in] a_nu The collision coefficients.
       * @param[in] a_num_interspecies_collisions The number of other
       *                                          interspecies collisions.
       */
      RosenbluthCommunicator(
         const KineticSpecies* a_species,
         int a_maxNumFields,
         const vector<double>& a_nu,
         int a_num_interspecies_collisions);

      /**
       * @brief Destructor.
       */
      virtual
      ~RosenbluthCommunicator();

      /**
       * @brief Send species interior data to the head processor of the species
       *        and construct full species domain data.
       *
       * @param[in] a_BR If true this communication is for the back reaction.
       */
      void
      CommunicateSpatialDataToSpeciesHead(
         bool a_BR);

      /**
       * @brief Send full species domain data to head processors of the other
       *        species.
       *
       * @param[in] a_BR If true this communication is for the back reaction.
       */
      void
      CommunicateSpatialDataToOtherSpecies(
         bool a_BR);

      /**
       * @brief Send the species interior data for the species with index
       *        a_otherSpecies to each non-head processor of this species.
       *
       * @param[in] a_otherSpecies Index of the other species in the collision.
       * @param[in] a_BR If true this communication is for the back reaction.
       */
      void
      CommunicateSpatialDataFromSpeciesHeadFromOtherSpecies(
         int a_otherSpecies,
         bool a_BR);

      /**
       * @brief Returns the interior data for this processor.
       */
      ParallelArray&
      MyInteriorData()
      {
         return m_interiorData[0];
      }


      /**
       * @brief Ensure that the head processor sends in
       *        CommunicateSpatialDataFromSpeciesHeadFromOtherSpecies have all
       *        completed.
       */
      void
      completeSends()
      {
         if (Loki_Utilities::s_my_id == m_speciesHead) {
            int numSends = m_proc_hi - m_speciesHead;
            MPI_Waitall(numSends, m_request, m_status);
         }
      }

   private:
      // This represents the interior data of the communicated data for each
      // unique partition of configuration space of the species associated with
      // the operator.  On the species head processor, there is an entry for
      // each unique partition of configuration space.  On non-head processors
      // is only one entry needed and that is for this processor's there
      // configuration space partition.
      vector<ParallelArray> m_interiorData;

      // This represents the range of each partition of configuration space.
      // The comments WRT m_interiorData about how many entries are needed for a
      // given processor apply to this data member.
     vector<vector<int> > m_interiorDataRange;

      // This represents the communicated data for the entire configuration
      // space.  It is only relevant on the species head processor where there
      // is 1 entry for each species.
      vector<ParallelArray> m_domainData;

      // This is the processor ID of the head processor of each species.
      vector<int> m_speciesHeads;

      // For each processor of the operator's species this indicates if the head
      // processor will receive data from it during a call to
      // CommunicateSpatialDataToSpeciesHead.  It is only relevant to the head
      // processor.  There are multiple processors with the same configuration
      // space partition.  As each of these should contain the same data to be
      // commiunicated within roundoff it is only necessary for one to send its
      // m_interiorData to the head process.  The processor containing (0,0) in
      // velocity space will send data, all others will not.
      deque<bool> m_receiveFrom;

      // For each true entry in m_receiveFrom, the ID of the processor being
      // received from.
      vector<int> m_receiveFromPID;

      // Pretty much the inverse of m_receiveFrom and, as such, it is only
      // relevant on the non-head processors.
      bool m_sendConfigData;

      // The highest rank processor for the KineticSpecies associated with this
      // operator.
      int m_proc_hi;

      // The species index of the KineticSpecies associated with this operator.
      int m_species_index;

      // The number of species in the problem.
      int m_num_species;

      // Collision coefficients for interactions with each species.  The self
      // collision coefficient for species i is the ith entry.  The other
      // entries are the collision coefficient with the species with that index.
      vector<double> m_nu;

      // The number of other species with which there is a collision.
      int m_num_interspecies_collisions;

      // The head processor for the KineticSpecies associated with this
      // operator.
      int m_speciesHead;

      // The maximum number of fields to be communicated.
      int m_max_num_communicated_fields;

      // Just the number of true entries in m_receiveFrom to avoid counting this
      // all the time.
      int m_num_species_proc_received_from;

      // An array of requests to use in the non-blocking communication.
      MPI_Request* m_request;

      // An array of statuses to use in the non-blocking communication,
      MPI_Status* m_status;

      // The index of the configuration space partition for each procesor for
      // the KineticSpecies associated with this operator.
      map<int, int> m_pid_to_index;

      // The processors for each index of the configuration space partition.
      // The inverse of m_pid_to_index.
      multimap<int, int> m_index_to_pids;

      // Tag for data range communication.
      static const int s_TAG_SETUP;

      // Tag for spatial data to head node communication.
      static const int s_TAG_SPATIAL_TO_HEAD;

      // Tag for spatial data from head node to other species' head nodes.
      static const int s_TAG_SPATIAL_TO_OTHER_HEAD;

      // Tag for other species' spatial data from head node to other nodes.
      static const int s_TAG_OTHER_SPATIAL_FROM_HEAD;
   };

   // Unimplemented default constructor.
   RosenbluthCollisionOperator();

   // Unimplemented copy constructor.
   RosenbluthCollisionOperator(
      const RosenbluthCollisionOperator& other);

   // Unimplemented assignment operator.
   RosenbluthCollisionOperator&
   operator = (
      const RosenbluthCollisionOperator& rhs);

   // Construct fundamental arrays used by operator.
   void
   constructArrays();

   // Construct the various momenta and kinetic energy array data members used
   // in the back reaction calculation.
   void
   constructBackReactionArrays();

   // Construct arrays needed only for collision diagnostics.
   void
   constructDiagnosticArrays();

   // Construct the momenta and kinetic energy reduction schedules.
   void
   constructReductionSchedules();

   // Computes momenta and kinetic energy and their 2D reductions used in the
   // back reaction calculation.
   void
   computeBackReactionMoments(
      const ParallelArray& a_d);

   // Computes the collision diagnostics.
   void
   computeDiagnosticFields(
      const ParallelArray& a_rhs,
      const ParallelArray& a_u,
      double a_dt);

   // Parse the input parameters.
   void
   parseParameters(
      LokiInputParser& a_pp,
      const KineticSpecies* a_species);

   // Perform the portions of evaluate necessary to compute Re(lambda).  This is
   // essentially the computation of all forward reaction diffusion tensors and
   // everything that they depend on.
   void
   evaluateReLam(
      const KineticSpecies& a_u);

   // Performs the forward collision with the species whose species index is
   // a_other_species.
   void
   evaluateInterspecies(
      ParallelArray& a_rhs,
      const ParallelArray& a_u,
      int a_other_species,
      bool a_last_interspecies_collision,
      double a_mass_ratio,
      double a_dt,
      bool a_last_rk_stage);

   // Computes the back reaction from the species whose species index is
   // a_other_species.
   void
   evaluateInterspeciesBackReaction(
      ParallelArray& a_rhs,
      const ParallelArray& a_u,
      bool a_last_interspecies_collision,
      double a_mass_ratio,
      double a_dt,
      bool a_last_rk_stage);

   // Fills m_rosenbluthComm->MyInteriorData() with information about this
   // species which will be communicated to the other species so that they may
   // perform the forward interspecies collisions.
   void
   getInterspeciesData() const;

   // Fills m_rosenbluthComm->MyInteriorData() with back reaction information
   // for this species with every other species which will be communicated to
   // the other species so that they may perform the back reaction.
   void
   getInterspeciesBRData() const;

   // Retrieve from m_rosenbluthComm->MyInteriorData() information about the
   // species with species index a_which_interspecies_collision prior to
   // evaluating that interspecies forward collision.
   void
   setOtherSpeciesData(
      int a_which_interspecies_collision,
      bool a_computeReLamOnly,
      double mboverma);

   // Retrieve from m_rosenbluthComm->MyInteriorData() information about the
   // species with species index a_which_interspecies_collision prior to
   // evaluating that interspecies back reaction.
   void
   setOtherSpeciesBRData(
      int a_which_interspecies_collision,
      double mboverma);

   // The name of this type of collision operator.
   static const string s_CLASS_NAME;

   // Indices of specific parameters in m_dparameters.
   enum {
      NU,
      NUM_DPARAMS
   };

   // Indices of specific parameters in m_iparameters.
   enum {
      SOLN_ORDER,
      NUM_GHOSTS,
      KERNEL_TYPE,
      DENSITY_DEP_NU,
      DO_RELATIVITY,
      NUM_IPARAMS
   };

   // Flag indicating if the diagnostic data should be computed and output.
   bool m_diagnostic;

   // Flag indicating if self-collisions should be computed.
   bool m_selfcoll;

   // Indicates whether interspecies collisions are to be performed.
   bool m_interspecies;

   // Collision coefficients for interactions with each species.  The self
   // collision coefficient for species i is the ith entry.  The other entries
   // are the collision coefficient with the species with that index.
   vector<double> m_nu;

   // Flag indicating if the back reactions should be computed.
   bool m_back_reaction;

   // Coordinates of the lower bounds of the range of the collision in velocity
   // space.
   vector<double> m_range_lo;

   // Coordinates of the upper bounds of the range of the collision in velocity
   // space.
   vector<double> m_range_hi;

   // 4D x momentum temps
   ParallelArray m_uMomx;

   // 4D y momentum temps
   ParallelArray m_uMomy;

   // 4D kinetic energy temps
   ParallelArray m_uKE;

   // 4D (unreduced) x momentum.
   ParallelArray m_cMomx;

   // 4D (unreduced) y momentum.
   ParallelArray m_cMomy;

   // 4D (unreduced) kinetic energy.
   ParallelArray m_cKE;

   // 2D (locally itegrated over velocity space) x momentum.
   ParallelArray m_rMomx;

   // 2D (locally itegrated over velocity space) y momentum.
   ParallelArray m_rMomy;

   // 2D (locally itegrated over velocity space) kinetic energy.
   ParallelArray m_rKE;

   // 2D x momentum denominator term.  Reduction of m_rMomx.
   ParallelArray m_IMomxD;

   // 2D y momentum denominator term.  Reduction of m_rMomy.
   ParallelArray m_IMomyD;

   // 2D kinetic energy denominator term.  Reduction of m_rKE.
   ParallelArray m_IKED;

   // 2D x momentum numerator term for each species collision (self and
   // interspecies).  Reduction of m_rMomx.
   vector<ParallelArray> m_IMomxN;

   // 2D x momentum numerator term communicated from another species.  This is
   // updated for each back reaction species pair.
   ParallelArray m_IMomxNb;

   // 2D y momentum numerator term for each species collision (self and
   // interspecies).  Reduction of m_rMomy.
   vector<ParallelArray> m_IMomyN;

   // 2D y momentum numerator term communicated from another species.  This is
   // updated for each back reaction species pair.
   ParallelArray m_IMomyNb;

   // 2D kinetic energy numerator term for each species collision (self and
   // interspecies).  Reduction of m_rKE.
   vector<ParallelArray> m_IKEN;

   // 2D kinetic energy numerator term communicated from another species.  This
   // is updated for each back reaction species pair.
   ParallelArray m_IKENb;

   // 2D diagnostic fields.  Reductions of m_rMomx, m_rMomy, m_rKE, and m_rN
   // respectively.
   ParallelArray m_dMomx, m_dMomy, m_dKE, m_dEnt;

   // 4D helper field for diagnostics.
   ParallelArray m_d4D;

   // x momentum reduction schedule to reduce m_rMomx.
   tbox::Pointer<ReductionSchedule4D> m_x_mom_reduction;

   // y momentum reduction schedule to reduce m_rMomy.
   tbox::Pointer<ReductionSchedule4D> m_y_mom_reduction;

   // kinetic energy reduction schedule to reduce m_rKE.
   tbox::Pointer<ReductionSchedule4D> m_ke_reduction;

   // 2D (locally integrated over velocity space) density.
   ParallelArray m_rN;

   // 2D (locally integrated over velocity space) x momentum.
   ParallelArray m_rGammax;

   // 2D (locally integrated over velocity space) y momentum.
   ParallelArray m_rGammay;

   // 2D (locally integrated over velocity space) kinetic energy.
   ParallelArray m_rsKEC;

   // 2D density.  Reduction of m_rN.
   ParallelArray m_IN;

   // 2D density communicated from another species.  This is updated for each
   // forward reaction species pair.
   ParallelArray m_INb;

   // 2D x momentum.  Reduction of m_rGammax.
   ParallelArray m_IGammax;

   // 2D y momentum.  Reduction of m_rGammay.
   ParallelArray m_IGammay;

   // 2D kinetic energy.  Reduction of m_rsKEC,
   ParallelArray m_IsKEC;

   // 2D x flow velocity.
   ParallelArray m_IVx;

   // 2D x flow velocity communicated from another species.  This is updated for
   // each forward reaction species pair.
   ParallelArray m_IVxb;

   // 2D y flow velocity.
   ParallelArray m_IVy;

   // 2D y flow velocity communicated from another species.  This is updated for
   // each forward reaction species pair.
   ParallelArray m_IVyb;

   // 2D thermal velocity.
   ParallelArray m_IVth;

   // 2D thermal velocity communicated from another species.  This is updated
   // for each forward reaction species pair.
   ParallelArray m_IVthb;

   // Temporary variable for the result of forward collision.  Needed to compute
   // back reaction numerators.
   ParallelArray m_coll;

   // Density reduction schedule
   tbox::Pointer<ReductionSchedule4D> m_N_reduction;

   // x momentum reduction schedule
   tbox::Pointer<ReductionSchedule4D> m_Gammax_reduction;

   // y momentum reduction schedule
   tbox::Pointer<ReductionSchedule4D> m_Gammay_reduction;

   // Kinetic energy reduction schedule
   tbox::Pointer<ReductionSchedule4D> m_skec_reduction;

   // Diffusion tensors for self collision forward reaction, self collision back
   // reaction, interspecies forward reaction, and interspecies back reaction
   // respectively.
   ParallelArray m_d_a, m_d_aa, m_d_ba, m_d_ab;

   // The real part of the forward collision eigenvalue at each point in phase
   // space.
   ParallelArray m_reLam;

   // The maximum of the real part of the forward collision eigenvalue over all
   // of phase space.
   double m_reLamMax;

   // The maximum number of fields communicated for interspecies collisions.
   int m_maxNumCommunicatedFields;

   // Index of species this operator is associated with,
   int m_species_index;

   // Number of species in problem.
   int m_num_species;

   // The number of other species with which there is a collision.
   int m_num_interspecies_collisions;

   // The mass of each species in the problem.
   vector<double> m_species_masses;

   // The entity that handles all interspecies communication for this operator.
   tbox::Pointer<RosenbluthCommunicator> m_rosenbluthComm;
};

} // end namespace Loki

#endif
