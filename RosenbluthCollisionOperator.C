/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "RosenbluthCollisionOperator.H"
#include "RosenbluthCollisionOperatorF.H"
#include "KineticSpecies.H"
#include "TimerManager.H"

namespace Loki {

const string
RosenbluthCollisionOperator::s_CLASS_NAME(
   "Rosenbluth Collision Operator");
  
const int RosenbluthCollisionOperator::RosenbluthCommunicator::s_TAG_SETUP = 4375;

const int RosenbluthCollisionOperator::RosenbluthCommunicator::s_TAG_SPATIAL_TO_HEAD = 4376;

const int RosenbluthCollisionOperator::RosenbluthCommunicator::s_TAG_SPATIAL_TO_OTHER_HEAD = 4377;

const int RosenbluthCollisionOperator::RosenbluthCommunicator::s_TAG_OTHER_SPATIAL_FROM_HEAD = 4378;


RosenbluthCollisionOperator::RosenbluthCollisionOperator(
   LokiInputParser& a_pp,
   const KineticSpecies* a_species)
   : CollisionOperator(a_species->domain(), NUM_DPARAMS, NUM_IPARAMS),
     m_nu(a_species->numSpecies()),
     m_range_lo(2),
     m_range_hi(2),
     m_maxNumCommunicatedFields(0),
     m_species_index(a_species->speciesIndex()),
     m_num_species(a_species->numSpecies()),
     m_num_interspecies_collisions(0)
{
   // Set what we know and read from user input what we do not know.
   int solution_order = a_species->spatialSolutionOrder();
   m_iparameters[SOLN_ORDER] = solution_order;
   m_iparameters[NUM_GHOSTS] = (solution_order == 4) ? 2 : 3;
   m_iparameters[DO_RELATIVITY] = Simulation::s_DO_RELATIVITY ? 1 : 0;
   parseParameters(a_pp, a_species);
   if (m_interspecies) {
      if (m_back_reaction) {
         m_maxNumCommunicatedFields = max(4, 3*m_num_interspecies_collisions);
      }
      else {
         m_maxNumCommunicatedFields = 4;
      }
   }
   if (m_back_reaction) {
      m_IMomxN.resize(m_num_species);
      m_IMomyN.resize(m_num_species);
      m_IKEN.resize(m_num_species);
   }
}


RosenbluthCollisionOperator::~RosenbluthCollisionOperator()
{
}


void
RosenbluthCollisionOperator::initialize(
   KineticSpecies* a_species)
{
   // Always get the information about how the species this operator is
   // associated with is partitioned.  Then if diagnostics are on, construct
   // the diagnostic arrays as all processors must have these.
   CollisionOperator::initialize(a_species);
   if (m_diagnostic) {
      constructDiagnosticArrays();
   }

   // The other arrays and schedules are only needed on processors that the
   // species is partitioned among.
   if (a_species->isInRange(Loki_Utilities::s_my_id)) {
      m_species_masses = a_species->speciesMasses();
      if (m_interspecies) {
         m_rosenbluthComm = new RosenbluthCommunicator(a_species,
            m_maxNumCommunicatedFields,
            m_nu,
            m_num_interspecies_collisions);
      }

      constructArrays();
      constructReductionSchedules();
      if (m_back_reaction) {
         constructBackReactionArrays();
      }
   }
}


void
RosenbluthCollisionOperator::evaluate(
   KineticSpecies& a_rhs,
   const KineticSpecies& a_u,
   double a_dt,
   bool a_last_rk_stage)
{
   // This routine will compute quantities needed throughout the operator.  If
   // self collisions are to be computed both the forward and back self
   // collision reactions are done here.  If interspecies collision are to be
   // done the forward and, if requested, the back interspecies reactions are
   // also done here.

   ParallelArray& rhs_distribution = a_rhs.distribution();
   const ParallelArray& current_distribution = a_u.distribution();

   // For diagnostics, set m_d4D to rhs_distribution, which is the initial rhs
   // value.
   if (m_diagnostic) {
      m_d4D = rhs_distribution;
   }

   // If neither of these is set the operator is a no-op so just return.
   if (!m_selfcoll && !m_interspecies) {
      return;
   }

   // We need several quantities computed whether we're doing self or
   // interspecies collisions.
   m_rN = 0.0;
   m_rGammax = 0.0;
   m_rGammay = 0.0;
   m_rsKEC = 0.0;

   // m_rsKEC is passed as a dummy and will be computed below.
   FORT_COMPUTE_ROSENBLUTH_SPECIES_MOMENTS(*m_rN.getData(),
      *m_rGammax.getData(),
      *m_rGammay.getData(),
      *m_rsKEC.getData(),
      *current_distribution.getData(),
      0,
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      *m_velocities.getData());

   m_IN = 0.0;
   m_IGammax = 0.0;
   m_IGammay = 0.0;

   m_N_reduction->execute(m_IN);
   m_Gammax_reduction->execute(m_IGammax);
   m_Gammay_reduction->execute(m_IGammay);
   FORT_COMPUTE_ROSENBLUTH_SPECIES_REDUCED_FIELDS(*m_IVx.getData(),
      *m_IVy.getData(),
      *m_IN.getData(),
      *m_IGammax.getData(),
      *m_IGammay.getData(),
      BOX4D_TO_FORT(m_data_box));

   m_rsKEC = 0.0;
   FORT_COMPUTE_ROSENBLUTH_SPECIES_KEC(*m_rsKEC.getData(),
      *m_IVx.getData(),
      *m_IVy.getData(),
      *current_distribution.getData(),
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      *m_velocities.getData());

   m_IsKEC = 0.0;
   m_skec_reduction->execute(m_IsKEC);

   FORT_COMPUTE_ROSENBLUTH_SPECIES_VTHERMAL(*m_IVth.getData(),
      *m_IsKEC.getData(),
      *m_IN.getData(),
      BOX4D_TO_FORT(m_data_box));

   // If we're doing self collisions then compute the self collision forward
   // reaction diffusion tensor.
   if (m_selfcoll) {
      m_dparameters[NU] = m_nu[m_species_index];
      FORT_COMPUTE_ROSENBLUTH_DIFFUSION_TENSOR(*m_d_a.getData(),
         *m_reLam.getData(),
         1.0,
         *m_IVx.getData(),
         *m_IVy.getData(),
         *m_IVth.getData(),
         *m_IN.getData(),
         BOX4D_TO_FORT(m_data_box),
         BOX4D_TO_FORT(m_interior_box),
         PROBLEMDOMAIN_TO_FORT(m_domain),
         *m_velocities.getData(),
         0,
         0,
         m_range_lo[0],
         m_range_hi[0],
         m_dparameters[0],
         m_iparameters[0]);
   }

   // If the back reaction is computed then there's a bunch of other things that
   // need to be done.
   if (m_back_reaction) {
      // The self collision back reaction moments must be computed at each RK
      // stage.
      if (m_selfcoll) {
         // Compute the self collision back reaction diffusion tensor.
         FORT_COMPUTE_ROSENBLUTH_DIFFUSION_TENSOR(*m_d_aa.getData(),
            *m_reLam.getData(),
            1.0,
            *m_IVx.getData(),
            *m_IVy.getData(),
            *m_IVth.getData(),
            *m_IN.getData(),
            BOX4D_TO_FORT(m_data_box),
            BOX4D_TO_FORT(m_interior_box),
            PROBLEMDOMAIN_TO_FORT(m_domain),
            *m_velocities.getData(),
            1,
            0,
            m_range_lo[0],
            m_range_hi[0],
            m_dparameters[0],
            m_iparameters[0]);

         computeBackReactionMoments(m_d_aa);
      }
   }

   // If there are self collisions, then compute the self collision forward
   // reaction, m_coll.
   if (m_selfcoll) {
      m_coll = 0.0;
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("collision kernel");
      FORT_APPEND_ROSENBLUTH_COLLISION(*m_coll.getData(),
         *current_distribution.getData(),
         *m_velocities.getData(),
         *m_IVx.getData(),
         *m_IVy.getData(),
         *m_d_a.getData(),
         *m_IVth.getData(),
         *m_IN.getData(),
         1.0,
         BOX4D_TO_FORT(m_data_box),
         BOX4D_TO_FORT(m_interior_box),
         m_domain.dx()[0],
         m_dparameters[0],
         m_iparameters[0]);
      timers->stopTimer("collision kernel");

      // Add the forward reaction to the RHS.
      FORT_ADD_TO_DOUBLEARRAY_4D(*rhs_distribution.getData(),
         *m_coll.getData(),
         1.0,
         1.0,
         BOX4D_TO_FORT(m_data_box));

      if (m_back_reaction) {
         // Use m_coll to compute the self collision back reaction terms
         // numerators.
         m_rMomx = 0.0;
         m_rMomy = 0.0;
         m_rKE = 0.0;
         FORT_COMPUTE_ROSENBLUTH_NUMERATORS(*m_rMomx.getData(),
            *m_rMomy.getData(),
            *m_rKE.getData(),
            *m_IVx.getData(),
            *m_IVy.getData(),
            *m_coll.getData(),
            BOX4D_TO_FORT(m_data_box),
            BOX4D_TO_FORT(m_interior_box),
            *m_velocities.getData());

         m_IMomxN[m_species_index] = 0.0;
         m_IMomyN[m_species_index] = 0.0;
         m_IKEN[m_species_index] = 0.0;
         m_x_mom_reduction->execute(m_IMomxN[m_species_index]);
         m_y_mom_reduction->execute(m_IMomyN[m_species_index]);
         m_ke_reduction->execute(m_IKEN[m_species_index]);

         // If there are self collisions then add the self collision back
         // reaction to the RHS.
         FORT_APPEND_ROSENBLUTH_BACK_REACTION(*rhs_distribution.getData(),
            *m_cMomx.getData(),
            *m_cMomy.getData(),
            *m_cKE.getData(),
            *m_IMomxN[m_species_index].getData(),
            *m_IMomyN[m_species_index].getData(),
            *m_IKEN[m_species_index].getData(),
            *m_IMomxD.getData(),
            *m_IMomyD.getData(),
            *m_IKED.getData(),
            BOX4D_TO_FORT(m_data_box),
            BOX4D_TO_FORT(m_interior_box),
            1.0);
      }
   }

   if (!m_interspecies) {
      // If diagnostics are requested and there's no interspecies collisions
      // then the operator is done so compute the diagnostics now.
      if (m_diagnostic) {
         computeDiagnosticFields(rhs_distribution, current_distribution, a_dt);
      }
   }
   else {
      // Since interspecies collisions are on there is more work to be done.

      // Put the information about this species needed for the interspecies
      // forward collisions into the communication buffer.
      getInterspeciesData();

      // Now communicate that information to the head processors of each other
      // species.
      m_rosenbluthComm->CommunicateSpatialDataToSpeciesHead(false);
      m_rosenbluthComm->CommunicateSpatialDataToOtherSpecies(false);

      // For each interspecies collision, communicate the other species' data
      // from the head processors of this species to the other processors of
      // this species.  Then get that data from the communication buffer and
      // perform the forward interspecies collision with the other species.
      for (int i = 0; i < m_num_species; ++i) {
         if (i != m_species_index && m_nu[i] != 0.0) {
            m_rosenbluthComm->CommunicateSpatialDataFromSpeciesHeadFromOtherSpecies(i, false);

            double mass_ratio =
               m_species_masses[i] / m_species_masses[m_species_index];
            setOtherSpeciesData(i, false, mass_ratio);

            bool last_interspecies_collision = false;
            if (m_species_index == m_num_species-1) {
               if (i == m_num_species-2) {
                  last_interspecies_collision = true;
               }
            }
            else if (i == m_num_species-1) {
               last_interspecies_collision = true;
            }
            evaluateInterspecies(rhs_distribution,
               current_distribution,
               i,
               last_interspecies_collision,
               mass_ratio,
               a_dt,
               a_last_rk_stage);

            // Ensure that all the head processor sends have completed.
            m_rosenbluthComm->completeSends();
         }
      }

      // If back reactions are on we must communicate the back reaction
      // numerators and perform the back reaction in much the same way that the
      // forward interspecies collisions were handled above.
      if (m_back_reaction) {
         // Put the information about this species back reaction numerators
         // needed for the interspecies back reactions into the communication
         // buffer.
         getInterspeciesBRData();

         // Now communicate that information to the head processors of each
         // other species.
         m_rosenbluthComm->CommunicateSpatialDataToSpeciesHead(true);
         m_rosenbluthComm->CommunicateSpatialDataToOtherSpecies(true);

         // For each interspecies collision, communicate the other species' data
         // from the head processors of this species to the other processors of
         // this species.  Then get that data from the communication buffer and
         // perform the back reaction with the other species.
         for (int i = 0; i < m_num_species; ++i) {
            if (i != m_species_index && m_nu[i] != 0.0) {
               m_rosenbluthComm->CommunicateSpatialDataFromSpeciesHeadFromOtherSpecies(i, true);

               double mass_ratio =
                  m_species_masses[i] / m_species_masses[m_species_index];
               setOtherSpeciesBRData(i, mass_ratio);

               bool last_interspecies_collision = false;
               if (m_species_index == m_num_species-1) {
                  if (i == m_num_species-2) {
                     last_interspecies_collision = true;
                  }
               }
               else if (i == m_num_species-1) {
                  last_interspecies_collision = true;
               }
               evaluateInterspeciesBackReaction(rhs_distribution,
                  current_distribution,
                  last_interspecies_collision,
                  mass_ratio,
                  a_dt,
                  a_last_rk_stage);

               // Ensure that all the head processor sends have completed.
               m_rosenbluthComm->completeSends();
            }
         }
      }
   }
}


void
RosenbluthCollisionOperator::evaluateReLam(
   const KineticSpecies& a_u)
{
   // This routine does the parts of evaluate that are necessary to compute the
   // Re(lambda).  Essentially all the forward reaction diffusion tensors and
   // everything that they depend on must be computed 

   const ParallelArray& current_distribution = a_u.distribution();

   m_reLam = 0.0;
   m_reLamMax = 0.0;

   // If neither of these is set the operator is a no-op so just return.
   if (!m_selfcoll && !m_interspecies) {
      return;
   }

   // We need several quantities computed whether we're doing self or
   // interspecies collisions.
   m_rN = 0.0;
   m_rGammax = 0.0;
   m_rGammay = 0.0;
   m_rsKEC = 0.0;

   // m_rsKEC is passed as a dummy and will be computed below.
   FORT_COMPUTE_ROSENBLUTH_SPECIES_MOMENTS(*m_rN.getData(),
      *m_rGammax.getData(),
      *m_rGammay.getData(),
      *m_rsKEC.getData(),
      *current_distribution.getData(),
      0,
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      *m_velocities.getData());

   m_IN = 0.0;
   m_IGammax = 0.0;
   m_IGammay = 0.0;

   m_N_reduction->execute(m_IN);
   m_Gammax_reduction->execute(m_IGammax);
   m_Gammay_reduction->execute(m_IGammay);
   FORT_COMPUTE_ROSENBLUTH_SPECIES_REDUCED_FIELDS(*m_IVx.getData(),
      *m_IVy.getData(),
      *m_IN.getData(),
      *m_IGammax.getData(),
      *m_IGammay.getData(),
      BOX4D_TO_FORT(m_data_box));

   m_rsKEC = 0.0;
   FORT_COMPUTE_ROSENBLUTH_SPECIES_KEC(*m_rsKEC.getData(),
      *m_IVx.getData(),
      *m_IVy.getData(),
      *current_distribution.getData(),
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      *m_velocities.getData());

   m_IsKEC = 0.0;
   m_skec_reduction->execute(m_IsKEC);

   FORT_COMPUTE_ROSENBLUTH_SPECIES_VTHERMAL(*m_IVth.getData(),
      *m_IsKEC.getData(),
      *m_IN.getData(),
      BOX4D_TO_FORT(m_data_box));

   // If we're doing self collisions then compute the self collision forward
   // reaction diffusion tensor.
   if (m_selfcoll) {
      m_dparameters[NU] = m_nu[m_species_index];
      FORT_COMPUTE_ROSENBLUTH_DIFFUSION_TENSOR(*m_d_a.getData(),
         *m_reLam.getData(),
         1.0,
         *m_IVx.getData(),
         *m_IVy.getData(),
         *m_IVth.getData(),
         *m_IN.getData(),
         BOX4D_TO_FORT(m_data_box),
         BOX4D_TO_FORT(m_interior_box),
         PROBLEMDOMAIN_TO_FORT(m_domain),
         *m_velocities.getData(),
         0,
         1,
         m_range_lo[0],
         m_range_hi[0],
         m_dparameters[0],
         m_iparameters[0]);

      // If only self collisions are being done then just compute the max of
      // reLam for the time step control.
      if (!m_interspecies) {
         FORT_GET_MAX_RELAM(m_reLamMax,
            *m_reLam.getData(),
            BOX4D_TO_FORT(m_interior_box),
            m_iparameters[0]);
         return;
      }
   }

   // Interspecies collisions are on so there is more work to be done.
   if (m_interspecies) {
      // Put the information about this species needed for the interspecies
      // forward collisions into the communication buffer.
      getInterspeciesData();

      // Now communicate that information to the head processors of each other
      // species.
      m_rosenbluthComm->CommunicateSpatialDataToSpeciesHead(false);
      m_rosenbluthComm->CommunicateSpatialDataToOtherSpecies(false);

      // For each interspecies collision, communicate the other species' data
      // from the head processors of this species to the other processors of
      // this species.  Then get that data from the communication buffer and
      // perform  the forward interspecies collision with the other species.
      for (int i = 0; i < m_num_species; ++i) {
         if (i != m_species_index && m_nu[i] != 0.0) {
            m_rosenbluthComm->CommunicateSpatialDataFromSpeciesHeadFromOtherSpecies(i, false);

            double mass_ratio =
               m_species_masses[i] / m_species_masses[m_species_index];
            setOtherSpeciesData(i, true, mass_ratio);

            bool last_interspecies_collision = false;
            if (m_species_index == m_num_species-1) {
               if (i == m_num_species-2) {
                  last_interspecies_collision = true;
               }
            }
            else if (i == m_num_species-1) {
               last_interspecies_collision = true;
            }
            if (last_interspecies_collision) {
               FORT_GET_MAX_RELAM(m_reLamMax,
                  *m_reLam.getData(),
                  BOX4D_TO_FORT(m_interior_box),
                  m_iparameters[0]);
            }

            // Ensure that all the head processor sends have completed.
            m_rosenbluthComm->completeSends();
         }
      }
   }
}


double
RosenbluthCollisionOperator::computeRealLam(
   KineticSpecies& a_u)
{
   evaluateReLam(a_u);
   double pi = 4.0*atan(1.0);
   return m_reLamMax*pi*pi*(1.0/pow(a_u.domain().dx(V1), 2.0)+1.0/pow(a_u.domain().dx(V2), 2.0));
}


void
RosenbluthCollisionOperator::printParameters() const
{
   // Print all operator parameters in play.
   Loki_Utilities::printF("  Using Rosenbluth collision:\n" );
   Loki_Utilities::printF("    collision velocity range = %e %e %e %e\n",
      m_range_lo[0], m_range_hi[0], m_range_lo[1], m_range_hi[1]);
   if (m_back_reaction) {
      Loki_Utilities::printF("    computing back reaction\n");
   }
   else {
      Loki_Utilities::printF("    not computing back reaction\n");
   }
   if (m_diagnostic) {
      Loki_Utilities::printF("    performing diagnostic computations and output\n");
   }
   if (m_interspecies) {
      Loki_Utilities::printF("    computing interspecies collisions\n");
      Loki_Utilities::printF("    nuCoeffInterspecies      = ");
      for (int i = 0; i < m_num_species; ++i) {
         if (i != m_species_index) {
            Loki_Utilities::printF("%e ", m_nu[i]);
         }
      }
      Loki_Utilities::printF("\n");
   }
   else {
      Loki_Utilities::printF("    not computing back reaction\n");
   }
   if (m_selfcoll) {
      Loki_Utilities::printF("    computing self collisions\n");
      Loki_Utilities::printF("    nuCoeffSelf              = %e\n",
         m_nu[m_species_index]);
   }
   else {
      Loki_Utilities::printF("    not computing self collisions\n");
   }
   if (m_iparameters[KERNEL_TYPE] == 0) {
      Loki_Utilities::printF("    using primitive kernel\n");
   } 
   else if (m_iparameters[KERNEL_TYPE] == 1) {
      Loki_Utilities::printF("    using log kernel\n");
   }
   else if (m_iparameters[KERNEL_TYPE] == 2) {
      Loki_Utilities::printF("    using original kernel\n");
   }
   else {
      Loki_Utilities::printF("    using asinh kernel\n");
   }
   if (m_iparameters[DENSITY_DEP_NU] == 0) {
      Loki_Utilities::printF("    no collision coefficient density modulation\n");
   }
   else {
      Loki_Utilities::printF("    with collision coefficient density modulation\n");
   }
}


void
RosenbluthCollisionOperator::copyDiagnosticFields(
   vector<ParallelArray>& a_d) const
{
   // If diagnostics were computed then copy them.  Otherwise we know nothing so
   // set the output to 0.
   if (m_diagnostic) {
      a_d[0] = m_dMomx;
      a_d[1] = m_dMomy;
      a_d[2] = m_dKE;
      a_d[3] = m_dEnt;
   }
   else {
      a_d[0] = 0.0;
      a_d[1] = 0.0;
      a_d[2] = 0.0;
      a_d[3] = 0.0;
   }
}


bool
RosenbluthCollisionOperator::isType(
   const string& a_name)
{
   if (a_name.compare(s_CLASS_NAME) == 0) {
      return true;
   }
   return false;
}


void
RosenbluthCollisionOperator::evaluateInterspecies(
   ParallelArray& a_rhs,
   const ParallelArray& a_u,
   int a_other_species,
   bool a_last_interspecies_collision,
   double a_mass_ratio,
   double a_dt,
   bool a_last_rk_stage)
{
   if (!m_interspecies) {
      return;
   }

   // If there are collisions with the other species then compute the
   // interspecies forward reaction, m_coll.
   m_coll = 0.0;
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("collision kernel");
   FORT_APPEND_ROSENBLUTH_COLLISION(*m_coll.getData(),
      *a_u.getData(),
      *m_velocities.getData(),
      *m_IVxb.getData(),
      *m_IVyb.getData(),
      *m_d_ba.getData(),
      *m_IVthb.getData(),
      *m_INb.getData(),
      a_mass_ratio,
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      m_domain.dx()[0],
      m_dparameters[0],
      m_iparameters[0]);
   timers->stopTimer("collision kernel");

   // Add the forward reaction to the RHS.
   FORT_ADD_TO_DOUBLEARRAY_4D(*a_rhs.getData(),
      *m_coll.getData(),
      1.0,
      1.0,
      BOX4D_TO_FORT(m_data_box));

   // If the interspecies back reaction is to be computed, use m_coll to compute
   // the interspecies back reaction numerators.
   if (m_back_reaction) {
      m_rMomx = 0.0;
      m_rMomy = 0.0;
      m_rKE = 0.0;
      FORT_COMPUTE_ROSENBLUTH_NUMERATORS(*m_rMomx.getData(),
         *m_rMomy.getData(),
         *m_rKE.getData(),
         *m_IVxb.getData(),
         *m_IVyb.getData(),
         *m_coll.getData(),
         BOX4D_TO_FORT(m_data_box),
         BOX4D_TO_FORT(m_interior_box),
         *m_velocities.getData());

      m_IMomxN[a_other_species] = 0.0;
      m_IMomyN[a_other_species] = 0.0;
      m_IKEN[a_other_species] = 0.0;
      m_x_mom_reduction->execute(m_IMomxN[a_other_species]);
      m_y_mom_reduction->execute(m_IMomyN[a_other_species]);
      m_ke_reduction->execute(m_IKEN[a_other_species]);
   }

   // If diagnostics are reqested, there's no back reaction, and this is the
   // last interspecies collision, then the operator is done so compute the
   // diagnostics now.
   if (a_last_rk_stage && m_diagnostic &&
       !m_back_reaction && a_last_interspecies_collision) {
      computeDiagnosticFields(a_rhs, a_u, a_dt);
   }
}


void
RosenbluthCollisionOperator::evaluateInterspeciesBackReaction(
   ParallelArray& a_rhs,
   const ParallelArray& a_u,
   bool a_last_interspecies_collision,
   double a_mass_ratio,
   double a_dt,
   bool a_last_rk_stage)
{
   if (!m_interspecies || !m_back_reaction) {
      return;
   }

   // If there collisions with the other species then add the interspecies back
   // reaction to the RHS.
   FORT_APPEND_ROSENBLUTH_BACK_REACTION(*a_rhs.getData(),
      *m_cMomx.getData(),
      *m_cMomy.getData(),
      *m_cKE.getData(),
      *m_IMomxNb.getData(),
      *m_IMomyNb.getData(),
      *m_IKENb.getData(),
      *m_IMomxD.getData(),
      *m_IMomyD.getData(),
      *m_IKED.getData(),
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      a_mass_ratio);

   // There's nothing else for the collision operator to do so if diagnostics
   // are reqested compute the diagnostics now.
   if (a_last_rk_stage && m_diagnostic && a_last_interspecies_collision) {
      computeDiagnosticFields(a_rhs, a_u, a_dt);
   }
}


void
RosenbluthCollisionOperator::computeDiagnosticFields(
   const ParallelArray& a_rhs,
   const ParallelArray& a_u,
   double a_dt)
{
   // Subtract the initial RHS, m_d4D, from the current RHS which give the
   // contribution of the collision operator to the RHS.
   FORT_ADD_TO_DOUBLEARRAY_4D(*m_d4D.getData(),
      *a_rhs.getData(),
      -1.0,
      1.0,
      BOX4D_TO_FORT(m_data_box));

   // Use FORT_COMPUTE_ROSENBLUTH_SPECIES_MOMENTS to compute change in momentum
   // and energy.  m_rN is passed as a dummy variable and will be computed
   // below.
   m_dEnt = 0.0;
   m_dMomx = 0.0;
   m_dMomy = 0.0;
   m_dKE = 0.0;
   FORT_COMPUTE_ROSENBLUTH_SPECIES_MOMENTS(*m_dEnt.getData(),
      *m_dMomx.getData(),
      *m_dMomy.getData(),
      *m_dKE.getData(),
      *m_d4D.getData(),
      1,
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      *m_velocities.getData());

   // Compute change in entropy: -dt * df/dt * (ln(|f_i| + floor) + 1)
   m_dEnt = 0.0;
   FORT_COMPUTE_ROSENBLUTH_ENTROPY_CHANGE(*m_dEnt.getData(),
      *a_u.getData(),
      *m_d4D.getData(),
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box));

   // Do not reduce the diagnostic fields.  When they are used they will be
   // summed over all processors which is a reduction without the multiplication
   // by dvx*dvy so to reduce here would result in double (or more) counting.
   // All we need to do here is to scale by dv.
   double dv = m_domain.dx(V1)*m_domain.dx(V2);
   m_dMomx *= dv*a_dt*m_species_masses[m_species_index];
   m_dMomy *= dv*a_dt*m_species_masses[m_species_index];
   m_dKE *= dv*a_dt*0.5*m_species_masses[m_species_index];
   m_dEnt *= dv*a_dt;
}


void
RosenbluthCollisionOperator::getInterspeciesData() const
{
   FORT_GET_COMMUNICATED_SPECIES_MOMENTS(*m_rosenbluthComm->MyInteriorData().getData(),
      *m_IN.getData(),
      *m_IVx.getData(),
      *m_IVy.getData(),
      *m_IVth.getData(),
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      m_maxNumCommunicatedFields);
}


void
RosenbluthCollisionOperator::getInterspeciesBRData() const
{
   int idx = 0;
   for (int i = 0; i < m_num_species; ++i) {
      if (i != m_species_index && m_nu[i] != 0.0) {
         FORT_GET_COMMUNICATED_BR_SPECIES_MOMENTS(*m_rosenbluthComm->MyInteriorData().getData(),
            *m_IMomxN[i].getData(),
            *m_IMomyN[i].getData(),
            *m_IKEN[i].getData(),
            BOX4D_TO_FORT(m_data_box),
            BOX4D_TO_FORT(m_interior_box),
            m_maxNumCommunicatedFields,
            idx);
         ++idx;
      }
   }
}


void
RosenbluthCollisionOperator::setOtherSpeciesData(
   int a_which_interspecies_collision,
   bool a_computeReLamOnly,
   double mboverma)
{
   m_dparameters[NU] = m_nu[a_which_interspecies_collision];
   FORT_SET_COMMUNICATED_SPECIES_MOMENTS(*m_INb.getData(),
      *m_IVxb.getData(),
      *m_IVyb.getData(),
      *m_IVthb.getData(),
      *m_rosenbluthComm->MyInteriorData().getData(),
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      m_maxNumCommunicatedFields);

   // Now that we have all the necessary information, compute the interspecies
   // collision diffusion tensor for the forward collision, m_d_ba.
   FORT_COMPUTE_ROSENBLUTH_DIFFUSION_TENSOR(*m_d_ba.getData(),
      *m_reLam.getData(),
      1.0 / mboverma,
      *m_IVxb.getData(),
      *m_IVyb.getData(),
      *m_IVthb.getData(),
      *m_INb.getData(),
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      PROBLEMDOMAIN_TO_FORT(m_domain),
      *m_velocities.getData(),
      0,
      a_computeReLamOnly ? 1 : 0,
      m_range_lo[0],
      m_range_hi[0],
      m_dparameters[0],
      m_iparameters[0]);
}


void
RosenbluthCollisionOperator::setOtherSpeciesBRData(
   int a_which_interspecies_collision,
   double mboverma)
{
   m_dparameters[NU] = m_nu[a_which_interspecies_collision];
   FORT_SET_COMMUNICATED_SPECIES_BR_MOMENTS(*m_IMomxNb.getData(),
      *m_IMomyNb.getData(),
      *m_IKENb.getData(),
      *m_rosenbluthComm->MyInteriorData().getData(),
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      m_maxNumCommunicatedFields);

   // Overwrite m_cKE, which currently stores the result for self-collisions and
   // was already used for the self back reaction, with
   // cKE_ba = RC(v_a^2*F_Ma, nu, alpha_ba, v_xa, v_ya, d_ab), to be
   // used by evaluateInterspecies for the interspecies back reaction. First,
   // compute m_d_ab, which is the same as m_d_a except for the mass scaling.
   // Same for m_cMomx, m_cMomy, m_IMomxD, m_IMomyD, and m_IKED.

   // Compute the interspecies collision diffusion tensor for the back reaction,
   // m_d_ab, and the interspecies collision back reaction moments.
   FORT_COMPUTE_ROSENBLUTH_DIFFUSION_TENSOR(*m_d_ab.getData(),
      *m_reLam.getData(),
      mboverma,
      *m_IVx.getData(),
      *m_IVy.getData(),
      *m_IVth.getData(),
      *m_IN.getData(),
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      PROBLEMDOMAIN_TO_FORT(m_domain),
      *m_velocities.getData(),
      1,
      0,
      m_range_lo[0],
      m_range_hi[0],
      m_dparameters[0],
      m_iparameters[0]);

   computeBackReactionMoments(m_d_ab);
}


void
RosenbluthCollisionOperator::parseParameters(
   LokiInputParser& a_pp,
   const KineticSpecies* a_species)
{
   // The upper and lower ranges of the collision in velocity space must be
   // supplied.  If this is a relativistic problem these ranges need to be
   // transformed into a momentum box.
   if (!a_pp.queryarr("collision_vel_range_lo", m_range_lo, 0, 2)) {
      LOKI_ABORT("Must supply collision_vel_range_lo.");
   }
   else if (a_pp.countval("collision_vel_range_lo") != 2) {
      LOKI_ABORT("collision_vel_range_lo must have 2 entries.");
   }
   if (Simulation::s_DO_RELATIVITY) {
      for (int i = 0; i < 2; ++i) {
         double v = m_range_lo[i];
         m_range_lo[i] =
            a_species->mass()*v/sqrt(1.0-pow(v/Simulation::s_LIGHT_SPEED, 2));
      }
   }
   if (!a_pp.queryarr("collision_vel_range_hi", m_range_hi, 0, 2)) {
      LOKI_ABORT("Must supply collision_vel_range_hi.");
   }
   else if (a_pp.countval("collision_vel_range_hi") != 2) {
      LOKI_ABORT("collision_vel_range_hi must have 2 entries.");
   }
   if (Simulation::s_DO_RELATIVITY) {
      for (int i = 0; i < 2; ++i) {
         double v = m_range_hi[i];
         m_range_hi[i] =
            a_species->mass()*v/sqrt(1.0-pow(v/Simulation::s_LIGHT_SPEED, 2));
      }
   }

   // Check sanity of upper and lower ranges of collisionality.
   double rolloff = m_iparameters[SOLN_ORDER] == 4 ? 3 : 4;
   double vxmin = m_domain.lower(V1) + rolloff * m_domain.dx(V1);
   double vxmax = m_domain.upper(V1) - rolloff * m_domain.dx(V1);
   double vymin = m_domain.lower(V2) + rolloff * m_domain.dx(V2);
   double vymax = m_domain.upper(V2) - rolloff * m_domain.dx(V2);
   if (m_range_lo[0] <= vxmin || m_range_hi[0] >= vxmax) {
      LOKI_ABORT("x collision_vel_range box too large");
   }
   if (m_range_lo[1] <= vymin || m_range_hi[1] >= vymax) {
      LOKI_ABORT("y collision_vel_range box too large");
   }
   if (m_range_lo[0] >= m_range_hi[0]) {
      LOKI_ABORT("x collision_vel_range_lo exceeds x collision_vel_range_hi");
   }
   if (m_range_lo[1] >= m_range_hi[1]) {
      LOKI_ABORT("y collision_vel_range_lo exceeds y collision_vel_range_hi");
   }

   // Read if the back reaction should be computed.
   string tmp("false");
   a_pp.query("collision_back_reaction", tmp);
   m_back_reaction = tmp.compare("false") == 0 ? false : true;

   // Read if interspecies collisions should be computed.
   string tmp2("false");
   a_pp.query("collision_interspecies", tmp2);
   m_interspecies = tmp2.compare("false") == 0 ? false : true;

   // If interspecies collisions are on the interspecies collision
   // coefficient(s) must be supplied.
   if (m_interspecies) {
      vector<double> nu_interspecies(m_num_species-1);
      if (!a_pp.queryarr("collision_nuCoeffInterspecies",
                         nu_interspecies, 0, m_num_species-1)) {
         LOKI_ABORT("Must supply collision_nuCoeffInterspecies when doing interspecies collisions.");
      }
      else if (a_pp.countval("collision_nuCoeffInterspecies") !=
               m_num_species-1) {
         LOKI_ABORT("collision_nuCoeffInterspecies must have num_species-1 entries.");
      }
      // If all the interspecies collision coefficients are zero then this
      // species is not really doing interspecies collisions so turn it off for
      // efficiency.
      for (int i = 0; i < m_num_species-1; ++i) {
         if (nu_interspecies[i] != 0.0) {
            ++m_num_interspecies_collisions;
         }
         if (i < m_species_index) {
            m_nu[i] = nu_interspecies[i];
         }
         else {
            m_nu[i+1] = nu_interspecies[i];
         }
      }
      if (m_num_interspecies_collisions == 0) {
         m_interspecies = false;
      }
   }

   // Read if diagnostic data should be computed and output.
   string tmp3("false");
   a_pp.query("collision_diagnostic", tmp3);
   m_diagnostic = tmp3.compare("false") == 0 ? false : true;

   // Read if self collisions should be computed.
   string tmp4("false");
   a_pp.query("collision_self", tmp4);
   m_selfcoll = tmp4.compare("false") == 0 ? false : true;

   // If self collisions are on the self collision coefficient must be supplied.
   // If the self collision coefficient is zero then self collisions are not
   // really being performed so turn them off.
   if (m_selfcoll) {
      double nu_self;
      if (!a_pp.query("collision_nuCoeffSelf", nu_self)) {
         LOKI_ABORT("Must supply collision_nuCoeffSelf when doing self collisions.");
      }
      else if (nu_self == 0.0) {
         m_selfcoll = false;
      }
      m_nu[m_species_index] = nu_self;
   }

   // Read the kernel algorithm.
   string tmp5("primitive");
   a_pp.query("kernel_alg", tmp5);
   if (tmp5.compare("primitive") == 0) {
      m_iparameters[KERNEL_TYPE] = 0;
   }
   else if (tmp5.compare("log") == 0) {
      m_iparameters[KERNEL_TYPE] = 1;
   }
   else if (tmp5.compare("original") == 0) {
      m_iparameters[KERNEL_TYPE] = 2;
   }
   else if (tmp5.compare("asinh") == 0) {
      m_iparameters[KERNEL_TYPE] = 3;
   }
   else {
      LOKI_ABORT("Unknown kernel type.");
   }

   // Read if the collision coefficient should be modulated by the local
   // density.
   string tmp6("false");
   a_pp.query("modulate_collision_rates", tmp6);
   if (tmp6.compare("false") == 0) {
      m_iparameters[DENSITY_DEP_NU] = 0;
   }
   else {
      m_iparameters[DENSITY_DEP_NU] = 1;
   }
}


void
RosenbluthCollisionOperator::constructDiagnosticArrays()
{
   int nGhosts = m_iparameters[NUM_GHOSTS];
   vector<int> num_global_cells4D(4);
   for (int i = 0; i < 4; ++i) {
      num_global_cells4D[i] = m_domain.numberOfCells(i);
   }

   // Dimension 4D helper array.
   m_d4D.partition(m_interior_box, 4, nGhosts, num_global_cells4D);

   // Configuration space box.
   ParallelArray::Box config_space(2);
   vector<int> num_global_cells2D(2);
   for (int i = 0; i < 2; ++i) {
      config_space.lower(i) = m_interior_box.lower(i);
      config_space.upper(i) = m_interior_box.upper(i);
      num_global_cells2D[i] = m_domain.numberOfCells(i);
   }

   // Dimension persistent 2D globally reduced momenta and kinetic energy.
   m_dMomx.partition(config_space, 2, nGhosts, num_global_cells2D);
   m_dMomy.partition(config_space, 2, nGhosts, num_global_cells2D);
   m_dKE.partition(config_space, 2, nGhosts, num_global_cells2D);
   m_dEnt.partition(config_space, 2, nGhosts, num_global_cells2D);
}


void
RosenbluthCollisionOperator::constructArrays()
{
   int nGhosts = m_iparameters[NUM_GHOSTS];

   // Configuration space box.
   ParallelArray::Box config_space(2);
   vector<int> num_global_cells2D(2);
   for (int i = 0; i < 2; ++i) {
      config_space.lower(i) = m_interior_box.lower(i);
      config_space.upper(i) = m_interior_box.upper(i);
      num_global_cells2D[i] = m_domain.numberOfCells(i);
   }

   // Dimension persistent arrays.
   vector<int> num_global_cells4D(4);
   for (int i = 0; i < 4; ++i) {
      num_global_cells4D[i] = m_domain.numberOfCells(i);
   }
   if (m_selfcoll) {
      ParallelArray::Box diff_ten_base_space(5);
      for (int i = 0; i < 2; ++i) {
         diff_ten_base_space.lower(i) = m_interior_box.lower(i);
         diff_ten_base_space.upper(i) = m_interior_box.upper(i);
      }
      for (int i = 2; i < 4; ++i) {
         diff_ten_base_space.lower(i) = m_interior_box.lower(i)-nGhosts;
         diff_ten_base_space.upper(i) = m_interior_box.upper(i)+nGhosts;
      }
      diff_ten_base_space.lower(4) = 0;
      diff_ten_base_space.upper(4) = 2;
      m_d_a.partition(diff_ten_base_space, 4, 0, num_global_cells4D);
   }
   ParallelArray::Box lambda_base_space(4);
   for (int i = 0; i < 2; ++i) {
      lambda_base_space.lower(i) = m_interior_box.lower(i);
      lambda_base_space.upper(i) = m_interior_box.upper(i);
   }
   for (int i = 2; i < 4; ++i) {
      lambda_base_space.lower(i) = m_interior_box.lower(i)-nGhosts;
      lambda_base_space.upper(i) = m_interior_box.upper(i)+nGhosts;
   }
   m_reLam.partition(lambda_base_space, 4, 0, num_global_cells4D);
   m_coll.partition(m_interior_box, 4, nGhosts, num_global_cells4D);
   m_IVx.partition(config_space, 2, nGhosts, num_global_cells2D);
   m_IVy.partition(config_space, 2, nGhosts, num_global_cells2D);
   m_IVth.partition(config_space, 2, nGhosts, num_global_cells2D);
   m_IN.partition(config_space, 2, nGhosts, num_global_cells2D);
   m_rN.partition(config_space, 2, nGhosts, num_global_cells2D);
   m_rGammax.partition(config_space, 2, nGhosts, num_global_cells2D);
   m_rGammay.partition(config_space, 2, nGhosts, num_global_cells2D);
   m_rsKEC.partition(config_space, 2, nGhosts, num_global_cells2D);
   m_IGammax.partition(config_space, 2, nGhosts, num_global_cells2D);
   m_IGammay.partition(config_space, 2, nGhosts, num_global_cells2D);
   m_IsKEC.partition(config_space, 2, nGhosts, num_global_cells2D);
   if (m_diagnostic || m_back_reaction) {
      m_rMomx.partition(config_space, 2, nGhosts, num_global_cells2D);
      m_rMomy.partition(config_space, 2, nGhosts, num_global_cells2D);
      m_rKE.partition(config_space, 2, nGhosts, num_global_cells2D);
   }
   if (m_interspecies) {
      ParallelArray::Box diff_ten_base_space(5);
      for (int i = 0; i < 2; ++i) {
         diff_ten_base_space.lower(i) = m_interior_box.lower(i);
         diff_ten_base_space.upper(i) = m_interior_box.upper(i);
      }
      for (int i = 2; i < 4; ++i) {
         diff_ten_base_space.lower(i) = m_interior_box.lower(i)-nGhosts;
         diff_ten_base_space.upper(i) = m_interior_box.upper(i)+nGhosts;
      }
      diff_ten_base_space.lower(4) = 0;
      diff_ten_base_space.upper(4) = 2;
      m_d_ba.partition(diff_ten_base_space, 4, 0, num_global_cells4D);
      m_INb.partition(config_space, 2, nGhosts, num_global_cells2D);
      m_IVthb.partition(config_space, 2, nGhosts, num_global_cells2D);
      m_IVxb.partition(config_space, 2, nGhosts, num_global_cells2D);
      m_IVyb.partition(config_space, 2, nGhosts, num_global_cells2D);
   }
}


void
RosenbluthCollisionOperator::constructReductionSchedules()
{
   // Build reduction schedules for momenta and kinetic energy denominators.
   deque<bool> collapse_dir(m_domain.dim(), false);
   collapse_dir[V1] = true;
   collapse_dir[V2] = true;
   double measure(m_domain.dx(2) * m_domain.dx(3));

   if (m_back_reaction) {
      // x-momentum reduction schedule
      m_x_mom_reduction = new ReductionSchedule4D(m_rMomx,
         m_interior_box,
         m_domain.box(),
         collapse_dir,
         measure,
         m_proc_lo,
         m_proc_hi,
         m_comm);

      // y-momentum reduction schedule
      m_y_mom_reduction = new ReductionSchedule4D(m_rMomy,
         m_interior_box,
         m_domain.box(),
         collapse_dir,
         measure,
         m_proc_lo,
         m_proc_hi,
         m_comm);

      // KE reduction schedule
      m_ke_reduction = new ReductionSchedule4D(m_rKE,
         m_interior_box,
         m_domain.box(),
         collapse_dir,
         measure,
         m_proc_lo,
         m_proc_hi,
         m_comm);
   }

   // Density reduction schedule
   m_N_reduction = new ReductionSchedule4D(m_rN,
      m_interior_box,
      m_domain.box(),
      collapse_dir,
      measure,
      m_proc_lo,
      m_proc_hi,
      m_comm);
   // Momentum reduction schedule
   m_Gammax_reduction = new ReductionSchedule4D(m_rGammax,
      m_interior_box,
      m_domain.box(),
      collapse_dir,
      measure,
      m_proc_lo,
      m_proc_hi,
      m_comm);
   // Momentum reduction schedule
   m_Gammay_reduction = new ReductionSchedule4D(m_rGammay,
      m_interior_box,
      m_domain.box(),
      collapse_dir,
      measure,
      m_proc_lo,
      m_proc_hi,
      m_comm);
   // Kinetic energy reduction schedule
   m_skec_reduction = new ReductionSchedule4D(m_rsKEC,
      m_interior_box,
      m_domain.box(),
      collapse_dir,
      measure,
      m_proc_lo,
      m_proc_hi,
      m_comm);
}


void
RosenbluthCollisionOperator::constructBackReactionArrays()
{
   int nGhosts = m_iparameters[NUM_GHOSTS];
   vector<int> num_global_cells2D(2);
   vector<int> num_global_cells4D(4);
   for (int i = 0; i < 2; ++i) {
      num_global_cells2D[i] = m_domain.numberOfCells(i);
      num_global_cells4D[i] = m_domain.numberOfCells(i);
   }
   for (int i = 2; i < 4; ++i) {
      num_global_cells4D[i] = m_domain.numberOfCells(i);
   }

   // Dimension persistent 4D momenta and kinetic energy.
   m_uMomx.partition(m_interior_box, 4, nGhosts, num_global_cells4D);
   m_uMomy.partition(m_interior_box, 4, nGhosts, num_global_cells4D);
   m_uKE.partition(m_interior_box, 4, nGhosts, num_global_cells4D);
   m_cMomx.partition(m_interior_box, 4, nGhosts, num_global_cells4D);
   m_cMomy.partition(m_interior_box, 4, nGhosts, num_global_cells4D);
   m_cKE.partition(m_interior_box, 4, nGhosts, num_global_cells4D);

   // m_IMomxN, m_IMomyN, and m_IKEN are used for back reactions of enabled
   // forward reactions.
   ParallelArray::Box config_space(2);
   for (int i = 0; i < 2; ++i) {
      config_space.lower(i) = m_interior_box.lower(i);
      config_space.upper(i) = m_interior_box.upper(i);
   }
   for (int i = 0; i < m_num_species; ++i) {
      if ((m_selfcoll && i == m_species_index) ||
          (m_interspecies && i != m_species_index && m_nu[i] != 0.0)) {
         m_IMomxN[i].partition(config_space, 2, nGhosts, num_global_cells2D);
         m_IMomyN[i].partition(config_space, 2, nGhosts, num_global_cells2D);
         m_IKEN[i].partition(config_space, 2, nGhosts, num_global_cells2D);
      }
   }

   // The self collision back reaction diffusion tensor.
   if (m_selfcoll) {
      ParallelArray::Box diff_ten_base_space(5);
      for (int i = 0; i < 2; ++i) {
         diff_ten_base_space.lower(i) = m_interior_box.lower(i);
         diff_ten_base_space.upper(i) = m_interior_box.upper(i);
      }
      for (int i = 2; i < 4; ++i) {
         diff_ten_base_space.lower(i) = m_interior_box.lower(i)-nGhosts;
         diff_ten_base_space.upper(i) = m_interior_box.upper(i)+nGhosts;
      }
      diff_ten_base_space.lower(4) = 0;
      diff_ten_base_space.upper(4) = 2;
      m_d_aa.partition(diff_ten_base_space, 4, 0, num_global_cells4D);
   }

   // The interspecies back reaction diffusion tensor and numerators.
   if (m_interspecies) {
      ParallelArray::Box diff_ten_base_space(5);
      for (int i = 0; i < 2; ++i) {
         diff_ten_base_space.lower(i) = m_interior_box.lower(i);
         diff_ten_base_space.upper(i) = m_interior_box.upper(i);
      }
      for (int i = 2; i < 4; ++i) {
         diff_ten_base_space.lower(i) = m_interior_box.lower(i)-nGhosts;
         diff_ten_base_space.upper(i) = m_interior_box.upper(i)+nGhosts;
      }
      diff_ten_base_space.lower(4) = 0;
      diff_ten_base_space.upper(4) = 2;
      m_d_ab.partition(diff_ten_base_space, 4, 0, num_global_cells4D);
      m_IMomxNb.partition(config_space, 2, nGhosts, num_global_cells2D);
      m_IMomyNb.partition(config_space, 2, nGhosts, num_global_cells2D);
      m_IKENb.partition(config_space, 2, nGhosts, num_global_cells2D);
   }

   // Dimension persistent 2D globally reduced momenta and kinetic energy
   // denominator terms.
   m_IMomxD.partition(config_space, 2, nGhosts, num_global_cells2D);
   m_IMomyD.partition(config_space, 2, nGhosts, num_global_cells2D);
   m_IKED.partition(config_space, 2, nGhosts, num_global_cells2D);
}


void
RosenbluthCollisionOperator::computeBackReactionMoments(
   const ParallelArray& a_d)
{
   // Compute temporary 4D momenta and kinetic energy.
   FORT_COMPUTE_ROSENBLUTH_TEMPS(*m_uMomx.getData(),
      *m_uMomy.getData(),
      *m_uKE.getData(),
      *m_IVx.getData(),
      *m_IVy.getData(),
      *m_IVth.getData(),
      *m_IN.getData(),
      BOX4D_TO_FORT(m_data_box),
      *m_velocities.getData());

   // Compute persistent 4D momenta and kinetic energy.
   m_cMomx = 0.0;
   m_cMomy = 0.0;
   m_cKE = 0.0;
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("collision kernel");
   FORT_APPEND_ROSENBLUTH_COLLISION(*m_cMomx.getData(),
      *m_uMomx.getData(),
      *m_velocities.getData(),
      *m_IVx.getData(),
      *m_IVy.getData(),
      *a_d.getData(),
      *m_IVth.getData(),
      *m_IN.getData(),
      1.0,
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      m_domain.dx()[0],
      m_dparameters[0],
      m_iparameters[0]);
   FORT_APPEND_ROSENBLUTH_COLLISION(*m_cMomy.getData(),
      *m_uMomy.getData(),
      *m_velocities.getData(),
      *m_IVx.getData(),
      *m_IVy.getData(),
      *a_d.getData(),
      *m_IVth.getData(),
      *m_IN.getData(),
      1.0,
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      m_domain.dx()[0],
      m_dparameters[0],
      m_iparameters[0]);
   FORT_APPEND_ROSENBLUTH_COLLISION(*m_cKE.getData(),
      *m_uKE.getData(),
      *m_velocities.getData(),
      *m_IVx.getData(),
      *m_IVy.getData(),
      *a_d.getData(),
      *m_IVth.getData(),
      *m_IN.getData(),
      1.0,
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      m_domain.dx()[0],
      m_dparameters[0],
      m_iparameters[0]);
   timers->stopTimer("collision kernel");

   // Fill m_rMomx, m_rMomy, and m_rKE with quantities to be reduced to form
   // the back reaction denominator terms.
   m_rMomx = 0.0;
   m_rMomy = 0.0;
   m_rKE = 0.0;
   FORT_COMPUTE_ROSENBLUTH_DENOMS(*m_rMomx.getData(),
      *m_rMomy.getData(),
      *m_rKE.getData(),
      *m_cMomx.getData(),
      *m_cMomy.getData(),
      *m_cKE.getData(),
      *m_IVx.getData(),
      *m_IVy.getData(),
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      *m_velocities.getData());

   // Reduce m_rMomx. m_rMomy, and m_rKE to get the back reaction denominator
   // terms.
   m_IMomxD = 0.0;
   m_IMomyD = 0.0;
   m_IKED = 0.0;
   m_x_mom_reduction->execute(m_IMomxD);
   m_y_mom_reduction->execute(m_IMomyD);
   m_ke_reduction->execute(m_IKED);
}


RosenbluthCollisionOperator::RosenbluthCommunicator::RosenbluthCommunicator(
   const KineticSpecies* a_species,
   int a_maxNumFields,
   const vector<double>& a_nu,
   int a_num_interspecies_collisions)
   : m_speciesHeads(a_species->speciesHeads()),
     m_proc_hi(a_species->procHi()),
     m_species_index(a_species->speciesIndex()),
     m_num_species(a_species->numSpecies()),
     m_nu(a_nu),
     m_num_interspecies_collisions(a_num_interspecies_collisions),
     m_speciesHead(a_species->headRank()),
     m_max_num_communicated_fields(a_maxNumFields)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("collision comm");

   // For each partition of this species, communicate an identifier for its
   // partition of configuration space to the head processor.
   const int numSpeciesProcs = m_proc_hi - m_speciesHead + 1;
   int* config_space_idents = new int [numSpeciesProcs];
   const ParallelArray::Box& interior_box = a_species->interiorBox();
   const tbox::Box& domain_box = a_species->domainBox();
   int config_space_ident = interior_box.lower(X2)*domain_box.numberCells(X1) +
      interior_box.lower(X1);
   MPI_Gather(&config_space_ident, 1, MPI_INT,
      config_space_idents, 1, MPI_INT, 0, a_species->communicator());

   // Now the head processor must determine which of this species' processors
   // are associated with a given configuration space partition as specified by
   // the just communicated identifiers.
   int config_space_index = 0;
   if (Loki_Utilities::s_my_id == m_speciesHead) {
      // Use a map to identify unique configuration space identifiers and assign
      // an index to each unique identifier.  Each processor will then be
      // associated with the indexes.
      map<int, int> config_space_ident_to_index;

      // Deal with the head processor first which is guaranteed to be a unique
      // identifier.
      config_space_ident_to_index[config_space_ident] = config_space_index;
      m_pid_to_index[Loki_Utilities::s_my_id] = config_space_index;
      m_index_to_pids.insert(make_pair(config_space_index,
         Loki_Utilities::s_my_id));
      ++config_space_index;

      // Now deal with the other processors.  Their identifiers will not all be
      // unique so first test by trial insertion if it has already been found.
      // If it was actually inserted and therefore had not yet been found then
      // assign the new index to this processor.  If it was not inserted because
      // it already exists then assign the index already associated with the
      // identifier to the processor,
      for (int srcId = m_speciesHead+1; srcId <= m_proc_hi; ++srcId) {
         int this_id = config_space_idents[srcId-m_speciesHead];
         pair<map<int, int>::iterator, bool> insert_test =
            config_space_ident_to_index.insert(make_pair(this_id, config_space_index));
         if (insert_test.second) {
            m_pid_to_index[srcId] = config_space_index;
            m_index_to_pids.insert(make_pair(config_space_index, srcId));
            ++config_space_index;
         }
         else {
            m_pid_to_index[srcId] = insert_test.first->second;
            m_index_to_pids.insert(make_pair(insert_test.first->second, srcId));
         }
      }
   }
   delete [] config_space_idents;

   // Construct an array holding the range of this processor's partition of
   // configuration space.
   vector<int> interiorRange(5);
   interiorRange[0] = interior_box.lower(X1);
   interiorRange[1] = interior_box.upper(X1);
   interiorRange[2] = interior_box.lower(X2);
   interiorRange[3] = interior_box.upper(X2);

   // Since the data is identical on all processors having the same
   // configuration space partition we only need to send the data from the
   // processor occupying (0,0) in velocity space.
   m_sendConfigData =
      interior_box.lower(V1) == 0 && interior_box.lower(V2) == 0;
   interiorRange[4] = m_sendConfigData ? 1 : 0;

   // Communicate interior data sizes from the other processors to the species
   // head process and allocate/resize arrays accordingly.
   ParallelArray::Box interior_base_space(3);
   vector<int> num_global_cells(2);
   interior_base_space.lower(0) = interior_box.lower(0);
   interior_base_space.upper(0) = interior_box.upper(0);
   num_global_cells[0] = domain_box.numberCells(0);
   interior_base_space.lower(1) = interior_box.lower(1);
   interior_base_space.upper(1) = interior_box.upper(1);
   num_global_cells[1] = domain_box.numberCells(1);
   interior_base_space.lower(2) = 0;
   interior_base_space.upper(2) = m_max_num_communicated_fields-1;
   if (Loki_Utilities::s_my_id == m_speciesHead) {
      // We need data storage and data ranges for each unique partition of
      // configuration space.
      m_interiorData.resize(config_space_index);
      m_interiorDataRange.resize(config_space_index);

      // Keep track of the data storage that has been dealt with.
      deque<bool> done_partition(config_space_index, false);

      // Dimension the data storage for this processor's configuration space
      // partition.
      m_interiorData[0].partition(interior_base_space, 2, 0, num_global_cells);
      done_partition[0] = true;

      // We need to know if the head processor will get data from each of this
      // species' processors.  The head processor does not send data to itself.
      m_receiveFrom.resize(numSpeciesProcs);
      m_receiveFrom[0] = false;

      // Fill in the data range for the head processor.
      m_interiorDataRange[0].resize(4);
      for (int i = 0; i < 4; ++i) {
         m_interiorDataRange[0][i] = interiorRange[i];
      }

      // Now receive the data range from all the other processors for this
      // species.
      for (int srcId = m_speciesHead+1; srcId <= m_proc_hi; ++srcId) {
         MPI_Status receiveStatus;
         MPI_Recv(&interiorRange[0],
            interiorRange.size(),
            MPI_INT,
            srcId,
            s_TAG_SETUP,
            MPI_COMM_WORLD,
            &receiveStatus);

         // See if this other proceesor will need to send data.
         m_receiveFrom[srcId - m_speciesHead] = (interiorRange[4] == 1);

         // Get the configuration space index for this processor and, if this is
         // the first processsor associated with this parition of configuration
         // space, dimension the data storage for this partition and fill in its
         // data range.
         int this_interior_data = m_pid_to_index[srcId];
         if (!done_partition[this_interior_data]) {
            ParallelArray::Box ir_base_space(3);
            ir_base_space.lower(0) = interiorRange[0];
            ir_base_space.upper(0) = interiorRange[1];
            ir_base_space.lower(1) = interiorRange[2];
            ir_base_space.upper(1) = interiorRange[3];
            ir_base_space.lower(2) = 0;
            ir_base_space.upper(2) = m_max_num_communicated_fields-1;
            m_interiorData[this_interior_data].partition(ir_base_space,
               2,
               0,
               num_global_cells);
            m_interiorDataRange[this_interior_data].resize(4);
            for (int i = 0; i < 4; ++i) {
               m_interiorDataRange[this_interior_data][i] = interiorRange[i];
            }
            done_partition[this_interior_data] = true;
         }
      }

      // Allocate storage for the communicated data from the entire
      // configuration space for each species.
      m_domainData.resize(m_num_species);
      ParallelArray::Box domain_base_space(3);
      domain_base_space.lower(0) = domain_box.lower(0);
      domain_base_space.upper(0) = domain_box.upper(0);
      domain_base_space.lower(1) = domain_box.lower(1);
      domain_base_space.upper(1) = domain_box.upper(1);
      domain_base_space.lower(2) = 0;
      for (int i = 0; i < m_num_species; ++i) {
         if (i == m_species_index) {
            domain_base_space.upper(2) = m_max_num_communicated_fields-1;
         }
         else if (i != m_species_index && m_nu[i] != 0.0) {
            domain_base_space.upper(2) = 3;
         }
         m_domainData[i].partition(domain_base_space, 2, 0, num_global_cells);
      }

      // Determine how many other processors of this species the head processor
      // receives data from.
      m_num_species_proc_received_from = 0;
      for (int srcId = m_speciesHead+1; srcId <= m_proc_hi; ++srcId) {
         if (m_receiveFrom[srcId - m_speciesHead]) {
            m_receiveFromPID.push_back(srcId);
            ++m_num_species_proc_received_from;
         }
      }

      // Allocate requests and statuses for the non-blocking communication.
      int max_num_status = max(m_proc_hi-m_speciesHead,
         2*m_num_interspecies_collisions);
      int max_num_requests =
         max(max_num_status, m_num_species_proc_received_from);
      m_request = new MPI_Request [max_num_requests];
      m_status = new MPI_Status [max_num_status];
   }
   else {
      MPI_Send(&interiorRange[0],
         interiorRange.size(),
         MPI_INT,
         m_speciesHead,
         s_TAG_SETUP,
         MPI_COMM_WORLD);

      // On non-head processors we only need to store data associated with its
      // configuration space partition.
      m_interiorData.resize(1);
      m_interiorDataRange.resize(1);
      m_interiorDataRange[0].resize(4);
      for (int i = 0; i < 4; ++i) {
         m_interiorDataRange[0][i] = interiorRange[i];
      }

      // Dimension the data storage for this processor's configuration space
      // partition.
      m_interiorData[0].partition(interior_base_space, 2, 0, num_global_cells);
   }
   timers->stopTimer("collision comm");
}


RosenbluthCollisionOperator::RosenbluthCommunicator::~RosenbluthCommunicator()
{
   if (Loki_Utilities::s_my_id == m_speciesHead) {
      delete [] m_request;
      delete [] m_status;
   }
}


void
RosenbluthCollisionOperator::RosenbluthCommunicator::CommunicateSpatialDataToSpeciesHead(
   bool a_BR)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("collision comm");
   // Determine number of fields communicated.  Forward and back reaction send
   // different amounts of data.
   int num_communicated_fields;
   if (a_BR) {
      num_communicated_fields = 3*m_num_interspecies_collisions;
   }
   else {
      num_communicated_fields = 4;
   }

   // Communicate interior data from the species other processors that send data
   // to the species head processor.
   if (Loki_Utilities::s_my_id == m_speciesHead) {
      // Issue the non-blocking receives.
      int recvIdx = 0;
      for (int srcId = m_speciesHead+1; srcId <= m_proc_hi; ++srcId) {
         // This processor sends data.
         if (m_receiveFrom[srcId - m_speciesHead]) {
            // Get the configuration space index for this processor and receive
            // its data into the array with that index.
            int this_interior_data = m_pid_to_index[srcId];
            int data_size = m_interiorData[this_interior_data].dataBox().size()*
               num_communicated_fields/m_max_num_communicated_fields;
            MPI_Irecv(m_interiorData[this_interior_data].getData(),
               data_size,
               MPI_DOUBLE,
               srcId,
               s_TAG_SPATIAL_TO_HEAD,
               MPI_COMM_WORLD,
               &m_request[recvIdx]);
            ++recvIdx;
         }
      }

      // Check for received messages and copy that message's interior data to
      // the domain array on species head process.  Include copy of head
      // processor's interior data which is not received.
      for (int i = 0; i < m_interiorData.size(); ++i) {
         int this_interior_data;
         if (i == 0) {
            // Head processor's interior data.
            this_interior_data = 0;
         }
         else {
            // Interior data sent to head processor.
            MPI_Status stat;
            MPI_Waitany(m_num_species_proc_received_from,
                  m_request,
                  &recvIdx,
                  &stat);
            this_interior_data = m_pid_to_index[m_receiveFromPID[recvIdx]];
         }

         // Locally copy this interior data to domain array on species head
         // process.
         for (int j = m_interiorDataRange[this_interior_data][0];
              j <= m_interiorDataRange[this_interior_data][1];
              ++j) {
            for (int k = m_interiorDataRange[this_interior_data][2];
                 k <= m_interiorDataRange[this_interior_data][3];
                 ++k) {
               for (int l = 0; l < num_communicated_fields; ++l) {
                  m_domainData[m_species_index](j,k,l) =
                     m_interiorData[this_interior_data](j,k,l);
               }
            }
         }
      }
   }
   else if (m_sendConfigData) {
      // Send the data to the head processor.
      int data_size = m_interiorData[0].dataBox().size()*
         num_communicated_fields/m_max_num_communicated_fields;
      MPI_Send(m_interiorData[0].getData(),
         data_size,
         MPI_DOUBLE,
         m_speciesHead,
         s_TAG_SPATIAL_TO_HEAD,
         MPI_COMM_WORLD);
   }
   timers->stopTimer("collision comm");
}


void
RosenbluthCollisionOperator::RosenbluthCommunicator::CommunicateSpatialDataToOtherSpecies(
   bool a_BR)
{
   // Only the species head processes participate in this communication.
   if (Loki_Utilities::s_my_id != m_speciesHead) {
      return;
   }

   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("collision comm");

   // Determine amount of communicated.  Forward and back reaction send
   // different amounts of data.
   int data_size;
   if (a_BR) {
      data_size = m_domainData[m_species_index].dataBox().size()*
         3/m_max_num_communicated_fields;
   }
   else {
      data_size = m_domainData[m_species_index].dataBox().size()*
         4/m_max_num_communicated_fields;
   }

   // Send this species full configuration space data to the other species head
   // processors.  For the back reaction only the data for the other species is
   // sent so we need to pick that data out of m_domainData.
   int commIdx = 0;
   int offset = 0;
   for (int i = 0; i < m_num_species; ++i) {
      if (i != m_species_index && m_nu[i] != 0.0) {
         MPI_Isend(m_domainData[m_species_index].getData()+offset,
            data_size,
            MPI_DOUBLE,
            m_speciesHeads[i],
            s_TAG_SPATIAL_TO_OTHER_HEAD,
            MPI_COMM_WORLD,
            &m_request[commIdx]);
         ++commIdx;
         if (a_BR) {
            offset += data_size;
         }
      }
   }

   // Receive other species full configuration space data from other species
   // head processors.
   for (int i = 0; i < m_num_species; ++i) {
      if (i != m_species_index && m_nu[i] != 0.0) {
         MPI_Irecv(m_domainData[i].getData(),
            data_size,
            MPI_DOUBLE,
            m_speciesHeads[i],
            s_TAG_SPATIAL_TO_OTHER_HEAD,
            MPI_COMM_WORLD,
            &m_request[commIdx]);
         ++commIdx;
      }
   }

   // Ensure all sends and receives have completed.
   int numComms = 2*m_num_interspecies_collisions;
   MPI_Waitall(numComms, m_request, m_status);
   timers->stopTimer("collision comm");
}


void
RosenbluthCollisionOperator::RosenbluthCommunicator::CommunicateSpatialDataFromSpeciesHeadFromOtherSpecies(
   int a_otherSpecies,
   bool a_BR)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("collision comm");
   // Determine number of fields communicated.  Forward and back reaction send
   // different amounts of data.
   int num_communicated_fields;
   if (a_BR) {
      num_communicated_fields = 3;
   }
   else {
      num_communicated_fields = 4;
   }

   // The head processor communicates the appropriate partition of configuration
   // space of the data for a_otherSpecies to the other processes for this
   // species.
   if (Loki_Utilities::s_my_id == m_speciesHead) {
      // Locally copy domain data on the species head process to the interior
      // data for each unique partition of configuration space.  Then
      // communicate that interior data to each processor which shared that
      // unique partition of conffiguration space.
      for (int i = 0; i < m_interiorData.size(); ++i) {
         for (int j = m_interiorDataRange[i][0];
              j <= m_interiorDataRange[i][1];
              ++j) {
            for (int k = m_interiorDataRange[i][2];
                 k <= m_interiorDataRange[i][3];
                 ++k) {
               for (int l = 0; l < num_communicated_fields; ++l) {
                  m_interiorData[i](j,k,l) =
                     m_domainData[a_otherSpecies](j,k,l);
               }
            }
         }
         int data_size = m_interiorData[i].dataBox().size()*
            num_communicated_fields/m_max_num_communicated_fields;
         pair<multimap<int, int>::iterator, multimap<int, int>::iterator> pids = m_index_to_pids.equal_range(i);
         for (multimap<int, int>::iterator it = pids.first;
              it != pids.second;
              ++it) {
            int destId = it->second;
            if (destId != Loki_Utilities::s_my_id) {
               MPI_Isend(m_interiorData[i].getData(),
                  data_size,
                  MPI_DOUBLE,
                  destId,
                  s_TAG_OTHER_SPATIAL_FROM_HEAD,
                  MPI_COMM_WORLD,
                  &m_request[destId - m_speciesHead - 1]);
            }
         }
      }
   }
   else {
      // Receive the data.
      int data_size = m_interiorData[0].dataBox().size()*
         num_communicated_fields/m_max_num_communicated_fields;
      MPI_Status receiveStatus;
      MPI_Recv(m_interiorData[0].getData(),
         data_size,
         MPI_DOUBLE,
         m_speciesHead,
         s_TAG_OTHER_SPATIAL_FROM_HEAD,
         MPI_COMM_WORLD,
         &receiveStatus);
   }
   timers->stopTimer("collision comm");
}

} // end namespace Loki
