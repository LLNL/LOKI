/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "PitchAngleCollisionOperator.H"
#include "PitchAngleCollisionOperatorF.H"
#include "KineticSpecies.H"
#include "Loki_Defines.H"
#include "TimerManager.H"

namespace Loki {

const string PitchAngleCollisionOperator::s_CLASS_NAME(
   "Pitch Angle Collision Operator");


PitchAngleCollisionOperator::PitchAngleCollisionOperator(
   LokiInputParser& a_pp,
   const KineticSpecies* a_species)
   : CollisionOperator(a_species->domain(), NUM_DPARAMS, NUM_IPARAMS),
     m_range_lo(2),
     m_range_hi(2)
{
   // Set a default for the conservative operator option.
   m_iparameters[ICONS] = 1;
   m_iparameters[SOLN_ORDER] = a_species->spatialSolutionOrder();
   m_iparameters[DO_RELATIVITY] = Simulation::s_DO_RELATIVITY ? 1 : 0;

   // Get what the user really wants.
   parseParameters(a_pp, a_species);
}


PitchAngleCollisionOperator::~PitchAngleCollisionOperator()
{
}


void
PitchAngleCollisionOperator::initialize(
   KineticSpecies* a_species)
{
   // Always get the information about how the species this operator is
   // associated with is partitioned.
   CollisionOperator::initialize(a_species);

   // The arrays and schedules are only needed on processors that the species is
   // partitioned among.
   if (a_species->isInRange(Loki_Utilities::s_my_id)) {
      constructArrays();
      constructReductionSchedules();
   }
}


void
PitchAngleCollisionOperator::evaluate(
   KineticSpecies& a_rhs,
   const KineticSpecies& a_u,
   double a_dt,
   bool a_last_rk_stage)
{
   NULL_USE(a_dt);
   NULL_USE(a_last_rk_stage);

   // We need the reduced density number, drift velocities, and thermal velocity
   // at each RK stage.
   m_rN = 0.0;
   m_rGammax = 0.0;
   m_rGammay = 0.0;

   const ParallelArray& current_distribution = a_u.distribution();
   FORT_COMPUTE_PITCH_ANGLE_SPECIES_MOMENTS(*m_rN.getData(),
      *m_rGammax.getData(),
      *m_rGammay.getData(),
      *current_distribution.getData(),
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      *m_velocities.getData());

   m_IN = 0.0;
   m_IGammax = 0.0;
   m_IGammay = 0.0;

   m_N_reduction->execute(m_IN);
   m_Gammax_reduction->execute(m_IGammax);
   m_Gammay_reduction->execute(m_IGammay);
   FORT_COMPUTE_PITCH_ANGLE_SPECIES_REDUCED_FIELDS(*m_IVx.getData(),
      *m_IVy.getData(),
      *m_IN.getData(),
      *m_IGammax.getData(),
      *m_IGammay.getData(),
      BOX4D_TO_FORT(m_data_box));

   m_rsKEC = 0.0;
   FORT_COMPUTE_PITCH_ANGLE_SPECIES_KEC(*m_rsKEC.getData(),
      *m_IVx.getData(),
      *m_IVy.getData(),
      *current_distribution.getData(),
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      *m_velocities.getData());

   m_IsKEC = 0.0;
   m_skec_reduction->execute(m_IsKEC);

   FORT_COMPUTE_PITCH_ANGLE_SPECIES_VTHERMAL(*m_IVth.getData(),
      *m_IsKEC.getData(),
      *m_IN.getData(),
      BOX4D_TO_FORT(m_data_box));

   // Delegate the operator to fortran.
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("collision kernel");
   FORT_APPEND_PITCH_ANGLE_COLLISION(*a_rhs.distribution().getData(),
      *current_distribution.getData(),
      *m_velocities.getData(),
      *m_IVx.getData(),
      *m_IVy.getData(),
      *m_IVth.getData(),
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      PROBLEMDOMAIN_TO_FORT(m_domain),
      m_range_lo[0],
      m_range_hi[0],
      m_dparameters[0],
      m_iparameters[0]);
   timers->stopTimer("collision kernel");
}


double
PitchAngleCollisionOperator::computeRealLam(
   KineticSpecies& a_u)
{
   double dv = min(a_u.domain().dx(V1), a_u.domain().dx(V2));
   double pi = 4.0*atan(1.0);
   return m_dparameters[NU]*pow(m_dparameters[VTHERMAL_DT], 3.0)*pi*pi/(dv*dv*max(dv, m_dparameters[VFLOOR]));
}


void
PitchAngleCollisionOperator::printParameters() const
{
   // Print the user selected input.
   Loki_Utilities::printF("  Using pitch angle collision:\n" );
   Loki_Utilities::printF("    collision velocity range = %e %e %e %e\n",
      m_range_lo[0],
      m_range_hi[0],
      m_range_lo[1],
      m_range_hi[1]);
   Loki_Utilities::printF("    vfloor                   = %e\n",
      m_dparameters[VFLOOR]);
   Loki_Utilities::printF("    vthermal_dt              = %e\n",
      m_dparameters[VTHERMAL_DT]);
   Loki_Utilities::printF("    nuCoeff                  = %e\n",
      m_dparameters[NU]);
   Loki_Utilities::printF("    conservative             = %d\n",
      m_iparameters[ICONS]);
}


void
PitchAngleCollisionOperator::copyDiagnosticFields(
   vector<ParallelArray>& a_d) const
{
   // Collision diagnostics not implemented for this operator.
   a_d[0] = 0.0;
   a_d[1] = 0.0;
   a_d[2] = 0.0;
   a_d[3] = 0.0;
}


bool
PitchAngleCollisionOperator::isType(
   const string& a_name)
{
   if (a_name.compare(s_CLASS_NAME) == 0) {
      return true;
   }
   return false;
}


void
PitchAngleCollisionOperator::parseParameters(
   LokiInputParser& a_pp,
   const KineticSpecies* a_species)
{
   // vfloor is required
   if (!a_pp.query("collision_vfloor", m_dparameters[VFLOOR])) {
      LOKI_ABORT("Must supply collision_vfloor");
   }

   // These inputs are required.
   if (!a_pp.query("collision_vthermal_dt", m_dparameters[VTHERMAL_DT])) {
      LOKI_ABORT("Must supply collision_vthermal_dt");
   }
   if (!a_pp.query("collision_nuCoeff", m_dparameters[NU])) {
      LOKI_ABORT("Must supply collision_nuCoeff");
   }

   // The default is the conservative operator.
   a_pp.query("collision_conservative", m_iparameters[ICONS]);
   if (m_iparameters[ICONS] != 1 && m_iparameters[SOLN_ORDER] == 6) {
      LOKI_ABORT("Non-conservative operator in 6th order not supported.");
   }

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
}


void
PitchAngleCollisionOperator::constructArrays()
{
   // Dimension arrays then delegate the computation of them to fortran.
   int nGhosts = m_iparameters[SOLN_ORDER] == 4 ? 2 : 3;
   ParallelArray::Box config_space(2);
   vector<int> num_global_cells(2);
   config_space.lower(0) = m_interior_box.lower(0);
   config_space.upper(0) = m_interior_box.upper(0);
   num_global_cells[0] = m_domain.numberOfCells(0);
   config_space.lower(1) = m_interior_box.lower(1);
   config_space.upper(1) = m_interior_box.upper(1);
   num_global_cells[1] = m_domain.numberOfCells(1);
   m_IVx.partition(config_space, 2, nGhosts, num_global_cells);
   m_IVy.partition(config_space, 2, nGhosts, num_global_cells);
   m_IVth.partition(config_space, 2, nGhosts, num_global_cells);
   m_IN.partition(config_space, 2, nGhosts, num_global_cells);
   m_rN.partition(config_space, 2, nGhosts, num_global_cells);
   m_rGammax.partition(config_space, 2, nGhosts, num_global_cells);
   m_rGammay.partition(config_space, 2, nGhosts, num_global_cells);
   m_rsKEC.partition(config_space, 2, nGhosts, num_global_cells);
   m_IGammax.partition(config_space, 2, nGhosts, num_global_cells);
   m_IGammay.partition(config_space, 2, nGhosts, num_global_cells);
   m_IsKEC.partition(config_space, 2, nGhosts, num_global_cells);
}


void
PitchAngleCollisionOperator::constructReductionSchedules()
{
   // Build reduction schedules for density number, momenta, and kinetic energy.
   deque<bool> collapse_dir(m_domain.dim(), false);
   collapse_dir[V1] = true;
   collapse_dir[V2] = true;
   double measure(m_domain.dx(2) * m_domain.dx(3));

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

} // end namespace Loki
