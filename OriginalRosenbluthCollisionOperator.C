/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "OriginalRosenbluthCollisionOperator.H"
#include "OriginalRosenbluthCollisionOperatorF.H"
#include "CollisionOperatorF.H"
#include "RosenbluthCollisionOperatorF.H"
#include "KineticSpecies.H"
#include "Loki_Defines.H"
#include "TimerManager.H"

namespace Loki {

const string OriginalRosenbluthCollisionOperator::s_CLASS_NAME(
   "Original Rosenbluth Collision Operator");


OriginalRosenbluthCollisionOperator::OriginalRosenbluthCollisionOperator(
   LokiInputParser& a_pp,
   const KineticSpecies* a_species)
   : CollisionOperator(a_species->domain(), NUM_DPARAMS, NUM_IPARAMS),
     m_dist_func_avg(0),
     m_vthermal_method(INPUT_VTHERMAL),
     m_range_lo(2),
     m_range_hi(2)
{
   // Set what we know and read from user input what we do not know.
   int solution_order = a_species->spatialSolutionOrder();
   m_iparameters[SOLN_ORDER] = solution_order;
   m_iparameters[NUM_GHOSTS] = (solution_order == 4) ? 2 : 3;
   m_iparameters[DO_RELATIVITY] = Simulation::s_DO_RELATIVITY ? 1 : 0;
   parseParameters(a_pp, a_species);
}


OriginalRosenbluthCollisionOperator::~OriginalRosenbluthCollisionOperator()
{
   if (m_dist_func_avg != 0) {
      delete [] m_dist_func_avg;
   }
}


void
OriginalRosenbluthCollisionOperator::initialize(
   KineticSpecies* a_species)
{
   // Always get the information about how the species this operator is
   // associated with is partitioned.
   CollisionOperator::initialize(a_species);

   // The arrays and schedules are only needed on processors that the species is
   // partitioned among.
   if (a_species->isInRange(Loki_Utilities::s_my_id)) {
      constructArrays();
      if (m_back_reaction) {
         constructBackReactionArrays();
         constructBackReactionReductionSchedules();
      }
   }
}


void
OriginalRosenbluthCollisionOperator::evaluate(
   KineticSpecies& a_rhs,
   const KineticSpecies& a_u,
   double a_dt,
   bool a_last_rk_stage)
{
   NULL_USE(a_dt);
   NULL_USE(a_last_rk_stage);

   ParallelArray& rhs_distribution = a_rhs.distribution();
   const ParallelArray& current_distribution = a_u.distribution();
   // Compute vthermal and alpha if the user specified a computation method.
   double vflowx = a_u.vflowinitx();
   double vflowy = a_u.vflowinity();
   m_dparameters[VFLOWX] = vflowx;
   m_dparameters[VFLOWY] = vflowy;
   if (m_vthermal_method == LOCAL_VTHERMAL) {
      double vth = computeVthLocal(current_distribution, vflowx, vflowy);
      m_dparameters[VTHERMAL] = vth;
      m_dparameters[ALPHA] = 1.0/(vth*vth);
   }
   else if (m_vthermal_method == GLOBAL_VTHERMAL) {
      double vth = computeVthGlobal(current_distribution, vflowx, vflowy);
      m_dparameters[VTHERMAL] = vth;
      m_dparameters[ALPHA] = 1.0/(vth*vth);
   }

   // Compute the diffusion tensor.
   FORT_COMPUTE_ORIGINAL_ROSENBLUTH_DIFFUSION_TENSOR(*m_d.getData(),
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      PROBLEMDOMAIN_TO_FORT(m_domain),
      *m_velocities.getData(),
      m_range_lo[0],
      m_range_hi[0],
      m_dparameters[0],
      m_iparameters[0]);

   // Compute collision operator
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("collision kernel");
   FORT_APPEND_ORIGINAL_ROSENBLUTH_COLLISION(*rhs_distribution.getData(),
      *current_distribution.getData(),
      *m_velocities.getData(),
      *m_d.getData(),
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      m_domain.dx()[0],
      m_dparameters[0],
      m_iparameters[0]);
   timers->stopTimer("collision kernel");

   // If the back reaction is on then add it in.
   if (m_back_reaction) {
      // Compute back reaction terms.
      computeBackReactionMoments();
      m_rMomx = 0.0;
      m_rMomy = 0.0;
      m_rKE = 0.0;
      FORT_COMPUTE_ORIGINAL_ROSENBLUTH_NUMERATORS(*m_rMomx.getData(),
         *m_rMomy.getData(),
         *m_rKE.getData(),
         *m_cMomx.getData(),
         *m_cMomy.getData(),
         *m_cKE.getData(),
         *m_fM.getData(),
         *current_distribution.getData(),
         BOX4D_TO_FORT(m_data_box),
         BOX4D_TO_FORT(m_interior_box));
      m_IMomxN = 0.0;
      m_IMomyN = 0.0;
      m_IKEN = 0.0;
      m_x_mom_reduction->execute(m_IMomxN);
      m_y_mom_reduction->execute(m_IMomyN);
      m_ke_reduction->execute(m_IKEN);

      // Fold back reaction terms into collision operator.
      FORT_APPEND_ORIGINAL_ROSENBLUTH_BACK_REACTION(*rhs_distribution.getData(),
         *m_cMomx.getData(),
         *m_cMomy.getData(),
         *m_cKE.getData(),
         *m_IMomxN.getData(),
         *m_IMomyN.getData(),
         *m_IKEN.getData(),
         *m_IMomxD.getData(),
         *m_IMomyD.getData(),
         *m_IKED.getData(),
         BOX4D_TO_FORT(m_data_box),
         BOX4D_TO_FORT(m_interior_box),
         m_dparameters[0]);
   }
}


double
OriginalRosenbluthCollisionOperator::computeRealLam(
   KineticSpecies& a_u)
{
   double pi = 4.0*atan(1.0);
   return -2.0*sqrt(2.0/pi)/3.0*m_dparameters[NU]*pi*pi*
      (1.0/pow(a_u.domain().dx(V1), 2.0)+1.0/pow(a_u.domain().dx(V2), 2.0));
}


void
OriginalRosenbluthCollisionOperator::printParameters() const
{
   // Print all operator parameters in play.
   Loki_Utilities::printF("  Using Original Rosenbluth collision:\n" );
   Loki_Utilities::printF("    collision velocity range = %e %e %e %e\n",
      m_range_lo[0], m_range_hi[0], m_range_lo[1], m_range_hi[1]);
   if (m_vthermal_method == INPUT_VTHERMAL) {
      Loki_Utilities::printF("    vthermal                 = %e\n",
         m_dparameters[VTHERMAL]);
   }
   else if (m_vthermal_method == LOCAL_VTHERMAL) {
      Loki_Utilities::printF("    computing vthermal locally\n");
   }
   else if (m_vthermal_method == GLOBAL_VTHERMAL) {
      Loki_Utilities::printF("    computing vthermal globally\n");
   }
   Loki_Utilities::printF("    nuCoeff                  = %e\n",
      m_dparameters[NU]);
   if (m_vthermal_method == INPUT_VTHERMAL) {
     Loki_Utilities::printF("    alpha                     = %e\n",
        m_dparameters[ALPHA]);
   }
   else if (m_vthermal_method == LOCAL_VTHERMAL) {
      Loki_Utilities::printF("    computing alpha locally\n");
   }
   else if (m_vthermal_method == GLOBAL_VTHERMAL) {
      Loki_Utilities::printF("    computing alpha globally\n");
   }
   if (m_back_reaction) {
      Loki_Utilities::printF("    massR                    = %e\n",
         m_dparameters[MASSR]);
      Loki_Utilities::printF("    including back reaction\n");
   }
   else {
      Loki_Utilities::printF("    not including back reaction\n");
   }
}


void
OriginalRosenbluthCollisionOperator::copyDiagnosticFields(
   vector<ParallelArray>& a_d) const
{
   a_d[0] = 0.0;
   a_d[1] = 0.0;
   a_d[2] = 0.0;
   a_d[3] = 0.0;
}


bool
OriginalRosenbluthCollisionOperator::isType(
   const string& a_name)
{
   if (a_name.compare(s_CLASS_NAME) == 0) {
      return true;
   }
   return false;
}


void
OriginalRosenbluthCollisionOperator::parseParameters(
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
   else if (Simulation::s_DO_RELATIVITY) {
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
   else if (Simulation::s_DO_RELATIVITY) {
      for (int i = 0; i < 2; ++i) {
         double v = m_range_hi[i];
         m_range_hi[i] =
            a_species->mass()*v/sqrt(1.0-pow(v/Simulation::s_LIGHT_SPEED, 2));
      }
   }

   // This is kind of messy.  There are 2 basic ways that vthermal can be
   // "computed":
   // 1) Through one of 3 calculation methods:
   //    a) The "local vthermal" calculation
   //       In this case the user should not specify vthermal as that is
   //       contradictory
   //    b) The "global vthermal" calculation
   //       In this case the user should not specify vthermal as that is
   //       contradictory
   //    c) The "input vthermal"
   //       In this case the user must specify vthermal otherwise we know
   //       nothing.
   // 2) The user just gives us vthermal which is a shortcut the 1c above.
   // There is also the alpha parameter which should be 1/vth**2.  There's
   // a bunch of logic here to check for user double talk about what they
   // want to do WRT alpha.  Truthfully, the ALPHA parameter should be
   // eliminated as it is really a derived quantity.
   bool has_vthermal_method = a_pp.contains("collision_vthermal_method");
   bool has_vthermal = a_pp.contains("collision_vthermal");
   bool has_alpha = a_pp.contains("collision_alpha");
   if (!has_vthermal_method) {
      if (!has_vthermal) {
         LOKI_ABORT("If inputting collision_vthermal then it must be supplied");
         
      }
      else if (!has_alpha) {
         LOKI_ABORT("If inputting collision_vthermal then collision_alpha must be supplied");
         
      }
      else {
         a_pp.query("collision_vthermal", m_dparameters[VTHERMAL]);
         a_pp.query("collision_alpha", m_dparameters[ALPHA]);
      }
      m_vthermal_method = INPUT_VTHERMAL;
   }
   else {
      string vthermal_method("input vthermal");
      a_pp.query("collision_vthermal_method", vthermal_method);
      if (vthermal_method.compare("input vthermal") == 0) {
         if (!has_vthermal) {
            LOKI_ABORT("If inputting collision_vthermal then it must be supplied");
         }
         else if (!has_alpha) {
            LOKI_ABORT("If inputting collision_vthermal then collision_alpha must be supplied");
         
         }
         else {
            a_pp.query("collision_vthermal", m_dparameters[VTHERMAL]);
            a_pp.query("collision_alpha", m_dparameters[ALPHA]);
         }
         m_vthermal_method = INPUT_VTHERMAL;
      }
      else if (vthermal_method.compare("local vthermal") == 0) {
         if (has_vthermal) {
            LOKI_ABORT("Either compute vthermal or specify vthermal but not both");
         }
         else if (has_alpha) {
            LOKI_ABORT("Either compute alpha or specify alpha but not both");
         }
         m_vthermal_method = LOCAL_VTHERMAL;
      }
      else if (vthermal_method.compare("global vthermal") == 0) {
         if (has_vthermal) {
            LOKI_ABORT("Either compute vthermal or specify vthermal but not both");
         }
         else if (has_alpha) {
            LOKI_ABORT("Either compute alpha or specify alpha but not both");
         }
         m_vthermal_method = GLOBAL_VTHERMAL;
      }
      else {
         LOKI_ABORT("Unknown input for collision_vthermal_method");
      }
   }

   // The collision_nuCoeff is required.
   if (!a_pp.query("collision_nuCoeff", m_dparameters[NU])) {
      LOKI_ABORT("Must supply collision_nuCoeff");
   }

   // Read if the back reaction should be computed.  If it is, then the mass
   // ratio must be supplied.
   string tmp("false");
   a_pp.query("collision_back_reaction", tmp);
   m_back_reaction = tmp.compare("false") == 0 ? false : true;
   if (m_back_reaction &&
       !a_pp.query("collision_massR", m_dparameters[MASSR])) {
      LOKI_ABORT("Must supply collision_massR");
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
}


double
OriginalRosenbluthCollisionOperator::computeVthLocal(
   const ParallelArray& a_u,
   double a_vflowx,
   double a_vflowy)
{
   // Compute the mean spatial distribution locally and use it in a global
   // calculation of vthermal.
   double vth;
   FORT_COMPUTE_VTH_LOCAL(BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      a_vflowx,
      a_vflowy,
      *m_velocities.getData(),
      m_comm,
      *a_u.getData(),
      vth);
   return vth;
}


double
OriginalRosenbluthCollisionOperator::computeVthGlobal(
   const ParallelArray& a_u,
   double a_vflowx,
   double a_vflowy)
{
   // Compute the mean spatial distribution globally.
   if (m_dist_func_avg == 0) {
      m_vel_space_size =
         m_interior_box.numberOfCells(2)*m_interior_box.numberOfCells(3);
      m_dist_func_avg = new double [m_vel_space_size];

      m_config_space_size =
         m_domain.box().numberCells(0)*m_domain.box().numberCells(1);

      int color = m_interior_box.lower(3)*m_domain.box().numberCells(2)+
                  m_interior_box.lower(2);

      int comm_id;
      MPI_Comm_rank(m_comm, &comm_id);
      const int status = MPI_Comm_split(m_comm,
         color,
         comm_id,
         &m_config_space_comm);
      if (status != MPI_SUCCESS) {
         LOKI_ABORT("Configuration space splitting of MPI communicator failed");
      }

      int config_space_comm_id;
      MPI_Comm_rank(m_config_space_comm, &config_space_comm_id);
      m_is_config_space_head_node = config_space_comm_id == 0 ? 1 : 0;
   }
   FORT_COMPUTE_DIST_FUNC_AVG(BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      m_vel_space_size,
      m_config_space_size,
      *a_u.getData(),
      *m_dist_func_avg);
   MPI_Allreduce(MPI_IN_PLACE,
      m_dist_func_avg,
      m_vel_space_size,
      MPI_DOUBLE,
      MPI_SUM,
      m_config_space_comm);

   // Use the mean spatial distribution in a global calculation of vthermal.
   double vth;
   FORT_COMPUTE_VTH_GLOBAL(BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      m_vel_space_size,
      a_vflowx,
      a_vflowy,
      *m_velocities.getData(),
      *m_dist_func_avg,
      m_is_config_space_head_node,
      m_comm,
      vth);
   return vth;
}


void
OriginalRosenbluthCollisionOperator::constructArrays()
{
   // Dimension arrays.
   int nGhosts = m_iparameters[NUM_GHOSTS];
   ParallelArray::Box base_space(5);
   vector<int> num_global_cells(4);
   base_space.lower(0) = m_interior_box.lower(0);
   base_space.upper(0) = m_interior_box.upper(0);
   num_global_cells[0] = m_domain.numberOfCells(0);
   base_space.lower(1) = m_interior_box.lower(1);
   base_space.upper(1) = m_interior_box.upper(1);
   num_global_cells[1] = m_domain.numberOfCells(1);
   base_space.lower(2) = m_interior_box.lower(2)-nGhosts;
   base_space.upper(2) = m_interior_box.upper(2)+nGhosts;
   num_global_cells[2] = m_domain.numberOfCells(2);
   base_space.lower(3) = m_interior_box.lower(3)-nGhosts;
   base_space.upper(3) = m_interior_box.upper(3)+nGhosts;
   num_global_cells[3] = m_domain.numberOfCells(3);
   base_space.lower(4) = 0;
   base_space.upper(4) = 2;
   m_d.partition(base_space, 4, 0, num_global_cells);
}


void
OriginalRosenbluthCollisionOperator::constructBackReactionArrays()
{
   int nGhosts = m_iparameters[NUM_GHOSTS];
   vector<int> num_global_cells4D(4);
   for (int dim = 0; dim < 4; ++dim) {
      num_global_cells4D[dim] = m_domain.numberOfCells(dim);
   }

   // Dimension persistent 4D momenta, kinetic energy, and maxwellian
   // distribution.
   m_cMomx.partition(m_interior_box, 4, nGhosts, num_global_cells4D);
   m_cMomy.partition(m_interior_box, 4, nGhosts, num_global_cells4D);
   m_cKE.partition(m_interior_box, 4, nGhosts, num_global_cells4D);
   m_fM.partition(m_interior_box, 4, nGhosts, num_global_cells4D);

   // Dimension intermediate 2D locally reduced momenta and kinetic energy.
   ParallelArray::Box config_space(2);
   vector<int> num_global_cells2D(2);
   config_space.lower(0) = m_interior_box.lower(0);
   config_space.upper(0) = m_interior_box.upper(0);
   num_global_cells2D[0] = m_domain.numberOfCells(0);
   config_space.lower(1) = m_interior_box.lower(1);
   config_space.upper(1) = m_interior_box.upper(1);
   num_global_cells2D[1] = m_domain.numberOfCells(1);
   m_rMomx.partition(config_space, 2, nGhosts, num_global_cells2D);
   m_rMomy.partition(config_space, 2, nGhosts, num_global_cells2D);
   m_rKE.partition(config_space, 2, nGhosts, num_global_cells2D);

   // Dimension persistent 2D globally reduced momenta and kinetic energy
   // denominator terms.
   m_IMomxD.partition(config_space, 2, nGhosts, num_global_cells2D);
   m_IMomyD.partition(config_space, 2, nGhosts, num_global_cells2D);
   m_IKED.partition(config_space, 2, nGhosts, num_global_cells2D);

   // Dimension persistent 2D globally reduced momenta and kinetic energy
   // numerator terms.
   m_IMomxN.partition(config_space, 2, nGhosts, num_global_cells2D);
   m_IMomyN.partition(config_space, 2, nGhosts, num_global_cells2D);
   m_IKEN.partition(config_space, 2, nGhosts, num_global_cells2D);
}


void
OriginalRosenbluthCollisionOperator::constructBackReactionReductionSchedules()
{
   // Build reduction schedules for momenta and kinetic energy denominators.
   deque<bool> collapse_dir(m_domain.dim(), false);
   collapse_dir[V1] = true;
   collapse_dir[V2] = true;
   double measure(m_domain.dx(2) * m_domain.dx(3));

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


void
OriginalRosenbluthCollisionOperator::computeBackReactionMoments()
{
   // Construct 4D momenta and kinetic energy temps.
   vector<int> num_global_cells(4);
   for (int i = 0; i < 4; ++i) {
      num_global_cells[i] = m_domain.numberOfCells(i);
   }
   ParallelArray uMomx(m_data_box, 4, 0, num_global_cells);
   ParallelArray uMomy(m_data_box, 4, 0, num_global_cells);
   ParallelArray uKE(m_data_box, 4, 0, num_global_cells);

   // Compute temporary 4D momenta and kinetic energy and persistent 4D
   // maxwellian distribution.
   FORT_COMPUTE_ORIGINAL_ROSENBLUTH_TEMPS(*m_fM.getData(),
      *uMomx.getData(),
      *uMomy.getData(),
      *uKE.getData(),
      BOX4D_TO_FORT(m_data_box),
      *m_velocities.getData(),
      m_dparameters[0]);

   // Compute persistent 4D momenta and kinetic energy.
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("collision kernel");
   m_cMomx = 0.0;
   FORT_APPEND_ORIGINAL_ROSENBLUTH_COLLISION(*m_cMomx.getData(),
      *uMomx.getData(),
      *m_velocities.getData(),
      *m_d.getData(),
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      m_domain.dx()[0],
      m_dparameters[0],
      m_iparameters[0]);
   m_cMomy = 0.0;
   FORT_APPEND_ORIGINAL_ROSENBLUTH_COLLISION(*m_cMomy.getData(),
      *uMomy.getData(),
      *m_velocities.getData(),
      *m_d.getData(),
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      m_domain.dx()[0],
      m_dparameters[0],
      m_iparameters[0]);
   m_cKE = 0.0;
   FORT_APPEND_ORIGINAL_ROSENBLUTH_COLLISION(*m_cKE.getData(),
      *uKE.getData(),
      *m_velocities.getData(),
      *m_d.getData(),
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
   FORT_COMPUTE_ORIGINAL_ROSENBLUTH_DENOMS(*m_rMomx.getData(),
      *m_rMomy.getData(),
      *m_rKE.getData(),
      *m_cMomx.getData(),
      *m_cMomy.getData(),
      *m_cKE.getData(),
      BOX4D_TO_FORT(m_data_box),
      BOX4D_TO_FORT(m_interior_box),
      *m_velocities.getData(),
      m_dparameters[0]);

   // Reduce m_rMomx. m_rMomy, and m_rKE to get the back reaction denominator
   // terms.
   m_x_mom_reduction->execute(m_IMomxD);
   m_y_mom_reduction->execute(m_IMomyD);
   m_ke_reduction->execute(m_IKED);
}

} // end namespace Loki
