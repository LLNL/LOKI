/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "KineticSpecies.H"

#include "tbox/IntVector.H"
#include "ICFactory.H"
#include "TZSourceFactory.H"
#include "ElectricFieldDriverFactory.H"
#include "CollisionOperatorFactory.H"
#include "Maxwell.H"
#include "Poisson.H"
#include "RestartManager.H"
#include "getVelocityF.H"

#include <deque>
#include <sstream>

namespace Loki {

KineticSpecies::KineticSpecies(
   RestartReader& a_reader,
   int a_species_num,
   int a_number_of_species,
   int a_spatial_solution_order,
   int a_temporal_solution_order,
   bool a_plot_ke_vel_bdy_flux,
   bool a_use_new_bcs,
   const string& a_name)
   : m_pdim(PDIM),
     m_cdim(CDIM),
     m_name("undefined"),
     m_species_index(a_species_num-1),
     m_number_of_species(a_number_of_species),
     m_mass(-1.0),
     m_charge(0.0),
     m_spatial_solution_order(a_spatial_solution_order),
     m_temporal_solution_order(a_temporal_solution_order),
     m_stencil_width(a_spatial_solution_order+1),
     m_global_box(m_pdim),
     m_number_of_procs(1),
     m_fixed_number_of_procs(false),
     m_partition_defined(false),
     m_comm(MPI_COMM_NULL),
     m_ef_drivers(0),
     m_num_external_drivers(0),
     m_old_driver_syntax(false),
     m_do_maxwell(false),
     m_collision_operators(0),
     m_num_collision_operators(0),
     m_problem_has_particles(false),
     m_use_new_bcs(a_use_new_bcs),
     m_plot_ke_vel_bdy_flux(a_plot_ke_vel_bdy_flux),
     m_bz_const(0.0), // IEO
     m_integrated_ke_e_dot(0.0)
{
   // Set number of ghosts based of order of solution.
   if (m_spatial_solution_order == 4) {
      m_n_ghosts = 2;
   }
   else {
      m_n_ghosts = 3;
   }

   // Read how many processors created this restart dump.
   int generating_processes;
   a_reader.readIntegerValue("generating processes", generating_processes);

   // find subdatabase with the name of this distribution
   m_name = a_name;
   a_reader.pushSubDir(m_name);

   // get the pdim, cdim from database
   int tmp_pdim, tmp_cdim;
   a_reader.readIntegerValue("pdim", tmp_pdim);
   a_reader.readIntegerValue("cdim", tmp_cdim);
   tbox::Dimension pdim(static_cast<unsigned short>(tmp_pdim));
   tbox::Dimension cdim(static_cast<unsigned short>(tmp_cdim));

   // get the mass and charge from database
   a_reader.readDoubleValue("mass", m_mass);
   a_reader.readDoubleValue("charge", m_charge);
   a_reader.readDoubleValue("bz_const", m_bz_const);

   cout << "In KineticSpecies db constructor, creating new ProblemDomain"
        << endl;
   m_domain = new ProblemDomain(pdim, a_reader);
   cout << "Returned from creating new ProblemDomain" << endl;

   // Dimension the distributed array.
   ParallelArray::Box base_space(m_pdim);
   vector<int> num_global_cells(m_pdim);
   for (int d = 0; d < m_pdim; ++d) {
      base_space.lower(d) = 0;
      base_space.upper(d) = m_domain->numberOfCells(d)-1;
      num_global_cells[d] = m_domain->numberOfCells(d);
   }
   m_dist_func.partition(base_space, m_pdim, m_n_ghosts, num_global_cells);

   // read restart distribution from database and now that the local data is
   // defined, get the local array
   a_reader.readParallelArray("distribution",
      m_dist_func,
      generating_processes);
   cout << "Returned from getDistributed (K.S. constructor)" << endl;
   a_reader.popSubDir();
}


KineticSpecies::KineticSpecies(
   const tbox::Pointer<ProblemDomain>& a_cfg_domain,
   LokiInputParser& a_pp,
   int a_species_num,
   int a_number_of_species,
   int a_spatial_solution_order,
   int a_temporal_solution_order,
   bool a_plot_ke_vel_bdy_flux,
   bool a_use_new_bcs,
   bool a_do_maxwell,
   double a_bz_const)
   : m_pdim(PDIM),
     m_cdim(CDIM),
     m_name("undefined"),
     m_species_index(a_species_num-1),
     m_number_of_species(a_number_of_species),
     m_mass(-1.0),
     m_charge(0.0),
     m_spatial_solution_order(a_spatial_solution_order),
     m_temporal_solution_order(a_temporal_solution_order),
     m_stencil_width(a_spatial_solution_order+1),
     m_global_box(m_pdim),
     m_number_of_procs(1),
     m_fixed_number_of_procs(false),
     m_partition_defined(false),
     m_comm(MPI_COMM_NULL),
     m_ef_drivers(0),
     m_num_external_drivers(0),
     m_old_driver_syntax(false),
     m_do_maxwell(a_do_maxwell),
     m_collision_operators(0),
     m_num_collision_operators(0),
     m_problem_has_particles(false),
     m_use_new_bcs(a_use_new_bcs),
     m_plot_ke_vel_bdy_flux(a_plot_ke_vel_bdy_flux),
     m_bz_const(a_bz_const), //IEO
     m_integrated_ke_e_dot(0.0)
{
   // Set number of ghosts based of order of solution.
   if (m_spatial_solution_order == 4) {
      m_n_ghosts = 2;
   }
   else {
      m_n_ghosts = 3;
   }

   // Read all the user input for this species and construct much of its data
   // members.
   parseParameters(a_pp, a_cfg_domain);

   // Create the global box to be the domain + ghosts.
   m_global_box = m_domain->box();
   m_global_box.grow(m_n_ghosts);

   // Now construct all the entities held by the species that must read the
   // input deck.

   // Get the sub-database for this species' initial condition and construct the
   // initial condition object.
   ostringstream input_string;
   input_string << "kinetic_species." << a_species_num << ".ic";
   LokiInputParser ic_pp(input_string.str().c_str());
   input_string.str("");
   m_initial_condition = ICFactory::create(ic_pp, this);

   // Get the sub-database for this species' twilight zone and construct the
   // twilight zone object if one is specified.
   input_string << "kinetic_species." << a_species_num << ".tz";
   LokiInputParser tz_pp(input_string.str().c_str());
   input_string.str("");
   m_tz_source = TZSourceFactory::create(tz_pp);

   // For each of this species' external E field drivers get its sub-database
   // and construct the field driver object.
   for (int i = 0; i < m_num_external_drivers; ++i) {
      if (m_old_driver_syntax) {
         input_string << "kinetic_species." << a_species_num
                      << ".external_driver";
      }
      else {
         input_string << "kinetic_species." << a_species_num
                      << ".external_driver." << i+1;
      }
      LokiInputParser efdf_pp(input_string.str().c_str());
      input_string.str("");
      m_ef_drivers[i] = ElectricFieldDriverFactory::create(efdf_pp,
         i+1,
         a_species_num);
   }

   // For each of this species' collision operators get its sub-database and
   // construct the collision operator object.
   for (int i = 0; i < m_num_collision_operators; ++i) {
      input_string << "kinetic_species." << a_species_num
                   << ".collision_operator." << i+1;
      LokiInputParser coll_op_pp(input_string.str().c_str());
      input_string.str("");
      m_collision_operators[i] = CollisionOperatorFactory::create(coll_op_pp,
         this);
   }

   // Get the sub-database for this species' Krook layer and construct it if
   // one is specified.
   input_string << "kinetic_species." << a_species_num << ".krook";
   LokiInputParser krook_pp(input_string.str().c_str());
   input_string.str("");
   m_krook_layer = new KrookLayer(m_cdim, krook_pp, *m_domain);

   // Get the sub-database for this species' ExternalDistKrookLayer and
   // construct it if one is specified.
   input_string << "kinetic_species." << a_species_num
                << ".external_dist_krook";
   LokiInputParser external_dist_krook_pp(input_string.str().c_str());
   input_string.str("");
   m_external_dist_krook = new ExternalDistKrookLayer(m_cdim,
      external_dist_krook_pp,
      *m_domain);

   // The KineticSpecies write restart data so register this object with the
   // RestartManager which will use the putToRestart/getFromRestart callbacks to
   // get the restart data written/read.
   RestartManager* restart_manager(RestartManager::getManager());
   restart_manager->registerRestart(this);
}


KineticSpecies::KineticSpecies(
   const KineticSpecies& a_other)
   : Load(),
     Serializable(),
     m_pdim(a_other.m_pdim),
     m_cdim(a_other.m_cdim),
     m_global_box(m_pdim)
{
   // Copy all the contents and build new schedules.
   copy(a_other);
   defineExtEfieldContractionSchedule();
   defineChargeDensityReductionSchedule();
   if (m_plot_ke_vel_bdy_flux) {
      defineKineticEnergySummationSchedule();
   }
   defineMomentReductionSchedule();
   defineDiagnosticSummationSchedule();
}


KineticSpecies::~KineticSpecies()
{
}


void
KineticSpecies::addData(
   const KineticSpecies& a_rhs,
   double a_factor)
{
   // Check that we're copying between similar species and add a_factor times
   // a_rhs' distribution function to this species' distribution function.  If
   // there are any external field drivers then do the same to the
   // integrated_ke_e_dot time history.
   if (conformsTo(a_rhs, false)) {
      FORT_XPBY_4D(*m_dist_func.getData(),
         *a_rhs.m_dist_func.getData(),
         a_factor,
         BOX4D_TO_FORT(dataBox()),
         BOX4D_TO_FORT(interiorBox()));
      if (m_num_external_drivers > 0) {
         m_integrated_ke_e_dot += a_rhs.m_integrated_ke_e_dot*a_factor;
      }
   }
}


bool
KineticSpecies::conformsTo(
   const KineticSpecies& a_rhs,
   bool a_include_ghost_cells) const
{
   // Check that the basic defining properties of the 2 species match.
   if (m_name != a_rhs.m_name ||
       m_species_index != a_rhs.m_species_index ||
       m_mass != a_rhs.m_mass ||
       m_charge != a_rhs.m_charge ||
       m_domain != a_rhs.m_domain ||
       m_proc_lo != a_rhs.m_proc_lo ||
       m_proc_hi != a_rhs.m_proc_hi ||
       m_partition_defined != a_rhs.m_partition_defined) {
      return false;
   }

   // Check that the boxes of these 2 species match.
   if (!a_include_ghost_cells) {
      return (m_domain->box() == (a_rhs.m_domain)->box());
   }
   return (m_global_box == (a_rhs.m_global_box));
}


void
KineticSpecies::copy(
   const KineticSpecies& a_rhs)
{
   if ((m_pdim != a_rhs.m_pdim) || (m_cdim != a_rhs.m_cdim)) {
      LOKI_ABORT("Attemtpt to copy incongruent species!");
   }

   // If the 2 species are different, then copy their internals.
   if (&a_rhs != this) {
      m_name = a_rhs.m_name;
      m_species_index = a_rhs.m_species_index;
      m_number_of_species = a_rhs.m_number_of_species;
      m_mass = a_rhs.m_mass;
      m_charge = a_rhs.m_charge;
      m_bz_const = a_rhs.m_bz_const; //IEO
      m_domain = a_rhs.m_domain;
      m_initial_condition = a_rhs.m_initial_condition;
      m_tz_source = a_rhs.m_tz_source;
      m_krook_layer = a_rhs.m_krook_layer;
      m_external_dist_krook = a_rhs.m_external_dist_krook;
      m_number_of_procs = a_rhs.m_number_of_procs;
      m_proc_lo = a_rhs.m_proc_lo;
      m_proc_hi = a_rhs.m_proc_hi;
      m_fixed_number_of_procs = a_rhs.m_fixed_number_of_procs;
      m_partition_defined = a_rhs.m_partition_defined;
      m_comm = a_rhs.m_comm;
      m_num_external_drivers = a_rhs.m_num_external_drivers;
      m_old_driver_syntax = a_rhs.m_old_driver_syntax;
      m_ef_drivers = a_rhs.m_ef_drivers;
      m_do_maxwell = a_rhs.m_do_maxwell;
      m_collision_operators = a_rhs.m_collision_operators;
      m_num_collision_operators = a_rhs.m_num_collision_operators;
      m_problem_has_particles = a_rhs.m_problem_has_particles;
      m_use_new_bcs = a_rhs.m_use_new_bcs;
      m_plot_ke_vel_bdy_flux = a_rhs.m_plot_ke_vel_bdy_flux;

      m_speciesHeads = a_rhs.m_speciesHeads;
      m_speciesMass = a_rhs.m_speciesMass;
      m_dist_func = a_rhs.m_dist_func;

      m_spatial_solution_order = a_rhs.m_spatial_solution_order;
      m_temporal_solution_order = a_rhs.m_temporal_solution_order;
      m_n_ghosts = a_rhs.m_n_ghosts;
      m_stencil_width = a_rhs.m_stencil_width;
      m_global_box = a_rhs.m_global_box;

      m_lambda_max.resize(a_rhs.m_lambda_max.size());
      for (int i(0); i < static_cast<int>(a_rhs.m_lambda_max.size()); ++i) {
         m_lambda_max[i] = a_rhs.m_lambda_max[i];
      }
      m_vel_face = a_rhs.m_vel_face;
      m_flux = a_rhs.m_flux;
      m_u_face = a_rhs.m_u_face;

      if (m_do_maxwell) {
         m_em_expansion_schedule = a_rhs.m_em_expansion_schedule;
         m_vz_expansion_schedule = a_rhs.m_vz_expansion_schedule;
         m_em_vars = a_rhs.m_em_vars;
         m_vz = a_rhs.m_vz;
      }
      else {
         m_efield_expansion_schedule = a_rhs.m_efield_expansion_schedule;
         m_accel = a_rhs.m_accel;
      }
      if (m_num_external_drivers) {
         m_ext_efield = a_rhs.m_ext_efield;
      }
      m_velocities = a_rhs.m_velocities;
      m_vxface_velocities = a_rhs.m_vxface_velocities;
      m_vyface_velocities = a_rhs.m_vyface_velocities;
      m_integrated_ke_e_dot = a_rhs.m_integrated_ke_e_dot;
   }
}


void
KineticSpecies::copySolnData(
   const KineticSpecies& a_rhs)
{
   if ((m_pdim != a_rhs.m_pdim) || (m_cdim != a_rhs.m_cdim)) {
      LOKI_ABORT("Attemtpt to copy incongruent species!");
   }

   // If the 2 species are different, then copy the distribution function.  If
   // there are any external field drivers then also copy the integrated
   // ke_e_dot time history.
   if (&a_rhs != this) {
      m_dist_func = a_rhs.m_dist_func;
      if (m_num_external_drivers > 0) {
         m_integrated_ke_e_dot = a_rhs.m_integrated_ke_e_dot;
      }
   }
}


void
KineticSpecies::printParameters() const
{
   // Print this species parameters and the parameters of the entities that it
   // holds and are not accessible elsewhere like the external E field drivers,
   // collision operators, initial conditions, and twilight zones.
   Loki_Utilities::printF("\n#*#*# Kinetic Species %s #*#*#\n", m_name.c_str());
   Loki_Utilities::printF("  species index:   %d\n", m_species_index);
   Loki_Utilities::printF("  mass:            %e\n", m_mass);
   Loki_Utilities::printF("  charge:          %e\n", m_charge);
   Loki_Utilities::printF("  constant Bz             = %e\n", m_bz_const); //IEO
   m_domain->printParameters();
   m_krook_layer->printParameters();
   m_external_dist_krook->printParameters();
   Loki_Utilities::printF("\n  num external drivers    = %i\n",
      m_num_external_drivers);
   for (int i = 0; i < m_num_external_drivers; ++i) {
      m_ef_drivers[i]->printParameters();
   }
   Loki_Utilities::printF("\n  num collision operators = %i\n",
      m_num_collision_operators);
   for (int i = 0; i < m_num_collision_operators; ++i) {
      m_collision_operators[i]->printParameters();
   }
   m_initial_condition->printParameters();
   if (m_tz_source) {
      m_tz_source->printParameters();
   }
}


float
KineticSpecies::netCost() const
{
   // The computational cost of a species is proportional to the size of its
   // domain.
   float cost_per_cell(1.0);
   return cost_per_cell *
      static_cast<float>((m_domain->numberOfCells()).getProduct());
}


int
KineticSpecies::numberOfProcessors() const
{
   return m_number_of_procs;
}


bool
KineticSpecies::fixedNumberOfProcessors() const
{
   return m_fixed_number_of_procs;
}


void
KineticSpecies::createPartition(
   int a_proc_lo,
   int a_proc_hi,
   const MPI_Comm& a_comm)
{
   m_comm = a_comm;
   m_proc_lo = a_proc_lo;
   m_proc_hi = a_proc_hi;
   m_number_of_procs = a_proc_hi-a_proc_lo+1;
   m_partition_defined = true;

   // Partition the distribution function among its processors.
   deque<bool> is_periodic(m_pdim);
   vector<int> num_cells(m_pdim);
   for (int dim = 0; dim < m_pdim; ++dim) {
      is_periodic[dim] = m_domain->isPeriodic(dim);
      num_cells[dim] = m_domain->numberOfCells(dim);
   }
   m_dist_func.partition(m_pdim,
      m_pdim,
      a_proc_lo,
      a_proc_hi,
      m_n_ghosts,
      is_periodic,
      num_cells);

   // Check to make sure decomposition makes sense; we need at least
   // m_stencil_width interior points in each direction
   if (isInRange(Loki_Utilities::s_my_id)) {
      for (int dir(0); dir < m_pdim; ++dir) {
         int npts = interiorBox().numberOfCells(dir);
         if (npts < m_stencil_width) {
            ostringstream msg;
            msg << "Too few interior points in decomposition in direction "
                << dir;
            LOKI_ABORT(msg.str());
         }
      }
   }

   // Now that we know the various boxes we can define the schedules and local
   // arrays needed by this species.
   defineChargeDensityReductionSchedule();
   if (m_plot_ke_vel_bdy_flux) {
      defineKineticEnergySummationSchedule();
   }
   defineMomentReductionSchedule();
   defineDiagnosticSummationSchedule();

   allocateAuxArrays();

   int config_space_id =
      interiorBox().lower(X2)*m_domain->box().numberCells(X1) +
      interiorBox().lower(X1);
   if (m_do_maxwell) {
      m_em_expansion_schedule =
         new ExpansionSchedule(m_dist_func, *m_domain, m_comm);
      m_vz_expansion_schedule =
         new ExpansionSchedule(m_dist_func, *m_domain, m_comm);
   }
   else {
      m_efield_expansion_schedule =
         new ExpansionSchedule(m_dist_func, *m_domain, m_comm);
   }

   // Now that we have the boxes and local arrays we can initialize the
   // velocities on each face.
   initializeVelocity();
}


void
KineticSpecies::SetSpeciesHeadsAndMass(
   const int numSpecies,
   const int *heads,
   const double *mass)
{
   const int thisHead = headRank();
   if (thisHead != heads[m_species_index]) {
      LOKI_ABORT("Inconsistent species head process ranks.");
   }

   m_speciesHeads.clear();
   m_speciesMass.clear();
   for (int i=0; i<numSpecies; ++i) {
      m_speciesHeads.push_back(heads[i]);
      m_speciesMass.push_back(mass[i]);
   }
   for (int i = 0; i < m_num_collision_operators; ++i) {
      m_collision_operators[i]->initialize(this);
   }
}


bool
KineticSpecies::isInRange(
   int a_proc_id) const
{
   // Returns true if this species is partitioned onto this processor.
   return ((a_proc_id >= m_dist_func.procLo()) &&
           (a_proc_id <= m_dist_func.procHi()));
}


void
KineticSpecies::printDecomposition() const
{
   // This function is only valid if we actually know the decomposition.
   if (m_partition_defined) {
      // Print some basic decomposition info.
      Loki_Utilities::printF("  Kinetic Species \"%s\" processor(s):  [%d,%d]\n",
         m_name.c_str(),
         m_dist_func.procLo(),
         m_dist_func.procHi());

      // Now do some sanity checking.
      // A prime number of processors is almost certainly a bad idea.
      int num_procs = m_dist_func.procHi() - m_dist_func.procLo() + 1;
      bool num_procs_prime;
      if (num_procs > 3) {
         int max_div = int(sqrt(num_procs));
         num_procs_prime = true;
         for (int i = max_div; i > 1; --i) {
            if (num_procs % i == 0) {
               num_procs_prime = false;
               break;
            }
         }
      }
      else if (num_procs == 1) {
         num_procs_prime = false;
      }
      else {
         num_procs_prime = true;
      }
      if (num_procs_prime) {
         Loki_Utilities::printF("  Kinetic Species \"%s\" is partitioned across %d processors\n"
            "  which is prime and most likely not what you want\n",
            m_name.c_str(),
            num_procs);
      }

      // Give the user some idea about how evenly the different dimensions are
      // subdivided.
      int loc_zones[4], glob_zones[4];
      if (isInRange(Loki_Utilities::s_my_id)) {
         const ParallelArray::Box& data_box = dataBox();
         loc_zones[0] = data_box.numberOfCells(0);
         loc_zones[1] = data_box.numberOfCells(1);
         loc_zones[2] = data_box.numberOfCells(2);
         loc_zones[3] = data_box.numberOfCells(3);
      }
      else {
         loc_zones[0] = INT_MIN;
         loc_zones[1] = INT_MIN;
         loc_zones[2] = INT_MIN;
         loc_zones[3] = INT_MIN;
      }
      MPI_Reduce(&loc_zones[0], &glob_zones[0], 4,
                 MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
      Loki_Utilities::printF("  Kinetic Species \"%s\" maximum dimensions:  [%d,%d,%d,%d]\n",
         m_name.c_str(),
         glob_zones[0], glob_zones[1],
         glob_zones[2], glob_zones[3]);
      if (!isInRange(Loki_Utilities::s_my_id)) {
         loc_zones[0] = INT_MAX;
         loc_zones[1] = INT_MAX;
         loc_zones[2] = INT_MAX;
         loc_zones[3] = INT_MAX;
      }
      MPI_Reduce(&loc_zones[0], &glob_zones[0], 4,
                 MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
      Loki_Utilities::printF("  Kinetic Species \"%s\" minimum dimensions:  [%d,%d,%d,%d]\n",
         m_name.c_str(),
         glob_zones[0], glob_zones[1], 
         glob_zones[2], glob_zones[3]);
   }
}


double
KineticSpecies::computeDt()
{
   // local max lambdas are found every time the velocity is recomputed
   vector<double> lambda_max(m_pdim);
   Loki_Utilities::getMaxValues(&m_lambda_max[0],
      &lambda_max[0],
      m_pdim,
      -1,
      m_comm);

   double reLam(0.0);
   double imLam(0.0);
   double pi = 4.0*atan(1.0);

   // eigenvalue of Vlasov operator
   for (int dir(X1); dir < m_pdim; ++dir) {
      imLam += pi * lambda_max[dir] / m_domain->dx(dir);
   }
   // eigenvalue of collision operator
   for (int i = 0; i < m_num_collision_operators; ++i) {
      double thisReLam = m_collision_operators[i]->computeRealLam(*this);
      if (abs(thisReLam) > reLam) {
         reLam = thisReLam;
      }
   }

   /*// return the maximal time step satisfying (reLam*dt)^2+(imLam*dt)^2=alpha^2
   double ddt;
   double alpha(2.6);
   ddt = sqrt(alpha*alpha/(reLam*reLam+imLam*imLam));*/

   // return the maximal time step satisfying (reLam*dt/alpha)^2+(imLam*dt/beta)^2=1
   double ddt, alpha, beta;
   if( m_temporal_solution_order == 4 ) {
     // RK4
     alpha = 2.6;
     beta  = 2.6;
   }
   else {
     // RK6
     alpha = 4.95;
     beta  = 3.168;
   }
   ddt = sqrt( 1.0/(reLam*reLam/(alpha*alpha)+imLam*imLam/(beta*beta)) );

   return(ddt);
}


void
KineticSpecies::computeAcceleration(
   Poisson& a_poisson,
   ParallelArray& a_ext_efield,
   double a_time,
   double a_dt,
   bool a_first_rk_stage)
{
   TimerManager* timers(TimerManager::getManager());

   // We need to know the self consistent E field computed by the Poisson
   // process.  Communicate the part of this 2D field corresponding to each
   // species' configuration space extent from the Poisson processor to each
   // species.
   m_accel = 0.0;
   m_efield_expansion_schedule->execute(a_poisson.getEMVars(), m_accel);

   // apply external driver if applicable
   //   Note that the driver is applied in all cells (including ghost cells)
   //   and periodicity is not thereafter enforced ... that is on the user.
   if (m_num_external_drivers > 0) {
      timers->startTimer("field driver");
      // The externally applied E field only needs to be tracked for the update
      // of particles and the integrated ke added by the field.  This is denoted
      // by the driver_sums_into flag.  A value of 1 indicates m_ext_efield
      // only, 2 indicates both m_ext_efield and m_accel.
      int driver_sums_into = 2;
      m_ext_efield = 0.0;

      for (int i = 0; i < m_num_external_drivers; ++i) {
         m_ef_drivers[i]->evaluate(m_accel,
            m_ext_efield,
            *m_domain,
            driver_sums_into,
            a_time,
            a_dt,
            a_first_rk_stage);
      }

      // If there are particles we must communicate this species' external
      // efield to a_ext_field which is defined on the Poisson processor(s)
      // which own the particles and compute their equations of motion.
      if (m_problem_has_particles) {
         m_ext_efield_schedule->execute(a_ext_efield);
      }
      timers->stopTimer("field driver");
   }

   timers->startTimer("blowout");
   // locally turn into acceleration (for non-relativistic case) or force (for
   // relativistic case)
   double normalization;
   if (Simulation::s_DO_RELATIVITY) {
      normalization = m_charge;
   }
   else {
      normalization = m_charge / m_mass;
   }
   m_accel *= normalization;

   // Compute the accelerations on the faces in the V1 and V2 directions and the
   // maximum acclerations which are necessary for the time step calculation.
   double axmax, aymax;
   FORT_SET_PHASE_SPACE_VEL_4D(*(m_vel_face[V1].getData()),
      *(m_vel_face[V2].getData()),
      BOX4D_TO_FORT(dataBox()),
      BOX4D_TO_FORT(interiorBox()),
      *(m_vxface_velocities.getData()),
      *(m_vyface_velocities.getData()),
      normalization,
      m_bz_const,
      *(m_accel.getData()),
      BOX2D_TO_FORT(m_accel.dataBox()),
      axmax, aymax);
   m_lambda_max[V1] = axmax;
   m_lambda_max[V2] = aymax;
   timers->stopTimer("blowout");
}


void
KineticSpecies::computeAcceleration(
   Maxwell& a_maxwell,
   ParallelArray& a_ext_efield,
   double a_time,
   double a_dt,
   bool a_first_rk_stage)
{
   TimerManager* timers(TimerManager::getManager());

   // We need to know the electromagnetic fields computed by the Maxwell
   // process.  Communicate the part of these 2D fields corresponding to each
   // species' configuration space extent from the Maxwell processor to each
   // species.a
   m_em_vars = 0.0;
   m_em_expansion_schedule->execute(a_maxwell.getEMVars(), m_em_vars);

   // apply external driver if applicable
   //   Note that the driver is applied in all cells (including ghost cells)
   //   and periodicity is not thereafter enforced ... that is on the user.
   if (m_num_external_drivers > 0) {
      timers->startTimer("field driver");
      // The externally applied E field only needs to be tracked for the update
      // of particles and the integrated ke added by the field.  This is denoted
      // by the driver_sums_into flag.  A value of 1 indicates m_ext_efield
      // only, 2 indicates both m_ext_efield and m_em_vars.
      int driver_sums_into = 2;
      m_ext_efield = 0.0;

      for (int i = 0; i < m_num_external_drivers; ++i) {
         m_ef_drivers[i]->evaluate(m_em_vars,
            m_ext_efield,
            *m_domain,
            driver_sums_into,
            a_time,
            a_dt,
            a_first_rk_stage);
      }

      // If there are particles we must communicate this species' external
      // efield to a_ext_field which is defined on the Maxwell processor(s)
      // which own the particles and compute their equations of motion.
      if (m_problem_has_particles) {
         m_ext_efield_schedule->execute(a_ext_efield);
      }
      timers->stopTimer("field driver");
   }

   timers->startTimer("blowout");
   // Compute the accelerations on the faces in the V1 and V2 directions and the
   // maximum acclerations which are necessary for the time step calculation.
   double axmax, aymax, normalization;
   if (Simulation::s_DO_RELATIVITY) {
      normalization = m_charge;
   }
   else {
      normalization = m_charge / m_mass;
   }
   FORT_SET_PHASE_SPACE_VEL_MAXWELL_4D(*(m_vel_face[V1].getData()),
      *(m_vel_face[V2].getData()),
      BOX4D_TO_FORT(dataBox()),
      BOX4D_TO_FORT(interiorBox()),
      *(m_vxface_velocities.getData()),
      *(m_vyface_velocities.getData()),
      normalization,
      m_bz_const,
      *m_em_vars.getData(),
      *m_vz.getData(),
      axmax,
      aymax);
   m_lambda_max[V1] = axmax;
   m_lambda_max[V2] = aymax;
   timers->stopTimer("blowout");
}


void
KineticSpecies::currentDensity(
   const Maxwell& a_maxwell,
   ParallelArray& a_Jx,
   ParallelArray& a_Jy,
   ParallelArray& a_Jz)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("current density");

   // First time through define the local current densities and the reduction
   // schedules needed to compute the current densities.
   if (!m_jx_schedule) {
      defineCurrentDensityReductionSchedules();
   }

   // Get the local z drift velocity from Maxwell.
   m_vz = 0.0;
   m_vz_expansion_schedule->execute(a_maxwell.getVZVar(m_species_index), m_vz);

   // Zero out the local current densities.
   m_Jx = 0.0;
   m_Jy = 0.0;
   m_Jz = 0.0;

   // Compute the 4D current densities.
   FORT_COMPUTE_CURRENTS_4D(BOX4D_TO_FORT(dataBox()),
      BOX4D_TO_FORT(interiorBox()),
      *m_velocities.getData(),
      *m_dist_func.getData(),
      *m_vz.getData(),
      *m_Jx.getData(),
      *m_Jy.getData(),
      *m_Jz.getData());

   // Now reduce the 4D current densities into 2D quantites defined on the
   // Maxwell processors by integrating over the velocities.
   m_jx_schedule->execute(a_Jx);
   m_jy_schedule->execute(a_Jy);
   m_jz_schedule->execute(a_Jz);

   timers->stopTimer("current density");
}


void
KineticSpecies::getFromRestart(
   RestartReader& a_reader)
{
   // Get the version.  One value may not exist but we may restart without it.
   int major_version, minor_version, patch_level;
   a_reader.readIntegerValue("major version", major_version);
   a_reader.readIntegerValue("minor version", minor_version);
   a_reader.readIntegerValue("patch level", patch_level);

   // Read the field driver state.
   for (int i = 0; i < m_num_external_drivers; ++i) {
      m_ef_drivers[i]->getFromDatabase(a_reader);
   }

   // find subdirectory with the name of this distribution
   a_reader.pushSubDir(m_name);

   m_domain->getFromDatabase(a_reader);

   // read restart distribution from database
   a_reader.readParallelArray("distribution", m_dist_func);

   // If there are any external drivers then read the integrated E dot J.
   // This quantity was saved to restart starting with version 3.0.1.
   if (m_num_external_drivers > 0 &&
       major_version >= 3 && minor_version >= 0 && patch_level >= 1) {
      a_reader.readBulkDoubleValue("integrated_e_dot_j", m_integrated_ke_e_dot);
   }

   // Zero out the ghost data.  Can this be removed now that ParallelArray is
   // used?
   FORT_ZERO_GHOST_4D(*m_dist_func.getData(),
      BOX4D_TO_FORT(interiorBox()),
      BOX4D_TO_FORT(dataBox()));

   // Read in any KrookLayer and ExternalDistKrookLayer.
   m_krook_layer->getFromDatabase(a_reader);
   m_external_dist_krook->getFromDatabase(a_reader);

   a_reader.popSubDir();
}


void
KineticSpecies::putToRestart(
   RestartWriter& a_writer,
   double a_time)
{
   // Write everything except Krook data.
   putToRestart_SkipKrook(a_writer, a_time, false);

   // Now write Krook data and ExternalDistKrookLayer.
   bool write_data = Loki_Utilities::s_my_id == m_dist_func.procLo();
   m_krook_layer->putToDatabase(a_writer, write_data);
   m_external_dist_krook->putToDatabase(a_writer, write_data);
   a_writer.popSubDir();
}


void
KineticSpecies::putToRestart_SkipKrook(
   RestartWriter& a_writer,
   double a_time,
   bool a_return_to_root)
{
   bool write_data = Loki_Utilities::s_my_id == m_dist_func.procLo();
   
   // Save any field drivers state.
   for (int i = 0; i < m_num_external_drivers; ++i) {
      m_ef_drivers[i]->putToDatabase(a_writer, write_data);
   }

   // make a subdirectory with the  name of this distribution
   a_writer.pushSubDir(m_name);

   // save the pdim, cdim dimensions
   a_writer.writeIntegerValue("pdim", m_pdim, write_data);
   a_writer.writeIntegerValue("cdim", m_cdim, write_data);

   // save the mass and charge
   a_writer.writeDoubleValue("mass", m_mass, write_data);
   a_writer.writeDoubleValue("charge", m_charge, write_data);
   a_writer.writeDoubleValue("bz_const", m_bz_const, write_data); //IEO

   // The ProblemDomain must be saved as well.
   m_domain->putToDatabase(a_writer, write_data);

   // save the distribution function
   if (m_tz_source) {
      // If a twilight zone is defined, we want to see the error in the
      // distribution functions, not the functions themselves.
      ParallelArray tz_error_array_test(m_dist_func);
      if (isInRange(Loki_Utilities::s_my_id)) {
         m_tz_source->computeError(tz_error_array_test,
            m_dist_func,
            *m_domain,
            a_time,
            m_velocities);
      }
      else {
         tz_error_array_test = 0.0;
      }
      a_writer.writeParallelArray("distribution",
         tz_error_array_test,
         write_data);
   }
   else {
      a_writer.writeParallelArray("distribution",
         m_dist_func,
         write_data);
   }

   // If there are any external drivers then write the integrated E dot J.
   if (m_num_external_drivers > 0) {
      a_writer.writeBulkDoubleValue("integrated_e_dot_j",
         m_integrated_ke_e_dot);
   }

   if (a_return_to_root) {
      a_writer.popSubDir();
   }
}


void
KineticSpecies::completeRHS(
   KineticSpecies& a_rhs,
   double a_time,
   double a_dt,
   bool a_use_new_alg,
   bool a_last_rk_stage) const
{
   // If using the flux based formulation sum in the flux divergence.
   if (!a_use_new_alg) {
      addFluxDivergence(a_rhs);
   }

   // If any collision operators are present apply each of them.
   TimerManager* timers(TimerManager::getManager());
   if (m_num_collision_operators > 0) {
      timers->startTimer("collisions");
   }
   for (int i = 0; i < m_num_collision_operators; ++i) {
      m_collision_operators[i]->evaluate(a_rhs, *this, a_dt, a_last_rk_stage);
   }
   if (m_num_collision_operators > 0) {
      timers->stopTimer("collisions");
   }

   // If any Krook layers are present apply them.
   if (m_krook_layer->hasKrookLayer() && m_krook_layer->overlaps()) {
      timers->startTimer("krook");
      FORT_APPEND_KROOK(
         BOX4D_TO_FORT(dataBox()),
         BOX4D_TO_FORT(interiorBox()),
         a_dt,
         int64_t(m_initial_condition.getPointer()),
         *m_krook_layer->nu().getData(),
         *m_dist_func.getData(),
         *(a_rhs.m_dist_func.getData()));
      timers->stopTimer("krook");
   }

   // If any ExternalDistKrookLayer is present apply it.
   if (m_external_dist_krook->hasKrookLayer() &&
       m_external_dist_krook->overlaps()) {
      timers->startTimer("krook");
      FORT_APPEND_KROOK(
         BOX4D_TO_FORT(dataBox()),
         BOX4D_TO_FORT(interiorBox()),
         a_dt,
         int64_t(m_external_dist_krook->externalDistIC()),
         *m_external_dist_krook->nu().getData(),
         *m_dist_func.getData(),
         *(a_rhs.m_dist_func.getData()));
      timers->stopTimer("krook");
   }

   // If using the Twilight Zone then set the rhs to the appropriate state.
   if (m_tz_source) {
      m_tz_source->set(a_rhs.m_dist_func, *m_domain, a_time, m_velocities);
   }

   // If there are any external drivers compute the instantanious rate at which
   // they add energy to the species.
   if (m_num_external_drivers > 0) {
      FORT_COMPUTE_KE_E_DOT(BOX4D_TO_FORT(dataBox()),
         BOX4D_TO_FORT(interiorBox()),
         PROBLEMDOMAIN_TO_FORT((*m_domain)),
         *m_dist_func.getData(),
         m_charge,
         *m_velocities.getData(),
         *m_ext_efield.getData(),
         a_rhs.m_integrated_ke_e_dot);
   }
}


void
KineticSpecies::copyCollisionDiagnostics(
   const int a_collOperIndex,
   vector<ParallelArray>& a_diags)
{
   // If the specified collision operator exists then get its diagnostic data.
   // Otherwise use 0 for the diagnostic data as initialized in
   // defineDiagnosticReductionSchedule.
   if (a_collOperIndex < m_num_collision_operators) {
      m_collision_operators[a_collOperIndex]->copyDiagnosticFields(m_diagnostics);
   }

   for (int i=0; i<CollisionOperator::s_DIAGNOSTIC_WORK_SIZE; ++i) {
      int offset = m_species_index*2*CollisionOperator::s_DIAGNOSTIC_WORK_SIZE;
      a_diags[offset + i] = 0.0;

      m_diagnostics[CollisionOperator::s_DIAGNOSTIC_WORK_SIZE] =
         m_diagnostics[i];
      m_diagnostic_schedule->execute(a_diags[offset + i]);
   }
}


void
KineticSpecies::copyMomentDiagnostics(
   vector<ParallelArray>& a_diags)
{
   // Species momentum components, kinetic energy, entropy.
   // Zero out the local 4D moments.
   m_momx = 0.0;
   m_momy = 0.0;
   m_ke = 0.0;
   m_ent = 0.0;

   // Compute the 4D moments.  This is only performed by the Vlasov processors.

   FORT_COMPUTE_MOM_4D(BOX4D_TO_FORT(dataBox()),
      BOX4D_TO_FORT(interiorBox()),
      *m_velocities.getData(),
      *m_dist_func.getData(),
      *m_momx.getData(),
      *m_momy.getData(),
      *m_ke.getData(),
      *m_ent.getData());

   // Now reduce the 4D moments into 2D quantities by integrating over the
   // velocities.  The 2D quantities are owned by the non-Vlasov processor(s).

   for (int i=CollisionOperator::s_DIAGNOSTIC_WORK_SIZE;
        i<2*CollisionOperator::s_DIAGNOSTIC_WORK_SIZE; ++i) {
      a_diags[(m_species_index*2*CollisionOperator::s_DIAGNOSTIC_WORK_SIZE) + i] = 0.0;
   }

   int offset = m_species_index*2*CollisionOperator::s_DIAGNOSTIC_WORK_SIZE;
   m_momx_schedule->execute(a_diags[offset + CollisionOperator::s_DIAGNOSTIC_WORK_SIZE]);
   m_momy_schedule->execute(a_diags[offset + CollisionOperator::s_DIAGNOSTIC_WORK_SIZE + 1]);
   m_ke_schedule->execute(a_diags[offset + CollisionOperator::s_DIAGNOSTIC_WORK_SIZE + 2]);
   m_ent_schedule->execute(a_diags[offset + CollisionOperator::s_DIAGNOSTIC_WORK_SIZE + 3]);
}


void
KineticSpecies::computeKEVelBdyFlux(
   EMSolverBase& a_em_solver)
{
   // We need the up to date acceleration boundary conditions and acceleration
   // fluxes to do this.
   setAccelerationBCs();
   evalAccelerationFluxes();

   // Loop over the 2 velocity boundaries in each velocity direction.  Then sum
   // the contribution to this flux from each processor and communicate the sum
   // to the Poisson processor.
   const tbox::Box& domain_box = m_domain->box();
   for (int dir = V1; dir <= V2; ++dir) {
      for (int side = LO; side <= HI; ++side) {
         m_ke_vel_bdry_flux = 0.0;
         FORT_COMPUTE_KE_VEL_SPACE_FLUX(BOX4D_TO_FORT(dataBox()),
            BOX4D_TO_FORT(interiorBox()),
            BOX4D_TO_FORT(domain_box),
            m_domain->dx()[0],
            *m_flux[V1].getData(),
            *m_flux[V2].getData(),
            *m_ke_vel_bdry_flux.getData(),
            m_mass,
            *m_vxface_velocities.getData(),
            *m_vyface_velocities.getData(),
            side,
            dir);
         m_ke_vel_bdry_flux_schedule->execute(a_em_solver.getKEFluxVar(m_species_index, dir, side));
      }
   }
}


void
KineticSpecies::accumulateSequences(
   vector<vector<double> >& a_sequences,
   int a_saved_seq,
   int& a_seq_idx,
   double a_time)
{
   // VP specific species Kinetic Energy.
   // Compute the local kinetic energy.
   double ke = 0.0;
   double ke_x = 0.0;
   double ke_y = 0.0;
   double px = 0.0;
   double py = 0.0;
   FORT_COMPUTE_KE(BOX4D_TO_FORT(dataBox()),
      BOX4D_TO_FORT(interiorBox()),
      PROBLEMDOMAIN_TO_FORT((*m_domain)),
      *m_dist_func.getData(),
      m_mass,
      *m_velocities.getData(),
      ke,
      ke_x,
      ke_y,
      px,
      py);

   // Sum the local kinetic energies and add them to the time history.
   a_sequences[a_seq_idx++][a_saved_seq] = Loki_Utilities::getSum(ke, -1);
   a_sequences[a_seq_idx++][a_saved_seq] = Loki_Utilities::getSum(ke_x, -1);
   a_sequences[a_seq_idx++][a_saved_seq] = Loki_Utilities::getSum(ke_y, -1);
   a_sequences[a_seq_idx++][a_saved_seq] = Loki_Utilities::getSum(px, -1);
   a_sequences[a_seq_idx++][a_saved_seq] = Loki_Utilities::getSum(py, -1);

   // Time histories common to VM and VP systems.
   accumulateSequencesCommon(a_sequences, a_saved_seq, a_seq_idx, a_time);
}


void
KineticSpecies::accumulateSequences(
   const Maxwell& a_maxwell,
   vector<vector<double> >& a_sequences,
   int a_saved_seq,
   int& a_seq_idx,
   double a_time)
{
   // VM specific species Kinetic Energy.
   // Compute the local kinetic energy.
   m_vz = 0.0;
   m_vz_expansion_schedule->execute(a_maxwell.getVZVar(m_species_index), m_vz);
   double ke = 0.0;
   double ke_x = 0.0;
   double ke_y = 0.0;
   double px = 0.0;
   double py = 0.0;
   FORT_COMPUTE_KE_MAXWELL(BOX4D_TO_FORT(dataBox()),
      BOX4D_TO_FORT(interiorBox()),
      PROBLEMDOMAIN_TO_FORT((*m_domain)),
      *m_dist_func.getData(),
      m_mass,
      *m_velocities.getData(),
      *m_vz.getData(),
      ke,
      ke_x,
      ke_y,
      px,
      py);

   // Sum the local kinetic energies and add them to the time history.
   a_sequences[a_seq_idx++][a_saved_seq] = Loki_Utilities::getSum(ke, -1);
   a_sequences[a_seq_idx++][a_saved_seq] = Loki_Utilities::getSum(ke_x, -1);
   a_sequences[a_seq_idx++][a_saved_seq] = Loki_Utilities::getSum(ke_y, -1);
   a_sequences[a_seq_idx++][a_saved_seq] = Loki_Utilities::getSum(px, -1);
   a_sequences[a_seq_idx++][a_saved_seq] = Loki_Utilities::getSum(py, -1);

   // Time histories common to VM and VP systems.
   accumulateSequencesCommon(a_sequences, a_saved_seq, a_seq_idx, a_time);
}


void
KineticSpecies::accumulateMomentSequences(
   bool a_is_2d_proc,
   vector<vector<double> >& a_sequences,
   int a_saved_seq,
   int& a_seq_idx,
   ParallelArray& a_momx_global,
   ParallelArray& a_momy_global,
   ParallelArray& a_ke_global,
   ParallelArray& a_ent_global)
{
   // Species momentum components, kinetic energy, entropy.
   // Zero out the local 4D moments.
   m_momx = 0.0;
   m_momy = 0.0;
   m_ke = 0.0;
   m_ent = 0.0;

   // Compute the 4D moments.  This is only performed by the Vlasov
   // processors.
   FORT_COMPUTE_MOM_4D(BOX4D_TO_FORT(dataBox()),
      BOX4D_TO_FORT(interiorBox()),
      *m_velocities.getData(),
      *m_dist_func.getData(),
      *m_momx.getData(),
      *m_momy.getData(),
      *m_ke.getData(),
      *m_ent.getData());

   // Now reduce the 4D moments into 2D quantities by integrating over the
   // velocities.  The 2D quantities are owned by the non-Vlasov processor(s).
   m_momx_schedule->execute(a_momx_global);
   m_momy_schedule->execute(a_momy_global);
   m_ke_schedule->execute(a_ke_global);
   m_ent_schedule->execute(a_ent_global);

   // Compute the local scalar moments from the local 2D moments
   // (integrate over the configuration space dimensions).  This is only
   // performed by the non-Vlasov processor(s).
   double momx = 0.0;
   double momy = 0.0;
   double ke = 0.0;
   double ent = 0.0;
   const tbox::Box& domain_box = m_domain->box();
   if (a_is_2d_proc) {
      FORT_COMPUTE_MOM_2D(BOX2D_TO_FORT(a_momx_global.dataBox()),
         BOX2D_TO_FORT(a_momx_global.interiorBox()),
         m_domain->dx()[0],
         *a_momx_global.getData(),
         *a_momy_global.getData(),
         *a_ke_global.getData(),
         *a_ent_global.getData(),
         momx, momy, ke, ent);
   }

   // Sum the local scalar kinetic energy and add it to the time history.
   a_sequences[a_seq_idx++][a_saved_seq] = Loki_Utilities::getSum(momx, -1);
   a_sequences[a_seq_idx++][a_saved_seq] = Loki_Utilities::getSum(momy, -1);
   a_sequences[a_seq_idx++][a_saved_seq] = Loki_Utilities::getSum(ke, -1);
   a_sequences[a_seq_idx++][a_saved_seq] = Loki_Utilities::getSum(ent, -1);
}


void
KineticSpecies::accumulateCollisionSequences(
   double a_dt,
   int a_coll_op_idx,
   vector<vector<double> >& a_sequences,
   int a_saved_seq,
   int& a_seq_idx)
{
   // The collision diagnostic fields have already been reduced to 2D arrays
   // in m_diagnostics by copyCollisionDiagnostics.

   // Compute the local scalar diagnostic fields from the local 2D fields
   // (integrate over the configuration space dimensions).
   // This is only performed by the non-Vlasov processor(s).

   // If the specified collision operator exists then get its diagnostic data.
   // Otherwise use 0 for the diagnostic data as initialized in
   // defineDiagnosticReductionSchedule.
   if (a_coll_op_idx < m_num_collision_operators) {
      m_collision_operators[a_coll_op_idx]->copyDiagnosticFields(m_diagnostics);
   }

   double scalar = 0.0;

   for (int i=0; i<CollisionOperator::s_DIAGNOSTIC_WORK_SIZE; ++i) {
      FORT_INTEGRATE_2D(BOX2D_TO_FORT(dataBox()),
         BOX2D_TO_FORT(interiorBox()),
         m_domain->dx()[0],
         *(m_diagnostics[i].getData()),
         scalar);

      // Sum the local scalar and add it to the time history.
      a_sequences[a_seq_idx++][a_saved_seq] =
         Loki_Utilities::getSum(scalar, -1);
   }
}

//// PRIVATE METHODS ////////////////////////////////////////////////////

void
KineticSpecies::parseParameters(
   LokiInputParser& a_pp,
   const tbox::Pointer<ProblemDomain> a_cfg_domain)
{
   // A species name is required.
   if (a_pp.contains("name")) {
      string tmp;
      a_pp.get("name", tmp);
      m_name = tmp;
   }
   else {
      LOKI_ABORT("Must supply name!");
   }

   // A species mass is required.
   if (a_pp.contains("mass")) {
      a_pp.get("mass", m_mass);
   }
   else {
      LOKI_ABORT("Must supply mass!");
   }

   // A species charge is required.
   if (a_pp.contains("charge")) {
      a_pp.get("charge", m_charge);
   }
   else {
      LOKI_ABORT("Must supply charge!");
   }

   // This is a bit of backward compatibility for old decks when there was only
   // one external driver possible for each species.
   bool new_driver_syntax = a_pp.contains("num_external_drivers");
   bool old_driver_syntax = a_pp.contains("apply_external_driver");
   if (new_driver_syntax) {
      if (old_driver_syntax) {
         LOKI_ABORT("Mixed old and new electric field driver syntax!");
      }
      else {
         a_pp.query("num_external_drivers", m_num_external_drivers);
      }
   }
   else if (old_driver_syntax) {
      m_old_driver_syntax = true;
      string driver_on("false");
      a_pp.query("apply_external_driver", driver_on);
      m_num_external_drivers = driver_on.compare("false") == 0 ? 0 : 1;
   }
   m_ef_drivers.resize(m_num_external_drivers,
      tbox::Pointer<ElectricFieldDriver>(0));

   // See if there's any collision operators.
   a_pp.query("num_collision_operators", m_num_collision_operators);
   m_collision_operators.resize(m_num_collision_operators,
      tbox::Pointer<CollisionOperator>(0));

   // Get the configuration space limits from the (2D) ProblemDomain.
   vector<double> limits(2 * m_pdim);
   for (int d(0); d < m_cdim; ++d) {
      limits[2*d] = a_cfg_domain->lower(d);
      limits[2*d+1] = a_cfg_domain->upper(d);
   }

   // Now get the velocity space limits which are required.
   if (a_pp.contains("velocity_limits")) {
      tbox::Dimension vdim(static_cast<unsigned short>(m_pdim - m_cdim));
      int num_vlimits = 2*vdim;
      vector<double> vlimits(num_vlimits);
      a_pp.getarr("velocity_limits", vlimits, 0, num_vlimits);
      for (int i = 0; i < vdim; ++i) {
         if (vlimits[2*i] >= vlimits[2*i+1]) {
            LOKI_ABORT("Velocity space lower bound >= upper bound.");
         }
      }
      if (Simulation::s_DO_RELATIVITY) {
         double vxmax = max(fabs(vlimits[0]), fabs(vlimits[1]));
         double vymax = max(fabs(vlimits[2]), fabs(vlimits[3]));
         if (vxmax >= Simulation::s_LIGHT_SPEED ||
             vymax >= Simulation::s_LIGHT_SPEED) {
            LOKI_ABORT("Max velocity exceeds light speed.");
         }
      }
      for (int i(0); i < num_vlimits; ++i) {
         int j(2 * m_cdim + i);
         if (Simulation::s_DO_RELATIVITY) {
            double v = vlimits[i];
            limits[j] = m_mass*v/sqrt(1.0-pow(v/Simulation::s_LIGHT_SPEED, 2));
         }
         else {
            limits[j] = vlimits[i];
         }
      }
   }
   else {
      LOKI_ABORT("Must supply velocity_limits!");
   }

   // Dump all the limits into x_lo and x_hi so we can create the species'
   // ProblemDomain.
   vector<double> x_lo(m_pdim);
   vector<double> x_hi(m_pdim);
   for (int d(0); d < m_pdim; ++d) {
      x_lo[d] = limits[2*d];
      x_hi[d] = limits[2*d+1];
   }

   // Get the configuration space periodicity from the (2D) ProblemDomain.
   // Velocity space periodicity is always off.
   deque<bool> is_periodic(m_pdim);
   for (int d(0); d < m_cdim; ++d) {
      is_periodic[d] = a_cfg_domain->isPeriodic(d);
   }
   for (int d(m_cdim); d < m_pdim; ++d) {
      is_periodic[d] = false;
   }

   // Get the configuration space number of cells from the (2D) ProblemDomain
   // and read the velocity space number of cells which is required.
   tbox::IntVector n_cells(m_pdim);
   for (int d(0); d < m_cdim; ++d) {
      n_cells[d] = a_cfg_domain->numberOfCells(d);
   }
   if (a_pp.contains("Nv")) {
      vector<int> tmp(m_pdim);
      a_pp.getarr("Nv", tmp, 0, m_pdim - m_cdim);
      for (int d(m_cdim); d < m_pdim; ++d) {
         n_cells[d] = tmp[d-m_cdim];
      }
   }
   else {
      LOKI_ABORT("Must supply Nv!");
   }

   // Now create this species' ProblemDomain.
   m_domain = new ProblemDomain(m_pdim, n_cells, x_lo, x_hi, is_periodic);

   // Check for invalid old initial condition flow velocity syntax.
   if (a_pp.contains("vflowinitx") || a_pp.contains("vflowinity")) {
      LOKI_ABORT("vflowinitx and vflowinity are now specified in the initial condition.");
   }

   // Now we know the velocity space info so set m_lambda_max for X1 and X2.
   // This is kind of hoky.  m_domain may be set up WRT momentum space where as
   // m_lambda_max depends on velocity space.  So we do a bunch of domain setup
   // like calculations here to figure out these quantities.
   // First zero out the lambdas
   m_lambda_max.resize(m_pdim);
   for (int dir(X1); dir < m_pdim; ++dir) {
      m_lambda_max[dir] = 0.0;
   }

   // Now set the lambdas from the upper and lower velocity bounds.
   double vlo, vhi;
   if (Simulation::s_DO_RELATIVITY) {
      double px, py, vx, vx1, vx2, vy, vy1, vy2;

      // Get the maximum x velocity.
      px = limits[4] + 0.5 * (limits[5] - limits[4])/n_cells[V1];
      py = 0.0;
      FORT_GET_VELOCITY(px, py, m_mass, Simulation::s_LIGHT_SPEED, vx1, vy);
      px = limits[5] + 0.5 * (limits[5] - limits[4])/n_cells[V1];
      FORT_GET_VELOCITY(px, py, m_mass, Simulation::s_LIGHT_SPEED, vx2, vy);
      m_lambda_max[X1] = max(fabs(vx1), fabs(vx2));

      // Get the maximum y velocity.
      px = 0.0;
      py = limits[6] + 0.5 * (limits[7] - limits[6])/n_cells[V2];
      FORT_GET_VELOCITY(px, py, m_mass, Simulation::s_LIGHT_SPEED, vx, vy1);
      py = limits[7] + 0.5 * (limits[7] - limits[6])/n_cells[V2];
      FORT_GET_VELOCITY(px, py, m_mass, Simulation::s_LIGHT_SPEED, vx, vy2);
      m_lambda_max[X2] = max(fabs(vy1), fabs(vy2));
   }
   else {
      double vlo, vhi;
      vlo = limits[4] + 0.5 * (limits[5] - limits[4])/n_cells[V1];
      vhi = limits[5] + 0.5 * (limits[5] - limits[4])/n_cells[V1];
      m_lambda_max[X1] = max(fabs(vlo), fabs(vhi));

      vlo = limits[6] + 0.5 * (limits[7] - limits[6])/n_cells[V2];
      vhi = limits[7] + 0.5 * (limits[7] - limits[6])/n_cells[V2];
      m_lambda_max[X2] = max(fabs(vlo), fabs(vhi));
   }
}


void
KineticSpecies::allocateAuxArrays()
{
   // We need the distribution function, the velocity/acceleration, and the
   // distribution function flux on each face.
   m_u_face.resize(m_pdim, ParallelArray(m_pdim, m_pdim));
   m_vel_face.resize(m_pdim, ParallelArray(m_pdim, m_pdim));
   m_flux.resize(m_pdim, ParallelArray(m_pdim, m_pdim));

   ParallelArray::Box base_space(m_pdim);
   vector<int> num_global_cells_pdim(m_pdim);
   for (int dir(0); dir < m_pdim; ++dir) {
      base_space.lower(0) = interiorBox().lower(dir);
      base_space.upper(0) = interiorBox().upper(dir) + 1;
      num_global_cells_pdim[0] = m_domain->numberOfCells(dir)+1;
      for (int i = 1; i < m_pdim; ++i) {
         int k = (dir+i)%m_pdim;
         base_space.lower(i) = interiorBox().lower(k);
         base_space.upper(i) = interiorBox().upper(k);
         num_global_cells_pdim[i] = m_domain->numberOfCells(k);
      }
      m_u_face[dir].partition(base_space, m_n_ghosts, num_global_cells_pdim);
      m_vel_face[dir].partition(base_space, m_n_ghosts, num_global_cells_pdim);
      m_flux[dir].partition(base_space, m_n_ghosts, num_global_cells_pdim);
   }

   vector<int> num_global_cells_cdim(m_cdim);
   for (int dim = 0; dim < m_cdim; ++dim) {
      num_global_cells_cdim[dim] = m_domain->numberOfCells(dim);
   }
   if (m_num_external_drivers != 0) {
      ParallelArray::Box base_space_ext_e(m_cdim+1);
      for (int dim = 0; dim < m_cdim; ++dim) {
         base_space_ext_e.lower(dim) = interiorBox().lower(dim);
         base_space_ext_e.upper(dim) = interiorBox().upper(dim);
      }
      base_space_ext_e.lower(m_cdim) = 0;
      base_space_ext_e.upper(m_cdim) = 1;
      m_ext_efield.partition(base_space_ext_e,
         m_cdim,
         m_n_ghosts,
         num_global_cells_cdim);
   }

   // Figure out the box corresponding to the configuration space part of this
   // species on this processor.
   ParallelArray::Box base_space_em(m_cdim+1);
   for (int dim(0); dim < m_cdim; ++dim) {
      base_space_em.lower(dim) = interiorBox().lower(dim);
      base_space_em.upper(dim) = interiorBox().upper(dim);
   }
   base_space_em.lower(m_cdim) = 0;

   // Create arrays that depend on the EM system requested.
   if (m_do_maxwell) {
      // Dimension the 2D electromagnetic field.
      base_space_em.upper(m_cdim) = Maxwell::NUM_EM_VARS-1;
      m_em_vars.partition(base_space_em,
         m_cdim,
         m_n_ghosts,
         num_global_cells_cdim);

      // Dimension the 2D z drift velocity of this species.
      ParallelArray::Box base_space_vz(m_cdim);
      for (int dim(0); dim < m_cdim; ++dim) {
         base_space_vz.lower(dim) = interiorBox().lower(dim);
         base_space_vz.upper(dim) = interiorBox().upper(dim);
      }
      m_vz.partition(base_space_vz,
         m_cdim,
         m_n_ghosts,
         num_global_cells_cdim);
   }
   else {
      // Dimension the 2D accleration.
      base_space_em.upper(m_cdim) = 1;
      m_accel.partition(base_space_em,
         m_cdim,
         m_n_ghosts,
         num_global_cells_cdim);
   }

   // Allocate 2D arrays on the "velocity" space and fill them with the
   // velocities.  We need both cell and face centered velocities.
   // If we're running a relativistic problem "velocity" space is actually
   // momentum space and these arrays will contain the relativistic
   // transformation of the momentum at each point.  In this case, storing the
   // velocities and looking them up when needed is much more efficient than
   // computing them each time.  If the problem is non-relativistic these arrays
   // are not strictly necessary as it is trivial to compute the velocity when
   // needed but it's still probably faster to do the look-up and the API of all
   // the code that needs velocities is independent of relativity.
   buildVelocityArrays();
}


void
KineticSpecies::initializeVelocity()
{
   // Zero out the velocities.
   for (int dir(X1); dir < m_pdim; ++dir) {
      m_vel_face[dir] = 0.0;
   }

   // Set the velocity on the X1 and X2 faces.
   for (int i3(dataBox().lower(V1));
        i3 <= dataBox().upper(V1); ++i3) {
     for (int i4(dataBox().lower(V2));
          i4 <= dataBox().upper(V2); ++i4) {
         double coord_x3 = m_velocities(i3, i4, 0);
         for (int i2(dataBox().lower(X2));
              i2 <= dataBox().upper(X2); ++i2) {
            for (int i1(dataBox().lower(X1));
                 i1 <= dataBox().upper(X1)+1; ++i1) {
               m_vel_face[X1](i1, i2, i3, i4) = coord_x3;
            }
         }
      }
   }

   for (int i4(dataBox().lower(V2));
        i4 <= dataBox().upper(V2); ++i4) {
      for (int i3(dataBox().lower(V1));
           i3 <= dataBox().upper(V1); ++i3) {
         double coord_x4 = m_velocities(i3, i4, 1);
         for (int i1(dataBox().lower(X1));
              i1 <= dataBox().upper(X1); ++i1) {
            for (int i2(dataBox().lower(X2));
                 i2 <= dataBox().upper(X2)+1; ++i2){
               m_vel_face[X2](i2, i3, i4, i1) = coord_x4;
            }
         }
      }
   }
}


void
KineticSpecies::setPeriodicBCs()
{
   m_dist_func.communicatePeriodicBoundaries();
}


void
KineticSpecies::addFluxDivergence(
   KineticSpecies& a_rhs) const
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("Vlasov");

   FORT_ACCUM_FLUX_DIV_4D(*(a_rhs.m_dist_func).getData(),
      BOX4D_TO_FORT(dataBox()),
      BOX4D_TO_FORT(interiorBox()),
      *(m_flux[X1].getData()),
      *(m_flux[X2].getData()),
      *(m_flux[V1].getData()),
      *(m_flux[V2].getData()),
      m_domain->dx()[0]);

   timers->stopTimer("Vlasov");
}


void
KineticSpecies::defineExtEfieldContractionSchedule()
{
   if (m_num_external_drivers != 0 && m_problem_has_particles) {
      m_ext_efield_schedule = new ContractionSchedule(m_dist_func,
         m_ext_efield,
         m_domain->box(),
         m_comm);
   }
}


void
KineticSpecies::defineChargeDensityReductionSchedule()
{
   // We want to integrate over velocity.
   deque<bool> collapse_dir(m_pdim, false);
   collapse_dir[V1] = true;
   collapse_dir[V2] = true;
   double measure(m_domain->dx(2) * m_domain->dx(3));

   // Now make the schedule for the charge density.
   m_moment_schedule = new ReductionSchedule(m_dist_func,
      *m_domain,
      collapse_dir,
      m_n_ghosts,
      measure,
      m_charge,
      m_comm);
}


void
KineticSpecies::defineMomentReductionSchedule()
{
   // Dimension this species momentum components.
   deque<bool> is_periodic(m_pdim);
   vector<int> num_cells(m_pdim);
   for (int dim = 0; dim < m_pdim; ++dim) {
      is_periodic[dim] = m_domain->isPeriodic(dim);
      num_cells[dim] = m_domain->numberOfCells(dim);
   }
   m_momx.partition(m_pdim,
      m_pdim,
      m_proc_lo,
      m_proc_hi,
      m_n_ghosts,
      is_periodic,
      num_cells);
   m_momy.partition(m_pdim,
      m_pdim,
      m_proc_lo,
      m_proc_hi,
      m_n_ghosts,
      is_periodic,
      num_cells);
   m_ke.partition(m_pdim,
      m_pdim,
      m_proc_lo,
      m_proc_hi,
      m_n_ghosts,
      is_periodic,
      num_cells);
   m_ent.partition(m_pdim,
      m_pdim,
      m_proc_lo,
      m_proc_hi,
      m_n_ghosts,
      is_periodic,
      num_cells);

   // We want to integrate over velocity.
   deque<bool> collapse_dir(m_pdim, false);
   collapse_dir[V1] = true;
   collapse_dir[V2] = true;
   double measure(m_domain->dx(2) * m_domain->dx(3));

   // Now make the schedule for Momx.
   m_momx_schedule = new ReductionSchedule(m_momx,
      *m_domain,
      collapse_dir,
      m_n_ghosts,
      measure,
      m_mass,
      m_comm);

   // Now make the schedule for Momy.
   m_momy_schedule = new ReductionSchedule(m_momy,
      *m_domain,
      collapse_dir,
      m_n_ghosts,
      measure,
      m_mass,
      m_comm);

   // Now make the schedule for KE.
   m_ke_schedule = new ReductionSchedule(m_ke,
      *m_domain,
      collapse_dir,
      m_n_ghosts,
      measure,
      m_mass,
      m_comm);

   // Now make the schedule for the entropy.
   m_ent_schedule = new ReductionSchedule(m_ent,
      *m_domain,
      collapse_dir,
      m_n_ghosts,
      measure,
      m_mass,
      m_comm);
}


void
KineticSpecies::defineKineticEnergySummationSchedule()
{
   // Dimension this species kinetic energy velocity boundary flux.
   const ParallelArray::Box& interior_box = interiorBox();
   ParallelArray::Box base_space(2);
   vector<int> num_global_cells(2);
   for (int i = 0; i < 2; ++i) {
      base_space.lower(i) = interior_box.lower(i);
      base_space.upper(i) = interior_box.upper(i);
      num_global_cells[i] = m_domain->numberOfCells(i);
   }
   m_ke_vel_bdry_flux.partition(base_space, 2, m_n_ghosts, num_global_cells);

   // Now make the schedule for the KE velocity boundary flux.
   m_ke_vel_bdry_flux_schedule = new SummationSchedule(m_ke_vel_bdry_flux,
      m_dist_func,
      m_domain->box(),
      m_comm);
}


void
KineticSpecies::defineDiagnosticSummationSchedule()
{
   // Dimension this species' diagnostic fields.
   m_diagnostics.resize(CollisionOperator::s_DIAGNOSTIC_WORK_SIZE+1);

   const ParallelArray::Box& interior_box = interiorBox();
   ParallelArray::Box base_space(2);
   vector<int> num_global_cells(2);
   for (int i = 0; i < 2; ++i) {
      base_space.lower(i) = interior_box.lower(i);
      base_space.upper(i) = interior_box.upper(i);
      num_global_cells[i] = m_domain->numberOfCells(i);
   }
   for (int i = 0; i < CollisionOperator::s_DIAGNOSTIC_WORK_SIZE+1; ++i) {
      m_diagnostics[i].partition(base_space,
         2,
         m_n_ghosts,
         num_global_cells);
   }

   // Now make the schedule for the diagnostic field.
   m_diagnostic_schedule =
      new SummationSchedule(m_diagnostics[CollisionOperator::s_DIAGNOSTIC_WORK_SIZE],
         m_dist_func,
         m_domain->box(),
         m_comm);
}


void
KineticSpecies::defineCurrentDensityReductionSchedules()
{
   // Dimension this species current densities.
   deque<bool> is_periodic(m_pdim);
   vector<int> num_cells(m_pdim);
   for (int dim = 0; dim < m_pdim; ++dim) {
      is_periodic[dim] = m_domain->isPeriodic(dim);
      num_cells[dim] = m_domain->numberOfCells(dim);
   }
   m_Jx.partition(m_pdim,
      m_pdim,
      m_proc_lo,
      m_proc_hi,
      m_n_ghosts,
      is_periodic,
      num_cells);
   m_Jy.partition(m_pdim,
      m_pdim,
      m_proc_lo,
      m_proc_hi,
      m_n_ghosts,
      is_periodic,
      num_cells);
   m_Jz.partition(m_pdim,
      m_pdim,
      m_proc_lo,
      m_proc_hi,
      m_n_ghosts,
      is_periodic,
      num_cells);

   // Now make the schedules for Jx, Jy, and Jz.
   // We want to integrate over velocity.
   deque<bool> collapse_dir(m_pdim, false);
   collapse_dir[V1] = true;
   collapse_dir[V2] = true;
   double measure(m_domain->dx(2) * m_domain->dx(3));
   m_jx_schedule = new ReductionSchedule(m_Jx,
      *m_domain,
      collapse_dir,
      m_n_ghosts,
      measure,
      m_charge,
      m_comm);
   m_jy_schedule = new ReductionSchedule(m_Jy,
      *m_domain,
      collapse_dir,
      m_n_ghosts,
      measure,
      m_charge,
      m_comm);
   m_jz_schedule = new ReductionSchedule(m_Jz,
      *m_domain,
      collapse_dir,
      m_n_ghosts,
      measure,
      m_charge,
      m_comm);
}


void
KineticSpecies::buildVelocityArrays()
{
   int n3lo = dataBox().lower(V1);
   int n3hi = dataBox().upper(V1);
   int n4lo = dataBox().lower(V2);
   int n4hi = dataBox().upper(V2);
   ParallelArray::Box base_space(PDIM-CDIM+1);
   vector<int> num_global_cells(PDIM-CDIM);
   for (int i = CDIM; i < PDIM; ++i) {
      base_space.lower(i-CDIM) = interiorBox().lower(i);
      base_space.upper(i-CDIM) = interiorBox().upper(i);
      num_global_cells[i-CDIM] = m_domain->numberOfCells(i);
   }
   base_space.lower(PDIM-CDIM) = 0;
   base_space.upper(PDIM-CDIM) = 1;
   m_velocities.partition(base_space,
      PDIM-CDIM,
      m_n_ghosts,
      num_global_cells);
   base_space.upper(0) = interiorBox().upper(V1)+1;
   num_global_cells[0] = m_domain->numberOfCells(V1)+1;
   m_vxface_velocities.partition(base_space,
      PDIM-CDIM,
      m_n_ghosts,
      num_global_cells);
   base_space.upper(0) = interiorBox().upper(V1);
   base_space.upper(1) = interiorBox().upper(V2)+1;
   num_global_cells[0] = m_domain->numberOfCells(V1);
   num_global_cells[1] = m_domain->numberOfCells(V2)+1;
   m_vyface_velocities.partition(base_space,
      PDIM-CDIM,
      m_n_ghosts,
      num_global_cells);
   if (Simulation::s_DO_RELATIVITY) {
      double pxlo = m_domain->lower(V1);
      double pylo = m_domain->lower(V2);
      double dpx = m_domain->dx(V1);
      double dpy = m_domain->dx(V2);
      double vx, vy;
      for (int i3 = n3lo; i3 <= n3hi; ++i3) {
         double px = pxlo + (i3+0.5)*dpx;
         for (int i4 = n4lo; i4 <= n4hi; ++i4) {
            double py = pylo + (i4+0.5)*dpy;
            FORT_GET_VELOCITY(px, py, m_mass, Simulation::s_LIGHT_SPEED, vx, vy);
            m_velocities(i3, i4, 0) = vx;
            m_velocities(i3, i4, 1) = vy;
         }
      }
      for (int i3 = n3lo; i3 <= n3hi+1; ++i3) {
         double px = pxlo + i3*dpx;
         for (int i4 = n4lo; i4 <= n4hi; ++i4) {
            double py = pylo + (i4+0.5)*dpy;
            FORT_GET_VELOCITY(px, py, m_mass, Simulation::s_LIGHT_SPEED, vx, vy);
            m_vxface_velocities(i3, i4, 0) = vx;
            m_vxface_velocities(i3, i4, 1) = vy;
         }
      }
      for (int i3 = n3lo; i3 <= n3hi; ++i3) {
         double px = pxlo + (i3+0.5)*dpx;
         for (int i4 = n4lo; i4 <= n4hi+1; ++i4) {
            double py = pylo + i4*dpy;
            FORT_GET_VELOCITY(px, py, m_mass, Simulation::s_LIGHT_SPEED, vx, vy);
            m_vyface_velocities(i3, i4, 0) = vx;
            m_vyface_velocities(i3, i4, 1) = vy;
         }
      }
   }
   else {
      double vxlo = m_domain->lower(V1);
      double vylo = m_domain->lower(V2);
      double dvx = m_domain->dx(V1);
      double dvy = m_domain->dx(V2);
      for (int i3 = n3lo; i3 <= n3hi; ++i3) {
         double vx = vxlo + (i3+0.5)*dvx;
         for (int i4 = n4lo; i4 <= n4hi; ++i4) {
            m_velocities(i3, i4, 0) = vx;
            m_velocities(i3, i4, 1) = vylo + (i4+0.5)*dvy;
         }
      }
      for (int i3 = n3lo; i3 <= n3hi+1; ++i3) {
         double vx = vxlo + i3*dvx;
         for (int i4 = n4lo; i4 <= n4hi; ++i4) {
            m_vxface_velocities(i3, i4, 0) = vx;
            m_vxface_velocities(i3, i4, 1) = vylo + (i4+0.5)*dvy;
         }
      }
      for (int i3 = n3lo; i3 <= n3hi; ++i3) {
         double vx = vxlo + (i3+0.5)*dvx;
         for (int i4 = n4lo; i4 <= n4hi+1; ++i4) {
            m_vyface_velocities(i3, i4, 0) = vx;
            m_vyface_velocities(i3, i4, 1) = vylo + i4*dvy;
         }
      }
   }
}


void
KineticSpecies::accumulateSequencesCommon(
   vector<vector<double> >& a_sequences,
   int a_saved_seq,
   int& a_seq_idx,
   double a_time)
{
   // The species kinetic energy flux is only calculated on the processors over
   // which that species is distribututed.  So don't do irrelevent work.
   if (isInRange(Loki_Utilities::s_my_id)) {
      // Species kinetic energy flux at all phase space boundaries.
      // We need the current advection fluxes and hence the current physical
      // boundary conditions to do this calculation.
      setPhysicalBCs();
      evalAdvectionFluxes();
      // We also need the up to date acceleration boundary conditions and
      // acceleration fluxes to do this.
      setAccelerationBCs();
      evalAccelerationFluxes();
   }

   // Loop over the 4 boundaries of configuation space.
   const tbox::Box& domain_box = m_domain->box();
   for (int dir = X1; dir <= V2; ++dir) {
      for (int side = LO; side <= HI; ++side) {
         double ke_flux = 0.0;
         FORT_COMPUTE_KE_FLUX(BOX4D_TO_FORT(dataBox()),
            BOX4D_TO_FORT(interiorBox()),
            BOX4D_TO_FORT(domain_box),
            m_domain->dx()[0],
            *m_flux[X1].getData(),
            *m_flux[X2].getData(),
            *m_flux[V1].getData(),
            *m_flux[V2].getData(),
            *m_velocities.getData(),
            *m_vxface_velocities.getData(),
            *m_vyface_velocities.getData(),
            dir,
            side,
            m_mass,
            ke_flux);

         // Sum the local kinetic energy flux and add it to the time history.
         a_sequences[a_seq_idx++][a_saved_seq] =
            Loki_Utilities::getSum(ke_flux, -1);
      }
   }

   // Compute the rate of kinetic energy added to the species by any external
   // electric field drivers and the sum of the time dependent part of these
   // drivers.
   double ke_e_dot = 0.0;
   double Et = 0.0;
   if (m_num_external_drivers > 0) {
      // The driver_sums_into flag indicated what the driver sums the field
      // into.  A value of 1 indicates m_ext_efield only, 2 indicates both
      // m_ext_efield and m_em_vars/m_accel.
      int driver_sums_into = 1;
      m_ext_efield = 0.0;
      for (int i = 0; i < m_num_external_drivers; ++i) {
         if (m_do_maxwell) {
            m_ef_drivers[i]->evaluate(m_em_vars,
               m_ext_efield,
               *m_domain,
               driver_sums_into,
               a_time,
               0.0,
               0);
         }
         else {
            m_ef_drivers[i]->evaluate(m_accel,
               m_ext_efield,
               *m_domain,
               driver_sums_into,
               a_time,
               0.0,
               0);
         }
         double envel;
         m_ef_drivers[i]->evaluateTimeEnvelope(envel, a_time);
         Et += envel;
      }
      FORT_COMPUTE_KE_E_DOT(BOX4D_TO_FORT(dataBox()),
         BOX4D_TO_FORT(interiorBox()),
         PROBLEMDOMAIN_TO_FORT((*m_domain)),
         *m_dist_func.getData(),
         m_charge,
         *m_velocities.getData(),
         *m_ext_efield.getData(),
         ke_e_dot);
   }

   // Sum the local rate of kinetic energy added by external electric field
   // drivers, the integrated kinetic energy added by external electric field
   // drivers, and the time depdendent part of the drivers and add them to the
   // time history.
   a_sequences[a_seq_idx++][a_saved_seq] = Loki_Utilities::getSum(ke_e_dot, -1);
   a_sequences[a_seq_idx++][a_saved_seq] =
      Loki_Utilities::getSum(m_integrated_ke_e_dot, -1);
   a_sequences[a_seq_idx++][a_saved_seq] = Et;
}

} // end namespace Loki
