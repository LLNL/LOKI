/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "EMSolverBase.H"
#include "PoissonF.H"
#include "CollisionOperator.H"
#include "RestartManager.H"
#include "TimeHistWriter.H"

namespace Loki {

EMSolverBase::EMSolverBase(
   const tbox::Pointer<ProblemDomain>& a_domain,
   int a_num_kinetic_species,
   int a_solution_order,
   int a_num_em_vars,
   int a_plot_times_per_file,
   bool a_plot_ke_vel_bdy_flux,
   bool a_coll_diag_on,
   double a_bz_const)
   : m_dim(a_domain->dim()),
     m_domain(a_domain),
     m_number_of_procs(1),
     m_partition_defined(false),
     m_is_emsolver_processor(false),
     m_comm(MPI_COMM_NULL),
     m_num_em_vars(a_num_em_vars),
     m_loki_solver(0),
     m_ke_flux_vx_hi(0),
     m_ke_flux_vx_lo(0),
     m_ke_flux_vy_hi(0),
     m_ke_flux_vy_lo(0),
     m_diag_plot(0),
     m_solution_order(a_solution_order),
     m_reflect_noise_source_particles(false),
     m_bz_const(a_bz_const),
     m_num_kinetic_species(a_num_kinetic_species)
{
   if (a_plot_ke_vel_bdy_flux) {
      m_ke_flux_vx_hi.resize(a_num_kinetic_species);
      m_ke_flux_vx_lo.resize(a_num_kinetic_species);
      m_ke_flux_vy_hi.resize(a_num_kinetic_species);
      m_ke_flux_vy_lo.resize(a_num_kinetic_species);
   }

   if (a_coll_diag_on) {
      m_diag_plot.resize(2*CollisionOperator::s_DIAGNOSTIC_WORK_SIZE * a_num_kinetic_species);
   }

   if (m_dim != 2) {
      LOKI_ABORT("Maxwell only implemented for D=2!");
   }

   // Set number of ghosts based of order of solution.
   if (m_solution_order == 4) {
      m_n_ghosts = 2;
   }
   else {
      m_n_ghosts = 3;
   }

   // Build the field writers.
   RestartManager* restart_manager(RestartManager::getManager());
   m_field_writer = new FieldWriter(restart_manager->restartWritePath(),
      *m_domain,
      m_solution_order,
      a_plot_times_per_file);
   ostringstream coll_field_file_name;
   coll_field_file_name << restart_manager->restartWritePath() << "_coll";
   m_field_writer_coll = new FieldWriter(coll_field_file_name.str(),
      *m_domain,
      m_solution_order,
      a_plot_times_per_file);
}

EMSolverBase::EMSolverBase(
   const EMSolverBase& a_other)
   : m_dim(a_other.m_dim),
     m_phi_solver(a_other.m_phi_solver),
     m_num_kinetic_species(a_other.m_num_kinetic_species)
{
   if (this != &a_other) {
      m_domain = a_other.m_domain;
      m_number_of_procs = a_other.m_number_of_procs;
      m_proc_lo = a_other.m_proc_lo;
      m_proc_hi = a_other.m_proc_hi;
      m_partition_defined = a_other.m_partition_defined;
      m_phi_global = a_other.m_phi_global;
      m_rho_global = a_other.m_rho_global;
      m_is_emsolver_processor = a_other.m_is_emsolver_processor;
      m_comm = a_other.m_comm;
      m_em_vars = a_other.m_em_vars;
      m_num_em_vars = a_other.m_num_em_vars;
      m_solution_order = a_other.m_solution_order;
      m_n_ghosts = a_other.m_n_ghosts;
      m_loki_solver = a_other.m_loki_solver;
      m_tracking_particle_file = a_other.m_tracking_particle_file;
      m_problem_num_tracking_particles =
         a_other.m_problem_num_tracking_particles;
      m_tracking_particles = a_other.m_tracking_particles;
      m_noise_source_particle_file = a_other.m_noise_source_particle_file;
      m_problem_num_noise_source_particles =
         a_other.m_problem_num_noise_source_particles;
      m_noise_source_particles = a_other.m_noise_source_particles;
      m_reflect_noise_source_particles =
         a_other.m_reflect_noise_source_particles;
      m_ke_flux_vx_hi = a_other.m_ke_flux_vx_hi;
      m_ke_flux_vx_lo = a_other.m_ke_flux_vx_lo;
      m_ke_flux_vy_hi = a_other.m_ke_flux_vy_hi;
      m_ke_flux_vy_lo = a_other.m_ke_flux_vy_lo;
      m_bz_const = a_other.m_bz_const; //IEO
      for (int i=0; i<2; ++i) {
         m_vmin[i] = a_other.m_vmin[i];
         m_vmax[i] = a_other.m_vmax[i];
      }
      m_field_writer = a_other.m_field_writer;
      m_field_writer_coll = a_other.m_field_writer_coll;
   }
}

EMSolverBase::~EMSolverBase()
{
}

float
EMSolverBase::netCost() const
{
   // The computational cost of an EMSolverBase is proportional to the size of
   // its domain.
   float cost_per_cell(1.0);
   return cost_per_cell *
      static_cast<float>((m_domain->numberOfCells()).getProduct());
}


int
EMSolverBase::numberOfProcessors() const
{
   return m_number_of_procs;
}


bool
EMSolverBase::fixedNumberOfProcessors() const
{
   return true;
}


bool
EMSolverBase::isInRange(
   int a_proc_id) const
{
   // Returns true if the EMSolver calculation is partitioned onto this
   // processor.
   return ((a_proc_id >= m_proc_lo) && (a_proc_id <= m_proc_hi));
}


void
EMSolverBase::printDecomposition() const
{
   // This function is only valid if we actually know the decomposition.
   if (m_partition_defined) {
      // Print some basic decomposition info.
      Loki_Utilities::printF("  EM processor(s):  [%d,%d]\n",
         m_proc_lo,
         m_proc_hi);
   }
}


void
EMSolverBase::newAuxVariable(
   ParallelArray& a_var,
   int a_depth,
   int a_n_ghosts) const
{
   // Create any "auxilliary" varible that is partitioned like this EMSolver.
   int num_dims = a_depth == 1 ? m_dim : m_dim+1;
   deque<bool> is_periodic;
   vector<int> num_cells;
   for (int i = 0; i < m_dim; ++i) {
      is_periodic.push_back(m_domain->isPeriodic(i));
      num_cells.push_back(m_domain->numberOfCells(i));
   }
   for (int i = m_dim; i < num_dims; ++i) {
      num_cells.push_back(a_depth);
   }
   a_var.partition(num_dims,
      m_dim,
      m_proc_lo,
      m_proc_hi,
      a_n_ghosts,
      is_periodic,
      num_cells);
}


void
EMSolverBase::readParticles(
   LokiInputParser& a_pp,
   const double* a_vmin,
   const double* a_vmax)
{
   for (int i=0; i<2; ++i) {
      m_vmin[i] = a_vmin[i];
      m_vmax[i] = a_vmax[i];
   }

   // If tracking particles are specified, read them in.
   // If this is an EMSolver processor then actually read the tracking
   // particles, otherwise just read how many particles are in the problem.
   if (a_pp.contains("tracking_particle_file")) {
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("Tracking Particles");
      string tmp_name;
      a_pp.get("tracking_particle_file", tmp_name);
      m_tracking_particle_file = tmp_name;
      m_problem_num_tracking_particles =
         Particle::readParticleData(m_tracking_particle_file,
                                    m_tracking_particles,
                                    m_domain,
                                    m_vmin,
                                    m_vmax,
                                    isEMSolverProcessor(),
                                    false);
      timers->stopTimer("Tracking Particles");
   }
   else {
      m_problem_num_tracking_particles = 0;
   }

   // If noise source particles are specified, read them in.
   // If this is an EMSolver processor then actually read the noise source
   // particles, otherwise just read how many particles are in the problem.
   if (a_pp.contains("noise_source_particle_file")) {
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("Noisy Particles");
      string tmp_name;
      a_pp.get("noise_source_particle_file", tmp_name);
      m_noise_source_particle_file = tmp_name;
      m_problem_num_noise_source_particles =
         Particle::readParticleData(m_noise_source_particle_file,
                                    m_noise_source_particles,
                                    m_domain,
                                    m_vmin,
                                    m_vmax,
                                    isEMSolverProcessor(),
                                    true);
      string tmp("false");
      a_pp.query("reflect_noise_source_particles", tmp);
      m_reflect_noise_source_particles =
         tmp.compare("false") == 0 ? false : true;
      timers->stopTimer("Noisy Particles");
   }
   else {
      m_problem_num_noise_source_particles = 0;
   }
}


void
EMSolverBase::electricField(
   ParallelArray& a_charge_density,
   double a_time)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("Poisson");

   // Add any noise to the charge distribution.
   addNoiseToChargeDist(timers, a_charge_density, a_time);

   // Solve the system for phi.
   // For serial EM, this is straight forward.  For parallel EM there are more
   // hoops to jump through as we currently solve Poisson's equation serially on
   // each EM processor and then siphon off the local potential from the solved
   // serial potential.
   if (m_number_of_procs == 1) {
      // Modify the charge density to have a net neutral charge.
      neutralizeCharge4D(BOX2D_TO_FORT(a_charge_density.dataBox()),
         BOX2D_TO_FORT(a_charge_density.interiorBox()),
         *a_charge_density.getData(),
         m_comm);

      // Solve for the potential
      m_loki_solver->solve(m_phi_solver, a_charge_density);

      // Set boundary conditions of phi.
      startBCTimer(timers);

      m_phi_solver.communicatePeriodicBoundaries();
      stopBCTimer(timers);
   }
   else {
      // Fill in this EM processor's contribution to the serialized charge
      // density.  Then allreduce each processor's contributions so that each
      // EM processor ends up with the serialized charge density.
      m_rho_global = 0.0;
      const ParallelArray::Box& lb = a_charge_density.localBox();
      for (int i2 = lb.lower(X2); i2 <= lb.upper(X2); ++i2) {
         for (int i1 = lb.lower(X1); i1 <= lb.upper(X1); ++i1) {
            m_rho_global(i1, i2) = a_charge_density(i1, i2);
         }
      }
      MPI_Allreduce(MPI_IN_PLACE,
         m_rho_global.getData(),
         m_rho_global.dataBox().size(),
         MPI_DOUBLE,
         MPI_SUM,
         m_comm);

      // Modify the charge density to have a net neutral charge.
      neutralizeCharge4D(BOX2D_TO_FORT(m_rho_global.dataBox()),
         BOX2D_TO_FORT(m_rho_global.interiorBox()),
         *m_rho_global.getData(),
         m_comm);

      // Solve for the serialized potential.
      m_loki_solver->solve(m_phi_global, m_rho_global);

      // Set boundary conditions of phi.
      startBCTimer(timers);

      m_phi_global.communicatePeriodicBoundaries();
      stopBCTimer(timers);

      // Fill m_phi_solver with this processor's part of the serialized
      // potential.
      const ParallelArray::Box& db = m_phi_solver.dataBox();
      for (int i2 = db.lower(X2); i2 <= db.upper(X2); ++i2) {
         for (int i1 = db.lower(X1); i1 <= db.upper(X1); ++i1) {
            m_phi_solver(i1, i2) = m_phi_global(i1, i2);
         }
      }
   }

   // Add potential from any external drivers.
   addExternalPotential(m_phi_solver, a_time);

   // Now we have all contributions to the local potential.
   // Compute E = -grad(phi) on the interior of each EM processor.
   m_em_vars = 0.0;
   computeEFieldFromPotential(BOX2D_TO_FORT(m_phi_solver.dataBox()),
      BOX2D_TO_FORT(m_phi_solver.interiorBox()),
      m_solution_order,
      m_num_em_vars,
      m_domain->dx()[0],
      *m_em_vars.getData(),
      *m_phi_solver.getData());

   // If this is parallel Poisson each processor needs to get the E field's
   // ghost zones.
   if (m_number_of_procs > 1) {
      m_em_vars.communicateGhostData();
   }

   // fix periodicity of the plottable electric field components
   startBCTimer(timers);
   m_em_vars.communicatePeriodicBoundaries();

   stopBCTimer(timers);
   timers->stopTimer("Poisson");
}


void
EMSolverBase::createPartitionCommon(
   bool a_enforce_periodicity,
   int a_proc_lo,
   int a_proc_hi,
   const MPI_Comm& a_comm)
{
   m_comm = a_comm;
   m_proc_lo = a_proc_lo;
   m_proc_hi = a_proc_hi;
   m_partition_defined = true;

   // Partition the electromagnetic fields among its processors.
   if (isInRange(Loki_Utilities::s_my_id)) {
      m_is_emsolver_processor = true;
   }
   deque<bool> is_periodic(m_dim);
   vector<int> num_cells(m_dim+1);
   for (int dim = 0; dim < m_dim; ++dim) {
      if (a_enforce_periodicity) {
         is_periodic[dim] = true;
      }
      else {
         is_periodic[dim] = m_domain->isPeriodic(dim);
      }
      num_cells[dim] = m_domain->numberOfCells(dim);
   }
   num_cells[m_dim] = m_num_em_vars;
   m_em_vars.partition(m_dim+1,
      m_dim,
      a_proc_lo,
      a_proc_hi,
      m_n_ghosts,
      is_periodic,
      num_cells);
}

  
void
EMSolverBase::initializeCommon(
   int a_num_species,
   bool a_plot_ke_vel_bdy_flux,
   bool a_coll_diag_on)
{
   // Set up arrays of optional plot variables.
   for (int i = 0; i < a_num_species; ++i) {
      if (a_plot_ke_vel_bdy_flux) {
         // We want to plot each species' velocity boundary KE flux so make
         // arrays for each.
         newAuxVariable(m_ke_flux_vx_lo[i], 1);
         newAuxVariable(m_ke_flux_vx_hi[i], 1);
         newAuxVariable(m_ke_flux_vy_lo[i], 1);
         newAuxVariable(m_ke_flux_vy_hi[i], 1);
      }

      if (a_coll_diag_on) {
         int offset = 2*CollisionOperator::s_DIAGNOSTIC_WORK_SIZE*i;
         for (int plot = 0;
              plot < 2*CollisionOperator::s_DIAGNOSTIC_WORK_SIZE;
              ++plot) {
            newAuxVariable(m_diag_plot[offset+plot], 1);
         }
      }
   }
}


void
EMSolverBase::plotCommon(
   const vector<vector<double> >& a_sequences,
   const vector<double>& a_time_seq,
   int a_num_probes,
   int a_saved_seq,
   int& a_saved_save,
   string& a_time_hist_file_name,
   vector<string>& a_frame_names,
   vector<string>& a_time_hist_names)
{
   // Write all components of the EM field.
   // Get the pointer to the first EM field and the size of each field.
   const ParallelArray::Box& data_box = m_em_vars.dataBox();
   const ParallelArray::Box& local_box = m_em_vars.localBox();
   double* this_em_var = m_em_vars.getData();
   int em_var_size = data_box.size()/m_num_em_vars;
   for (int i = 0; i < m_num_em_vars; ++i) {
      // Write this component of the EM field to its frame.
      m_field_writer->writeField(a_frame_names[i],
         this_em_var,
         data_box,
         local_box,
         m_n_ghosts);
      this_em_var += em_var_size;
   }

   // Write the time histories.
   if (Loki_Utilities::s_my_id == m_proc_lo) {
      writeTimeHistories(a_sequences,
         a_time_hist_names,
         a_time_seq,
         a_num_probes,
         a_saved_seq,
         a_saved_save,
         a_time_hist_file_name);
   }
}


void
EMSolverBase::writeTimeHistories(
   const vector<vector<double> >& a_sequences,
   const vector<string>& a_name,
   const vector<double>& a_time_seq,
   int a_saved_seq,
   int& a_saved_save,
   string& a_time_hist_file_name)
{
   if (!isEMSolverProcessor()) {
      return;
   }

   // Form the name of this time history file and a writer for it.
   ostringstream time_hist_file;
   time_hist_file << a_time_hist_file_name << "_" << a_saved_save << ".hdf";
   TimeHistWriter time_hist_writer(time_hist_file.str(),
      a_saved_seq,
      a_time_seq);

   // Write each time history.
   int num_seq = static_cast<int>(a_name.size());
   for (int sequence = 0; sequence < num_seq; ++sequence) {
      time_hist_writer.writeTimeHistory(a_name[sequence],
         a_saved_seq,
         a_sequences[sequence]);
   }
}


void
EMSolverBase::writeTimeHistories(
   const vector<vector<double> >& a_sequences,
   const vector<string>& a_name,
   const vector<double>& a_time_seq,
   int a_num_probes,
   int a_saved_seq,
   int& a_saved_save,
   string& a_time_hist_file_name)
{
   // Form the name of this time history file and open it.
   ostringstream time_hist_file;
   time_hist_file << a_time_hist_file_name << "_" << a_saved_save << ".hdf";
   TimeHistWriter time_hist_writer(time_hist_file.str(),
      a_num_probes,
      numProblemTrackingParticles(),
      a_saved_seq,
      a_time_seq);

   // Write each time history.
   int num_seq = static_cast<int>(a_name.size());
   for (int sequence = 0; sequence < num_seq; ++sequence) {
      time_hist_writer.writeTimeHistory(a_name[sequence],
         a_saved_seq,
         a_sequences[sequence]);
   }
}


void
EMSolverBase::parseParametersCommon(
   LokiInputParser& a_pp)
{
   // See if the user has specified how many processors to use for the EM
   // solve.
   a_pp.query("number_of_processors", m_number_of_procs);
}

} // end namespace Loki
