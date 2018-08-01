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
#include "Poisson.H"
#include "PoissonF.H"

#include "MappedGridOperators.h"
#include "OgesParameters.h"
#include "TridiagonalSolver.h"
#include "PlotStuff.h"
#include "PlotStuffParameters.h"

#include "BoxOps.H"
#include "ElectricPotentialDriver.H"
#include "ElectricPotentialDriverFactory.H"
#include "Loki_Utilities.H"
#include "RestartManager.H"
#include "Interpolator4.H"
#include "Interpolator6.H"

#include "hdf5.h"

namespace Loki {

const int Poisson::TIME_HISTS_PER_PROBE = 2;
const int Poisson::GLOBAL_TIME_HISTS = 5;
int Poisson::NUM_FRAME_SERIES;

Poisson::Poisson(ParmParse& a_pp,
   const tbox::Pointer<ProblemDomain>& a_domain,
   int a_num_kinetic_species,
   int a_solution_order)
   : m_dim(a_domain->dim()),
     m_domain(a_domain),
     m_number_of_procs(1),
     m_partition_defined(false),
     m_n_ghosts(m_dim),
     m_is_poisson_processor(false),
     m_comm(MPI_COMM_NULL),
     m_solver(0),
     m_ke_flux_vx_hi(a_num_kinetic_species),
     m_ke_flux_vx_lo(a_num_kinetic_species),
     m_ke_flux_vy_hi(a_num_kinetic_species),
     m_ke_flux_vy_lo(a_num_kinetic_species),
#ifdef USE_SUPERLU
     m_solver_method(SUPERLU_DIRECT),
#else
     m_solver_method(OVERTURE_BEST_ITERATIVE),
#endif
     m_stencil_width(a_solution_order+1),
     m_max_iterations(10000),
     m_apply_external_potential(false),
     m_solution_order(a_solution_order),
     m_ex_interp(0),
     m_ey_interp(0),
     m_rho_interp(0)
{
   NUM_FRAME_SERIES = 2 + 4*a_num_kinetic_species;

   if (m_dim != 2) {
      OV_ABORT("Poisson only implemented for D=2!");
   }

   // Set number of ghosts based of order of solution.
   int num_ghosts;
   if (m_solution_order == 4) {
      num_ghosts = 2;
   }
   else {
      num_ghosts = 3;
   }
   for (int i = 0; i < m_dim; ++i) {
      m_n_ghosts[i] = num_ghosts;
   }

   // Set the boundary condition stencil size.
   m_stencil_size = int(pow(m_stencil_width, m_dim));

   // Read all input related to this object.
   parseParameters(a_pp);

   // The following code could probably be placed into parseParameters but
   // isn't because it involves parsing separate sub-databases of this Poisson.
   if (m_apply_external_potential) {
      // Get the sub-database for any external potential driver and create it.
      ParmParse epdf_pp("poisson.external_potential");
      m_ep_driver =
         ElectricPotentialDriverFactory::create(epdf_pp);
   }

   // The Poisson object writes restart data so register this object with the
   // RestartManager which will use the putToRestart/getFromRestart callbacks to
   // get the restart data written/read.
   RestartManager* restart_manager(RestartManager::getManager());
   restart_manager->registerRestart(this);

   printParameters();
}


Poisson::Poisson(
   const Poisson& a_other)
   : m_dim(a_other.m_dim),
     m_n_ghosts(m_dim),
     m_ex_interp(0),
     m_ey_interp(0),
     m_rho_interp(0)
{
   // If the 2 Poissons are different copy their internals.
   if (this != &a_other) {
      // DO NOT copy m_solver_composite_grid and m_composite_grid.  There are 2
      // reasons for this.  First it generates an inexplicable error from the
      // align method in bsparti.c.  The error makes no sense and should be
      // investigated further.  Secondly, it's simply not necessary to copy
      // them.  m_solver_composite_grid is set up, handed to m_phi_solver and
      // m_rho_solver and never used again.  Copying m_phi_solver and
      // m_rho_solver is sufficient.  Similarly, m_composite_grid is set up and
      // handed to m_phi, m_electric_field_*, and m_ke_flux_* and never used
      // again.  Copying these things is sufficient.
      m_domain = a_other.m_domain;
      m_number_of_procs = a_other.m_number_of_procs;
      m_processor_range = a_other.m_processor_range;
      m_partition_defined = a_other.m_partition_defined;
      m_n_ghosts = a_other.m_n_ghosts;
      m_is_poisson_processor = a_other.m_is_poisson_processor;
      m_comm = a_other.m_comm;
      m_partition = a_other.m_partition;
      m_solver = a_other.m_solver;
      m_rho_solver = a_other.m_rho_solver;
      m_phi_solver = a_other.m_phi_solver;
      m_phi = a_other.m_phi;
      m_electric_field_x = a_other.m_electric_field_x;
      m_electric_field_y = a_other.m_electric_field_y;
      m_ke_flux_vx_hi = a_other.m_ke_flux_vx_hi;
      m_ke_flux_vx_lo = a_other.m_ke_flux_vx_lo;
      m_ke_flux_vy_hi = a_other.m_ke_flux_vy_hi;
      m_ke_flux_vy_lo = a_other.m_ke_flux_vy_lo;
      m_solver_method = a_other.m_solver_method;
      m_tol = a_other.m_tol;
      m_stencil_size = a_other.m_stencil_size;
      m_stencil_width = a_other.m_stencil_width;
      m_max_iterations = a_other.m_max_iterations;
      m_ep_driver = a_other.m_ep_driver;
      m_apply_external_potential = a_other.m_apply_external_potential;
      m_solution_order = a_other.m_solution_order;
      m_problem_num_tracking_particles =
         a_other.m_problem_num_tracking_particles;
      m_tracking_particles = a_other.m_tracking_particles;
      m_problem_num_noise_source_particles =
         a_other.m_problem_num_noise_source_particles;
      m_noise_source_particles = a_other.m_noise_source_particles;
      buildInterpolators();
   }
}


Poisson::~Poisson()
{
   if (m_ex_interp) {
      for (int i = 0; i < m_solution_order; ++i) {
         delete [] m_ex_interp[i];
      }
      delete [] m_ex_interp;
   }
   if (m_ey_interp) {
      for (int i = 0; i < m_solution_order; ++i) {
         delete [] m_ey_interp[i];
      }
      delete [] m_ey_interp;
   }
   if (m_rho_interp) {
      for (int i = 0; i < m_solution_order; ++i) {
         delete [] m_rho_interp[i];
      }
      delete [] m_rho_interp;
   }
}


void
Poisson::initialize(
   const std::vector<string> a_species_names)
{
   // Build some of the mappings, composite grids, etc. needed for the Poisson
   // solve.
   defineSolver();

   // Now that we have the composite grid make the grid function containing the
   // potential.  This is used in the actual solve.
   m_phi_solver = realCompositeGridFunction(m_solver_composite_grid);
   m_phi_solver = 0.0;

   // Make another composite grid.
   Mapping* physical_mapping = createMapping(-1);
   MappedGrid physical_mg(*physical_mapping);
   setupMappedGrid(physical_mg, 0);
   m_composite_grid.add(physical_mg);
   m_composite_grid.updateReferences();
   m_composite_grid.update(MappedGrid::THEvertex | MappedGrid::THEcenter);

   // Make another grid function containing the potential.  This version of the
   // potential is copied from m_phi_solver.  Don't know why we need 2 versions
   // of these.
   m_phi = realCompositeGridFunction(m_composite_grid);
   m_phi.setName("physical function");
   m_phi.setName("potential", 0);
   m_phi = 0.0;

   // We want to plot the x and y components of the E field so make grid
   // functions for each.
   m_electric_field_x = realCompositeGridFunction(m_composite_grid);
   m_electric_field_x.setName("Ex");
   m_electric_field_x.setName("Ex", 0);
   m_electric_field_x = 0.0;

   m_electric_field_y = realCompositeGridFunction(m_composite_grid);
   m_electric_field_y.setName("Ey");
   m_electric_field_y.setName("Ey", 0);
   m_electric_field_y = 0.0;

   // We want to plot each species' velocity boundary KE flux so make grid
   // functions for each.
   int num_kinetic_species = static_cast<int>(a_species_names.size());
   for (int i = 0; i < num_kinetic_species; ++i) {
      char buffer[80];
      sprintf(buffer, "%s ke flux vx hi", a_species_names[i].c_str());
      m_ke_flux_vx_hi[i] = realCompositeGridFunction(m_composite_grid);
      m_ke_flux_vx_hi[i].setName(buffer);
      m_ke_flux_vx_hi[i].setName(buffer, 0);
      m_ke_flux_vx_hi[i] = 0.0;

      sprintf(buffer, "%s ke flux vx lo", a_species_names[i].c_str());
      m_ke_flux_vx_lo[i] = realCompositeGridFunction(m_composite_grid);
      m_ke_flux_vx_lo[i].setName(buffer);
      m_ke_flux_vx_lo[i].setName(buffer, 0);
      m_ke_flux_vx_lo[i] = 0.0;

      sprintf(buffer, "%s ke flux vy hi", a_species_names[i].c_str());
      m_ke_flux_vy_hi[i] = realCompositeGridFunction(m_composite_grid);
      m_ke_flux_vy_hi[i].setName(buffer);
      m_ke_flux_vy_hi[i].setName(buffer, 0);
      m_ke_flux_vy_hi[i] = 0.0;

      sprintf(buffer, "%s ke flux vy lo", a_species_names[i].c_str());
      m_ke_flux_vy_lo[i] = realCompositeGridFunction(m_composite_grid);
      m_ke_flux_vy_lo[i].setName(buffer);
      m_ke_flux_vy_lo[i].setName(buffer, 0);
      m_ke_flux_vy_lo[i] = 0.0;
   }

   // Now we know enough to build the appropriate field interpolators if needed.
   buildInterpolators();
}


void
Poisson::buildInterpolators()
{
   // If this object is on the Poisson processor and there are either tracking
   // or noise source particles then we need to construct the interpolator of
   // the appropriate order.
   if (isPoissonProcessor() &&
       (m_problem_num_tracking_particles > 0 ||
        m_problem_num_noise_source_particles > 0)) {
#ifdef USE_PPP
      RealArray Ex;
      getLocalArrayWithGhostBoundaries(m_electric_field_x[0], Ex);
      RealArray rho;
      getLocalArrayWithGhostBoundaries(m_rho_solver[0], rho);
#else
      RealArray& Ex = m_electric_field_x[0];
      RealArray& rho = m_rho_solver[0];
#endif
      const tbox::Box e_local_box(BoxOps::getLocalBox(Ex));
      const tbox::Box rho_local_box(BoxOps::getLocalBox(rho));
      if (m_solution_order == 4) {
         m_e_interpolator = new Interpolator4(e_local_box, *m_domain);
         m_rho_interpolator = new Interpolator4(rho_local_box, *m_domain);
      }
      else {
         m_e_interpolator = new Interpolator6(e_local_box, *m_domain);
         m_rho_interpolator = new Interpolator6(rho_local_box, *m_domain);
      }
      if (m_problem_num_tracking_particles > 0 ||
          m_problem_num_noise_source_particles > 0) {
         m_ex_interp = new double* [m_solution_order];
         m_ey_interp = new double* [m_solution_order];
         for (int i = 0; i < m_solution_order; ++i) {
            m_ex_interp[i] = new double [m_solution_order];
            m_ey_interp[i] = new double [m_solution_order];
         }
      }
      if (m_problem_num_noise_source_particles > 0) {
         m_rho_interp = new double* [m_solution_order];
         for (int i = 0; i < m_solution_order; ++i) {
            m_rho_interp[i] = new double [m_solution_order];
         }
      }
   }
}


void
Poisson::computeFactorization()
{
   // Do an initial solve to obtain the static factorization.
   m_solver->solve(m_phi_solver, m_rho_solver);
   m_solver->setRefactor(FALSE);
   m_solver->setReorder(FALSE);
}


bool
Poisson::isInRange(
   int a_proc_id) const
{
   // Returns true if the Poisson calculation is partitioned onto this
   // processor.
   return ((a_proc_id >= m_processor_range.getBase()) &&
           (a_proc_id <= m_processor_range.getBound()));
}


void
Poisson::printDecomposition() const
{
   // This function is only valid if we actually know the decomposition.
   if (m_partition_defined) {
      // Print some basic decomposition info.
      printF("  Poisson processor(s):  [%d,%d]\n",
             m_processor_range.getBase(),
             m_processor_range.getBound());
   }
}


bool
Poisson::conformsTo(
   const Poisson& a_other) const
{
   // If the 2 Poissons are different, check that they are defined on the same
   // piece of space and are partitioned the same.
   if (m_domain != a_other.m_domain ||
       m_processor_range != a_other.m_processor_range ||
       m_partition_defined != a_other.m_partition_defined ||
       numTrackingParticles() != a_other.numTrackingParticles() ||
       numNoiseSourceParticles() != a_other.numNoiseSourceParticles()) {
      return false;
   }
   else {
      return true;
   }
}


void
Poisson::addData(
   const Poisson& a_increment_poisson,
   real a_factor,
   bool a_final_rk)
{
   // Check that we're copying between similar Poissons and to each particle's
   // x, y, vx, and vy add a_factor times x, y, vx, and vy of the corresponding
   // particle in a_increment_poisson.
   if (conformsTo(a_increment_poisson)) {
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("Tracking Particles");
      addParticleData(numTrackingParticles(),
         m_tracking_particles,
         a_increment_poisson.m_tracking_particles,
         a_factor,
         a_final_rk);
      timers->stopTimer("Tracking Particles");
      timers->startTimer("Noisy Particles");
      addParticleData(numNoiseSourceParticles(),
         m_noise_source_particles,
         a_increment_poisson.m_noise_source_particles,
         a_factor,
         a_final_rk);
      timers->stopTimer("Noisy Particles");
   }
}


void
Poisson::addParticleData(
   int a_num_particles,
   std::vector<Particle>& a_particles,
   const std::vector<Particle>& a_increment_particles,
   real a_factor,
   bool a_final_rk)
{
   for (int i = 0; i < a_num_particles; ++i) {
      a_particles[i].x() += a_increment_particles[i].x()*a_factor;
      a_particles[i].y() += a_increment_particles[i].y()*a_factor;
      a_particles[i].vx() += a_increment_particles[i].vx()*a_factor;
      a_particles[i].vy() += a_increment_particles[i].vy()*a_factor;
      // This is the final RK update so the particle positions need to be
      // limited.  If a particle is outside the physical domain in a periodic
      // direction then wrap it back around to its periodic location.  If a
      // particle is outside the physical domain in a non-periodic direction
      // then limit its position to the edge of the physical domain in that
      // direction.
      if (a_final_rk) {
         if (a_particles[i].x() < m_domain->lower(X1)) {
            if (m_domain->isPeriodic(X1)) {
               a_particles[i].x() += m_domain->upper(X1) - m_domain->lower(X1);
            }
            else {
               a_particles[i].x() = m_domain->lower(X1);
            }
         }
         else if (a_particles[i].x() > m_domain->upper(X1)) {
            if (m_domain->isPeriodic(X1)) {
               a_particles[i].x() += m_domain->lower(X1) - m_domain->upper(X1);
            }
            else {
               a_particles[i].x() = m_domain->upper(X1);
            }
         }
         if (a_particles[i].y() < m_domain->lower(X2)) {
            if (m_domain->isPeriodic(X2)) {
               a_particles[i].y() += m_domain->upper(X2) - m_domain->lower(X2);
            }
            else {
               a_particles[i].y() = m_domain->lower(X2);
            }
         }
         else if (a_particles[i].y() > m_domain->upper(X2)) {
            if (m_domain->isPeriodic(X2)) {
               a_particles[i].y() += m_domain->lower(X2) - m_domain->upper(X2);
            }
            else {
               a_particles[i].y() = m_domain->upper(X2);
            }
         }
         if (a_particles[i].vx() < m_domain->lower(V1)) {
            a_particles[i].vx() = m_domain->lower(V1);
         }
         else if (a_particles[i].vx() > m_domain->upper(V1)) {
            a_particles[i].vx() = m_domain->upper(V1);
         }
         if (a_particles[i].vy() < m_domain->lower(V2)) {
            a_particles[i].vy() = m_domain->lower(V2);
         }
         else if (a_particles[i].vy() > m_domain->upper(V2)) {
            a_particles[i].vy() = m_domain->upper(V2);
         }
      }
   }
}


void
Poisson::copySolnData(
   const Poisson& a_rhs)
{
   // If the 2 Poissons are different copy each tracking particle's x, y, vx,
   // and vy.
   if (&a_rhs != this) {
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("Tracking Particles");
      for (int i = 0; i < numTrackingParticles(); ++i) {
         m_tracking_particles[i].x() = a_rhs.m_tracking_particles[i].x();
         m_tracking_particles[i].y() = a_rhs.m_tracking_particles[i].y();
         m_tracking_particles[i].vx() = a_rhs.m_tracking_particles[i].vx();
         m_tracking_particles[i].vy() = a_rhs.m_tracking_particles[i].vy();
      }
      timers->stopTimer("Tracking Particles");
      timers->startTimer("Noisy Particles");
      for (int i = 0; i < numNoiseSourceParticles(); ++i) {
         m_noise_source_particles[i].x() =
            a_rhs.m_noise_source_particles[i].x();
         m_noise_source_particles[i].y() =
            a_rhs.m_noise_source_particles[i].y();
         m_noise_source_particles[i].vx() =
            a_rhs.m_noise_source_particles[i].vx();
         m_noise_source_particles[i].vy() =
            a_rhs.m_noise_source_particles[i].vy();
      }
      timers->stopTimer("Noisy Particles");
   }
}


void
Poisson::defineSolver()
{
   Mapping* poisson_mapping = createMapping(1);

   MappedGrid::setMinimumNumberOfDistributedGhostLines(m_n_ghosts.max());
   MappedGrid poisson_mg(*poisson_mapping);
   setupMappedGrid(poisson_mg, m_stencil_width);

   // Set solver parameters given requested solver type and other user input.
   setupSolver(poisson_mg);
}


void
Poisson::parseParameters(
   ParmParse& a_pp)
{
   // See if there is an external potential.
   aString potential_str("false");
   a_pp.query("apply_external_potential", potential_str);
   m_apply_external_potential = potential_str.matches("false") ? false : true;

   // See if the user has specified how many processors to use for the Poisson
   // solve.  This shouldn't currently be used because we don't have a scalable
   // Poisson solver and if we did, there's plenty of assumptions that would
   // probably break.
   a_pp.query("number_of_processors", m_number_of_procs);

   // See if the user wants a specific solver method and that it's valid.
   if (a_pp.contains("solver_method")) {
      aString solver_method;
      a_pp.get("solver_method", solver_method);
      if (solver_method.matches("original iterative")) {
         m_solver_method = ORIGINAL_ITERATIVE;
      }
      else if (solver_method.matches("overture best iterative")) {
         m_solver_method = OVERTURE_BEST_ITERATIVE;
      }
      else if (solver_method.matches("superlu direct")) {
#ifndef USE_SUPERLU
         OV_ABORT("Loki not compiled to use superlu solvers.");
#endif
         m_solver_method = SUPERLU_DIRECT;
      }
      else {
         aString msg("Unknown solver_method \"" + solver_method + "\" ... quitting");
         OV_ABORT(msg);
      }
   }

   // Get the requested tolerance.  If none, then set default according to
   // solver method.
   if (a_pp.contains("solver_tolerance")) {
      a_pp.get("solver_tolerance", m_tol);
   }
   else if (m_solver_method == ORIGINAL_ITERATIVE) {
      m_tol = 1.0e-6;
   }
   else if (m_solver_method == OVERTURE_BEST_ITERATIVE) {
      m_tol = 1.0e-12;
   }
   else {
      m_tol = -1.0;
   }
}


void
Poisson::readParticles(
   ParmParse& a_pp)
{
   // If tracking particles are specified, read them in.
   // If this is a Poisson processor then actually read the tracking particles,
   // otherwise just read how many particles are in the problem.
   if (a_pp.contains("tracking_particle_file")) {
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("Tracking Particles");
      aString tracking_particle_file;
      a_pp.get("tracking_particle_file", tracking_particle_file);
      m_problem_num_tracking_particles =
         readParticleData(tracking_particle_file, m_tracking_particles, false);
      timers->stopTimer("Tracking Particles");
   }
   else {
      m_problem_num_tracking_particles = 0;
   }

   // If noise source particles are specified, read them in.
   // If this is a Poisson processor then actually read the noise source
   // particles, otherwise just read how many particles are in the problem.
   if (a_pp.contains("noise_source_particle_file")) {
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("Noisy Particles");
      aString noise_source_particle_file;
      a_pp.get("noise_source_particle_file", noise_source_particle_file);
      m_problem_num_noise_source_particles =
         readParticleData(noise_source_particle_file, m_noise_source_particles, true);
      timers->stopTimer("Noisy Particles");
   }
   else {
      m_problem_num_noise_source_particles = 0;
   }
}


void
Poisson::printParameters() const
{
   // Print this Poisson's parameters and the parameters of the entities that it
   // holds and are not accessible elsewhere like the external potential driver.
   printF("\n#*#*# Poisson #*#*#\n");
   printF("  Number of Processsors = %i\n", m_number_of_procs);
   printF("  apply external potential = %s\n",
      m_apply_external_potential ? "true" : "false");
   if (!m_ep_driver.isNull() && m_apply_external_potential) {
      m_ep_driver->printParameters();
   }
}


float
Poisson::netCost() const
{
   // The computational cost of a Poisson is proportional to the size of its
   // domain.
   float cost_per_cell(1.0);
   return cost_per_cell * static_cast<float>((m_domain->numberOfCells()).getProduct());
}


int
Poisson::numberOfProcessors() const
{
   return m_number_of_procs;
}


bool
Poisson::fixedNumberOfProcessors() const
{
   return true;
}


void
Poisson::createPartition(
   const Range& a_range,
   const MPI_Comm& a_comm)
{
   m_comm = a_comm;
   m_processor_range = a_range;
   m_partition_defined = true;

   const int my_id(std::max(0, Communication_Manager::My_Process_Number));
   if (isInRange(my_id)) {
      m_is_poisson_processor = true;
   }
   
   // Overture is compiled such that each parallel array has MAX_ARRAY_DIMENSION
   // dimensions.  So we need to tell the partition how to partition all of
   // these dimensions.  Naturally we only care about the first m_dim of them.
   // The others are ignored.
   m_partition.SpecifyProcessorRange(m_processor_range);
   m_partition.SpecifyDecompositionAxes(m_dim);

   for (int dir(0); dir < m_dim; ++dir) {
      m_partition.partitionAlongAxis(dir, true, m_n_ghosts[dir]);
   }
   for (int dir(m_dim); dir < MAX_ARRAY_DIMENSION; ++dir) {
      m_partition.partitionAlongAxis(dir, false, 0);
   }
}


void
Poisson::newAuxVariable(
   realArray& a_var,
   int a_depth,
   const tbox::IntVector& a_n_ghosts) const
{
   // Create any "auxilliary" varible that is partitioned like this Poisson.
   if (m_partition_defined) {

      tbox::Box box = m_domain->box();
      box.grow(a_n_ghosts);

      if (m_dim == tbox::Dimension(2)) {
         a_var.redim(BoxOps::range(box, X1),
            BoxOps::range(box, X2),
            Range(0, a_depth - 1));
      }
      else {
         OV_ABORT("Not implemented for phase D!=2!");
      }
      a_var.partition(m_partition);
      RealArray a_var_local;
      getLocalArrayWithGhostBoundaries(a_var, a_var_local);
      a_var_local = 0.0;
   }
   else {
      OV_ABORT("Attempt to define variables without partition!");
   }
}


void
Poisson::evalRHS(
   Poisson& a_rhs,
   const realArray& a_net_ext_efield,
   real a_time)
{
   // The RHS evaluation only applies to the evaluation of the equations of
   // motion of the particles.
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("Poisson RHS");

   timers->startTimer("Tracking Particles");
   int num_particles = numTrackingParticles();
   if (num_particles > 0) {
      evalRHSParticles(num_particles,
         m_tracking_particles,
         a_rhs.m_tracking_particles,
         a_net_ext_efield,
         a_time);
   }
   timers->stopTimer("Tracking Particles");

   timers->startTimer("Noisy Particles");
   num_particles = numNoiseSourceParticles();
   if (num_particles > 0) {
      evalRHSParticles(num_particles,
         m_noise_source_particles,
         a_rhs.m_noise_source_particles,
         a_net_ext_efield,
         a_time);
   }
   timers->stopTimer("Noisy Particles");

   timers->stopTimer("Poisson RHS");
}


void
Poisson::evalRHSParticles(
   int a_num_particles,
   const std::vector<Particle>& a_particles,
   std::vector<Particle>& a_rhs_particles,
   const realArray& a_net_ext_efield,
   real a_time)
{
   // Get the local x and y components of the self consistent E field as well
   // as the net external applied E field.
#ifdef USE_PPP
   RealArray Ex;
   getLocalArrayWithGhostBoundaries(m_electric_field_x[0], Ex);
   RealArray Ey;
   getLocalArrayWithGhostBoundaries(m_electric_field_y[0], Ey);
   RealArray E_ext;
   getLocalArrayWithGhostBoundaries(a_net_ext_efield, E_ext);
#else
   RealArray& Ex  = m_electric_field_x[0];
   RealArray& Ey  = m_electric_field_y[0];
   RealArray& E_ext = a_net_ext_efield;
#endif
   const double* Ex_ptr = Ex.getDataPointer();
   const double* Ey_ptr = Ey.getDataPointer();
   const double* E_ext_ptr = E_ext.getDataPointer();

   // For each particle use an interpolating polynomial to compute the E field
   // that it sees and update x, y, vx, and vy.
   for (int i = 0; i < a_num_particles; ++i) {
      const Particle& this_particle = a_particles[i];
      if (a_time >= this_particle.startingTime()) {
         // Set the evaluation point to be the location of this particle.
         m_e_interpolator->setEvaluationPoint(this_particle.x(),
                                            this_particle.y());

         // Figure out Ex and Ey at each interpolating point.
         int nx = m_e_interpolator->nx();
         int nxny = m_e_interpolator->nxny();
         for (int xidx = 0; xidx < m_solution_order; ++xidx) {
            int xGridIdx = m_e_interpolator->xInterpolatingGridIndex(xidx);
            for (int yidx = 0; yidx < m_solution_order; ++yidx) {
               int idx1 =
                  m_e_interpolator->yInterpolatingGridIndex(yidx)*nx + xGridIdx;
               int idx2 = idx1 + nxny;
               m_ex_interp[xidx][yidx] = Ex_ptr[idx1] + E_ext_ptr[idx1];
               m_ey_interp[xidx][yidx] = Ey_ptr[idx1] + E_ext_ptr[idx2];
            }
         }

         // Now compute Ex and Ey for this particle.
         double Ex_particle = m_e_interpolator->interpolate(m_ex_interp);
         double Ey_particle = m_e_interpolator->interpolate(m_ey_interp);

         // The new x/y is the old vx/vy.  The new vx/vy is E*q/m.
         double qoverm = this_particle.charge()/this_particle.mass();
         a_rhs_particles[i].x() = this_particle.vx();
         a_rhs_particles[i].y() = this_particle.vy();
         a_rhs_particles[i].vx() = Ex_particle*qoverm;
         a_rhs_particles[i].vy() = Ey_particle*qoverm;
      }
   }
}


void
Poisson::electricField(
   realArray& a_efield,
   const realArray& a_charge_density,
   real a_time)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("Poisson");

   const tbox::Box& domain_box(m_domain->box());

   // FIXME: a_charge_density and m_rho_solver should be defined on the
   // same set of configuration space processors, i.e., we should have to
   // use a parallel copy here
   // In any event, the solver's global representation of rho needs to get its
   // data from a_charge_density.
   Index dest[4], src[4];
   dest[0] = Range(BoxOps::range(domain_box, X1));
   dest[1] = Range(BoxOps::range(domain_box, X2));
   dest[2] = Range(0, 0);
   dest[3] = Range(0, 0);
   src[0] = Range(BoxOps::range(domain_box, X1));
   src[1] = Range(BoxOps::range(domain_box, X2));
   src[2] = Range(0, 0);
   src[3] = Range(0, 0);
   ParallelUtility::copy(m_rho_solver[0], dest, a_charge_density, src, 4);

   // Now that the solver's global representation of rho is set, we can get the
   // appropriate local array.
#ifdef USE_PPP
   RealArray rhs_local;
   getLocalArrayWithGhostBoundaries(m_rho_solver[0], rhs_local);
#else
   RealArray& rhs_local = m_rho_solver[0];
#endif

   // If there are noise source particles, add the contribution of each to the
   // net charge density.
   timers->startTimer("Noisy Particles");
   int num_noise_source_particles = numNoiseSourceParticles();
   double* rhs_local_ptr = rhs_local.getDataPointer();
   for (int i = 0; i < num_noise_source_particles; ++i) {
      // Set the evaluation point to be the location of this particle.
      Particle& this_particle = m_noise_source_particles[i];
      m_rho_interpolator->setEvaluationPoint(this_particle.x(),
                                             this_particle.y());

      // Now compute the charge at the interpolating points that is required
      // to give the charge of this noise source particle.
      int nx = m_rho_interpolator->nx();
      double noise_charge =
         this_particle.charge()*this_particle.noiseSourceWeight();
      m_rho_interpolator->deterpolate(noise_charge, m_rho_interp);
      for (int xidx = 0; xidx < 4; ++xidx) {
         int xGridIdx = m_rho_interpolator->xInterpolatingGridIndex(xidx);
         for (int yidx = 0; yidx < 4; ++yidx) {
            int rhs_idx =
               m_rho_interpolator->yInterpolatingGridIndex(yidx)*nx+xGridIdx;
            rhs_local_ptr[rhs_idx] += m_rho_interp[xidx][yidx];
         }
      }
   }
   timers->stopTimer("Noisy Particles");

   // Modify the charge density to have a net neutral charge.
   const tbox::Box rhs_box(BoxOps::getLocalBox(rhs_local));
   const tbox::Box global_box(tbox::Box::grow(domain_box, m_n_ghosts));
   const tbox::Box interior_box = BoxOps::getOperationalBox(rhs_local,
      m_rho_solver[0],
      domain_box,
      global_box);
   neutralizeCharge4D(BOX2D_TO_FORT(rhs_box),
      BOX2D_TO_FORT(interior_box),
      *rhs_local_ptr,
      m_comm);
   // IS THIS NECESSARY FOR PARALLEL?
   m_rho_solver[0].updateGhostBoundaries();

   // Solve the system for phi.
   m_solver->solve(m_phi_solver, m_rho_solver);
   m_solver->setRefactor(FALSE);
   m_solver->setReorder(FALSE);

   timers->startTimer("BC (Poisson)");
   // copy back to physical domain
   ParallelUtility::copy(m_phi[0], dest, m_phi_solver[0], src, 4);

   Loki_Utilities::fixPeriodicity(m_phi[0], *m_domain, m_n_ghosts, 1);
   timers->stopTimer("BC (Poisson)");

   // this needs to be fixed once we distribute the Poisson solve?
   //   be particularly careful with the external field so that positions
   //   are computed correctly
   // In any event, we need to do perform a local calculation of Ex ane Ey from
   // phi.
#ifdef USE_PPP
   RealArray Ex_local;
   getLocalArrayWithGhostBoundaries(m_electric_field_x[0], Ex_local);

   RealArray Ey_local;
   getLocalArrayWithGhostBoundaries(m_electric_field_y[0], Ey_local);

   RealArray phi_local;
   getLocalArrayWithGhostBoundaries(m_phi[0], phi_local);
#else
   RealArray& Ex_local  = m_electric_field_x[0];
   RealArray& Ey_local  = m_electric_field_y[0];
   RealArray& phi_local = m_phi[0];
#endif

   // Get externally defined phi and add to phi_local.
   tbox::Box fill_box(BoxOps::getLocalBox(phi_local));
   if (m_ep_driver && m_apply_external_potential) {
      Range Rx(phi_local.getBase(0), phi_local.getBound(0));
      Range Ry(phi_local.getBase(1), phi_local.getBound(1));
      RealArray phi_extern(Rx, Ry);
      phi_extern = 0.0;
      m_ep_driver->evaluate(phi_extern, fill_box, *m_domain, a_time);
      phi_local += phi_extern;
   }

   // Now we have all contributions to phi.  Compute E = -grad(phi).
   computeEFieldFromPotential(BOX2D_TO_FORT(fill_box),
      BOX2D_TO_FORT(interior_box),
      *(m_domain->dx()).dataPtr(),
      *Ex_local.getDataPointer(),
      *Ey_local.getDataPointer(),
      *phi_local.getDataPointer(),
      m_solution_order);

   // fix periodicity of the plottable electric field components
   timers->startTimer("BC (Poisson)");
   Loki_Utilities::fixPeriodicity(m_electric_field_x[0],
      *m_domain,
      m_n_ghosts,
      1);
   Loki_Utilities::fixPeriodicity(m_electric_field_y[0],
      *m_domain,
      m_n_ghosts,
      1);

   // Now copy the plottable E field components to a_efield, a plain old
   // parallel array.
   Index all;
   dest[0] = all;
   dest[1] = all;
   dest[2] = Range(0, 0);
   dest[3] = Range(0, 0);
   src[0]  = all;
   src[1]  = all;
   src[2]  = Range(0, 0);
   src[3]  = Range(0, 0);
   ParallelUtility::copy(a_efield, dest, m_electric_field_x[0], src, 4);
   dest[2] = Range(1, 1);
   ParallelUtility::copy(a_efield, dest, m_electric_field_y[0], src, 4);
   timers->stopTimer("BC (Poisson)");
   timers->stopTimer("Poisson");
}


void
Poisson::copyPlotFields(
   Poisson& a_other)
{
   // Hack for 6th order RK.  Read description of this function in the header.
   a_other.m_electric_field_x = m_electric_field_x;
   a_other.m_electric_field_y = m_electric_field_y;
}


void
Poisson::plot(
   real a_time,
   real a_dt,
   const RealArray& a_sequences,
   const std::vector<string>& a_species_names,
   const RealArray& a_time_seq,
   const RealArray& a_probes, 
   int a_num_probes,
   int a_num_seq,
   int a_saved_seq,
   int& a_saved_save,
   Ogshow& a_show)
{
   const tbox::IntVector num_cells(m_domain->numberOfCells());
   int num_species = static_cast<int>(a_species_names.size());
   int num_tracking_particles = numProblemTrackingParticles();
   char buffer[80];

   {
      // Set up the frame for Ex.
      a_show.setCurrentFrameSeries("Ex");
      a_show.startFrame();
      HDF_DataBase *dbp(a_show.getFrame());
      if (dbp == 0) {
         OV_ABORT("Can not get frame from a_show");
      }

      // Write some general plot info for the frame.
      dbp->put(a_time, "time");
      dbp->put(a_dt, "dt");
      dbp->put(a_num_probes, "numProbes");
      dbp->put(num_tracking_particles, "numTrackingParticles");

      // Fill in the time history names and probe positions.
      aString name[a_num_seq];
      int idx = 0;
      name[idx++] = "E_max";
      name[idx++] = "norm E";
      name[idx++] = "Ex_max";
      name[idx++] = "Ey_max";
      name[idx++] = "field_energy";
      for (int ip(0); ip < a_num_probes; ++ip) {
         sprintf(buffer, "ix_pos%i", ip);
         dbp->put(int(floor(a_probes(0, ip) * num_cells[X1])), buffer);
         name[idx++] = buffer;

         sprintf(buffer, "iy_pos%i", ip);
         dbp->put(int(floor(a_probes(1, ip) * num_cells[X2])), buffer);
         name[idx++] = buffer;
      }
      for (int ipart = 0; ipart < num_tracking_particles; ++ipart) {
         sprintf(buffer, "particle%i_x", ipart);
         name[idx++] = buffer;
         sprintf(buffer, "particle%i_y", ipart);
         name[idx++] = buffer;
         sprintf(buffer, "particle%i_vx", ipart);
         name[idx++] = buffer;
         sprintf(buffer, "particle%i_vy", ipart);
         name[idx++] = buffer;
      }
      for (int is = 0; is < num_species; ++is) {
         aString tmp(a_species_names[is] + "_ke");
         name[idx++] = tmp;
         tmp = a_species_names[is] + "_xlo_flux";
         name[idx++] = tmp;
         tmp = a_species_names[is] + "_xhi_flux";
         name[idx++] = tmp;
         tmp = a_species_names[is] + "_ylo_flux";
         name[idx++] = tmp;
         tmp = a_species_names[is] + "_yhi_flux";
         name[idx++] = tmp;
      }

      // Write a general comment, the time and time step for the frame.
      a_show.saveComment(0, sPrintF(buffer, "VlasovPoisson4D"));
      a_show.saveComment(1,
         sPrintF(buffer, "t=%5.2f, dt=%5.2e", a_time, a_dt));

      // Write the Ex plot.
      a_show.saveSolution(m_electric_field_x);

      // Write the time histories.
      Range N(0, a_saved_seq-1), S = a_num_seq;
      a_show.saveSequence("time histories",
         a_time_seq(N),
         a_sequences(N, S),
         name);
      a_show.endFrame();
   }

   {
      // Set up the frame for Ey.
      a_show.setCurrentFrameSeries("Ey");
      a_show.startFrame();
      HDF_DataBase *dbp(a_show.getFrame());
      if (dbp == 0) {
         OV_ABORT("Can not get frame from a_show");
      }

      // Write some general plot info for the frame.
      dbp->put(a_time, "time");
      dbp->put(a_dt, "dt");

      // Write a general comment, the time and time step for the frame.
      a_show.saveComment(0, sPrintF(buffer, "VlasovPoisson4D"));
      a_show.saveComment(1,
         sPrintF(buffer, "t=%5.2f, dt=%5.2e", a_time, a_dt));

      // Write the Ey plot.
      a_show.saveSolution(m_electric_field_y);
      a_show.endFrame();
   }

   // Write the plot data for each species following the same pattern as above.
   for (int is = 0; is < num_species; ++is) {
      // Vx lo KE flux
      sprintf(buffer, "%s vx lo ke flux", a_species_names[is].c_str());
      a_show.setCurrentFrameSeries(buffer);
      a_show.startFrame();
      HDF_DataBase *dbp(a_show.getFrame());
      if (dbp == 0) {
         OV_ABORT("Can not get frame from a_show");
      }

      dbp->put(a_time, "time");
      dbp->put(a_dt, "dt");
      a_show.saveComment(0, sPrintF(buffer, "VlasovPoisson4D"));
      a_show.saveComment(1,
         sPrintF(buffer, "t=%5.2f, dt=%5.2e", a_time, a_dt));
      a_show.saveSolution(m_ke_flux_vx_lo[is]);
      a_show.endFrame();

      // Vx hi KE flux
      sprintf(buffer, "%s vx hi ke flux", a_species_names[is].c_str());
      a_show.setCurrentFrameSeries(buffer);
      a_show.startFrame();
      dbp = a_show.getFrame();
      if (dbp == 0) {
         OV_ABORT("Can not get frame from a_show");
      }

      dbp->put(a_time, "time");
      dbp->put(a_dt, "dt");
      a_show.saveComment(0, sPrintF(buffer, "VlasovPoisson4D"));
      a_show.saveComment(1,
         sPrintF(buffer, "t=%5.2f, dt=%5.2e", a_time, a_dt));
      a_show.saveSolution(m_ke_flux_vx_hi[is]);
      a_show.endFrame();

      // Vy lo KE flux
      sprintf(buffer, "%s vy lo ke flux", a_species_names[is].c_str());
      a_show.setCurrentFrameSeries(buffer);
      a_show.startFrame();
      dbp = a_show.getFrame();
      if (dbp == 0) {
         OV_ABORT("Can not get frame from a_show");
      }

      dbp->put(a_time, "time");
      dbp->put(a_dt, "dt");
      a_show.saveComment(0, sPrintF(buffer, "VlasovPoisson4D"));
      a_show.saveComment(1,
         sPrintF(buffer, "t=%5.2f, dt=%5.2e", a_time, a_dt));
      a_show.saveSolution(m_ke_flux_vy_lo[is]);
      a_show.endFrame();

      // Vy hi KE flux
      sprintf(buffer, "%s vy hi ke flux", a_species_names[is].c_str());
      a_show.setCurrentFrameSeries(buffer);
      a_show.startFrame();
      dbp = a_show.getFrame();
      if (dbp == 0) {
         OV_ABORT("Can not get frame from a_show");
      }

      dbp->put(a_time, "time");
      dbp->put(a_dt, "dt");
      a_show.saveComment(0, sPrintF(buffer, "VlasovPoisson4D"));
      a_show.saveComment(1,
         sPrintF(buffer, "t=%5.2f, dt=%5.2e", a_time, a_dt));
      a_show.saveSolution(m_ke_flux_vy_hi[is]);
      a_show.endFrame();
   }

   // We finished saving another plot cycle.
   ++a_saved_save;
}


void
Poisson::accumulateSequences(
   RealArray& a_sequences,
   const RealArray& a_probes,
   int a_num_probes,
   int a_saved_seq,
   int& a_seq_idx)
{
   // Initialize the time history to 0.
   // EM time histories
   real e_sum_tot(0.0);
   real e_max(0.0);
   real e_tot(0.0);
   real ex_max(0.0);
   real ey_max(0.0);

   // Probe time histories
   RealArray eProbe(TIME_HISTS_PER_PROBE, a_num_probes);
   RealArray eProbeTotal(TIME_HISTS_PER_PROBE, a_num_probes);
   for (int ix(0); ix < TIME_HISTS_PER_PROBE; ++ix) {
      for (int ip(0); ip < a_num_probes; ++ip) {
         eProbe(ix, ip) = 0.0;
      }
   }

   // Only processors onto which the Poisson calculation is distributed need to
   // do this.
   if (isPoissonProcessor()) {
      // Get the local x and y components of the E field.
#ifdef USE_PPP
      RealArray Ex;
      getLocalArrayWithGhostBoundaries(m_electric_field_x[0], Ex);

      RealArray Ey;
      getLocalArrayWithGhostBoundaries(m_electric_field_y[0], Ey);
#else
      RealArray& Ex = m_electric_field_x[0];
      RealArray& Ey = m_electric_field_y[0];
#endif

      // Get a bunch of box info needed for indexing into the E field arrays.
      const tbox::Box local_box(BoxOps::getLocalBox(Ex));
      const tbox::Box& domain_box(m_domain->box());
      const tbox::Box global_box(tbox::Box::grow(domain_box, m_n_ghosts));
      const tbox::Box interior_box = BoxOps::getOperationalBox(Ex,
         m_electric_field_x[0],
         domain_box,
         global_box);
      const int x_lo = interior_box.lower(0) - local_box.lower(0);
      const int x_hi = interior_box.upper(0) - local_box.lower(0);
      const int y_lo = interior_box.lower(1) - local_box.lower(1);
      const int y_hi = interior_box.upper(1) - local_box.lower(1);
      const int nx = local_box.upper(0) - local_box.lower(0) + 1;

      // Get some other misc quantities we'll need along the way.
      const tbox::IntVector N(m_domain->numberOfCells());
      const Array<double>& dx(m_domain->dx());
      const real area(dx[X1] * dx[X2]);
      const double* Ex_ptr = Ex.getDataPointer();
      const double* Ey_ptr = Ey.getDataPointer();

      // Compute the local parts of the E field time histories.
      for (int i2(y_lo); i2 <= y_hi; ++i2) {
         int idx = i2*nx+x_lo;
         for (int i1(x_lo); i1 <= x_hi; ++i1) {
            real ex = Ex_ptr[idx];
            real ey = Ey_ptr[idx];
            real tmp = ex * ex + ey * ey;
            real e_loc = sqrt(tmp);
            e_sum_tot += 0.5*tmp;
            e_max = max(e_max, e_loc);
            e_tot += e_loc;
            ex_max = max(ex_max, abs(ex));
            ey_max = max(ey_max, abs(ey));
            ++idx;
         }
      }
      // These are integrated over configuration space.
      e_sum_tot *= area;
      e_tot     *= area;

      // Compute the local part of the probe time histories.
      for (int np(0); np < a_num_probes; ++np) {
         const int ip(int(floor(a_probes(0, np) * N[X1])));
         const int jp(int(floor(a_probes(1, np) * N[X2])));
         if (ip >= interior_box.lower(0) && ip <= interior_box.upper(0) &&
             jp >= interior_box.lower(1) && jp <= interior_box.upper(1)) {
            int idx = (jp-local_box.lower(1))*nx+(ip-local_box.lower(0));
            eProbe(0, np) = Ex_ptr[idx];
            eProbe(1, np) = Ey_ptr[idx];
         }
      }
   }

   // Now get the appropriate max/sum of the E field time histories from all the
   // processors and store to the time history.
   a_sequences(a_saved_seq, a_seq_idx++) =
      ParallelUtility::getMaxValue(e_max, -1);
   a_sequences(a_saved_seq, a_seq_idx++) =
      ParallelUtility::getSum(e_tot, -1);
   a_sequences(a_saved_seq, a_seq_idx++) =
      ParallelUtility::getMaxValue(ex_max, -1);
   a_sequences(a_saved_seq, a_seq_idx++) =
      ParallelUtility::getMaxValue(ey_max, -1);
   a_sequences(a_saved_seq, a_seq_idx++) =
      ParallelUtility::getSum(e_sum_tot, -1);

   // Get the sum of the probe time histories from all the processors and store
   // to the time history.
   ParallelUtility::getSums(&(*eProbe.getDataPointer()),
      &(*eProbeTotal.getDataPointer()),
      TIME_HISTS_PER_PROBE * a_num_probes,
      -1);
   for (int ip(0); ip < a_num_probes; ++ip) {
      a_sequences(a_saved_seq, a_seq_idx++) = eProbeTotal(0, ip);
      a_sequences(a_saved_seq, a_seq_idx++) = eProbeTotal(1, ip);
   }

   // Store the tracking particle positions and velocities to the time history.
   for (int ip = 0; ip < numTrackingParticles(); ++ip) {
      a_sequences(a_saved_seq, a_seq_idx++) = m_tracking_particles[ip].x();
      a_sequences(a_saved_seq, a_seq_idx++) = m_tracking_particles[ip].y();
      a_sequences(a_saved_seq, a_seq_idx++) = m_tracking_particles[ip].vx();
      a_sequences(a_saved_seq, a_seq_idx++) = m_tracking_particles[ip].vy();
   }
}


void
Poisson::getFromRestart(
   const HDF_DataBase& a_db)
{
   // No restart data needs to be read.
   NULL_USE(a_db);
}


void
Poisson::putToRestart(
   HDF_DataBase& a_db,
   real a_time)
{
   // We write the number of ghosts only so that the post processor knows it.
   NULL_USE(a_time);
   a_db.put(m_n_ghosts[0], "nGhost");
}

int
Poisson::readParticleData(
   aString& a_particle_file,
   std::vector<Particle>& a_particles,
   bool a_read_noise_source_weight)
{
   herr_t errf;

   // Open particle file.
   hid_t file_id =
      H5Fopen(a_particle_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   if (file_id < 0) {
      OV_ABORT("Unable to open particle file.");
   }

   // Open and read the dataset containing the number of particles.
   int num_particles;
   hid_t dataset_id;
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5DOpen_vers) && H5DOpen_vers == 1
   dataset_id = H5Dopen(file_id, "num_particles");
#else
   dataset_id = H5Dopen(file_id, "num_particles", H5P_DEFAULT);
#endif
#else
   dataset_id = H5Dopen(file_id, "num_particles");
#endif
   if (dataset_id < 0) {
      OV_ABORT("Can not open dataset \"num_particles\".");
   }
   else {
      hid_t dspace = H5Dget_space(dataset_id);
      if (dspace < 0) {
         OV_ABORT("Can not get dataspace for dataset \"num_particles\".");
      }
      errf = H5Dread(dataset_id,
                     H5T_NATIVE_INT,
                     H5S_ALL,
                     H5S_ALL,
                     H5P_DEFAULT,
                     &num_particles);
      if (errf < 0) {
         OV_ABORT("Can not read dataset \"num_particles\".");
      }

      // Close the dataspace.
      errf = H5Sclose(dspace);
      if (errf < 0) {
         OV_ABORT("Can not close dataspace for dataset \"num_particles\".");
      }

      // Close the dataset.
      errf = H5Dclose(dataset_id);
      if (errf < 0) {
         OV_ABORT("Can not close dataset \"num_particles\".");
      }
   }

   // If this Poisson object is not on a Poisson processor we don't read the
   // particles themselves and we're done.
   if (!isPoissonProcessor()) {
      // Close the file.
      errf = H5Fclose(file_id);
      if (errf < 0) {
         OV_ABORT("Can not close particle file.");
      }

      return num_particles;
   }

   // Open and read the dataset containing the initial phase space location
   // of all the particles.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
   dataset_id = H5Dopen(file_id, "particle_data");
#else
   dataset_id = H5Dopen(file_id, "particle_data", H5P_DEFAULT);
#endif
#else
   dataset_id = H5Dopen(file_id, "particle_data");
#endif
   if (dataset_id < 0) {
      OV_ABORT("Can not open dataset \"particle_data\".");
   }
   else {
      // Check that the dataset is the right size.
      int data_size;
      if (a_read_noise_source_weight) {
         data_size = 8*num_particles;
      }
      else {
         data_size = 7*num_particles;
      }
      hid_t dspace = H5Dget_space(dataset_id);
      if (dspace < 0) {
         OV_ABORT("Can not get dataspace for dataset \"particle_data\".");
      }
      hsize_t nsel = H5Sget_select_npoints(dspace);
      if (static_cast<int>(nsel) != data_size) {
         OV_ABORT("Incorrect amount of particle data.");
      }

      // Now that everything checks out, read the dataset.
      double* particle_data = new double [data_size];
      errf = H5Dread(dataset_id,
                     H5T_NATIVE_DOUBLE,
                     H5S_ALL,
                     H5S_ALL,
                     H5P_DEFAULT,
                     particle_data);
      if (errf < 0) {
         OV_ABORT("Can not read dataset \"particle_data\".");
      }

      // Close the dataspace.
      errf = H5Sclose(dspace);
      if (errf < 0) {
         OV_ABORT("Can not close dataspace for dataset \"particle_data\".");
      }

      // Close the dataset.
      errf = H5Dclose(dataset_id);
      if (errf < 0) {
         OV_ABORT("Can not close dataset \"particle_data\".");
      }

      int idx = 0;
      for (int i = 0; i < num_particles; ++i) {
         double start = particle_data[idx++];
         double charge = particle_data[idx++];
         double mass = particle_data[idx++];
         double xinit = particle_data[idx++];
         double yinit = particle_data[idx++];
         double vxinit = particle_data[idx++];
         double vyinit = particle_data[idx++];
         double noise_source_weight;
         if (a_read_noise_source_weight) {
            noise_source_weight = particle_data[idx++];
         }
         else {
            noise_source_weight = 0.0;
         }
         // Let's make sure that the initial location of the particle is inside
         // the physical domain.
         if (xinit < m_domain->lower(X1) ) {
            if (m_domain->isPeriodic(X1)) {
               int domain_blocks =
                  int(1.0 + (m_domain->lower(X1) - xinit)/
                      (m_domain->upper(X1) - m_domain->lower(X1)));
               xinit +=
                  domain_blocks*(m_domain->upper(X1) - m_domain->lower(X1));
            }
            else {
               xinit = m_domain->lower(X1);
            }
         }
         else if (xinit > m_domain->upper(X1)) {
            if (m_domain->isPeriodic(X1)) {
               int domain_blocks =
                  int(1.0 + (xinit - m_domain->upper(X1))/
                      (m_domain->upper(X1) - m_domain->lower(X1)));
               xinit -=
                  domain_blocks*(m_domain->upper(X1) - m_domain->lower(X1));
            }
            else {
               xinit = m_domain->upper(X1);
            }
         }
         if (yinit < m_domain->lower(X2) ) {
            if (m_domain->isPeriodic(X2)) {
               int domain_blocks =
                  int(1.0 + (m_domain->lower(X2) - yinit)/
                      (m_domain->upper(X2) - m_domain->lower(X2)));
               yinit +=
                  domain_blocks*(m_domain->upper(X2) - m_domain->lower(X2));
            }
            else {
               yinit = m_domain->lower(X2);
            }
         }
         else if (yinit > m_domain->upper(X2)) {
            if (m_domain->isPeriodic(X2)) {
               int domain_blocks =
                  int(1.0 + (yinit - m_domain->upper(X2))/
                      (m_domain->upper(X2) - m_domain->lower(X2)));
               yinit -=
                  domain_blocks*(m_domain->upper(X2) - m_domain->lower(X2));
            }
            else {
               yinit = m_domain->upper(X2);
            }
         }
         if (vxinit < m_domain->lower(V1)) {
            vxinit = m_domain->lower(V1);
         }
         else if (vxinit > m_domain->upper(V1)) {
            vxinit = m_domain->upper(V1);
         }
         if (vyinit < m_domain->lower(V2)) {
            vyinit = m_domain->lower(V2);
         }
         else if (vyinit > m_domain->upper(V2)) {
            vyinit = m_domain->upper(V2);
         }
         Particle particle(start,
                           charge,
                           mass,
                           xinit,
                           yinit,
                           vxinit,
                           vyinit,
                           noise_source_weight);
         a_particles.push_back(particle);
      }
      delete [] particle_data;
   }

   // Close the file.
   errf = H5Fclose(file_id);
   if (errf < 0) {
      OV_ABORT("Can not close particle file.");
   }

   return num_particles;
}

} // end namespace Loki
