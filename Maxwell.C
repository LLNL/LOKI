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
#include "Maxwell.H"

#include "MappedGridOperators.h"
#include "OgesParameters.h"
#include "TridiagonalSolver.h"
#include "PlotStuff.h"
#include "PlotStuffParameters.h"

#include "BoxOps.H"
#include "CurrentDriverFactory.H"
#include "EMICFactory.H"
#include "VELICFactory.H"
#include "PoissonF.H"
#include "RestartManager.H"
#include "Loki_Utilities.H"
#include "Interpolator4.H"
#include "Interpolator6.H"

#include "hdf5.h"

#include <sstream>

namespace Loki {

int Maxwell::TIME_HISTS_PER_PROBE;
const int Maxwell::GLOBAL_TIME_HISTS = 12;
int Maxwell::NUM_FRAME_SERIES;

Maxwell::Maxwell(ParmParse& a_pp,
   tbox::Pointer<ProblemDomain> a_domain,
   int a_num_kinetic_species,
   int a_solution_order)
   : m_dim(a_domain->dim()),
     m_domain(a_domain),
     m_number_of_procs(1),
     m_partition_defined(false),
     m_n_ghosts(m_dim),
     m_local_box(m_dim),
     m_interior_box(m_dim),
     m_global_box(m_dim),
     m_is_maxwell_processor(false),
     m_comm(MPI_COMM_NULL),
     m_num_kinetic_species(a_num_kinetic_species),
     m_vz_global(a_num_kinetic_species),
#ifdef USE_PPP
     m_vz_local(a_num_kinetic_species),
#else
     m_em_vars_local(m_em_vars_global),
     m_vz_local(m_vz_global),
#endif
     m_c(1.0),
     m_avWeak(0.1),
     m_avStrong(0.0),
     m_solver(0),
     m_vz_plot(a_num_kinetic_species),
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
     m_current_drivers(0),
     m_num_current_drivers(0),
     m_solution_order(a_solution_order),
     m_em_ics(0),
     m_vel_ics(0),
     m_ex_interp(0),
     m_ey_interp(0),
     m_bz_interp(0)
{
   // 6 components of EM + vz for each species.
   TIME_HISTS_PER_PROBE = 6 + a_num_kinetic_species;
   NUM_FRAME_SERIES = 6 + 5*a_num_kinetic_species;

   if (m_dim != 2) {
      OV_ABORT("Maxwell only implemented for D=2!");
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

   // The following code could all probably be placed into parseParameters but
   // isn't because it involves parsing separate sub-databases of this Maxwell.

   // Get the sub-database for any current drivers and create them.
   char buffer[100];
   for (int i = 0; i < m_num_current_drivers; ++i) {
      sprintf(buffer, "maxwell.current_driver.%i", i+1);
      ParmParse cdf_pp(buffer);
      m_current_drivers[i] = CurrentDriverFactory::create(cdf_pp, i+1);
   }

   // Get the sub-database for each EM initial condition and create them.
   // There must be 2 and only 2--one for E and one for B.
   if (a_pp.contains("maxwell.emic.3.name")) {
      OV_ABORT("There may be only 2 electromagnetic initial conditions.");
   }
   m_em_ics.resize(2, tbox::Pointer<EMICInterface>(0));
   int num_E_initializers = 0;
   int num_B_initializers = 0;
   for (int i = 0; i < 2; ++i) {
      sprintf(buffer, "maxwell.em_ic.%i", i+1);
      ParmParse emic_pp(buffer);
      m_em_ics[i] = EMICFactory::create(emic_pp);
      if (m_em_ics[i]->initializesE()) {
         ++num_E_initializers;
      }
      else {
         ++num_B_initializers;
      }
   }
   if (num_E_initializers != 1 || num_B_initializers != 1) {
      OV_ABORT("There must be 1 initial condition for each of E and B");
   }

   // Get the sub-database for each species transverse drift velocity initial
   // condition and create them.
   m_vel_ics.resize(m_num_kinetic_species, tbox::Pointer<VELICInterface>(0));
   for (int i = 0; i < m_num_kinetic_species; ++i) {
      sprintf(buffer, "maxwell.vel_ic.%i", i+1);
      ParmParse vel_ic_pp(buffer);
      m_vel_ics[i] = VELICFactory::create(vel_ic_pp);
   }

   // The Maxwell object writes restart data so register this object with the
   // RestartManager which will use the putToRestart/getFromRestart callbacks to
   // get the restart data written/read.
   RestartManager* restart_manager(RestartManager::getManager());
   restart_manager->registerRestart(this);

   printParameters();
}


Maxwell::Maxwell(
   const Maxwell& a_other,
   bool a_deep_copy)
   : m_dim(a_other.m_dim),
     m_n_ghosts(m_dim),
     m_local_box(m_dim),
     m_interior_box(m_dim),
     m_global_box(m_dim),
     m_vz_global(a_other.m_num_kinetic_species),
#ifdef USE_PPP
     m_vz_local(a_other.m_num_kinetic_species),
#else
     m_em_vars_local(m_em_vars_global),
     m_vz_local(m_vz_global),
#endif
     m_ex_interp(0),
     m_ey_interp(0),
     m_bz_interp(0)
{
   copy(a_other, a_deep_copy);
}


Maxwell::~Maxwell()
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
   if (m_bz_interp) {
      for (int i = 0; i < m_solution_order; ++i) {
         delete [] m_bz_interp[i];
      }
      delete [] m_bz_interp;
   }
}


float
Maxwell::netCost() const
{
   // The computational cost of a Maxwell is proportional to the size of its
   // domain.
   float cost_per_cell(1.0);
   return cost_per_cell *
      static_cast<float>((m_domain->numberOfCells()).getProduct());
}


int
Maxwell::numberOfProcessors() const
{
   return m_number_of_procs;
}


bool
Maxwell::fixedNumberOfProcessors() const
{
   return true;
}


void
Maxwell::createPartition(
   const Range& a_range,
   const MPI_Comm& a_comm)
{
   m_comm = a_comm;
   m_processor_range = a_range;
   m_partition_defined = true;

   const int my_id(std::max(0, Communication_Manager::My_Process_Number));
   if (isInRange(my_id)) {
      m_is_maxwell_processor = true;
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

   // Partition the global array and add the ghosts to the global box.  This
   // gives the range of each dimension which is needed in order to redimension
   // it.  I think that Overture's term "partition" is not accurate.  I think
   // that the partition call just defines which processors the global array
   // lives on and how it CAN be partitioned among those processors.  HOW it is
   // partitioned is determined by the redim call.  It can't very well be
   // partitioned unless its extent is known.
   m_em_vars_global.partition(m_partition);
   m_global_box = m_domain->box();
   m_global_box.grow(m_n_ghosts);
   if (m_dim == tbox::Dimension(2)) {
      m_em_vars_global.redim(BoxOps::range(m_global_box, X1),
         BoxOps::range(m_global_box, X2),
         Range(0, NUM_EM_VARS-1));
   }
   else {
      OV_ABORT("Not implemented for phase D!=2!");
   }

   // Now that the global array has been divvied up we can get the local array
   // for the EM vars.
#ifdef USE_PPP
   getLocalArrayWithGhostBoundaries(m_em_vars_global, m_em_vars_local);
#endif

   // Repeat the above process for the global array of each species' transverse
   // drift velocity.
   if (m_dim == tbox::Dimension(2)) {
      for (int i = 0; i < m_num_kinetic_species; ++i) {
         m_vz_global[i].partition(m_partition);
         m_vz_global[i].redim(BoxOps::range(m_global_box, X1),
            BoxOps::range(m_global_box, X2),
            Range(0, 0));
      }
   }
   else {
      OV_ABORT("Not implemented for phase D!=2!");
   }

   // Now that the global array has been divvied up we can get the local array
   // for the transverse drift velocity of each species.
#ifdef USE_PPP
   for (int i = 0; i < m_num_kinetic_species; ++i) {
      getLocalArrayWithGhostBoundaries(m_vz_global[i], m_vz_local[i]);
   }
#endif

   // Figure out the local and interior boxes.  Recall that all Maxwells exist
   // on all processors that only 1 (or possibly a few) processors have a
   // Maxwell that does work.
   if (isInRange(my_id)) {
      tbox::Box tmp = BoxOps::getLocalBox(m_em_vars_local);
      for (int dir(0); dir < m_dim; ++dir) {
         m_local_box.lower(dir) = tmp.lower(dir);
         m_local_box.upper(dir) = tmp.upper(dir);
      }
      m_interior_box = m_domain->box();
      tbox::IntVector diff(m_interior_box.upper() - m_interior_box.lower());
      for (int dir(0); dir < m_dim; ++dir) {
         if (diff[dir] < m_n_ghosts[dir]) {
            OV_ABORT("Too few interior points in decomposition!");
         }
      }
   }
   else {
      for (int dir(0); dir < m_dim; ++dir) {
         m_local_box.lower(dir) = 0;
         m_local_box.upper(dir) = -1;
         m_interior_box.lower(dir) = 0;
         m_interior_box.upper(dir) = -1;
      }
   }
}


bool
Maxwell::isInRange(
   int a_proc_id) const
{
   // Returns true if the Maxwell calculation is partitioned onto this
   // processor.
   return ((a_proc_id >= m_processor_range.getBase()) &&
           (a_proc_id <= m_processor_range.getBound()));
}


void
Maxwell::printDecomposition() const
{
   // This function is only valid if we actually know the decomposition.
   if (m_partition_defined) {
      // Print some basic decomposition info.
      printF("  Maxwell processor(s):  [%d,%d]\n",
             m_processor_range.getBase(),
             m_processor_range.getBound());
   }
}


bool
Maxwell::conformsTo(
   const Maxwell& a_other,
   bool a_include_ghost_cells) const
{
   // If the 2 Maxwells are different, check that they are defined on the same
   // piece of space and are partitioned the same.
   if (m_domain != a_other.m_domain ||
       m_processor_range != a_other.m_processor_range ||
       m_partition_defined != a_other.m_partition_defined ||
       numTrackingParticles() != a_other.numTrackingParticles()) {
      return false;
   }
   else if (!a_include_ghost_cells) {
      return (m_domain->box() == a_other.m_domain->box());
   }
   else {
      return m_global_box == a_other.m_global_box;
   }
}


void
Maxwell::addData(
   const Maxwell& a_increment_maxwell,
   real a_factor,
   bool a_final_rk)
{
   // Check that we're copying between similar Maxwells.  To each tracking
   // particle's x, y, vx, and vy add a_factor times x, y, vx, and vy of the
   // corresponding tracking particle in a_increment_maxwell.  Add a_factor
   // times a_increment_maxwell to this Maxwell's EM fields.
   if (conformsTo(a_increment_maxwell, false)) {
      // Handle the EM fields and the transverse drift velocities of each
      // species.
      tbox::Box intersect_box(
         m_interior_box * a_increment_maxwell.m_interior_box);
      FORT_XPBY_2D(*m_em_vars_local.getDataPointer(),
         BOX2D_TO_FORT(m_local_box),
         *a_increment_maxwell.m_em_vars_local.getDataPointer(),
         BOX2D_TO_FORT(intersect_box),
         a_factor,
         NUM_EM_VARS);
      for (int i = 0; i < m_num_kinetic_species; ++i) {
         FORT_XPBY_2D(*m_vz_local[i].getDataPointer(),
            BOX2D_TO_FORT(m_local_box),
            *a_increment_maxwell.m_vz_local[i].getDataPointer(),
            BOX2D_TO_FORT(intersect_box),
            a_factor,
            1);
      }

      // Handle the tracking particles.
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("Tracking Particles");
      for (int i = 0; i < numTrackingParticles(); ++i) {
         m_tracking_particles[i].x() +=
            a_increment_maxwell.m_tracking_particles[i].x()*a_factor;
         m_tracking_particles[i].y() +=
            a_increment_maxwell.m_tracking_particles[i].y()*a_factor;
         m_tracking_particles[i].vx() +=
            a_increment_maxwell.m_tracking_particles[i].vx()*a_factor;
         m_tracking_particles[i].vy() +=
            a_increment_maxwell.m_tracking_particles[i].vy()*a_factor;
         // This is the final RK update so the tracking particle positions
         // need to be limited.  If a tracking particle is outside the physical
         // domain in a periodic direction then wrap it back around to its
         // periodic location.  If a tracking particle is outside the physical
         // domain in a non-periodic direction then limit its position to the
         // edge of the physical domain in that direction.
         if (a_final_rk) {
            if (m_tracking_particles[i].x() < m_domain->lower(X1)) {
               if (m_domain->isPeriodic(X1)) {
                  m_tracking_particles[i].x() += m_domain->upper(X1) -
                     m_domain->lower(X1);
               }
               else {
                  m_tracking_particles[i].x() = m_domain->lower(X1);
               }
            }
            else if (m_tracking_particles[i].x() > m_domain->upper(X1)) {
               if (m_domain->isPeriodic(X1)) {
                  m_tracking_particles[i].x() += m_domain->lower(X1) -
                     m_domain->upper(X1);
               }
               else {
                  m_tracking_particles[i].x() = m_domain->upper(X1);
               }
            }
            if (m_tracking_particles[i].y() < m_domain->lower(X2)) {
               if (m_domain->isPeriodic(X2)) {
                  m_tracking_particles[i].y() += m_domain->upper(X2) -
                     m_domain->lower(X2);
               }
               else {
                  m_tracking_particles[i].y() = m_domain->lower(X2);
               }
            }
            else if (m_tracking_particles[i].y() > m_domain->upper(X2)) {
               if (m_domain->isPeriodic(X2)) {
                  m_tracking_particles[i].y() += m_domain->lower(X2) -
                     m_domain->upper(X2);
               }
               else {
                  m_tracking_particles[i].y() = m_domain->upper(X2);
               }
            }
            if (m_tracking_particles[i].vx() < m_domain->lower(V1)) {
               m_tracking_particles[i].vx() = m_domain->lower(V1);
            }
            else if (m_tracking_particles[i].vx() > m_domain->upper(V1)) {
               m_tracking_particles[i].vx() = m_domain->upper(V1);
            }
            if (m_tracking_particles[i].vy() < m_domain->lower(V2)) {
               m_tracking_particles[i].vy() = m_domain->lower(V2);
            }
            else if (m_tracking_particles[i].vy() > m_domain->upper(V2)) {
               m_tracking_particles[i].vy() = m_domain->upper(V2);
            }
         }
      }
      timers->stopTimer("Tracking Particles");
   }
}


void
Maxwell::copySolnData(
   const Maxwell& a_rhs)
{
   if (m_dim != a_rhs.m_dim) {
      OV_ABORT("Attemtpt to copy incongruent Maxwells!");
   }

   // If the 2 Maxwells are different copy the EM fields, each species
   // transverse drift velocity, and the tracking particles' x, y, vx, and vy.
   if (&a_rhs != this) {
      // If there is > 1 Maxwell process this check is wrong.
      // In any event, deal with the EM and transverse drift velocities.  If the
      // box configuation of the 2 Maxwells is different then we need parallel
      // copy operations.  Otherwise it's a simple assignment.
      if (m_global_box.lower(0) != a_rhs.m_global_box.lower(0) ||
          m_global_box.upper(0) != a_rhs.m_global_box.upper(0) ||
          m_global_box.lower(1) != a_rhs.m_global_box.lower(1) ||
          m_global_box.upper(1) != a_rhs.m_global_box.upper(1)) {
         m_em_vars_global.redim(0);
         m_em_vars_global.partition(a_rhs.m_em_vars_global.getPartition());
         m_em_vars_global.redim(a_rhs.m_em_vars_global);
         for (int i = 0; i < m_num_kinetic_species; ++i) {
            m_vz_global[i].redim(0);
            m_vz_global[i].partition(a_rhs.m_vz_global[i].getPartition());
            m_vz_global[i].redim(a_rhs.m_vz_global[i]);
         }

         Index dst[4];
         Index* src = dst;
         for (int i(0); i < m_dim+1; ++i) {
            dst[i] =
               Range(a_rhs.m_em_vars_global.getBase(i),
                     a_rhs.m_em_vars_global.getBound(i));
         }
         CopyArray::copyArray(m_em_vars_global, dst, a_rhs.m_em_vars_global, src);

         dst[m_dim] = Range(0, 0);
         for (int i = 0; i < m_num_kinetic_species; ++i) {
            CopyArray::copyArray(m_vz_global[i], dst, a_rhs.m_vz_global[i], src);
         }
#ifdef USE_PPP
         getLocalArrayWithGhostBoundaries(m_em_vars_global, m_em_vars_local);
         for (int i = 0; i < m_num_kinetic_species; ++i) {
            getLocalArrayWithGhostBoundaries(m_vz_global[i], m_vz_local[i]);
         }
#endif
      }
      else {
         m_em_vars_local = a_rhs.m_em_vars_local;
         for (int i = 0; i < m_num_kinetic_species; ++i) {
            m_vz_local[i] = a_rhs.m_vz_local[i];
         }
      }

      // Copy the tracking particles' data.
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("Tracking Particles");
      for (int i = 0; i < numTrackingParticles(); ++i) {
         m_tracking_particles[i].x() = a_rhs.m_tracking_particles[i].x();
         m_tracking_particles[i].y() = a_rhs.m_tracking_particles[i].y();
         m_tracking_particles[i].vx() = a_rhs.m_tracking_particles[i].vx();
         m_tracking_particles[i].vy() = a_rhs.m_tracking_particles[i].vy();
      }
      timers->stopTimer("Tracking Particles");
   }
}


real
Maxwell::computeDt()
{
   return 1.0/(m_c*(1.0/m_domain->dx(0) + 1.0/m_domain->dx(1)));
}


void
Maxwell::initialize(
   const std::vector<string> a_species_names)
{
   // This is lifted straight from Poisson so you can read it for any details.
   // Currently this is pretty irrelevant as we don't do a Poisson solve in the
   // electrodynamic case.
   defineSolver();

   m_phi_solver = realCompositeGridFunction(m_solver_composite_grid);
   m_phi_solver = 0.0;

   Mapping* physical_mapping = createMapping(-1);
   MappedGrid physical_mg(*physical_mapping);
   setupMappedGrid(physical_mg, 0);
   m_composite_grid.add(physical_mg);
   m_composite_grid.updateReferences();
   m_composite_grid.update(MappedGrid::THEvertex | MappedGrid::THEcenter);

   m_phi = realCompositeGridFunction(m_composite_grid);
   m_phi.setName("physical function");
   m_phi.setName("potential", 0);
   m_phi = 0.0;

   m_ex_plot = realCompositeGridFunction(m_composite_grid);
   m_ex_plot.setName("Ex");
   m_ex_plot.setName("Ex", 0);
   m_ex_plot = 0.0;

   m_ey_plot = realCompositeGridFunction(m_composite_grid);
   m_ey_plot.setName("Ey");
   m_ey_plot.setName("Ey", 0);
   m_ey_plot = 0.0;

   m_ez_plot = realCompositeGridFunction(m_composite_grid);
   m_ez_plot.setName("Ez");
   m_ez_plot.setName("Ez", 0);
   m_ez_plot = 0.0;

   m_bx_plot = realCompositeGridFunction(m_composite_grid);
   m_bx_plot.setName("Bx");
   m_bx_plot.setName("Bx", 0);
   m_bx_plot = 0.0;

   m_by_plot = realCompositeGridFunction(m_composite_grid);
   m_by_plot.setName("By");
   m_by_plot.setName("By", 0);
   m_by_plot = 0.0;

   m_bz_plot = realCompositeGridFunction(m_composite_grid);
   m_bz_plot.setName("Bz");
   m_bz_plot.setName("Bz", 0);
   m_bz_plot = 0.0;

   for (int i = 0; i < m_num_kinetic_species; ++i) {
      char buffer[10];
      sprintf(buffer, "vz%i", i+1);
      m_vz_plot[i] = realCompositeGridFunction(m_composite_grid);
      m_vz_plot[i].setName(buffer);
      m_vz_plot[i].setName(buffer, 0);
      m_vz_plot[i] = 0.0;

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

   if (isMaxwellProcessor()) {
      m_em_vars_local = 0.0;
      for (int i = 0; i < m_num_kinetic_species; ++i) {
         m_vz_local[i] = 0.0;
      }
   }

   // Now we know enough to build the appropriate field interpolator if needed.
   buildInterpolator();
}


void
Maxwell::buildInterpolator()
{
   // If this object is on the Maxwell processor and there are tracking
   // particles then we need to construct the interpolator of the appropriate
   // order.
   if (isMaxwellProcessor() && m_problem_num_tracking_particles > 0) {
      const tbox::Box local_box(BoxOps::getLocalBox(m_em_vars_local));
      if (m_solution_order == 4) {
         m_interpolator = new Interpolator4(local_box, *m_domain);
      }
      else {
         m_interpolator = new Interpolator6(local_box, *m_domain);
      }
      m_ex_interp = new double* [m_solution_order];
      for (int i = 0; i < m_solution_order; ++i) {
         m_ex_interp[i] = new double [m_solution_order];
      }
      m_ey_interp = new double* [m_solution_order];
      for (int i = 0; i < m_solution_order; ++i) {
         m_ey_interp[i] = new double [m_solution_order];
      }
      m_bz_interp = new double* [m_solution_order];
      for (int i = 0; i < m_solution_order; ++i) {
         m_bz_interp[i] = new double [m_solution_order];
      }
   }
}


void
Maxwell::initializeEM()
{
   // The Maxwell processor(s) must initialize the EM fields.
   if (isMaxwellProcessor()) {
      for (int i = 0; i < 2; ++i) {
         m_em_ics[i]->set(m_em_vars_local, *m_domain);
      }
   }
}


void
Maxwell::initializeVZ()
{
   // The Maxwell processor(s) must initialize the transverse drift velocities.
   if (isMaxwellProcessor()) {
      for (int i = 0; i < m_num_kinetic_species; ++i) {
         m_vel_ics[i]->set(m_vz_local[i], *m_domain);
      }
   }
}


void
Maxwell::computeFactorization()
{
   // Do an initial solve to obtain the static factorization.
   // Again, not relevant as we do not currently do a Poisson solve here.
   m_solver->solve(m_phi_solver, m_rho_solver);
   m_solver->setRefactor(FALSE);
   m_solver->setReorder(FALSE);
}


void
Maxwell::applyExternalCurrent(
   realArray& a_Jx,
   realArray& a_Jy,
   realArray& a_Jz,
   real a_time)
{
   // Get the local current components.
#ifdef USE_PPP
   RealArray Jx_local, Jy_local, Jz_local;
   getLocalArrayWithGhostBoundaries(a_Jx, Jx_local);
   getLocalArrayWithGhostBoundaries(a_Jy, Jy_local);
   getLocalArrayWithGhostBoundaries(a_Jz, Jz_local);
#else
   RealArray& Jx_local = a_Jx;
   RealArray& Jy_local = a_Jy;
   RealArray& Jz_local = a_Jz;
#endif
   // Add the external current driver to the existing current.
   for (int i = 0; i < m_num_current_drivers; ++i) {
      m_current_drivers[i]->evaluate(Jx_local,
         Jy_local,
         Jz_local,
         m_local_box,
         *m_domain,
         a_time);
   }
}


void
Maxwell::electricField(
   const realArray& a_net_charge_density)
{
   // Lifted from Poisson.C.  Again, currently irrelevant as we do not solve for
   // the electrostatic E field in the electrodynamic case.
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("Poisson");

   const tbox::Box& domain_box(m_domain->box());

   Index dest[4], src[4];
   dest[0] = Range(BoxOps::range(domain_box, X1));
   dest[1] = Range(BoxOps::range(domain_box, X2));
   dest[2] = Range(0, 0);
   dest[3] = Range(0, 0);
   src[0] = Range(BoxOps::range(domain_box, X1));
   src[1] = Range(BoxOps::range(domain_box, X2));
   src[2] = Range(0, 0);
   src[3] = Range(0, 0);
   ParallelUtility::copy(m_rho_solver[0], dest, a_net_charge_density, src, 4);

#ifdef USE_PPP
   RealArray rhs_local;
   getLocalArrayWithGhostBoundaries(m_rho_solver[0], rhs_local);
#else
   RealArray& rhs_local = m_rho_solver[0];
#endif

   const tbox::Box rhs_box(BoxOps::getLocalBox(rhs_local));
   const tbox::Box global_box(tbox::Box::grow(domain_box, m_n_ghosts));
   const tbox::Box interior_box = BoxOps::getOperationalBox(rhs_local,
      m_rho_solver[0],
      domain_box,
      global_box);
   neutralizeCharge4D(BOX2D_TO_FORT(rhs_box),
      BOX2D_TO_FORT(interior_box),
      *rhs_local.getDataPointer(),
      m_comm);
   // IS THIS NECESSARY FOR PARALLEL?
   m_rho_solver[0].updateGhostBoundaries();

   m_solver->solve(m_phi_solver, m_rho_solver);
   m_solver->setRefactor(FALSE);
   m_solver->setReorder(FALSE);

   timers->startTimer("BC (Maxwell)");
   // copy back to physical domain
   ParallelUtility::copy(m_phi[0], dest, m_phi_solver[0], src, 4);

   Loki_Utilities::fixPeriodicity(m_phi[0], *m_domain, m_n_ghosts, 1);
   timers->stopTimer("BC (Maxwell)");

   // this needs to be fixed once we distribute the Poisson solve?
   //   be particularly careful with the external field so that positions
   //   are computed correctly
#ifdef USE_PPP
   RealArray phi_local;
   getLocalArrayWithGhostBoundaries(m_phi[0], phi_local);
#else
   RealArray& phi_local = m_phi[0];
#endif

   computeEFieldFromPotentialMaxwell(BOX2D_TO_FORT(m_local_box),
      BOX2D_TO_FORT(interior_box),
      m_solution_order,
      *(m_domain->dx()).dataPtr(),
      *m_em_vars_local.getDataPointer(),
      *phi_local.getDataPointer());

   // fix periodicity of the electromagnetic fields
   timers->startTimer("BC (Maxwell)");
   Loki_Utilities::fixPeriodicity(m_em_vars_global,
      *m_domain,
      m_n_ghosts,
      NUM_EM_VARS);

   timers->stopTimer("BC (Maxwell)");
   timers->stopTimer("Poisson");
}


void
Maxwell::newAuxVariable(
   realArray& a_var,
   int a_depth,
   const tbox::IntVector& a_n_ghosts) const
{
   // Create any "auxilliary" variable that is partitioned like this Maxwell.
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
   }
   else {
      OV_ABORT("Attempt to define variables without partition!");
   }
}


void
Maxwell::evalRHS(
   Maxwell& a_rhs,
   const KineticSpeciesPtrVect& a_kspv,
   const realArray& a_Jx,
   const realArray& a_Jy,
   const realArray& a_Jz,
   const realArray& a_net_ext_efield,
   real a_time)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("Maxwell RHS");

   // fix periodicity of the electromagnetic fields
   timers->startTimer("BC (Maxwell)");
   Loki_Utilities::fixPeriodicity(m_em_vars_global,
      *m_domain,
      m_n_ghosts,
      NUM_EM_VARS);

   timers->stopTimer("BC (Maxwell)");

   // Get the local components of the current densities.
#ifdef USE_PPP
   RealArray Jx_local, Jy_local, Jz_local, E_ext;
   getLocalArrayWithGhostBoundaries(a_Jx, Jx_local);
   getLocalArrayWithGhostBoundaries(a_Jy, Jy_local);
   getLocalArrayWithGhostBoundaries(a_Jz, Jz_local);
   getLocalArrayWithGhostBoundaries(a_net_ext_efield, E_ext);
#else
   RealArray& Jx_local = a_Jx;
   RealArray& Jy_local = a_Jy;
   RealArray& Jz_local = a_Jz;
   RealArray& E_ext = a_net_ext_efield;
#endif

   // Evaluate the RHS of Maxwell's equations.
   FORT_MAXWELL_EVAL_RHS(BOX2D_TO_FORT(m_local_box),
      BOX2D_TO_FORT(m_interior_box),
      PROBLEMDOMAIN_TO_FORT((*m_domain)),
      m_c,
      m_avWeak,
      m_avStrong,
      m_solution_order,
      *(m_em_vars_local.getDataPointer()),
      *(Jx_local.getDataPointer()),
      *(Jy_local.getDataPointer()),
      *(Jz_local.getDataPointer()),
      *(a_rhs.m_em_vars_local.getDataPointer()));

   // For each species, evaluate the RHS of the transverse drift velocity
   // equations, dvz/dt.
   for (int i = 0; i < m_num_kinetic_species; ++i) {
      const KineticSpeciesPtr& ksp = a_kspv[i];
      FORT_MAXWELL_EVAL_VZ_RHS(BOX2D_TO_FORT(m_local_box),
         BOX2D_TO_FORT(m_interior_box),
         ksp->charge()/ksp->mass(),
         *(m_em_vars_local.getDataPointer()),
         *(a_rhs.m_vz_local[i].getDataPointer()));
   }

   timers->startTimer("Tracking Particles");
   int num_tracking_particles = numTrackingParticles();
   if (num_tracking_particles > 0) {
      // Get the local net external applied E field.
#ifdef USE_PPP
      RealArray E_ext;
      getLocalArrayWithGhostBoundaries(a_net_ext_efield, E_ext);
#else
      RealArray& E_ext = a_net_ext_efield;
#endif
      const double* em_vars_ptr = m_em_vars_local.getDataPointer();
      const double* E_ext_ptr = E_ext.getDataPointer();

      // For each tracking particle use an interpolating polynomial to compute
      // the E and B fields that it sees and update x, y, vx, and vy.
      for (int i = 0; i < num_tracking_particles; ++i) {
         const Particle& this_particle = m_tracking_particles[i];
         if (a_time >= this_particle.startingTime()) {
            // Set the evaluation point to be the location of this particle.
            m_interpolator->setEvaluationPoint(this_particle.x(),
                                               this_particle.y());

            // Figure out Ex, Ey, and Bz at each interpolating point.
            int nx = m_interpolator->nx();
            int nxny = m_interpolator->nxny();
            for (int xidx = 0; xidx < 4; ++xidx) {
               int xGridIdx = m_interpolator->xInterpolatingGridIndex(xidx);
               for (int yidx = 0; yidx < 4; ++yidx) {
                  int yGridIdx = m_interpolator->yInterpolatingGridIndex(yidx);
                  int ex_idx = yGridIdx*nx+xGridIdx+EX*nxny;
                  int ey_idx = yGridIdx*nx+xGridIdx+EY*nxny;
                  int bz_idx = yGridIdx*nx+xGridIdx+BZ*nxny;
                  m_ex_interp[xidx][yidx] = em_vars_ptr[ex_idx] +
                     E_ext_ptr[ex_idx];
                  m_ey_interp[xidx][yidx] = em_vars_ptr[ey_idx] +
                     E_ext_ptr[ey_idx];
                  m_bz_interp[xidx][yidx] = em_vars_ptr[bz_idx];
               }
            }

            // Now compute Ex, Ey, Bx, for this tracking particle.
            double Ex_particle = m_interpolator->interpolate(m_ex_interp);
            double Ey_particle = m_interpolator->interpolate(m_ey_interp);
            double Bz_particle = m_interpolator->interpolate(m_bz_interp);

            // The new x/y is the old vx/vy.  The new vx/vy is determined by
            // the Lorentz force.
            double qoverm = this_particle.charge()/this_particle.mass();
            double vx = this_particle.vx();
            double vy = this_particle.vy();
            a_rhs.m_tracking_particles[i].x() = vx;
            a_rhs.m_tracking_particles[i].y() = vy;
            a_rhs.m_tracking_particles[i].vx() =
               (Ex_particle+vy*Bz_particle)*qoverm;
            a_rhs.m_tracking_particles[i].vy() =
               (Ey_particle-vx*Bz_particle)*qoverm;
         }
      }
   }
   timers->stopTimer("Tracking Particles");
   timers->stopTimer("Maxwell RHS");
}


void
Maxwell::plot(
   real                       a_time,
   real                       a_dt,
   const RealArray&           a_sequences,
   const std::vector<string>& a_species_names,
   const RealArray&           a_time_seq,
   const RealArray&           a_probes, 
   int                        a_num_probes,
   int                        a_num_seq,
   int                        a_saved_seq,
   int&                       a_saved_save,
   Ogshow&                    a_show)
{
   const tbox::IntVector num_cells(m_domain->numberOfCells());
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
      name[idx++] = "Ez_max";
      name[idx++] = "E_tot";
      name[idx++] = "B_max";
      name[idx++] = "norm B";
      name[idx++] = "Bx_max";
      name[idx++] = "By_max";
      name[idx++] = "Bz_max";
      name[idx++] = "B_tot";
      for (int ip(0); ip < a_num_probes; ++ip) {
         sprintf(buffer, "ix_pos%i", ip);
         dbp->put(int(floor(a_probes(0, ip) * num_cells[X1])), buffer);

         sprintf(buffer, "iy_pos%i", ip);
         dbp->put(int(floor(a_probes(1, ip) * num_cells[X2])), buffer);

         sprintf(buffer, "Ex_pos%i", ip);
         name[idx++] = buffer;

         sprintf(buffer, "Ey_pos%i", ip);
         name[idx++] = buffer;

         sprintf(buffer, "Ez_pos%i", ip);
         name[idx++] = buffer;

         sprintf(buffer, "Bx_pos%i", ip);
         name[idx++] = buffer;

         sprintf(buffer, "By_pos%i", ip);
         name[idx++] = buffer;

         sprintf(buffer, "Bz_pos%i", ip);
         name[idx++] = buffer;

         for (int is = 0; is < m_num_kinetic_species; ++is) {
            sprintf(buffer, "vz_pos%i_species%i", ip, is);
            name[idx++] = buffer;
         }
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
      for (int is = 0; is < m_num_kinetic_species; ++is) {
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
      a_show.saveComment(0, sPrintF(buffer, "Maxwell"));
      a_show.saveComment(1,
         sPrintF(buffer, "t=%5.2f, dt=%5.2e", a_time, a_dt));

      // Transfer the x component of the E field from m_em_vars_local into
      // m_ex_plot.
#ifdef USE_PPP
      RealArray ex_plot_local;
      getLocalArrayWithGhostBoundaries(m_ex_plot[0], ex_plot_local);
#else
      RealArray& ex_plot_local = m_ex_plot[0];
#endif
      Range R0(m_em_vars_local.getBase(0), m_em_vars_local.getBound(0));
      Range R1(m_em_vars_local.getBase(1), m_em_vars_local.getBound(1));
      Range R2(EX, EX);
      ex_plot_local(R0, R1, R2) = m_em_vars_local(R0, R1, R2);

      // Write the Ex plot.
      a_show.saveSolution(m_ex_plot);

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
      a_show.saveComment(0, sPrintF(buffer, "Maxwell"));
      a_show.saveComment(1,
         sPrintF(buffer, "t=%5.2f, dt=%5.2e", a_time, a_dt));

      // Transfer the y component of the E field from m_em_vars_local into
      // m_ex_plot.
#ifdef USE_PPP
      RealArray ey_plot_local;
      getLocalArrayWithGhostBoundaries(m_ey_plot[0], ey_plot_local);
#else
      RealArray& ey_plot_local = m_ey_plot[0];
#endif
      Range R0(m_em_vars_local.getBase(0), m_em_vars_local.getBound(0));
      Range R1(m_em_vars_local.getBase(1), m_em_vars_local.getBound(1));
      Range R2(EY, EY);
      Range R3(0, 0);
      ey_plot_local(R0, R1, R3) = m_em_vars_local(R0, R1, R2);

      // Write the Ey plot.
      a_show.saveSolution(m_ey_plot);
      a_show.endFrame();
   }

   {
      // Set up the frame for Ez.
      a_show.setCurrentFrameSeries("Ez");
      a_show.startFrame();
      HDF_DataBase *dbp(a_show.getFrame());
      if (dbp == 0) {
         OV_ABORT("Can not get frame from a_show");
      }

      // Write some general plot info for the frame.
      dbp->put(a_time, "time");
      dbp->put(a_dt, "dt");

      // Write a general comment, the time and time step for the frame.
      a_show.saveComment(0, sPrintF(buffer, "Maxwell"));
      a_show.saveComment(1,
         sPrintF(buffer, "t=%5.2f, dt=%5.2e", a_time, a_dt));

      // Transfer the z component of the E field from m_em_vars_local into
      // m_ex_plot.
#ifdef USE_PPP
      RealArray ez_plot_local;
      getLocalArrayWithGhostBoundaries(m_ez_plot[0], ez_plot_local);
#else
      RealArray& ez_plot_local = m_ez_plot[0];
#endif
      Range R0(m_em_vars_local.getBase(0), m_em_vars_local.getBound(0));
      Range R1(m_em_vars_local.getBase(1), m_em_vars_local.getBound(1));
      Range R2(EZ, EZ);
      Range R3(0, 0);
      ez_plot_local(R0, R1, R3) = m_em_vars_local(R0, R1, R2);

      // Write the Ez plot.
      a_show.saveSolution(m_ez_plot);
      a_show.endFrame();
   }

   {
      // Set up the frame for Bx.
      a_show.setCurrentFrameSeries("Bx");
      a_show.startFrame();
      HDF_DataBase *dbp(a_show.getFrame());
      if (dbp == 0) {
         OV_ABORT("Can not get frame from a_show");
      }

      // Write some general plot info for the frame.
      dbp->put(a_time, "time");
      dbp->put(a_dt, "dt");

      // Write a general comment, the time and time step for the frame.
      a_show.saveComment(0, sPrintF(buffer, "Maxwell"));
      a_show.saveComment(1,
         sPrintF(buffer, "t=%5.2f, dt=%5.2e", a_time, a_dt));

      // Transfer the x component of the B field from m_em_vars_local into
      // m_ex_plot.
#ifdef USE_PPP
      RealArray bx_plot_local;
      getLocalArrayWithGhostBoundaries(m_bx_plot[0], bx_plot_local);
#else
      RealArray& bx_plot_local = m_bx_plot[0];
#endif
      Range R0(m_em_vars_local.getBase(0), m_em_vars_local.getBound(0));
      Range R1(m_em_vars_local.getBase(1), m_em_vars_local.getBound(1));
      Range R2(BX, BX);
      Range R3(0, 0);
      bx_plot_local(R0, R1, R3) = m_em_vars_local(R0, R1, R2);

      // Write the Bx plot.
      a_show.saveSolution(m_bx_plot);
      a_show.endFrame();
   }

   {
      // Set up the frame for By.
      a_show.setCurrentFrameSeries("By");
      a_show.startFrame();
      HDF_DataBase *dbp(a_show.getFrame());
      if (dbp == 0) {
         OV_ABORT("Can not get frame from a_show");
      }

      // Write some general plot info for the frame.
      dbp->put(a_time, "time");
      dbp->put(a_dt, "dt");

      // Write a general comment, the time and time step for the frame.
      a_show.saveComment(0, sPrintF(buffer, "Maxwell"));
      a_show.saveComment(1,
         sPrintF(buffer, "t=%5.2f, dt=%5.2e", a_time, a_dt));

      // Transfer the y component of the B field from m_em_vars_local into
      // m_ex_plot.
#ifdef USE_PPP
      RealArray by_plot_local;
      getLocalArrayWithGhostBoundaries(m_by_plot[0], by_plot_local);
#else
      RealArray& by_plot_local = m_by_plot[0];
#endif
      Range R0(m_em_vars_local.getBase(0), m_em_vars_local.getBound(0));
      Range R1(m_em_vars_local.getBase(1), m_em_vars_local.getBound(1));
      Range R2(BY, BY);
      Range R3(0, 0);
      by_plot_local(R0, R1, R3) = m_em_vars_local(R0, R1, R2);

      // Write the By plot.
      a_show.saveSolution(m_by_plot);
      a_show.endFrame();
   }

   {
      // Set up the frame for Bz.
      a_show.setCurrentFrameSeries("Bz");
      a_show.startFrame();
      HDF_DataBase *dbp(a_show.getFrame());
      if (dbp == 0) {
         OV_ABORT("Can not get frame from a_show");
      }

      // Write some general plot info for the frame.
      dbp->put(a_time, "time");
      dbp->put(a_dt, "dt");

      // Write a general comment, the time and time step for the frame.
      a_show.saveComment(0, sPrintF(buffer, "Maxwell"));
      a_show.saveComment(1,
         sPrintF(buffer, "t=%5.2f, dt=%5.2e", a_time, a_dt));

      // Transfer the z component of the B field from m_em_vars_local into
      // m_ex_plot.
#ifdef USE_PPP
      RealArray bz_plot_local;
      getLocalArrayWithGhostBoundaries(m_bz_plot[0], bz_plot_local);
#else
      RealArray& bz_plot_local = m_bz_plot[0];
#endif
      Range R0(m_em_vars_local.getBase(0), m_em_vars_local.getBound(0));
      Range R1(m_em_vars_local.getBase(1), m_em_vars_local.getBound(1));
      Range R2(BZ, BZ);
      Range R3(0, 0);
      bz_plot_local(R0, R1, R3) = m_em_vars_local(R0, R1, R2);

      // Write the Bz plot.
      a_show.saveSolution(m_bz_plot);
      a_show.endFrame();
   }


   {
      // Write the plot data for each species following the same pattern as
      // above.
      for (int is(0); is < m_num_kinetic_species; ++is) {
         sprintf(buffer, "%s vz", a_species_names[is].c_str());
         a_show.setCurrentFrameSeries(buffer);
         a_show.startFrame();
         HDF_DataBase *dbp(a_show.getFrame());
         if (dbp == 0) {
            OV_ABORT("Can not get frame from a_show");
         }

         dbp->put(a_time, "time");
         dbp->put(a_dt, "dt");
         a_show.saveComment(0, sPrintF(buffer, "Maxwell"));
         a_show.saveComment(1,
            sPrintF(buffer, "t=%5.2f, dt=%5.2e", a_time, a_dt));

#ifdef USE_PPP
         RealArray vz_plot_local;
         getLocalArrayWithGhostBoundaries(m_vz_plot[is][0], vz_plot_local);
#else
         RealArray& vz_plot_local = m_vz_plot[is][0];
#endif
         Range R0(m_vz_local[is].getBase(0), m_vz_local[is].getBound(0));
         Range R1(m_vz_local[is].getBase(1), m_vz_local[is].getBound(1));
         Range R2(0, 0);
         vz_plot_local(R0, R1, R2) = m_vz_local[is](R0, R1, R2);
         a_show.saveSolution(m_vz_plot[is]);
         a_show.endFrame();

         sprintf(buffer, "%s vx lo ke flux", a_species_names[is].c_str());
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
         a_show.saveSolution(m_ke_flux_vx_lo[is]);
         a_show.endFrame();

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
   }

   // We finished saving another plot cycle.
   ++a_saved_save;
}


void
Maxwell::accumulateSequences(
   RealArray&       a_sequences,
   const RealArray& a_probes,
   int              a_num_probes,
   int              a_saved_seq,
   int&             a_seq_idx)
{
   // These are the scalar time history quantities.  Initialize them to 0.
   real e_sum_tot(0.0);
   real e_max(0.0);
   real e_tot(0.0);
   real ex_max(0.0);
   real ey_max(0.0);
   real ez_max(0.0);

   real b_sum_tot(0.0);
   real b_max(0.0);
   real b_tot(0.0);
   real bx_max(0.0);
   real by_max(0.0);
   real bz_max(0.0);

   // Do the same to the probes' time history data.
   RealArray emProbe(TIME_HISTS_PER_PROBE, a_num_probes);
   RealArray emProbeTotal(TIME_HISTS_PER_PROBE, a_num_probes);
   for (int ix(0); ix < TIME_HISTS_PER_PROBE; ++ix) {
      for (int ip(0); ip < a_num_probes; ++ip) {
         emProbe(ix, ip) = 0.0;
      }
   }

   // Only processors onto which the Maxwell calculation is distributed nedd to
   // do this.
   if (isMaxwellProcessor()) {
      // Get a bunch of box info needed for indexing into the EM field arrays.
      const tbox::Box local_box(BoxOps::getLocalBox(m_em_vars_local));
      const tbox::Box& domain_box(m_domain->box());
      const tbox::Box global_box(tbox::Box::grow(domain_box, m_n_ghosts));
      const tbox::Box interior_box = BoxOps::getOperationalBox(m_em_vars_local,
         m_em_vars_global,
         domain_box,
         global_box);
      const int x_lo = interior_box.lower(0) - local_box.lower(0);
      const int x_hi = interior_box.upper(0) - local_box.lower(0);
      const int y_lo = interior_box.lower(1) - local_box.lower(1);
      const int y_hi = interior_box.upper(1) - local_box.lower(1);
      const int nx = local_box.upper(0) - local_box.lower(0) + 1;
      const int ny = local_box.upper(1) - local_box.lower(1) + 1;
      const int nxy = nx*ny;

      // Get some other misc quantities we'll need along the way.
      const tbox::IntVector N(m_domain->numberOfCells());
      const Array<double>& dx(m_domain->dx());
      const real area(dx[X1] * dx[X2]);
      const double* em_vars_ptr = m_em_vars_local.getDataPointer();

      // Compute the local parts of the EM field time histories.
      for (int i1(x_lo); i1 <= x_hi; ++i1) {
         for (int i2(y_lo); i2 <= y_hi; ++i2) {
            int idx = i2*nx+i1;
            const real ex = em_vars_ptr[idx+EX*nxy];
            const real ey = em_vars_ptr[idx+EY*nxy];
            const real ez = em_vars_ptr[idx+EZ*nxy];
            const real bx = em_vars_ptr[idx+BX*nxy];
            const real by = em_vars_ptr[idx+BY*nxy];
            const real bz = em_vars_ptr[idx+BZ*nxy];
            real tmp = ex*ex + ey*ey + ez*ez;
            real e_loc = sqrt(tmp);
            e_sum_tot += tmp;
            e_max = max(e_max, e_loc);
            e_tot += e_loc;
            ex_max = max(ex_max, abs(ex));
            ey_max = max(ey_max, abs(ey));
            ez_max = max(ez_max, abs(ez));

            tmp = bx*bx + by*by + bz*bz;
            real b_loc = sqrt(tmp);
            b_sum_tot += tmp;
            b_max = max(b_max, b_loc);
            b_tot += b_loc;
            bx_max = max(bx_max, abs(bx));
            by_max = max(by_max, abs(by));
            bz_max = max(bz_max, abs(bz));
         }
      }
      // These are integrated over configuration space.
      e_sum_tot *= area;
      e_tot     *= area;
      b_sum_tot *= area;
      b_tot     *= area;

      // Compute the local part of the probe time histories.
      for (int np(0); np < a_num_probes; ++np) {
         const int ip(int(floor(a_probes(0, np) * N[X1])));
         const int jp(int(floor(a_probes(1, np) * N[X2])));
         if (ip >= interior_box.lower(0) && ip <= interior_box.upper(0) &&
             jp >= interior_box.lower(1) && jp <= interior_box.upper(1)) {
            int idx = (jp-local_box.lower(1))*nx+(ip-local_box.lower(0));
            emProbe(0, np) = em_vars_ptr[idx+EX*nxy];
            emProbe(1, np) = em_vars_ptr[idx+EY*nxy];
            emProbe(2, np) = em_vars_ptr[idx+EZ*nxy];
            emProbe(3, np) = em_vars_ptr[idx+BX*nxy];
            emProbe(4, np) = em_vars_ptr[idx+BY*nxy];
            emProbe(5, np) = em_vars_ptr[idx+BZ*nxy];
            for (int ns(0); ns < m_num_kinetic_species; ++ns) {
               const double* vz_ptr = m_vz_local[ns].getDataPointer();
               emProbe(6+ns, np) = vz_ptr[idx];
            }                
         }
      }
   }

   // Now get the appropriate max/sum of the EM field time histories from all
   // the processors and store to the time history.
   a_sequences(a_saved_seq, a_seq_idx++) =
      ParallelUtility::getMaxValue(e_max, -1);
   a_sequences(a_saved_seq, a_seq_idx++) =
      ParallelUtility::getSum(e_tot, -1);
   a_sequences(a_saved_seq, a_seq_idx++) =
      ParallelUtility::getMaxValue(ex_max, -1);
   a_sequences(a_saved_seq, a_seq_idx++) =
      ParallelUtility::getMaxValue(ey_max, -1);
   a_sequences(a_saved_seq, a_seq_idx++) =
      ParallelUtility::getMaxValue(ez_max, -1);
   a_sequences(a_saved_seq, a_seq_idx++) =
      ParallelUtility::getSum(e_sum_tot, -1);

   a_sequences(a_saved_seq, a_seq_idx++) =
      ParallelUtility::getMaxValue(b_max, -1);
   a_sequences(a_saved_seq, a_seq_idx++) =
      ParallelUtility::getSum(b_tot, -1);
   a_sequences(a_saved_seq, a_seq_idx++) =
      ParallelUtility::getMaxValue(bx_max, -1);
   a_sequences(a_saved_seq, a_seq_idx++) =
      ParallelUtility::getMaxValue(by_max, -1);
   a_sequences(a_saved_seq, a_seq_idx++) =
      ParallelUtility::getMaxValue(bz_max, -1);
   a_sequences(a_saved_seq, a_seq_idx++) =
      ParallelUtility::getSum(b_sum_tot, -1);

   // Get the sum of the probe time histories from all the processors and store
   // to the time history.
   ParallelUtility::getSums(&(*emProbe.getDataPointer()),
      &(*emProbeTotal.getDataPointer()),
      TIME_HISTS_PER_PROBE * a_num_probes,
      -1);
   for (int ip(0); ip < a_num_probes; ++ip) {
      a_sequences(a_saved_seq, a_seq_idx++) = emProbeTotal(0, ip);
      a_sequences(a_saved_seq, a_seq_idx++) = emProbeTotal(1, ip);
      a_sequences(a_saved_seq, a_seq_idx++) = emProbeTotal(2, ip);
      a_sequences(a_saved_seq, a_seq_idx++) = emProbeTotal(3, ip);
      a_sequences(a_saved_seq, a_seq_idx++) = emProbeTotal(4, ip);
      a_sequences(a_saved_seq, a_seq_idx++) = emProbeTotal(5, ip);
      for (int is(0); is < m_num_kinetic_species;  ++is) {
         a_sequences(a_saved_seq, a_seq_idx++) = emProbeTotal(6+is, ip);
      }
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
Maxwell::getFromRestart(
   const HDF_DataBase& a_db)
{
   // Get the subdatabase containing this objects data.
   HDF_DataBase sub_db;
   a_db.locate(sub_db, "Maxwell");

   // Get the electromagnetic variables from the subdatabase.
   char buffer[80];
   sub_db.getDistributed(m_em_vars_global, "EMVars");
   for (int i = 0; i < m_num_kinetic_species; ++i) {
      sprintf(buffer, "vz%i", i);
      sub_db.getDistributed(m_vz_global[i], buffer);
   }

#ifdef USE_PPP
   getLocalArrayWithGhostBoundaries(m_em_vars_global, m_em_vars_local);
   for (int i = 0; i < m_num_kinetic_species; ++i) {
      getLocalArrayWithGhostBoundaries(m_vz_global[i], m_vz_local[i]);
   }
#endif

   // As I recall, Overture does not read either any of the ghosts or the extra
   // ghosts on a physical boundary so we need to zero out that data.
   FORT_ZERO_GHOST_2D(*m_em_vars_local.getDataPointer(),
      BOX2D_TO_FORT(m_interior_box),
      BOX2D_TO_FORT(m_local_box),
      NUM_EM_VARS);
   for (int i = 0; i < m_num_kinetic_species; ++i) {
      FORT_ZERO_GHOST_2D(*m_vz_global[i].getDataPointer(),
         BOX2D_TO_FORT(m_interior_box),
         BOX2D_TO_FORT(m_local_box),
         1);
   }
}


void
Maxwell::putToRestart(
   HDF_DataBase& a_db,
   real a_time)
{
   NULL_USE(a_time);

   // We write the number of ghosts only so that the post processor knows it.
   a_db.put(m_n_ghosts[0], "nGhost");

   // Make a subdatabase.
   HDF_DataBase sub_db;
   a_db.create(sub_db, "Maxwell", "directory");

   // Put the electromagnetic variables in the subdatabase.
   char buffer[80];
   sub_db.putDistributed(m_em_vars_global, "EMVars");
   for (int i = 0; i < m_num_kinetic_species; ++i) {
      sprintf(buffer, "vz%i", i);
      sub_db.putDistributed(m_vz_global[i], buffer);
   }
}


void
Maxwell::readParticles(
   ParmParse& a_pp)
{
   // If tracking particles are specified, read them in.
   if (a_pp.contains("tracking_particle_file")) {
      herr_t errf;
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("Tracking Particles");
      aString tracking_particle_file;
      a_pp.get("tracking_particle_file", tracking_particle_file);

      // Open tracking particle file.
      hid_t file_id =
         H5Fopen(tracking_particle_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      if (file_id < 0) {
         OV_ABORT("Unable to open tracking particle file.");
      }

      // Open and read the dataset containing the number of tracking particles.
      int num_tracking_particles;
      hid_t dataset_id;
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
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
                        &num_tracking_particles);
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
      m_problem_num_tracking_particles = num_tracking_particles;

      // If this Maxwell object is not on a Maxwell processor we don't read the
      // particles themselves and we're done.
      if (!isMaxwellProcessor()) {
         // Close the file.
         errf = H5Fclose(file_id);
         if (errf < 0) {
            OV_ABORT("Can not close particle file.");
         }
         timers->stopTimer("Tracking Particles");
         return;
      }

      // Open and read the dataset containing the initial phase space location
      // of all the tracking particles.
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
         int data_size = 7*num_tracking_particles;
         hid_t dspace = H5Dget_space(dataset_id);
         if (dspace < 0) {
            OV_ABORT("Can not get dataspace for dataset \"particle_data\".");
         }
         hsize_t nsel = H5Sget_select_npoints(dspace);
         if (static_cast<int>(nsel) != data_size) {
            OV_ABORT("Incorrect amount of tracking particle data.");
         }

         // Now that everything checks out, read the dataset.
         double* tracking_particle_data = new double [data_size];
         errf = H5Dread(dataset_id,
                        H5T_NATIVE_DOUBLE,
                        H5S_ALL,
                        H5S_ALL,
                        H5P_DEFAULT,
                        tracking_particle_data);
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
         for (int i = 0; i < num_tracking_particles; ++i) {
            double start_tracking = tracking_particle_data[idx++];
            double charge = tracking_particle_data[idx++];
            double mass = tracking_particle_data[idx++];
            double xinit = tracking_particle_data[idx++];
            double yinit = tracking_particle_data[idx++];
            double vxinit = tracking_particle_data[idx++];
            double vyinit = tracking_particle_data[idx++];
            double noise_source_weight = 0.0;
            // Let's make sure that the initial location of the particle is
            // inside the physical domain.
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
            Particle particle(start_tracking,
                              charge,
                              mass,
                              xinit,
                              yinit,
                              vxinit,
                              vyinit,
                              noise_source_weight);
            m_tracking_particles.push_back(particle);
         }
      }

      // Close the file.
      errf = H5Fclose(file_id);
      if (errf < 0) {
         OV_ABORT("Can not close particle tracking file.");
      }
      timers->stopTimer("Tracking Particles");
   }
   else {
      m_problem_num_tracking_particles = 0;
   }
}


// Private member functions.


void
Maxwell::parseParameters(
   ParmParse& a_pp)
{
   // See if the user has specified how many processors to use for the Maxwell
   // solve.  This shouldn't currently be used because we just don't need > 1
   // Maxwell processor and if we did, there's plenty of assumptions that would
   // probably break.
   a_pp.query("number_of_processors", m_number_of_procs);

   // See if these optional inputs have been specified.
   a_pp.query("light_speed", m_c);
   a_pp.query("avWeak", m_avWeak);
   a_pp.query("avStrong", m_avStrong);

   // See if there are any current drivers.
   a_pp.query("num_current_drivers", m_num_current_drivers);
   m_current_drivers.resize(m_num_current_drivers,
      tbox::Pointer<CurrentDriver>(0));

   // All of this is currently meaningless as we do not solve for the
   // electrostatic E field.
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
Maxwell::printParameters() const
{
   // Print this Maxwell's parameters and the parameters of the entities that it
   // holds and are not accessible elsewhere like the current drivers and EM
   // initial conditions.
   printF("\n#*#*# Maxwell #*#*#\n");
   printF("  Number of Processsors = %i\n", m_number_of_procs);
   printF("  Light Speed = %f\n", m_c);
   printF("  avWeak = %f\n", m_avWeak);
   printF("  avStrong = %f\n", m_avStrong);
   printF("\n  num current drivers = %i\n", m_num_current_drivers);
   for (int i = 0; i < m_num_current_drivers; ++i) {
      m_current_drivers[i]->printParameters();
   }
   for (int i = 0; i < 2; ++i) {
      m_em_ics[i]->printParameters();
   }
   for (int i = 0; i < m_num_kinetic_species; ++i) {
      m_vel_ics[i]->printParameters();
   }
}


void
Maxwell::copy(
   const Maxwell& a_rhs,
   bool a_deep_copy)
{
   if (m_dim != a_rhs.m_dim) {
      OV_ABORT("Attempt to copy incongruent Maxwells");
   }

   // If these 2 Maxwells are different, then copy their internals.
   if (&a_rhs != this) {
      m_domain = a_rhs.m_domain;
      m_number_of_procs = a_rhs.m_number_of_procs;
      m_processor_range = a_rhs.m_processor_range;
      m_partition_defined = a_rhs.m_partition_defined;
      m_is_maxwell_processor = a_rhs.m_is_maxwell_processor;
      m_comm = a_rhs.m_comm;
      m_partition = a_rhs.m_partition;
      m_c = a_rhs.m_c;
      m_avWeak = a_rhs.m_avWeak;
      m_avStrong = a_rhs.m_avStrong;
      m_solution_order = a_rhs.m_solution_order;
      m_num_kinetic_species = a_rhs.m_num_kinetic_species;
      m_problem_num_tracking_particles = a_rhs.m_problem_num_tracking_particles;
      m_tracking_particles = a_rhs.m_tracking_particles;

      if (m_local_box.lower(0) != a_rhs.m_local_box.lower(0) ||
          m_local_box.upper(0) != a_rhs.m_local_box.upper(0) ||
          m_local_box.lower(1) != a_rhs.m_local_box.lower(1) ||
          m_local_box.upper(1) != a_rhs.m_local_box.upper(1)) {
         m_em_vars_global.redim(0);
         m_em_vars_global.partition(a_rhs.m_em_vars_global.getPartition());
         m_em_vars_global.redim(a_rhs.m_em_vars_global);
         for (int i = 0; i < m_num_kinetic_species; ++i) {
            m_vz_global[i].redim(0);
            m_vz_global[i].partition(a_rhs.m_vz_global[i].getPartition());
            m_vz_global[i].redim(a_rhs.m_vz_global[i]);
         }
      }
      if (a_deep_copy && isMaxwellProcessor()) {
         Index dst[4];
         Index* src = dst;
         for (int i = 0; i < m_dim+1; ++i) {
            dst[i] = Range(a_rhs.m_em_vars_global.getBase(i),
                           a_rhs.m_em_vars_global.getBound(i));
         }
         CopyArray::copyArray(m_em_vars_global,
            dst,
            a_rhs.m_em_vars_global,
            src);

         dst[m_dim] = Range(0, 0);
         for (int i = 0; i < m_num_kinetic_species; ++i) {
            CopyArray::copyArray(m_vz_global[i],
               dst,
               a_rhs.m_vz_global[i],
               src);
         }
      }

#ifdef USE_PPP
      getLocalArrayWithGhostBoundaries(m_em_vars_global, m_em_vars_local);
      for (int i = 0; i < m_num_kinetic_species; ++i) {
         getLocalArrayWithGhostBoundaries(m_vz_global[i], m_vz_local[i]);
      }
#endif

      m_n_ghosts = a_rhs.m_n_ghosts;
      m_local_box = a_rhs.m_local_box;
      m_interior_box = a_rhs.m_interior_box;
      m_global_box = a_rhs.m_global_box;
      m_solver = a_rhs.m_solver;
      m_rho_solver = a_rhs.m_rho_solver;
      m_phi_solver = a_rhs.m_phi_solver;
      m_phi = a_rhs.m_phi;
      m_ex_plot = a_rhs.m_ex_plot;
      m_ey_plot = a_rhs.m_ey_plot;
      m_ez_plot = a_rhs.m_ez_plot;
      m_bx_plot = a_rhs.m_bx_plot;
      m_by_plot = a_rhs.m_by_plot;
      m_bz_plot = a_rhs.m_bz_plot;
      m_vz_plot = a_rhs.m_vz_plot;
      m_solver_method = a_rhs.m_solver_method;
      m_tol = a_rhs.m_tol;
      m_stencil_size = a_rhs.m_stencil_size;
      m_stencil_width = a_rhs.m_stencil_width;
      m_max_iterations = a_rhs.m_max_iterations;
      m_current_drivers = a_rhs.m_current_drivers;
      m_num_current_drivers = a_rhs.m_num_current_drivers;
      m_em_ics = a_rhs.m_em_ics;
      m_vel_ics = a_rhs.m_vel_ics;
      buildInterpolator();
   }
}


void
Maxwell::defineSolver()
{
   // Lifted directly from Poisson and currently irrelevant.
   Mapping* poisson_mapping = createMapping(1);

   MappedGrid::setMinimumNumberOfDistributedGhostLines(m_n_ghosts.max());
   MappedGrid poisson_mg(*poisson_mapping);
   setupMappedGrid(poisson_mg, m_stencil_width);

   setupSolver(poisson_mg);
}

} // end namespace Loki
