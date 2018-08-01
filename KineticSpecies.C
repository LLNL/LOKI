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
#include "KineticSpecies.H"

#include "tbox/IntVector.H"
#include "ICFactory.H"
#include "TZSourceFactory.H"
#include "BoxOps.H"
#include "ElectricFieldDriverFactory.H"
#include "CollisionOperatorFactory.H"
#include "Maxwell.H"
#include "Poisson.H"
#include "RestartManager.H"

#include "IntegralOp.H"

#include <sstream>

namespace Loki {

static
void
allocateFaceArray(
   RealArray& a_array,
   const tbox::Box& a_box,
   int a_dir);


KineticSpecies::KineticSpecies(
   HDF_DataBase& a_db,
   int a_species_num,
   int a_number_of_species,
   int a_spatial_solution_order,
   int a_temporal_solution_order,
   bool a_use_new_bcs,
   const aString& a_name)
   : m_pdim(PDIM),
     m_cdim(CDIM),
     m_name("undefined"),
     m_species_index(a_species_num-1),
     m_number_of_species(a_number_of_species),
     m_mass(-1.0),
     m_charge(0.0),
#ifndef USE_PPP
     m_dist_func_local(m_dist_func_global),
#endif
     m_spatial_solution_order(a_spatial_solution_order),
     m_temporal_solution_order(a_temporal_solution_order),
     m_n_ghosts(m_pdim),
     m_stencil_width(a_spatial_solution_order+1),
     m_local_box(m_pdim),
     m_interior_box(m_pdim),
     m_global_box(m_pdim),
     m_accel_box(m_pdim),
     m_em_vars_box(m_pdim),
     m_vz_box(m_pdim),
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
     m_vflowx(0.0),
     m_vflowy(0.0),
     m_problem_has_particles(false),
     m_use_new_bcs(a_use_new_bcs)
{
   // Set number of ghosts based of order of solution.
   int num_ghosts;
   if (m_spatial_solution_order == 4) {
      num_ghosts = 2;
   }
   else {
      num_ghosts = 3;
   }
   for (int i = 0; i < m_pdim; ++i) {
      m_n_ghosts[i] = num_ghosts;
   }

   // find subdatabase with the name of this distribution
   m_name = a_name;
   HDF_DataBase sub_db;
   a_db.locate(sub_db, m_name);

   // get the pdim, cdim from database
   int tmp_pdim, tmp_cdim;
   sub_db.get(tmp_pdim, "pdim");
   sub_db.get(tmp_cdim, "cdim");
   tbox::Dimension pdim(static_cast<unsigned short>(tmp_pdim));
   tbox::Dimension cdim(static_cast<unsigned short>(tmp_cdim));

   // get the mass and charge from database
   sub_db.get(m_mass, "mass");
   sub_db.get(m_charge, "charge");

   std::cout << "In KineticSpecies db constructor, creating new ProblemDomain"
             << std::endl;
   m_domain = new ProblemDomain(pdim, sub_db);
   std::cout << "Returned from creating new ProblemDomain" << std::endl;

   // read restart distribution from database and now that the local data is
   // defined, get the local array
   sub_db.getDistributed(m_dist_func_global, "distribution");
   std::cout << "Returned from getDistributed (K.S. constructor)" << std::endl;
#ifdef USE_PPP
   getLocalArrayWithGhostBoundaries(m_dist_func_global, m_dist_func_local);
#endif
}


KineticSpecies::KineticSpecies(
   const tbox::Pointer<ProblemDomain>& a_cfg_domain,
   ParmParse& a_pp,
   int a_species_num,
   int a_number_of_species,
   int a_spatial_solution_order,
   int a_temporal_solution_order,
   bool a_use_new_bcs,
   bool a_do_maxwell)
   : m_pdim(PDIM),
     m_cdim(CDIM),
     m_name("undefined"),
     m_species_index(a_species_num-1),
     m_number_of_species(a_number_of_species),
     m_mass(-1.0),
     m_charge(0.0),
#ifndef USE_PPP
     m_dist_func_local(m_dist_func_global),
#endif
     m_spatial_solution_order(a_spatial_solution_order),
     m_temporal_solution_order(a_temporal_solution_order),
     m_n_ghosts(m_pdim),
     m_stencil_width(a_spatial_solution_order+1),
     m_local_box(m_pdim),
     m_interior_box(m_pdim),
     m_global_box(m_pdim),
     m_accel_box(m_pdim),
     m_em_vars_box(m_pdim),
     m_vz_box(m_pdim),
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
     m_vflowx(0.0),
     m_vflowy(0.0),
     m_problem_has_particles(false),
     m_use_new_bcs(a_use_new_bcs)
{
   // Set number of ghosts based of order of solution.
   int num_ghosts;
   if (m_spatial_solution_order == 4) {
      num_ghosts = 2;
   }
   else {
      num_ghosts = 3;
   }
   for (int i = 0; i < m_pdim; ++i) {
      m_n_ghosts[i] = num_ghosts;
   }

   // Read all the user input for this species and construct much of its data
   // members.
   parseParameters(a_pp, a_cfg_domain);

   // The following code could probably all be placed into parseParameters but
   // isn't because it involves parsing separate sub-databases of this species.

   // Get the sub-database for this species' initial condition and construct the
   // initial condition object.
   char buffer[100];
   sprintf(buffer, "kinetic_species.%i.ic", a_species_num);
   ParmParse ic_pp(buffer);
   m_initial_condition = ICFactory::create(ic_pp, m_vflowx, m_vflowy);

   // Get the sub-database for this species' twilight zone and construct the
   // twilight zone object if one is specified.
   sprintf(buffer, "kinetic_species.%i.tz", a_species_num);
   ParmParse tz_pp(buffer);
   m_tz_source = TZSourceFactory::create(tz_pp);

   // For each of this species' external E field drivers get its sub-database
   // and construct the field driver object.
   for (int i = 0; i < m_num_external_drivers; ++i) {
      if (m_old_driver_syntax) {
         sprintf(buffer, "kinetic_species.%i.external_driver", a_species_num);
      }
      else {
         sprintf(buffer,
            "kinetic_species.%i.external_driver.%i",
            a_species_num,
            i+1);
      }
      ParmParse efdf_pp(buffer);
      m_ef_drivers[i] = ElectricFieldDriverFactory::create(efdf_pp, i+1);
   }

   // For each of this species' collision operators get its sub-database and
   // construct the collision operator object.
   for (int i = 0; i < m_num_collision_operators; ++i) {
      sprintf(buffer,
         "kinetic_species.%i.collision_operator.%i",
         a_species_num,
         i+1);
      ParmParse coll_op_pp(buffer);
      m_collision_operators[i] = CollisionOperatorFactory::create(coll_op_pp,
         m_spatial_solution_order);
   }

   // Get the sub-database for this species' Krook layer and construct it if
   // one is specified.
   sprintf(buffer, "kinetic_species.%i.krook", a_species_num);
   ParmParse krook_pp(buffer);
   m_krook_layer = new KrookLayer(m_cdim, krook_pp, *m_domain);

   // The KineticSpecies write restart data so register this object with the
   // RestartManager which will use the putToRestart/getFromRestart callbacks to
   // get the restart data written/read.
   RestartManager* restart_manager(RestartManager::getManager());
   restart_manager->registerRestart(this);

   printParameters();
}


KineticSpecies::KineticSpecies(
   const KineticSpecies& a_other,
   bool a_deep_copy)
   : Load(),
     Serializable(),
     m_pdim(a_other.m_pdim),
     m_cdim(a_other.m_cdim),
#ifndef USE_PPP
     m_dist_func_local(m_dist_func_global),
#endif
     m_n_ghosts(m_pdim, 3),
     m_local_box(m_pdim),
     m_interior_box(m_pdim),
     m_global_box(m_pdim),
     m_accel_box(m_pdim),
     m_em_vars_box(m_pdim),
     m_vz_box(m_pdim)
{
   // Copy all the contents and build new schedules.
   copy(a_other, a_deep_copy);
   defineExtEfieldContractionSchedule();
   defineChargeDensityReductionSchedule();
   defineKineticEnergyReductionSchedule();
}


KineticSpecies::~KineticSpecies()
{
}


void
KineticSpecies::addData(
   const KineticSpecies& a_rhs,
   real a_factor)
{
   // Check that we're copying between similar species and add a_factor times
   // a_rhs' distribution function to this species' distribution function.
   if (conformsTo(a_rhs, false)) {
      tbox::Box intersect_box(m_interior_box * a_rhs.m_interior_box);
      FORT_XPBY_4D(*m_dist_func_local.getDataPointer(),
         BOX4D_TO_FORT(m_local_box),
         *a_rhs.m_dist_func_local.getDataPointer(),
         BOX4D_TO_FORT(a_rhs.m_local_box),
         BOX4D_TO_FORT(intersect_box),
         a_factor);
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
       m_processor_range != a_rhs.m_processor_range ||
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
   const KineticSpecies& a_rhs,
   bool a_deep_copy)
{
   if ((m_pdim != a_rhs.m_pdim) || (m_cdim != a_rhs.m_cdim)) {
      OV_ABORT("Attemtpt to copy incongruent species!");
   }

   // If the 2 species are different, then copy their internals.
   if (&a_rhs != this) {
      m_name = a_rhs.m_name;
      m_species_index = a_rhs.m_species_index;
      m_number_of_species = a_rhs.m_number_of_species;
      m_mass = a_rhs.m_mass;
      m_charge = a_rhs.m_charge;
      m_domain = a_rhs.m_domain;
      m_initial_condition = a_rhs.m_initial_condition;
      m_tz_source = a_rhs.m_tz_source;
      m_krook_layer = a_rhs.m_krook_layer;
      m_number_of_procs = a_rhs.m_number_of_procs;
      m_processor_range = a_rhs.m_processor_range;
      m_fixed_number_of_procs = a_rhs.m_fixed_number_of_procs;
      m_partition_defined = a_rhs.m_partition_defined;
      m_comm = a_rhs.m_comm;
      m_num_external_drivers = a_rhs.m_num_external_drivers;
      m_old_driver_syntax = a_rhs.m_old_driver_syntax;
      m_ef_drivers = a_rhs.m_ef_drivers;
      m_do_maxwell = a_rhs.m_do_maxwell;
      m_collision_operators = a_rhs.m_collision_operators;
      m_num_collision_operators = a_rhs.m_num_collision_operators;
      m_vflowx = a_rhs.m_vflowx;
      m_vflowy = a_rhs.m_vflowy;
      m_problem_has_particles = a_rhs.m_problem_has_particles;
      m_use_new_bcs = a_rhs.m_use_new_bcs;

      if (m_local_box.lower(0) != a_rhs.m_local_box.lower(0) ||
          m_local_box.upper(0) != a_rhs.m_local_box.upper(0) ||
          m_local_box.lower(1) != a_rhs.m_local_box.lower(1) ||
          m_local_box.upper(1) != a_rhs.m_local_box.upper(1) ||
          m_local_box.lower(2) != a_rhs.m_local_box.lower(2) ||
          m_local_box.upper(2) != a_rhs.m_local_box.upper(2) ||
          m_local_box.lower(3) != a_rhs.m_local_box.lower(3) ||
          m_local_box.upper(3) != a_rhs.m_local_box.upper(3)) {
         // this is probably super costly and is only OK if this is done only
         // once the first time this code is run through.
         m_dist_func_global.redim(0);
         m_dist_func_global.partition(a_rhs.m_dist_func_global.getPartition());
         m_dist_func_global.redim(a_rhs.m_dist_func_global);
      }
      if (a_deep_copy) {
         Index dst[4];
         Index* src = dst;
         for (int i(0); i < m_pdim; ++i) {
            dst[i] = Range(a_rhs.m_dist_func_global.getBase(i),
                           a_rhs.m_dist_func_global.getBound(i));
         }
         CopyArray::copyArray(m_dist_func_global,
            dst,
            a_rhs.m_dist_func_global,
            src);
      }

#ifdef USE_PPP
      getLocalArrayWithGhostBoundaries(m_dist_func_global, m_dist_func_local);
#endif

      m_spatial_solution_order = a_rhs.m_spatial_solution_order;
      m_temporal_solution_order = a_rhs.m_temporal_solution_order;
      m_n_ghosts = a_rhs.m_n_ghosts;
      m_stencil_width = a_rhs.m_stencil_width;
      m_local_box = a_rhs.m_local_box;
      m_interior_box = a_rhs.m_interior_box;
      m_global_box = a_rhs.m_global_box;
      m_accel_box = a_rhs.m_accel_box;
      m_em_vars_box = a_rhs.m_em_vars_box;
      m_vz_box = a_rhs.m_vz_box;

      m_efield_expansion_schedule = a_rhs.m_efield_expansion_schedule;
      m_em_expansion_schedule = a_rhs.m_em_expansion_schedule;
      m_vz_expansion_schedule = a_rhs.m_vz_expansion_schedule;

      m_lambda_max.resize(a_rhs.m_lambda_max.size());
      m_u_face.resize(a_rhs.m_u_face.size());
      m_vel_face.resize(a_rhs.m_vel_face.size());
      m_flux.resize(a_rhs.m_flux.size());
      for (int i(0); i < static_cast<int>(a_rhs.m_lambda_max.size()); ++i) {
         m_lambda_max[i] = a_rhs.m_lambda_max[i];
         m_u_face[i] = a_rhs.m_u_face[i];
         m_vel_face[i] = a_rhs.m_vel_face[i];
         m_flux[i] = a_rhs.m_flux[i];
      }

      m_accel = a_rhs.m_accel;
      if (m_problem_has_particles) {
         m_ext_efield_local = a_rhs.m_ext_efield_local;
      }
      m_em_vars = a_rhs.m_em_vars;
      m_vz = a_rhs.m_vz;
   }
}


void
KineticSpecies::copySolnData(
   const KineticSpecies& a_rhs)
{
   if ((m_pdim != a_rhs.m_pdim) || (m_cdim != a_rhs.m_cdim)) {
      OV_ABORT("Attemtpt to copy incongruent species!");
   }

   // If the 2 species are different, then copy the distribution function.
   if (&a_rhs != this) {
      // I would like to replace these with a check that the shape and
      // partition are the same
      if (m_local_box.lower(0) != a_rhs.m_local_box.lower(0) ||
          m_local_box.upper(0) != a_rhs.m_local_box.upper(0) ||
          m_local_box.lower(1) != a_rhs.m_local_box.lower(1) ||
          m_local_box.upper(1) != a_rhs.m_local_box.upper(1) ||
          m_local_box.lower(2) != a_rhs.m_local_box.lower(2) ||
          m_local_box.upper(2) != a_rhs.m_local_box.upper(2) ||
          m_local_box.lower(3) != a_rhs.m_local_box.lower(3) ||
          m_local_box.upper(3) != a_rhs.m_local_box.upper(3)) {
         m_dist_func_global.redim(0);
         m_dist_func_global.partition(a_rhs.m_dist_func_global.getPartition());
         m_dist_func_global.redim(a_rhs.m_dist_func_global);

         Index dst[4];
         Index* src = dst;
         for (int i(0); i < m_pdim; ++i) {
            dst[i] = Range(a_rhs.m_dist_func_global.getBase(i),
                           a_rhs.m_dist_func_global.getBound(i));
         }
         CopyArray::copyArray(m_dist_func_global,
            dst,
            a_rhs.m_dist_func_global,
            src);

         // not sure if this is needed after the reshape
#ifdef USE_PPP
         getLocalArrayWithGhostBoundaries(m_dist_func_global, m_dist_func_local);
#endif
      }
      else {
         m_dist_func_local = a_rhs.m_dist_func_local;
      }
   }
}


void
KineticSpecies::printParameters() const
{
   // Print this species parameters and the parameters of the entities that it
   // holds and are not accessible elsewhere like the external E field drivers,
   // collision operators, initial conditions, and twilight zones.
   printF("\n#*#*# Kinetic Species %s #*#*#\n", m_name.c_str());
   printF("  species index:   %d\n", m_species_index);
   printF("  mass:            %e\n", m_mass);
   printF("  charge:          %e\n", m_charge);
   printF("  initial x vflow: %e\n", m_vflowx);
   printF("  initial y vflow: %e\n", m_vflowy);
   m_domain->printParameters();
   m_krook_layer->printParameters();
   printF("\n  num external drivers = %i\n", m_num_external_drivers);
   for (int i = 0; i < m_num_external_drivers; ++i) {
      m_ef_drivers[i]->printParameters();
   }
   printF("\n  num collision operators = %i\n", m_num_collision_operators);
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
   const Range& a_range,
   const MPI_Comm& a_comm)
{
   m_comm = a_comm;
   m_processor_range = a_range;
   m_number_of_procs = a_range.length();
   m_partition_defined = true;

   // Overture is compiled such that each parallel array has MAX_ARRAY_DIMENSION
   // dimensions.  So we need to tell the partition how to partition all of
   // these dimensions.  Naturally we only care about the first m_pdim of them.
   // The others are ignored.
   Partitioning_Type partition;
   partition.SpecifyProcessorRange(m_processor_range);
   partition.SpecifyDecompositionAxes(m_pdim);

   for (int dir(0); dir < m_pdim; ++dir) {
      partition.partitionAlongAxis(dir, true, m_n_ghosts[dir]);
   }
   for (int dir(m_pdim); dir < MAX_ARRAY_DIMENSION; ++dir) {
      partition.partitionAlongAxis(dir, false, 0);
   }

   // Partition the global array and add the ghosts to the global box.  This
   // gives the range of each dimension which is needed in order to redimension
   // it.  I think that Overture's term "partition" is not accurate.  I think
   // that the partition call just defines which processors the global array
   // lives on and how it CAN be partitioned among those processors.  HOW it is
   // partitioned is determined by the redim call.  It can't very well be
   // partitioned unless its extent is known.
   m_dist_func_global.partition(partition);
   m_global_box = m_domain->box();
   m_global_box.grow(m_n_ghosts);
   if (m_pdim == tbox::Dimension(4)) {
      m_dist_func_global.redim(BoxOps::range(m_global_box, X1),
         BoxOps::range(m_global_box, X2),
         BoxOps::range(m_global_box, V1),
         BoxOps::range(m_global_box, V2));
   }
   else {
      OV_ABORT("Not implemented for phase D!=4!");
   }

   // Now that the global array has been divvied up we can get the local array.
#ifdef USE_PPP
   getLocalArrayWithGhostBoundaries(m_dist_func_global, m_dist_func_local);
#endif

   // Figure out the local and interior boxes.  Recall that all KineticSpecies
   // exist on all processors but in order to load balance things each species
   // is only defined and therefore "active" on a range of processors.
   int my_id(Communication_Manager::My_Process_Number);
   my_id = std::max(0, my_id);
   if (isInRange(my_id)) {
      m_local_box = BoxOps::getLocalBox(m_dist_func_local);
      m_interior_box = BoxOps::getOperationalBox(m_dist_func_local,
                                                 m_dist_func_global,
                                                 m_domain->box(),
                                                 m_global_box);
   }
   else {
      for (int dir(0); dir < m_pdim; ++dir) {
         m_local_box.lower(dir) = 0;
         m_local_box.upper(dir) = -1;
         m_interior_box.lower(dir) = 0;
         m_interior_box.upper(dir) = -1;
      }
   }

   // Check to make sure decomposition makes sense; we need at least
   // m_stencil_width interior points in each direction
   if (isInRange(my_id)) {
      tbox::IntVector npts(m_interior_box.upper() - m_interior_box.lower() + 1);
      for (int dir(0); dir < m_pdim; ++dir) {
         if (npts[dir] < m_stencil_width) {
            char msg[80];
            sprintf(msg,
                    "Too few interior points in decomposition in direction %d",
                    dir);
            OV_ABORT(msg);
         }
      }
   }

   // Now that we know the various boxes we can define the schedules and local
   // arrays needed by this species.
   defineExtEfieldContractionSchedule();
   defineChargeDensityReductionSchedule();
   defineKineticEnergyReductionSchedule();

   allocateLocalAuxArrays();

   int config_space_id =
      m_interior_box.lower(X2)*m_domain->box().numberCells(X1) +
      m_interior_box.lower(X1);
   if (m_do_maxwell) {
      m_em_expansion_schedule = new ExpansionSchedule(m_em_vars_box,
         m_number_of_species,
         config_space_id,
         m_comm);
      m_vz_expansion_schedule = new ExpansionSchedule(m_vz_box,
         m_number_of_species,
         config_space_id,
         m_comm);
   }
   else {
      m_efield_expansion_schedule = new ExpansionSchedule(m_accel_box,
         m_number_of_species,
         config_space_id,
         m_comm);
   }

   // Now that we have a communicator and processor range give that to the
   // collision operators.
   for (int i = 0; i < m_num_collision_operators; ++i) {
      m_collision_operators[i]->setCommunicationInfo(m_processor_range,
         m_comm);
   }

   // Now that we have the boxes and local arrays we can initialize the
   // velocities on each face.
   initializeVelocity();
}


bool
KineticSpecies::isInRange(
   int a_proc_id) const
{
   // Returns true if this species is partitioned onto this processor.
   return ((a_proc_id >= m_processor_range.getBase()) &&
           (a_proc_id <= m_processor_range.getBound()));
}


void
KineticSpecies::printDecomposition() const
{
   // This function is only valid if we actually know the decomposition.
   if (m_partition_defined) {
      // Print some basic decomposition info.
      printF("  Kinetic Species \"%s\" processor(s):  [%d,%d]\n",
             m_name.c_str(),
             m_processor_range.getBase(),
             m_processor_range.getBound());

      // Now do some sanity checking.
      // A prime number of processors is almost certainly a bad idea.
      int num_procs =
         m_processor_range.getBound() - m_processor_range.getBase() + 1;
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
      else {
         num_procs_prime = true;
      }
      if (num_procs_prime) {
         printF("  Kinetic Species \"%s\" is partitioned across %d processors\n"
                "  which is prime and most likely not what you want\n",
                m_name.c_str(), num_procs);
      }

      // Give the user some idea about how evenly the different dimensions are
      // subdivided.
      int loc_zones[4], glob_zones[4];
      int my_id(Communication_Manager::My_Process_Number);
      my_id = std::max(0, my_id);
      if (isInRange(my_id)) {
         loc_zones[0] = m_local_box.upper(0) - m_local_box.lower(0) + 1;
         loc_zones[1] = m_local_box.upper(1) - m_local_box.lower(1) + 1;
         loc_zones[2] = m_local_box.upper(2) - m_local_box.lower(2) + 1;
         loc_zones[3] = m_local_box.upper(3) - m_local_box.lower(3) + 1;
      }
      else {
         loc_zones[0] = INT_MIN;
         loc_zones[1] = INT_MIN;
         loc_zones[2] = INT_MIN;
         loc_zones[3] = INT_MIN;
      }
      MPI_Reduce(&loc_zones[0], &glob_zones[0], 4,
                 MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
      printF("  Kinetic Species \"%s\" maximum dimensions:  [%d,%d,%d,%d]\n",
             m_name.c_str(),
             glob_zones[0], glob_zones[1],
             glob_zones[2], glob_zones[3]);
      if (!isInRange(my_id)) {
         loc_zones[0] = INT_MAX;
         loc_zones[1] = INT_MAX;
         loc_zones[2] = INT_MAX;
         loc_zones[3] = INT_MAX;
      }
      MPI_Reduce(&loc_zones[0], &glob_zones[0], 4,
                 MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
      printF("  Kinetic Species \"%s\" minimum dimensions:  [%d,%d,%d,%d]\n",
             m_name.c_str(),
             glob_zones[0], glob_zones[1], 
             glob_zones[2], glob_zones[3]);
   }
}


real
KineticSpecies::computeDt()
{
   // local max lambdas are found every time the velocity is recomputed
   std::vector<real> lambda_max(m_pdim);
   ParallelUtility::getMaxValues(m_lambda_max, lambda_max, m_comm);

   real reLam(0.0);
   real imLam(0.0);
   real pi = 4.0*atan(1.0);

   // eigenvalue of Vlasov operator
   for (int dir(X1); dir < m_pdim; ++dir) {
      imLam += pi * lambda_max[dir] / m_domain->dx(dir);
   }
   // eigenvalue of collision operator
   real dv = std::min(m_domain->dx(V1), m_domain->dx(V2));
   for (int i = 0; i < m_num_collision_operators; ++i) {
      real thisReLam = m_collision_operators[i]->computeRealLam(dv);
      if (thisReLam > reLam) {
         reLam = thisReLam;
      }
   }

   /*// return the maximal time step satisfying (reLam*dt)^2+(imLam*dt)^2=alpha^2
   real ddt;
   real alpha(2.6);
   ddt = sqrt(alpha*alpha/(reLam*reLam+imLam*imLam));*/

   // return the maximal time step satisfying (reLam*dt/alpha)^2+(imLam*dt/beta)^2=1
   real ddt, alpha, beta;
   if( m_temporal_solution_order == 4 ) {
     // RK4
     alpha = 2.6;
     beta  = 2.6;
   } else {
     // RK6
     alpha = 4.95;
     beta  = 3.168;
   }
   ddt = sqrt( 1.0/(reLam*reLam/(alpha*alpha)+imLam*imLam/(beta*beta)) );

   return(ddt);
}


void
KineticSpecies::computeAcceleration(
   const realArray& a_efield,
   realArray&       a_ext_efield_global,
   real             a_time,
   real             a_dt,
   int              a_stage)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("phys to phase");

   // We need to know the self consistent E field computed by the Poisson
   // process.  Communicate the part of this 2D field corresponding to each
   // species' configuration space extent from the Poisson processor to each
   // species.
   m_accel = 0.0;
   m_efield_expansion_schedule->execute(a_efield, m_accel);

   timers->stopTimer("phys to phase");

   // apply external driver if applicable
   //   Note that the driver is applied in all cells (including ghost cells)
   //   and periodicity is not thereafter enforced ... that is on the user.
   timers->startTimer("driver");
   if (m_num_external_drivers > 0) {
      // The externally applied E field only needs to be tracked for the update
      // of particles.
      if (m_problem_has_particles) {
         m_ext_efield_local = 0.0;
      }

      // HMM, if I was a little smarter here we could get rid of this swizzling
      // and just pass m_accel to evaluate.
      Range Rx(m_accel.getBase(0), m_accel.getBound(0));
      Range Ry(m_accel.getBase(1), m_accel.getBound(1));
      RealArray Ex_extern(Rx, Ry);
      RealArray Ey_extern(Rx, Ry);
      for (int i = 0; i < m_num_external_drivers; ++i) {
         Ex_extern = 0.0;
         Ey_extern = 0.0;

         m_ef_drivers[i]->evaluate(Ex_extern,
            Ey_extern,
            m_accel_box,
            *m_domain,
            a_time,
            a_dt,
            a_stage);
         m_accel(Rx, Ry, Maxwell::EX) += Ex_extern(Rx, Ry);
         m_accel(Rx, Ry, Maxwell::EY) += Ey_extern(Rx, Ry);

         // Add this drivers' field into the total external E field.
         if (m_problem_has_particles) {
            m_ext_efield_local(Rx, Ry, Maxwell::EX) += Ex_extern(Rx, Ry);
            m_ext_efield_local(Rx, Ry, Maxwell::EY) += Ey_extern(Rx, Ry);
         }
      }

      // Communicate this species' external efield to a_ext_field which is
      // defined on the Poisson processor(s).
      if (m_problem_has_particles) {
         m_ext_efield_schedule->execute(a_ext_efield_global);
      }
   }
   timers->stopTimer("driver");

   timers->startTimer("blowout");
   // locally turn into acceleration
   m_accel *= (m_charge / m_mass);

   // Compute the accelerations on the faces in the V1 and V2 directions and the
   // maximum acclerations which are necessary for the time step calculation.
   real axmax, aymax;
   FORT_SET_PHASE_SPACE_VEL_4D(*(m_vel_face[V1].getDataPointer()),
      *(m_vel_face[V2].getDataPointer()),
      BOX4D_TO_FORT(m_local_box),
      *(m_accel.getDataPointer()),
      BOX4D_TO_FORT(m_accel_box),
      axmax, aymax);

   m_lambda_max[V1] = axmax;
   m_lambda_max[V2] = aymax;
   timers->stopTimer("blowout");
}


void
KineticSpecies::computeAcceleration(
   Maxwell&   a_maxwell,
   realArray& a_ext_efield_global,
   real       a_time,
   real       a_dt,
   int        a_stage)
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("phys to phase");

   // We need to know the electromagnetic fields computed by the Maxwell
   // process.  Communicate the part of these 2D fields corresponding to each
   // species' configuration space extent from the Maxwell processor to each
   // species.
   m_em_vars = 0.0;
   m_em_expansion_schedule->execute(a_maxwell.getGlobalEMVars(), m_em_vars);

   timers->stopTimer("phys to phase");

   // apply external driver if applicable
   //   Note that the driver is applied in all cells (including ghost cells)
   //   and periodicity is not thereafter enforced ... that is on the user.
   timers->startTimer("driver");
   if (m_num_external_drivers > 0) {
      // The externally applied E field only needs to be tracked for the update
      // of particles.
      if (m_problem_has_particles) {
         m_ext_efield_local = 0.0;
      }

      // HMM, if I was a little smarter here we could get rid of this swizzling
      // and just pass m_accel to evaluate.
      Range Rx(m_accel_box.lower(X1), m_accel_box.upper(X1));
      Range Ry(m_accel_box.lower(X2), m_accel_box.upper(X2));
      RealArray Ex_extern(Rx, Ry);
      RealArray Ey_extern(Rx, Ry);
      for (int i = 0; i < m_num_external_drivers; ++i) {
         Ex_extern = 0.0;
         Ey_extern = 0.0;

         m_ef_drivers[i]->evaluate(Ex_extern,
            Ey_extern,
            m_accel_box,
            *m_domain,
            a_time,
            a_dt,
            a_stage);

         m_em_vars(Rx, Ry, Maxwell::EX) += Ex_extern(Rx, Ry);
         m_em_vars(Rx, Ry, Maxwell::EY) += Ey_extern(Rx, Ry);

         // Add this drivers' field into the total external E field.
         if (m_problem_has_particles) {
            m_ext_efield_local(Rx, Ry, Maxwell::EX) += Ex_extern(Rx, Ry);
            m_ext_efield_local(Rx, Ry, Maxwell::EY) += Ey_extern(Rx, Ry);
         }
      }

      // Communicate this species' external efield to a_ext_field which is
      // defined on the Poisson processor(s).
      if (m_problem_has_particles) {
         m_ext_efield_schedule->execute(a_ext_efield_global);
      }
   }
   timers->stopTimer("driver");

   timers->startTimer("blowout");
   // Compute the accelerations on the faces in the V1 and V2 directions and the
   // maximum acclerations which are necessary for the time step calculation.
   real axmax, aymax, charge_per_mass = m_charge / m_mass;
   FORT_SET_PHASE_SPACE_VEL_MAXWELL_4D(*(m_vel_face[V1].getDataPointer()),
      *(m_vel_face[V2].getDataPointer()),
      BOX4D_TO_FORT(m_local_box),
      PROBLEMDOMAIN_TO_FORT((*m_domain)),
      charge_per_mass,
      *m_em_vars.getDataPointer(),
      *m_vz.getDataPointer(),
      axmax,
      aymax);

   m_lambda_max[V1] = axmax;
   m_lambda_max[V2] = aymax;
   timers->stopTimer("blowout");
}


void
KineticSpecies::currentDensity(
   const Maxwell& a_maxwell,
   realArray& a_Jx,
   realArray& a_Jy,
   realArray& a_Jz)
{
   // First time through define the local current densities and the reduction
   // schedules needed to compute the current densities.
   if (!m_jx_schedule) {
      defineCurrentDensityReductionSchedules();
   }

   // Get the local z drift velocity from Maxwell.
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("phys to phase");

   m_vz = 0.0;
   m_vz_expansion_schedule->execute(a_maxwell.getGlobalVZVar(m_species_index),
      m_vz);

   timers->stopTimer("phys to phase");

   // Zero out the local current densities.
   m_Jx_local = 0.0;
   m_Jy_local = 0.0;
   m_Jz_local = 0.0;

   // Compute the 4D current densities.
   FORT_COMPUTE_CURRENTS_4D(BOX4D_TO_FORT(m_local_box),
      BOX4D_TO_FORT(m_interior_box),
      PROBLEMDOMAIN_TO_FORT((*m_domain)),
      *m_dist_func_local.getDataPointer(),
      *m_vz.getDataPointer(),
      *m_Jx_local.getDataPointer(),
      *m_Jy_local.getDataPointer(),
      *m_Jz_local.getDataPointer());

   // Now reduce the 4D current densities into 2D quantites defined on the
   // Maxwell processors by integrating over the velocities.
   m_jx_schedule->execute(a_Jx, DensityKernel(m_charge));
   m_jy_schedule->execute(a_Jy, DensityKernel(m_charge));
   m_jz_schedule->execute(a_Jz, DensityKernel(m_charge));
}


void
KineticSpecies::getFromRestart(
   const HDF_DataBase& a_db)
{
   // find subdatabase with the name of this distribution
   HDF_DataBase sub_db;
   a_db.locate(sub_db, m_name);

   m_domain->getFromDatabase(sub_db);

   // read restart distribution from database
   sub_db.getDistributed(m_dist_func_global, "distribution");
#ifdef USE_PPP
   getLocalArrayWithGhostBoundaries(m_dist_func_global, m_dist_func_local);
#endif

   // As I recall, Overture does not read either any of the ghosts or the extra
   // ghosts on a physical boundary so we need to zero out that data.
   FORT_ZERO_GHOST_4D(*m_dist_func_local.getDataPointer(),
      BOX4D_TO_FORT(m_interior_box),
      BOX4D_TO_FORT(m_local_box));
   m_krook_layer->getFromDatabase(sub_db);
}


void
KineticSpecies::putToRestart(
   HDF_DataBase& a_db,
   real a_time)
{
   // Write everthing except Krook data.
   putToRestart_SkipKrook(a_db, a_time);

   // Now write Krook data.
   HDF_DataBase sub_db;
   a_db.locate(sub_db, m_name, "directory");
   m_krook_layer->putToDatabase(sub_db);
}


void
KineticSpecies::putToRestart_SkipKrook(
   HDF_DataBase& a_db,
   real a_time)
{
   // make a subdatabase with the  name of this distribution
   HDF_DataBase sub_db;
   a_db.create(sub_db, m_name, "directory");

   // save the pdim, cdim dimensions
   sub_db.put(m_pdim, "pdim");
   sub_db.put(m_cdim, "cdim");

   // save the mass and charge
   sub_db.put(m_mass, "mass");
   sub_db.put(m_charge, "charge");

   // save the distribution function
   if (m_tz_source) {
      // If a twilight zone is defined, we want to see the error in the
      // distribution functions, not the functions themselves.
      realArray tz_error_global;
      tz_error_global.partition(m_dist_func_global.getPartition());
      tz_error_global.redim(BoxOps::range(m_global_box, X1),
         BoxOps::range(m_global_box, X2),
         BoxOps::range(m_global_box, V1),
         BoxOps::range(m_global_box, V2));
#ifdef USE_PPP
      RealArray tz_error_local;
      getLocalArrayWithGhostBoundaries(tz_error_global, tz_error_local);
#else
      RealArray& tz_error_local = tz_error_global;
#endif
      int my_id(Communication_Manager::My_Process_Number);
      my_id = std::max(0, my_id);
      if (isInRange(my_id)) {
         m_tz_source->computeError(tz_error_local,
            m_dist_func_local,
            *m_domain,
            a_time);
      }
      else {
         tz_error_local = 0.0;
      }
      sub_db.putDistributed(tz_error_global, "distribution");
   }
   else {
      sub_db.putDistributed(m_dist_func_global, "distribution");
   }

   // The ProblemDomain must be saved as well.
   m_domain->putToDatabase(sub_db);
}


void
KineticSpecies::computeKEVelBdyFlux(
   Poisson& a_poisson,
   realArray& a_species_vel_bdry_flux)
{
   // We need the up to date acceleration ghosts and acceleration fluxes to do
   // this.
   fillAccelerationGhostCells();
   evalAccelerationFluxes();

   // Loop over the 2 velocity boundaries in each velocity direction.  In each
   // case, if this processor is on that boundary compute the KE flux across it.
   // Then sum the contribution to this flux from each processor and communicate
   // the sum to the Poisson processor.
   const tbox::Box& domain_box = m_domain->box();
   for (int dir = V1; dir <= V2; ++dir) {
      for (int side = LO; side <= HI; ++side) {
         m_ke_vel_bdry_flux_local = 0.0;
         bool on_vel_bdy;
         if (side == LO) {
            on_vel_bdy = m_interior_box.lower(dir) == domain_box.lower(dir);
         }
         else {
            on_vel_bdy = m_interior_box.upper(dir) == domain_box.upper(dir);
         }
         if (on_vel_bdy) {
            FORT_COMPUTE_KE_VEL_SPACE_FLUX(BOX4D_TO_FORT(m_local_box),
               BOX4D_TO_FORT(m_interior_box),
               BOX4D_TO_FORT(domain_box),
               PROBLEMDOMAIN_TO_FORT((*m_domain)),
               *m_flux[V1].getDataPointer(),
               *m_flux[V2].getDataPointer(),
               *m_ke_vel_bdry_flux_local.getDataPointer(),
               m_mass,
               side,
               dir);
         }
         m_ke_vel_bdry_flux_schedule->execute(a_species_vel_bdry_flux);

         Index dest[4], src[4];
         dest[0] = Range(BoxOps::range(domain_box, X1));
         dest[1] = Range(BoxOps::range(domain_box, X2));
         dest[2] = Range(0, 0);
         dest[3] = Range(0, 0);
         src[0] = Range(BoxOps::range(domain_box, X1));
         src[1] = Range(BoxOps::range(domain_box, X2));
         src[2] = Range(0, 0);
         src[3] = Range(0, 0);
         ParallelUtility::copy(a_poisson.getKEFluxVar(m_species_index, dir, side),
                               dest,
                               a_species_vel_bdry_flux,
                               src,
                               4);
      }
   }
}


void
KineticSpecies::computeKEVelBdyFlux(
   Maxwell& a_maxwell,
   realArray& a_species_vel_bdry_flux)
{
   // We need the up to date acceleration ghosts and acceleration fluxes to do
   // this.
   fillAccelerationGhostCells();
   evalAccelerationFluxes();

   // Loop over the 2 velocity boundaries in each velocity direction.  In each
   // case, if this processor is on that boundary compute the KE flux across it.
   // Then sum the contribution to this flux from each processor and communicate
   // the sum to the Maxwell processor.
   const tbox::Box& domain_box = m_domain->box();
   for (int dir = V1; dir <= V2; ++dir) {
      for (int side = LO; side <= HI; ++side) {
         m_ke_vel_bdry_flux_local = 0.0;
         bool on_vel_bdy;
         if (side == LO) {
            on_vel_bdy = m_interior_box.lower(dir) == domain_box.lower(dir);
         }
         else {
            on_vel_bdy = m_interior_box.upper(dir) == domain_box.upper(dir);
         }
         if (on_vel_bdy) {
            FORT_COMPUTE_KE_VEL_SPACE_FLUX(BOX4D_TO_FORT(m_local_box),
               BOX4D_TO_FORT(m_interior_box),
               BOX4D_TO_FORT(domain_box),
               PROBLEMDOMAIN_TO_FORT((*m_domain)),
               *m_flux[V1].getDataPointer(),
               *m_flux[V2].getDataPointer(),
               *m_ke_vel_bdry_flux_local.getDataPointer(),
               m_mass,
               side,
               dir);
         }
         m_ke_vel_bdry_flux_schedule->execute(a_species_vel_bdry_flux);

         Index dest[4], src[4];
         dest[0] = Range(BoxOps::range(domain_box, X1));
         dest[1] = Range(BoxOps::range(domain_box, X2));
         dest[2] = Range(0, 0);
         dest[3] = Range(0, 0);
         src[0] = Range(BoxOps::range(domain_box, X1));
         src[1] = Range(BoxOps::range(domain_box, X2));
         src[2] = Range(0, 0);
         src[3] = Range(0, 0);
         ParallelUtility::copy(a_maxwell.getKEFluxVar(m_species_index, dir, side),
                               dest,
                               a_species_vel_bdry_flux,
                               src,
                               4);
      }
   }
}


void
KineticSpecies::accumulateSequences(
   bool a_is_2d_proc,
   RealArray& a_sequences,
   int a_saved_seq,
   int& a_seq_idx,
   realArray& a_ke_global,
   realArray& a_ke_flux_global)
{
   // Species Kinetic Energy.
   // Zero out the local 4D kinetic energy.
   m_ke_local = 0.0;

   // Compute the 4D kinetic energy.  This is only performed by the Vlasov
   // processors.
   FORT_COMPUTE_KE_4D(BOX4D_TO_FORT(m_local_box),
      BOX4D_TO_FORT(m_interior_box),
      PROBLEMDOMAIN_TO_FORT((*m_domain)),
      *m_dist_func_local.getDataPointer(),
      *m_ke_local.getDataPointer());

   // Now reduce the 4D kinetic energy into a 2D quantity by integrating over
   // the velocities.  The 2D quantity is owned by the non-Vlasov processor(s).
   m_ke_schedule->execute(a_ke_global, DensityKernel(m_mass));

   // Compute the local scalar kinetic energy from the local 2D kinetic energy
   // (integrate over the configuration space dimensions).  This is only
   // performed by the non-Vlasov processor(s).
   real ke = 0.0;
   const tbox::Box& domain_box = m_domain->box();
   if (a_is_2d_proc) {
#ifdef USE_PPP
      RealArray ke_local;
      getLocalArrayWithGhostBoundaries(a_ke_global, ke_local);
#else
      RealArray& ke_local = a_ke_global;
#endif
      const tbox::Box sequence_box(BoxOps::getLocalBox(ke_local));
      const tbox::Box interior_box = BoxOps::getOperationalBox(ke_local,
         a_ke_global,
         domain_box,
         m_global_box);
      FORT_COMPUTE_KE_2D(BOX2D_TO_FORT(sequence_box),
         BOX2D_TO_FORT(interior_box),
         PROBLEMDOMAIN_TO_FORT((*m_domain)),
         *ke_local.getDataPointer(),
         ke);
   }

   // Sum the local scalar kinetic energy and add it to the time history.
   a_sequences(a_saved_seq, a_seq_idx++) = ParallelUtility::getSum(ke, -1);

   // Species kinetic energy flux at left, right, top, and bottom edges.
   // We need the current advection fluxes and hence the current physical
   // boundary conditions to do this calculation.
   setPhysicalBCs();
   evalAdvectionFluxes();

   // Loop over the 2 boundaries in each configuation space dimension.
   for (int dir = X1; dir <= X2; ++dir) {
      // Zero out the local 4D kinetic energy flux.
      m_ke_phys_bdry_flux_local = 0.0;

      // Compute the 4D kinetic energy flux.  This is only performed by the
      // Vlasov processors.
      FORT_COMPUTE_KE_FLUX_4D(BOX4D_TO_FORT(m_local_box),
         BOX4D_TO_FORT(m_interior_box),
         BOX4D_TO_FORT(domain_box),
         PROBLEMDOMAIN_TO_FORT((*m_domain)),
         *m_flux[X1].getDataPointer(),
         *m_flux[X2].getDataPointer(),
         *m_vel_face[X1].getDataPointer(),
         *m_vel_face[X2].getDataPointer(),
         *m_ke_phys_bdry_flux_local.getDataPointer(),
         dir);

      // Now reduce the 4D kinetic energy flux into a 2D quantity by
      // integrating over the velocities.  The 2D quantity is owned by the
      // non-Vlasov processor(s).
      m_ke_phys_bdry_flux_schedule->execute(a_ke_flux_global,
         DensityKernel(m_mass));

      for (int side = LO; side <= HI; ++side) {
         // Compute the local scalar kinetic energy flux from the local 2D
         // energy flux (integrate over the configuration space dimensions).
         // This is only performed by the non-Vlasov processor(s).
         real ke_flux = 0.0;
         if (a_is_2d_proc) {
#ifdef USE_PPP
            RealArray ke_flux_local;
            getLocalArrayWithGhostBoundaries(a_ke_flux_global, ke_flux_local);
#else
            RealArray& ke_flux_local = a_ke_flux_global;
#endif
            const tbox::Box sequence_box(BoxOps::getLocalBox(ke_flux_local));
            const tbox::Box interior_box =
               BoxOps::getOperationalBox(ke_flux_local,
                  a_ke_flux_global,
                  domain_box,
                  m_global_box);
            FORT_COMPUTE_KE_FLUX_2D(BOX2D_TO_FORT(sequence_box),
               BOX2D_TO_FORT(interior_box),
               BOX2D_TO_FORT(domain_box),
               PROBLEMDOMAIN_TO_FORT((*m_domain)),
               *ke_flux_local.getDataPointer(),
               ke_flux,
               dir, side);
         }

         // Sum the local scalar kinetic energy flux and add it to the time
         // history.
         a_sequences(a_saved_seq, a_seq_idx++) =
            ParallelUtility::getSum(ke_flux, -1);
      }
   }
}


//// PRIVATE METHODS ////////////////////////////////////////////////////

void
KineticSpecies::parseParameters(
   ParmParse& a_pp,
   const tbox::Pointer<ProblemDomain> a_cfg_domain)
{
   // A species name is required.
   aString tmp;
   if (a_pp.contains("name")) {
      a_pp.get("name", tmp);
      m_name = tmp;
   }
   else {
      OV_ABORT("Must supply name!");
   }

   // A species maxx is required.
   if (a_pp.contains("mass")) {
      a_pp.get("mass", m_mass);
   }
   else {
      OV_ABORT("Must supply mass!");
   }

   // A species charge is required.
   if (a_pp.contains("charge")) {
      a_pp.get("charge", m_charge);
   }
   else {
      OV_ABORT("Must supply charge!");
   }

   // The user may dictate the number of processors for this species to use.
   // Don't think this has ever been tested so it's not a very good idea.
   if (a_pp.contains("number_of_processors")) {
      a_pp.get("number_of_processors", m_number_of_procs);
      m_fixed_number_of_procs = true;
   }

   // This is a bit of backward compatibility for old decks when there was only
   // one external driver possible for each species.
   bool new_driver_syntax = a_pp.contains("num_external_drivers");
   bool old_driver_syntax = a_pp.contains("apply_external_driver");
   if (new_driver_syntax) {
      if (old_driver_syntax) {
         OV_ABORT("Mixed old and new electric field driver syntax!");
      }
      else {
         a_pp.query("num_external_drivers", m_num_external_drivers);
      }
   }
   else if (old_driver_syntax) {
      m_old_driver_syntax = true;
      aString driver_on("false");
      a_pp.query("apply_external_driver", driver_on);
      m_num_external_drivers = driver_on.matches("false") ? 0 : 1;
   }
   m_ef_drivers.resize(m_num_external_drivers,
      tbox::Pointer<ElectricFieldDriver>(0));

   // See if there's any collision operators.
   a_pp.query("num_collision_operators", m_num_collision_operators);
   m_collision_operators.resize(m_num_collision_operators,
      tbox::Pointer<CollisionOperator>(0));

   // Get the configuration space limits from the (2D) ProblemDomain.
   Array<double> limits(2 * m_pdim);
   for (int d(0); d < m_cdim; ++d) {
      limits[2*d] = a_cfg_domain->lower(d);
      limits[2*d+1] = a_cfg_domain->upper(d);
   }

   // Now get the velocity space limits which are required.
   if (a_pp.contains("velocity_limits")) {
      tbox::Dimension vdim(static_cast<unsigned short>(m_pdim - m_cdim));
      Array<double> vlimits(2 * vdim);
      a_pp.getarr("velocity_limits",
         vlimits,
         0,
         static_cast<int>(vlimits.length()));
      for (int i(0); i < vlimits.length(); ++i) {
         int j(2 * m_cdim + i);
         limits[j] = vlimits[i];
      }
   }
   else {
      OV_ABORT("Must supply velocity_limits!");
   }

   // Dump all the limits into x_lo and x_hi so we can create the species'
   // ProblemDomain.
   Array<double> x_lo(m_pdim);
   Array<double> x_hi(m_pdim);
   for (int d(0); d < m_pdim; ++d) {
      x_lo[d] = limits[2*d];
      x_hi[d] = limits[2*d+1];
   }

   // Get the configuration space periodicity from the (2D) ProblemDomain.
   // Velocity space periodicity is always off.
   Array<bool> is_periodic(m_pdim);
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
      Array<int> tmp(m_pdim);
      a_pp.getarr("Nv", tmp, 0, m_pdim - m_cdim);
      for (int d(m_cdim); d < m_pdim; ++d) {
         n_cells[d] = tmp[d-m_cdim];
      }
   }
   else {
      OV_ABORT("Must supply Nv!");
   }

   // Now create this species' ProblemDomain.
   m_domain = new ProblemDomain(m_pdim, n_cells, x_lo, x_hi, is_periodic);

   // Read any initial condition flow velocity.
   a_pp.query("vflowinitx", m_vflowx);
   a_pp.query("vflowinity", m_vflowy);
}


void
KineticSpecies::allocateLocalAuxArrays()
{
   // We need the distribution function, the velocity/acceleration, and the
   // distribution function flux on each face.
   m_u_face.resize(m_pdim);
   m_vel_face.resize(m_pdim);
   m_flux.resize(m_pdim);

   for (int dir(0); dir < m_pdim; ++dir) {
      allocateFaceArray(m_u_face[dir],   m_local_box, dir);
      allocateFaceArray(m_vel_face[dir], m_local_box, dir);
      allocateFaceArray(m_flux[dir],     m_local_box, dir);
   }

   // Figure out the box corresponding to the configuration space part of this
   // species on this processor.
   m_accel_box = m_local_box;
   m_accel_box.lower(CDIM) = 0;
   m_accel_box.upper(CDIM) = 1;
   for (int d(CDIM+1); d < PDIM; ++d) {
      m_accel_box.lower(d) = 0;
      m_accel_box.upper(d) = 0;
   }

   // Dimension the 2D accleration.
   m_accel.redim(Range(m_accel_box.lower(X1), m_accel_box.upper(X1)),
      Range(m_accel_box.lower(X2), m_accel_box.upper(X2)),
      Range(0, 1));

   // If this is an electrodynamic problem then we need more stuff.
   if (m_do_maxwell) {
      // Figure out the box corresponding to the electromagnetic fields on the
      // configuarion space part of this species.
      m_em_vars_box = m_local_box;
      m_em_vars_box.lower(CDIM) = 0;
      m_em_vars_box.upper(CDIM) = Maxwell::NUM_EM_VARS-1;
      for (int d(CDIM+1); d < PDIM; ++d) {
         m_em_vars_box.lower(d) = 0;
         m_em_vars_box.upper(d) = 0;
      }

      // Dimension the 2D electromagnetic field.
      m_em_vars.redim(Range(m_em_vars_box.lower(X1), m_em_vars_box.upper(X1)),
         Range(m_em_vars_box.lower(X2), m_em_vars_box.upper(X2)),
         Range(0, Maxwell::NUM_EM_VARS-1));

      // Figure out the box corresponding to the z drift velocity on the
      // configuarion space part of this species.
      m_vz_box = m_local_box;
      m_vz_box.lower(CDIM) = 0;
      m_vz_box.upper(CDIM) = 0;
      for (int d(CDIM+1); d < PDIM; ++d) {
         m_vz_box.lower(d) = 0;
         m_vz_box.upper(d) = 0;
      }

      // Dimension the 2D z drift velocity of this species.
      m_vz.redim(Range(m_vz_box.lower(X1), m_vz_box.upper(X1)),
         Range(m_vz_box.lower(X2), m_vz_box.upper(X2)),
         Range(0, 0));
   }
}


void
KineticSpecies::initializeVelocity()
{
   // Zero out the velocities.
   m_vel_face.resize(m_pdim);
   for (int dir(X1); dir < m_pdim; ++dir) {
      m_vel_face[dir] = 0.0;
   }

   // Set the velocity on the X1 and X2 faces.
   const Array<double>& x_lo(m_domain->lower());
   const Array<double>& dx(m_domain->dx());
   const tbox::IntVector& ncells(m_domain->numberOfCells());
   for (int i3(m_local_box.lower(V1)); i3 <= m_local_box.upper(V1); ++i3) {
      real coord_x3(x_lo[V1] + (i3 + 0.5) * dx[V1]);
      for (int i4(m_local_box.lower(V2)); i4 <= m_local_box.upper(V2); ++i4) {
         for (int i2(m_local_box.lower(X2));
              i2 <= m_local_box.upper(X2); ++i2) {
            for (int i1(m_local_box.lower(X1));
                 i1 <= m_local_box.upper(X1)+1; ++i1) {
               m_vel_face[X1](i1, i2, i3, i4) = coord_x3;
            }
         }
      }
   }

   for (int i4(m_local_box.lower(V2)); i4 <= m_local_box.upper(V2); ++i4) {
      real coord_x4(x_lo[V2] + (i4 + 0.5) * dx[V2]);
      for (int i3(m_local_box.lower(V1)); i3 <= m_local_box.upper(V1); ++i3) {
         for (int i1(m_local_box.lower(X1));
              i1 <= m_local_box.upper(X1); ++i1) {
            for (int i2(m_local_box.lower(X2));
                 i2 <= m_local_box.upper(X2)+1; ++i2){
               m_vel_face[X2](i2, i3, i4, i1) = coord_x4;
            }
         }
      }
   }

   // Zero out the lambdas
   m_lambda_max.resize(m_pdim);
   for (int dir(X1); dir < m_pdim; ++dir) {
      m_lambda_max[dir] = 0.0;
   }

   // Set the lambdas from the upper and lower velocity bounds.
   m_lambda_max[X1] = std::max(abs(x_lo[V1] + 0.5 * dx[V1]),
                               abs(x_lo[V1] + (0.5 + ncells[V1]) * dx[V1]));

   m_lambda_max[X2] = std::max(abs(x_lo[V2] + 0.5 * dx[V2]),
                               abs(x_lo[V2] + (0.5 + ncells[V2]) * dx[V2]));
}


void
KineticSpecies::setPeriodicBCs()
{
   const tbox::Box& domain_box(m_domain->box());

   // For each configuration space dimension set the upper ghosts to the lower
   // interior and lower ghosts to the upper interior.
   for (int dir(X1); dir <= X2; ++dir) {
      if (m_domain->isPeriodic(dir)) {

         Index dst[4];
         Index src[4];
         for (int i(0); i < m_pdim; ++i) {
            dst[i] = Range(BoxOps::range(m_global_box, i));
            src[i] = Range(BoxOps::range(m_global_box, i));
         }

         tbox::Box inner_box(domain_box);
         inner_box.grow(-m_n_ghosts);

         dst[dir] = Range(m_global_box.lower(dir), domain_box.lower(dir) - 1);
         src[dir] = Range(inner_box.upper(dir) + 1, domain_box.upper(dir));
         ParallelUtility::copy(m_dist_func_global,
            dst,
            m_dist_func_global,
            src,
            4);

         dst[dir] = Range(domain_box.upper(dir) + 1, m_global_box.upper(dir));
         src[dir] = Range(domain_box.lower(dir), inner_box.lower(dir) - 1);
         ParallelUtility::copy(m_dist_func_global,
            dst,
            m_dist_func_global,
            src,
            4);
      }
      // There must be a barrier here as the copies above are all non-blocking.
      // We can not interlace the x and y communications.
      MPI_Barrier(m_comm);
   }
}


void
KineticSpecies::addFluxDivergence(
   KineticSpecies& a_rhs) const
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("Vlasov");

   FORT_ACCUM_FLUX_DIV_4D(*(a_rhs.m_dist_func_local).getDataPointer(),
      BOX4D_TO_FORT(m_local_box),
      BOX4D_TO_FORT(m_interior_box),
      *(m_flux[X1].getDataPointer()),
      *(m_flux[X2].getDataPointer()),
      *(m_flux[V1].getDataPointer()),
      *(m_flux[V2].getDataPointer()),
      PROBLEMDOMAIN_TO_FORT((*m_domain)));

   timers->stopTimer("Vlasov");
}


void
KineticSpecies::defineExtEfieldContractionSchedule()
{
   if (m_problem_has_particles) {
      m_ext_efield_schedule = new ContractionSchedule(m_ext_efield_local,
         m_interior_box,
         m_global_box,
         *m_domain,
         m_processor_range,
         m_comm);
   }
}


void
KineticSpecies::defineChargeDensityReductionSchedule()
{
   // We want to integrate over velocity.
   Array<bool> collapse_dir(m_pdim, false);
   collapse_dir[V1] = true;
   collapse_dir[V2] = true;
   double measure(m_domain->dx(2) * m_domain->dx(3));

   // Now make the schedule for the charge density.
   tbox::Pointer<ReductionOp> int_op(new IntegralOp(measure));
   tbox::Pointer<Reduction> sum_reduction(
      new Reduction(m_pdim, collapse_dir, int_op));
   m_moment_schedule = new ReductionSchedule(m_dist_func_local,
      m_interior_box,
      *m_domain,
      sum_reduction,
      m_processor_range,
      m_comm);
}


void
KineticSpecies::defineKineticEnergyReductionSchedule()
{
   // Dimension this species kinetic energy.
   m_ke_local.redim(Range(m_local_box.lower(X1), m_local_box.upper(X1)),
      Range(m_local_box.lower(X2), m_local_box.upper(X2)),
      Range(m_local_box.lower(V1), m_local_box.upper(V1)),
      Range(m_local_box.lower(V2), m_local_box.upper(V2)),
      Range(0, 0));

   // We want to integrate over velocity.
   Array<bool> collapse_dir(m_pdim, false);
   collapse_dir[V1] = true;
   collapse_dir[V2] = true;
   double measure(m_domain->dx(2) * m_domain->dx(3));

   // Now make the schedule for the KE.
   tbox::Pointer<ReductionOp> int_op_ke(new IntegralOp(measure));
   tbox::Pointer<Reduction> sum_reduction_ke(
      new Reduction(m_pdim, collapse_dir, int_op_ke));
   m_ke_schedule = new ReductionSchedule(m_ke_local,
      m_interior_box,
      *m_domain,
      sum_reduction_ke,
      m_processor_range,
      m_comm);

   // Dimension this species kinetic energy physical boundary flux.
   m_ke_phys_bdry_flux_local.redim(Range(m_local_box.lower(X1), m_local_box.upper(X1)),
      Range(m_local_box.lower(X2), m_local_box.upper(X2)),
      Range(m_local_box.lower(V1), m_local_box.upper(V1)),
      Range(m_local_box.lower(V2), m_local_box.upper(V2)),
      Range(0, 0));

   // Now make the schedule for the KE physical boundary flux.
   tbox::Pointer<ReductionOp> int_op_ke_flux(new IntegralOp(measure));
   tbox::Pointer<Reduction> sum_reduction_ke_flux(
      new Reduction(m_pdim, collapse_dir, int_op_ke_flux));
   m_ke_phys_bdry_flux_schedule = new ReductionSchedule(m_ke_phys_bdry_flux_local,
      m_interior_box,
      *m_domain,
      sum_reduction_ke_flux,
      m_processor_range,
      m_comm);

   // Dimension this species kinetic energy velocity boundary flux.
   m_ke_vel_bdry_flux_local.redim(Range(m_local_box.lower(X1), m_local_box.upper(X1)),
      Range(m_local_box.lower(X2), m_local_box.upper(X2)),
      Range(0, 0));
   
   // Now make the schedule for the KE velocity boundary flux.
   m_ke_vel_bdry_flux_schedule = new SummationSchedule(m_ke_vel_bdry_flux_local,
      m_interior_box,
      m_domain->box(),
      m_processor_range,
      m_comm);
}


void
KineticSpecies::defineCurrentDensityReductionSchedules()
{
   // Dimension this species current densities.
   m_Jx_local.redim(Range(m_local_box.lower(X1), m_local_box.upper(X1)),
      Range(m_local_box.lower(X2), m_local_box.upper(X2)),
      Range(m_local_box.lower(V1), m_local_box.upper(V1)),
      Range(m_local_box.lower(V2), m_local_box.upper(V2)),
      Range(0, 0));
   m_Jy_local.redim(Range(m_local_box.lower(X1), m_local_box.upper(X1)),
      Range(m_local_box.lower(X2), m_local_box.upper(X2)),
      Range(m_local_box.lower(V1), m_local_box.upper(V1)),
      Range(m_local_box.lower(V2), m_local_box.upper(V2)),
      Range(0, 0));
   m_Jz_local.redim(Range(m_local_box.lower(X1), m_local_box.upper(X1)),
      Range(m_local_box.lower(X2), m_local_box.upper(X2)),
      Range(m_local_box.lower(V1), m_local_box.upper(V1)),
      Range(m_local_box.lower(V2), m_local_box.upper(V2)),
      Range(0, 0));

   // We want to integrate over velocity.
   Array<bool> collapse_dir(m_pdim, false);
   collapse_dir[V1] = true;
   collapse_dir[V2] = true;
   double measure(m_domain->dx(2) * m_domain->dx(3));

   // Now make the schedules for Jx, Jy, and Jz.
   tbox::Pointer<ReductionOp> int_op_x(new IntegralOp(measure));
   tbox::Pointer<Reduction> sum_reduction_x(
      new Reduction(m_pdim, collapse_dir, int_op_x));
   m_jx_schedule = new ReductionSchedule(m_Jx_local,
      m_interior_box,
      *m_domain,
      sum_reduction_x,
      m_processor_range,
      m_comm);
   tbox::Pointer<ReductionOp> int_op_y(new IntegralOp(measure));
   tbox::Pointer<Reduction> sum_reduction_y(
      new Reduction(m_pdim, collapse_dir, int_op_y));
   m_jy_schedule = new ReductionSchedule(m_Jy_local,
      m_interior_box,
      *m_domain,
      sum_reduction_y,
      m_processor_range,
      m_comm);
   tbox::Pointer<ReductionOp> int_op_z(new IntegralOp(measure));
   tbox::Pointer<Reduction> sum_reduction_z(
      new Reduction(m_pdim, collapse_dir, int_op_z));
   m_jz_schedule = new ReductionSchedule(m_Jz_local,
      m_interior_box,
      *m_domain,
      sum_reduction_z,
      m_processor_range,
      m_comm);
}


void
allocateFaceArray(
   RealArray& a_array,
   const tbox::Box& a_box,
   int a_dir)
{
   const unsigned int dim(a_box.getDim());
   std::vector<Range> r(dim);
   r[0] = Range(a_box.lower(a_dir), a_box.upper(a_dir) + 1);
   for (int j(1); j < static_cast<int>(dim); ++j) {
      int k = (a_dir+j)%dim;
      r[j] = Range(a_box.lower(k), a_box.upper(k));
   }

   if (dim == tbox::Dimension(4)) {
      a_array.redim(r[0], r[1], r[2], r[3]);
   }
   else {
      OV_ABORT("Not implemented for phase D!=4!");
   }

}

} // end namespace Loki
