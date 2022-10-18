/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "Particle.H"
#include "Directions.H"
#include "Loki_Defines.H"
#include "Loki_Utilities.H"
#include "RestartManager.H"

#include "hdf5.h"

namespace Loki {

Particle::Particle()
   : m_starting_time(0.0),
     m_charge(0.0),
     m_mass(0.0),
     m_x(0.0),
     m_y(0.0),
     m_vx(0.0),
     m_vy(0.0),
     m_noise_source_weight(0.0)
{
}


Particle::Particle(
   double a_starting_time,
   double a_charge,
   double a_mass,
   double a_xinit,
   double a_yinit,
   double a_vxinit,
   double a_vyinit,
   double a_noise_source_weight)
   : m_starting_time(a_starting_time),
     m_charge(a_charge),
     m_mass(a_mass),
     m_x(a_xinit),
     m_y(a_yinit),
     m_vx(a_vxinit),
     m_vy(a_vyinit),
     m_noise_source_weight(a_noise_source_weight)
{
}


Particle::Particle(
   const Particle& a_other)
   : m_starting_time(a_other.m_starting_time),
     m_charge(a_other.m_charge),
     m_mass(a_other.m_mass),
     m_x(a_other.m_x),
     m_y(a_other.m_y),
     m_vx(a_other.m_vx),
     m_vy(a_other.m_vy),
     m_noise_source_weight(a_other.m_noise_source_weight)
{
}


Particle::~Particle()
{
}


void
Particle::addParticleData(
   vector<Particle>& a_particles,
   const vector<Particle>& a_increment_particles,
   double a_factor,
   bool a_sum_reduce_inc,
   MPI_Comm a_comm)
{
   // It would be good to restructure the particles so that it is possible to
   // do the sum reduction of x, y, vx, and vy without these copies.

   // The particle RHS is computed by each EM processor only for those particles
   // on that EM processor.  Therefore, for parallel EM and when updating an
   // intermediate or final RK solution state the particle RHS must be sum
   // reduced across all EM processors so that the particle state is correctly
   // updated on each EM processor.  In either parallel or serial EM add
   // a_factor times a_increment_particles' x, y, vx, vy to a_particles.
   int num_particles = static_cast<int>(a_particles.size());
   if (a_sum_reduce_inc) {
      int num_locs = 4*num_particles;
      vector<double> locations(num_locs);
      int idx = 0;
      for (int i = 0; i < num_particles; ++i) {
         const Particle& inc_particle = a_increment_particles[i];
         locations[idx++] = inc_particle.x();
         locations[idx++] = inc_particle.y();
         locations[idx++] = inc_particle.vx();
         locations[idx++] = inc_particle.vy();
      }
      MPI_Allreduce(MPI_IN_PLACE,
         &locations[0],
         num_locs,
         MPI_DOUBLE,
         MPI_SUM,
         a_comm);
      idx = 0;
      for (int i = 0; i < num_particles; ++i) {
         Particle& this_particle = a_particles[i];
         this_particle.x() += locations[idx++]*a_factor;
         this_particle.y() += locations[idx++]*a_factor;
         this_particle.vx() += locations[idx++]*a_factor;
         this_particle.vy() += locations[idx++]*a_factor;
      }
   }
   else {
      for (int i = 0; i < num_particles; ++i) {
         Particle& this_particle = a_particles[i];
         const Particle& inc_particle = a_increment_particles[i];
         this_particle.x() += inc_particle.x()*a_factor;
         this_particle.y() += inc_particle.y()*a_factor;
         this_particle.vx() += inc_particle.vx()*a_factor;
         this_particle.vy() += inc_particle.vy()*a_factor;
      }
   }
}


void
Particle::copyParticleData(
   vector<Particle>& a_particles,
   const vector<Particle>& a_rhs_particles)
{
   // Copy the information from a_rhs_particles into a_particles.
   int num_particles = a_particles.size();
   for (int i = 0; i < num_particles; ++i) {
      Particle& this_particle = a_particles[i];
      const Particle& rhs_particle = a_rhs_particles[i];
      this_particle.x() = rhs_particle.x();
      this_particle.y() = rhs_particle.y();
      this_particle.vx() = rhs_particle.vx();
      this_particle.vy() = rhs_particle.vy();
   }
}


void
Particle::setBC(
   vector<Particle>& a_particles,
   const tbox::Pointer<ProblemDomain>& a_domain,
   const double* a_vmin,
   const double* a_vmax,
   bool a_reflect)
{
   // If a particle is outside the physical domain in a periodic direction then
   // wrap it back around to its periodic location.  If a particle is outside
   // the physical domain in a non-periodic direction then either limit its
   // position to the edge of the physical domain in that direction or, if
   // reflection has been requested, reflect the particle back inside the
   // domain.
   int num_particles = static_cast<int>(a_particles.size());
   for (int i = 0; i < num_particles; ++i) {
      Particle& this_particle = a_particles[i];
      if (this_particle.x() < a_domain->lower(X1)) {
         if (a_domain->isPeriodic(X1)) {
            this_particle.x() += a_domain->upper(X1) - a_domain->lower(X1);
         }
         else if (a_reflect) {
            this_particle.x() = 2.0*a_domain->lower(X1) - this_particle.x();
            this_particle.vx() = -this_particle.vx();
         }
         else {
            this_particle.x() = a_domain->lower(X1);
         }
      }
      else if (this_particle.x() > a_domain->upper(X1)) {
         if (a_domain->isPeriodic(X1)) {
            this_particle.x() += a_domain->lower(X1) - a_domain->upper(X1);
         }
         else if (a_reflect) {
            this_particle.x() = 2.0*a_domain->upper(X1) - this_particle.x();
            this_particle.vx() = -this_particle.vx();
         }
         else {
            this_particle.x() = a_domain->upper(X1);
         }
      }
      if (this_particle.y() < a_domain->lower(X2)) {
         if (a_domain->isPeriodic(X2)) {
            this_particle.y() += a_domain->upper(X2) - a_domain->lower(X2);
         }
         else if (a_reflect) {
            this_particle.y() =  2.0*a_domain->lower(X2) - this_particle.y();
            this_particle.vy() = -this_particle.vy();
         }
         else {
            this_particle.y() = a_domain->lower(X2);
         }
      }
      else if (this_particle.y() > a_domain->upper(X2)) {
         if (a_domain->isPeriodic(X2)) {
            this_particle.y() += a_domain->lower(X2) - a_domain->upper(X2);
         }
         else if (a_reflect) {
            this_particle.y() = 2.0*a_domain->upper(X2) - this_particle.y();
            this_particle.vy() = -this_particle.vy();
         }
         else {
            this_particle.y() = a_domain->upper(X2);
         }
      }

      this_particle.vx() = max(this_particle.vx(), a_vmin[0]);
      this_particle.vx() = min(this_particle.vx(), a_vmax[0]);

      this_particle.vy() = max(this_particle.vy(), a_vmin[1]);
      this_particle.vy() = min(this_particle.vy(), a_vmax[1]);
   }
}


void
Particle::putParticlesToRestart(const string& a_particle_file_name,
                                const vector<Particle>& a_particles,
                                bool a_write_noise_source_weight)
{
   // Glue the restart index onto the supplied particle file name to get
   // the particle restart file.
   ostringstream full_file_name;
   RestartManager* restart_manager(RestartManager::getManager());
   // Note that the RestartManager has already incremented the restart index at
   // the start of the restart write process so we need to subtract one from
   // the index here to get the right index.
   full_file_name << restart_manager->restartWritePath() << "/"
                  << a_particle_file_name << "_"
                  << restart_manager->restartIndex() - 1;

   // Create particle restart file.
   herr_t errf;
   hid_t file_id =
      H5Fcreate(full_file_name.str().c_str(),
         H5F_ACC_TRUNC,
         H5P_DEFAULT,
         H5P_DEFAULT);
   if (file_id < 0) {
      LOKI_ABORT("Unable to create particle restart file.");
   }

   // Create the dataspace for the number of particles.
   hsize_t dims1 = 1;
   hid_t dspace1 = H5Screate_simple(1, &dims1, NULL);
   if (dspace1 < 0) {
      LOKI_ABORT("Unable to create dataspace for number of particles.");
   }

   // Create the dataset for the number of particles.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dcreate_vers) && H5Dcreate_vers == 1
   hid_t dset1 = H5Dcreate(file_id,
                           "/num_particles",
                           H5T_STD_I32BE,
                           dspace1,
                           H5P_DEFAULT);
#else
   hid_t dset1 = H5Dcreate(file_id,
                           "/num_particles",
                           H5T_STD_I32BE,
                           dspace1,
                           H5P_DEFAULT,
                           H5P_DEFAULT,
                           H5P_DEFAULT);
#endif
#else
   hid_t dset1 = H5Dcreate(file_id,
                           "/num_particles",
                           H5T_STD_I32BE,
                           dspace1,
                           H5P_DEFAULT);
#endif
   if (dset1 < 0) {
      LOKI_ABORT("Unable to create dataset for number of particles.");
   }

   // Write the number of particles to the dataset.
   int num_particles = static_cast<int>(a_particles.size());
   errf = H5Dwrite(dset1,
                   H5T_NATIVE_INT,
                   H5S_ALL,
                   H5S_ALL,
                   H5P_DEFAULT,
                   &num_particles);
   if (errf < 0) {
      LOKI_ABORT("Unable to write to dataset for number of particles.");
   }

   // Close the dataset and dataspace for the number of particles.
   errf = H5Dclose(dset1);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataset for number of particles.");
   }
   errf = H5Sclose(dspace1);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataspace for number of particles.");
   }

   // Create the dataspace for the particle data.
   hsize_t dims2[2];
   dims2[0] = num_particles;
   if (a_write_noise_source_weight) {
      dims2[1] = 8;
   }
   else {
      dims2[1] = 7;
   }
   hid_t dspace2 = H5Screate_simple(2, dims2, NULL);
   if (dspace2 < 0) {
      LOKI_ABORT("Unable to create dataspace for particle data.");
   }

   // Create the dataset for the particle data but don't write anything
   // to it yet.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dcreate_vers) && H5Dcreate_vers == 1
   hid_t dset2 = H5Dcreate(file_id,
                           "/particle_data",
                           H5T_NATIVE_DOUBLE,
                           dspace2,
                           H5P_DEFAULT);
#else
   hid_t dset2 = H5Dcreate(file_id,
                           "/particle_data",
                           H5T_NATIVE_DOUBLE,
                           dspace2,
                           H5P_DEFAULT,
                           H5P_DEFAULT,
                           H5P_DEFAULT);
#endif
#else
   hid_t dset2 = H5Dcreate(file_id,
                           "/particle_data",
                           H5T_NATIVE_DOUBLE,
                           dspace2,
                           H5P_DEFAULT);
#endif
   if (dset2 < 0) {
      LOKI_ABORT("Unable to create dataset for particle data.");
   }

   // Close the dataset and dataspace for the particle data.
   errf = H5Dclose(dset2);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataset for particle data.");
   }
   errf = H5Sclose(dspace2);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataspace for particle data.");
   }
   errf = H5Fclose(file_id);
   if (errf < 0) {
      LOKI_ABORT("Unable to close particle restart file.");
   }

   // Reopen the particle restart file and the dataset for the particle data.
   file_id = H5Fopen(full_file_name.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
   if (file_id < 0) {
      LOKI_ABORT("Unable to reopen particle restart file.");
   }
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
   dset2 = H5Dopen(file_id, "/particle_data");
#else
   dset2 = H5Dopen(file_id, "/particle_data", H5P_DEFAULT);
#endif
#else
   dset2 = H5Dopen(file_id, "/particle_data");
#endif
   if (dset2 < 0) {
      LOKI_ABORT("Unable to reopen dataset for particle data.");
   }

   // Specify the size and shape of each particle's data.
   hsize_t particle_dim[2];
   particle_dim[0] = 1;
   particle_dim[1] = dims2[1];

   hsize_t stride[2];
   stride[0] = 1;
   stride[1] = 1;

   hsize_t block[2];
   block[0] = 1;
   block[1] = 1;

   hsize_t offset[2];
   offset[1] = 0;

   // Create memory space with size of a particle's data.  Get the file
   // dataspace and select subset from the file dataspace.
   hid_t memspace_id = H5Screate_simple(2, particle_dim, NULL);
   if (memspace_id < 0) {
      LOKI_ABORT("Unable to create particle data memory space.");
   }
   hid_t dataspace_id = H5Dget_space(dset2);
   if (dataspace_id < 0) {
      LOKI_ABORT("Unable to get particle data dataspace.");
   }

   // For each particle select the subset from the file dataspace and write
   // the particle data.
   double particle_data[8];
   for (int i = 0; i < num_particles; ++i) {
      offset[0] = i;
      errf = H5Sselect_hyperslab(dataspace_id,
                                 H5S_SELECT_SET,
                                 offset,
                                 stride,
                                 particle_dim,
                                 block);
      if (errf < 0) {
         LOKI_ABORT("Unable to select particle data hyperslab.");
      }
      const Particle& this_particle = a_particles[i];
      particle_data[0] = this_particle.startingTime();
      particle_data[1] = this_particle.charge();
      particle_data[2] = this_particle.mass();
      particle_data[3] = this_particle.x();
      particle_data[4] = this_particle.y();
      particle_data[5] = this_particle.vx();
      particle_data[6] = this_particle.vy();
      if (a_write_noise_source_weight) {
         particle_data[7] = this_particle.noiseSourceWeight();
      }
      errf = H5Dwrite(dset2,
                      H5T_NATIVE_DOUBLE,
                      memspace_id,
                      dataspace_id,
                      H5P_DEFAULT,
                      particle_data);
      if (errf < 0) {
         LOKI_ABORT("Unable to write to dataset for particle data.");
      }
   }

   // Close the dataset and dataspace for the particle data.
   errf = H5Sclose(memspace_id);
   if (errf < 0) {
      LOKI_ABORT("Unable to close particle data memory space.");
   }
   errf = H5Sclose(dataspace_id);
   if (errf < 0) {
      LOKI_ABORT("Unable to close particle data dataspace.");
   }
   errf = H5Dclose(dset2);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataset for particle data.");
   }
   errf = H5Fclose(file_id);
   if (errf < 0) {
      LOKI_ABORT("Unable to close particle restart file.");
   }
}


int
Particle::readParticleData(string& a_particle_file,
   vector<Particle>& a_particles,
   const tbox::Pointer<ProblemDomain>& a_domain,
   const double* a_vmin,
   const double* a_vmax,
   bool a_is_em_process,
   bool a_read_noise_source_weight)
{
   // If we're running from a restart then construct the name of the current
   // particle data restart file.
   ostringstream full_file_name;
   RestartManager* restart_manager(RestartManager::getManager());
   if (restart_manager->isFromRestart()) {
      // Glue the restart index onto the supplied particle file name to get
      // the particle restart file.
      full_file_name << restart_manager->restartReadPath() << "/"
                     << a_particle_file << "_"
                     << restart_manager->restartIndex();
   }
   else {
      full_file_name << a_particle_file;
   }

   herr_t errf;

   // Open particle file.
   hid_t file_id =
      H5Fopen(full_file_name.str().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   if (file_id < 0) {
      LOKI_ABORT("Unable to open particle file.");
   }

   // Open and read the dataset containing the number of particles.
   int num_particles;
   hid_t dataset_id, dataspace_id, memspace_id;
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
      LOKI_ABORT("Can not open dataset \"num_particles\".");
   }
   else {
      dataspace_id = H5Dget_space(dataset_id);
      if (dataspace_id < 0) {
         LOKI_ABORT("Can not get dataspace for dataset \"num_particles\".");
      }
      errf = H5Dread(dataset_id,
                     H5T_NATIVE_INT,
                     H5S_ALL,
                     H5S_ALL,
                     H5P_DEFAULT,
                     &num_particles);
      if (errf < 0) {
         LOKI_ABORT("Can not read dataset \"num_particles\".");
      }

      // Close the dataspace.
      errf = H5Sclose(dataspace_id);
      if (errf < 0) {
         LOKI_ABORT("Can not close dataspace for dataset \"num_particles\".");
      }

      // Close the dataset.
      errf = H5Dclose(dataset_id);
      if (errf < 0) {
         LOKI_ABORT("Can not close dataset \"num_particles\".");
      }
   }

   // If this not on a processor doing electromagnetics we don't read the
   // particles themselves and we're done.
   if (!a_is_em_process) {
      // Close the file.
      errf = H5Fclose(file_id);
      if (errf < 0) {
         LOKI_ABORT("Can not close particle file.");
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
      LOKI_ABORT("Can not open dataset \"particle_data\".");
   }
   else {
      // Check that the dataset is the right size.
      int vals_per_particle, data_size;
      if (a_read_noise_source_weight) {
         vals_per_particle = 8;
      }
      else {
         vals_per_particle = 7;
      }
      data_size = vals_per_particle*num_particles;
      dataspace_id = H5Dget_space(dataset_id);
      if (dataspace_id < 0) {
         LOKI_ABORT("Can not get dataspace for dataset \"particle_data\".");
      }
      hsize_t nsel = H5Sget_select_npoints(dataspace_id);
      if (static_cast<int>(nsel) != data_size) {
         LOKI_ABORT("Incorrect amount of particle data.");
      }

      // Now, specify the size and shape of each particle's data.
      hsize_t particle_dim[2];
      particle_dim[0] = 1;
      particle_dim[1] = vals_per_particle;

      hsize_t stride[2];
      stride[0] = 1;
      stride[1] = 1;

      hsize_t block[2];
      block[0] = 1;
      block[1] = 1;

      hsize_t offset[2];
      offset[0] = 0;
      offset[1] = 0;

      // Define the memory dataspace.
      memspace_id = H5Screate_simple(2, particle_dim, NULL);

      // Define the memory hyperslab.
      errf = H5Sselect_hyperslab(memspace_id,
                                 H5S_SELECT_SET,
                                 offset,
                                 NULL,
                                 particle_dim,
                                 NULL);
      if (errf < 0) {
         LOKI_ABORT("Cannot define particle hyperslab.");
      }

      // For each particle select the subset from the file dataspace and read
      // the particle data.
      double particle_data[8];
      for (int i = 0; i < num_particles; ++i) {
         offset[0] = i;
         errf = H5Sselect_hyperslab(dataspace_id,
                                    H5S_SELECT_SET,
                                    offset,
                                    NULL,
                                    particle_dim,
                                    NULL);
         if (errf < 0) {
            LOKI_ABORT("Unable to select particle data hyperslab");
         }
         errf = H5Dread(dataset_id,
                        H5T_NATIVE_DOUBLE,
                        memspace_id,
                        dataspace_id,
                        H5P_DEFAULT,
                        particle_data);
         if (errf < 0) {
            LOKI_ABORT("Can not read dataset \"particle_data\".");
         }

         double start = particle_data[0];
         double charge = particle_data[1];
         double mass = particle_data[2];
         double xinit = particle_data[3];
         double yinit = particle_data[4];
         double vxinit = particle_data[5];
         double vyinit = particle_data[6];
         double noise_source_weight;
         if (a_read_noise_source_weight) {
            noise_source_weight = particle_data[7];
         }
         else {
            noise_source_weight = 0.0;
         }
         // Let's make sure that the initial location of the particle is inside
         // the physical domain.
         if (xinit < a_domain->lower(X1) ) {
            if (a_domain->isPeriodic(X1)) {
               int domain_blocks =
                  int(1.0 + (a_domain->lower(X1) - xinit)/
                      (a_domain->upper(X1) - a_domain->lower(X1)));
               xinit +=
                  domain_blocks*(a_domain->upper(X1) - a_domain->lower(X1));
            }
            else {
               xinit = a_domain->lower(X1);
            }
         }
         else if (xinit > a_domain->upper(X1)) {
            if (a_domain->isPeriodic(X1)) {
               int domain_blocks =
                  int(1.0 + (xinit - a_domain->upper(X1))/
                      (a_domain->upper(X1) - a_domain->lower(X1)));
               xinit -=
                  domain_blocks*(a_domain->upper(X1) - a_domain->lower(X1));
            }
            else {
               xinit = a_domain->upper(X1);
            }
         }
         if (yinit < a_domain->lower(X2) ) {
            if (a_domain->isPeriodic(X2)) {
               int domain_blocks =
                  int(1.0 + (a_domain->lower(X2) - yinit)/
                      (a_domain->upper(X2) - a_domain->lower(X2)));
               yinit +=
                  domain_blocks*(a_domain->upper(X2) - a_domain->lower(X2));
            }
            else {
               yinit = a_domain->lower(X2);
            }
         }
         else if (yinit > a_domain->upper(X2)) {
            if (a_domain->isPeriodic(X2)) {
               int domain_blocks =
                  int(1.0 + (yinit - a_domain->upper(X2))/
                      (a_domain->upper(X2) - a_domain->lower(X2)));
               yinit -=
                  domain_blocks*(a_domain->upper(X2) - a_domain->lower(X2));
            }
            else {
               yinit = a_domain->upper(X2);
            }
         }

         vxinit = max(vxinit, a_vmin[0]);
         vxinit = min(vxinit, a_vmax[0]);

         vyinit = max(vyinit, a_vmin[1]);
         vyinit = min(vyinit, a_vmax[1]);

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
   }

   // Close the dataset.
   errf = H5Dclose(dataset_id);
   if (errf < 0) {
      LOKI_ABORT("Can not close dataset \"particle_data\".");
   }

   // Close the dataspace.
   errf = H5Sclose(dataspace_id);
   if (errf < 0) {
      LOKI_ABORT("Can not close dataspace for dataset \"particle_data\".");
   }

   // Close the memory dataspace.
   errf = H5Sclose(memspace_id);
   if (errf < 0) {
      LOKI_ABORT("Can not close memory dataspace.");
   }

   // Close the file.
   errf = H5Fclose(file_id);
   if (errf < 0) {
      LOKI_ABORT("Can not close particle file.");
   }

   return num_particles;
}

} // end namespace Loki
