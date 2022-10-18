/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include "tbox/Pointer.H"
#include "ProblemDomain.H"

namespace Loki {

/**
 * Glorified struct holding the information needed to describe a charged
 * particle.  We use these as tracking and noise source particles.  This is not
 * a PIC code so these do not represent the particles making up the plasma.
 */
class Particle
{
public:
   /**
    * @brief Default constructor.
    */
   Particle();

   /**
    * @brief Constructor.
    *
    * @param[in] a_starting_time For tracking particles the time at which
    *                            particles start to move.  For noise source
    *                            particles the time over which their strength
    *                            ramps up.
    * @param[in] a_charge Charge
    * @param[in] a_mass Mass
    * @param[in] a_xinit Initial x position
    * @param[in] a_yinit Initial y position
    * @param[in] a_vxinit Initial vx position
    * @param[in] a_vyinit Initial vy position
    * @param[in] a_noise_source_weight Weighting of noise source particles.
    */
   Particle(
      double a_starting_time,
      double a_charge,
      double a_mass,
      double a_xinit,
      double a_yinit,
      double a_vxinit,
      double a_vyinit,
      double a_noise_source_weight);

   /**
    * @brief Copy constructor.
    *
    * @param[in] a_other The thing to copy.
    */
   Particle(
      const Particle& a_other);

   /**
    * @brief Descructor.
    */
   ~Particle();

   // Assignment operator.
   Particle&
   operator = (
      const Particle& a_rhs)
   {
      if (this != &a_rhs) {
         m_starting_time = a_rhs.m_starting_time;
         m_charge = a_rhs.m_charge;
         m_mass = a_rhs.m_mass;
         m_x = a_rhs.m_x;
         m_y = a_rhs.m_y;
         m_vx = a_rhs.m_vx;
         m_vy = a_rhs.m_vy;
         m_noise_source_weight = a_rhs.m_noise_source_weight;
      }
      return *this;
   }

   /**
    * @brief Accessor for starting time.
    */
   const double&
   startingTime() const
   {
      return m_starting_time;
   }

   /**
    * @brief Accessor for charge.
    */
   const double&
   charge() const
   {
      return m_charge;
   }

   /**
    * @brief Accessor for mass.
    */
   const double&
   mass() const
   {
      return m_mass;
   }

   /**
    * @brief Accessor for x position.
    */
   double&
   x()
   {
      return m_x;
   }

   /**
    * @brief Accessor for x position.
    */
   const double&
   x() const
   {
      return m_x;
   }

   /**
    * @brief Accessor for y position.
    */
   double&
   y()
   {
      return m_y;
   }

   /**
    * @brief Accessor for y position.
    */
   const double&
   y() const
   {
      return m_y;
   }

   /**
    * @brief Accessor for x velocity.
    */
   double&
   vx()
   {
      return m_vx;
   }

   /**
    * @brief Accessor for x velocity.
    */
   const double&
   vx() const
   {
      return m_vx;
   }

   /**
    * @brief Accessor for y velocity.
    */
   double&
   vy()
   {
      return m_vy;
   }

   /**
    * @brief Accessor for y velocity.
    */
   const double&
   vy() const
   {
      return m_vy;
   }

   /**
    * @brief Accessor for noise source weight.
    */
   double&
   noiseSourceWeight()
   {
      return m_noise_source_weight;
   }

   /**
    * @brief Accessor for noise source weight.
    */
   const double&
   noiseSourceWeight() const
   {
      return m_noise_source_weight;
   }

   /**
    * @brief Adds a multiple of a vector of particles to to another.
    *
    * @param[in,out] a_particles The particles to be incremented.
    * @param[in] a_increment_particles The particle increment.
    * @param[in] a_factor Factor applied to a_increment_particles.
    * @param[in] a_sum_reduce_inc If true, sum reduce the particle increment.
    *                             This is necessary at different points in the
    *                             RK schemes.
    * @param[in] a_comm The communicator the sum reduction is performed on.
    */
   static
   void
   addParticleData(
      vector<Particle>& a_particles,
      const vector<Particle>& a_increment_particles,
      double a_factor,
      bool a_sum_reduce_inc,
      MPI_Comm a_comm);

   /**
    * @brief Copy the information from a_rhs_particles into a_particles.
    *
    * @param a_particles Particle information destination.
    * @param a_rhs_particles Particle information source.
    */
   static
   void
   copyParticleData(
      vector<Particle>& a_particles,
      const vector<Particle>& a_rhs_particles);

   /**
    * @brief Set the boundary conditions appropriate to the particles.
    *
    * @param a_particles The particles.
    * @param a_domain Information about the configuration space extents of the
    *                 problem so that the particles are properly bounded.
    * @param a_vmin Lower bound of the velocity extents of the problem.
    * @param a_vmax Upper bound of the velocity extents of the problem.
    * @param a_reflect If true reflect particles crossing physical boundary.
   */
   static
   void
   setBC(
      vector<Particle>& a_particles,
      const tbox::Pointer<ProblemDomain>& a_domain,
      const double* a_vmin,
      const double* a_vmax,
      bool a_reflect);

   /**
    * @brief Write the supplied particles to a particle restart file.
    *
    * @param a_particle_file_name The base part of the particle restart file.
    * @param a_particles The particles whose data is to be written.
    * @param a_write_noise_source_weight If true these are noise source
    *                                    particles.
    */
   static
   void
   putParticlesToRestart(const string& a_particle_file_name,
      const vector<Particle>& a_particles,
      bool a_write_noise_source_weight);

   /**
    * @brief Read the particle data into the supplied particles.
    *
    * @param a_particle_file The base part of the particle restart file.
    * @param a_particles The container for the particles whose data is to be
    *                    read.
    * @param a_domain Information about the configuration space extents of the
    *                 problem so that the particles are properly bounded.
    * @param a_vmin Lower bound of the velocity extents of the problem.
    * @param a_vmax Upper bound of the velocity extents of the problem.
    * @param a_is_em_process True if this is a Poisson/Maxwell process.
    * @param a_read_noise_source_weight If true these are noise source
    *                                   particles.
    *
    * @return The number of particles read.
    */
   static
   int
   readParticleData(string& a_particle_file,
      vector<Particle>& a_particles,
      const tbox::Pointer<ProblemDomain>& a_domain,
      const double* a_vmin,
      const double* a_vmax,
      bool a_is_em_process,
      bool a_read_noise_source_weight);

private:
   // Starting time.  For tracking particles the time at which particles start
   // to move.  For noise source particles the time over which their strength
   // ramps up.
   double m_starting_time;

   // Charge.
   double m_charge;

   // Mass.
   double m_mass;

   // x position.
   double m_x;

   // y position.
   double m_y;

   // x velocity.
   double m_vx;

   // y velocity.
   double m_vy;

   // Noise source weight.  Effective contribution of this particle to the
   // net charge density is m_noise_source_weight*m_charge.
   double m_noise_source_weight;
};

} // end namespace Loki

#endif
