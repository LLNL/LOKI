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
#ifndef _PARTICLE_H_
#define _PARTICLE_H_

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
    * @param[in] a_starting_time Time at which particles start to move.
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

private:
   // Starting time.
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
