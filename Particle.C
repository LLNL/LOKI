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
#include "Particle.H"

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
   : m_starting_time(0.0),
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

} // end namespace Loki
