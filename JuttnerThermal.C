/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "JuttnerThermal.H"
#include "Simulation.H"
#include "Loki_Defines.H"
#include <math.h>

namespace Loki {

JuttnerThermal::JuttnerThermal(
   double a_vflowinitx,
   double a_vflowinity,
   double a_vx0,
   double a_vy0,
   double a_mass,
   double a_tx,
   double a_ty)
   : ThermalBehavior(a_vflowinitx,
        a_vflowinity,
        a_vx0,
        a_vy0,
        a_mass,
        a_tx,
        a_ty)
{
   if (!Simulation::s_DO_RELATIVITY) {
     LOKI_ABORT("Juttner initial condition requires a relativistic problem.");
   }

   double pi = 4.0*atan(1.0);
   double c = Simulation::s_LIGHT_SPEED;
   if (c < pow(10, 10)) {
      m_fnorm = pow(c, 2)/(2.0*pi*(1.0+pow(c, 2)));
   }
   else {
      m_fnorm = 1.0/(2.0*pi);
   }
}


JuttnerThermal::~JuttnerThermal()
{
}


double
JuttnerThermal::thermalFactor(
   double a_spatial_factor,
   double a_x3,
   double a_x4) const
{
   double f;
   double x3 = m_x0*a_spatial_factor+m_flowinitx;
   double x4 = m_y0*a_spatial_factor+m_flowinity;
   double c2 = pow(Simulation::s_LIGHT_SPEED, 2);
   double gam = sqrt(1.0+(pow((a_x3-x3), 2)+pow((a_x4-x4), 2))/c2);
   f = exp(c2-c2*gam);
   return f;
}
  
} // end namespace Loki
