/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "MaxwellianThermal.H"
#include "Simulation.H"
#include <math.h>

namespace Loki {

MaxwellianThermal::MaxwellianThermal(
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
   double pi = 4.0*atan(1.0);
   if (Simulation::s_DO_RELATIVITY) {
      m_fnorm = 1.0/(2.0*pi*m_thermalx*m_thermaly);
   }
   else {
      m_fnorm = a_mass/(2.0*pi*sqrt(a_tx*a_ty));
   }
}


MaxwellianThermal::~MaxwellianThermal()
{
}


double
MaxwellianThermal::thermalFactor(
   double a_spatial_factor,
   double a_x3,
   double a_x4) const
{
   double f;
   double x3 = m_x0*a_spatial_factor+m_flowinitx;
   double x4 = m_y0*a_spatial_factor+m_flowinity;
   if (Simulation::s_DO_RELATIVITY) {
      f = exp(-0.5*(pow((a_x3-x3)/m_thermalx, 2)+pow((a_x4-x4)/m_thermaly, 2)));
   }
   else {
      f = exp(-0.5*(pow(a_x3-x3, 2)/m_thermalx+pow(a_x4-x4, 2)/m_thermaly));
   }
   return f;
}
  
} // end namespace Loki
