/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "ThermalBehavior.H"
#include "Simulation.H"
#include "getMomentumF.H"
#include <math.h>

namespace Loki {

ThermalBehavior::ThermalBehavior(
   double a_vflowinitx,
   double a_vflowinity,
   double a_vx0,
   double a_vy0,
   double a_mass,
   double a_tx,
   double a_ty)
{
   if (Simulation::s_DO_RELATIVITY) {
      FORT_GET_MOMENTUM(a_vflowinitx,
         a_vflowinity,
         a_mass,
         Simulation::s_LIGHT_SPEED,
         m_flowinitx,
         m_flowinity);
      FORT_GET_MOMENTUM(a_vx0,
         a_vy0,
         a_mass,
         Simulation::s_LIGHT_SPEED,
         m_x0,
         m_y0);
      double vthermalx = sqrt(a_tx/a_mass);
      double vthermaly = sqrt(a_ty/a_mass);
      double tmp;
      FORT_GET_MOMENTUM(vthermalx,
         0.0,
         a_mass,
         Simulation::s_LIGHT_SPEED,
         m_thermalx,
         tmp);
      FORT_GET_MOMENTUM(0.0,
         vthermaly,
         a_mass,
         Simulation::s_LIGHT_SPEED,
         tmp,
         m_thermaly);
   }
   else {
      m_flowinitx = a_vflowinitx;
      m_flowinity = a_vflowinity;
      m_x0 = a_vx0;
      m_y0 = a_vy0;
      m_thermalx = a_tx/a_mass;
      m_thermaly = a_ty/a_mass;
   }
}


ThermalBehavior::~ThermalBehavior()
{
}

} // end namespace Loki
