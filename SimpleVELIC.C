/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "SimpleVELIC.H"
#include "SimpleVELICF.H"
#include "Loki_Utilities.H"

namespace Loki {

const string SimpleVELIC::s_CLASS_NAME("SimpleVELIC");

bool
SimpleVELIC::isType(
   const string& a_name)
{
   if (a_name.compare(s_CLASS_NAME) == 0) {
      return true;
   }
   return false;
}

SimpleVELIC::SimpleVELIC(
   LokiInputParser& a_pp)
{
   // All this input is required as there's no obvious defaults.
   int num_amp, num_kx, num_ky, num_phi;
   if (!a_pp.contains("amp")) {
      LOKI_ABORT("Must supply amp");
   }
   else {
      num_amp = a_pp.countval("amp");
      m_amp.resize(num_amp);
      a_pp.queryarr("amp", m_amp, 0, num_amp);
   }

   if (!a_pp.contains("x_wave_number")) {
      LOKI_ABORT("Must supply x_wave_number");
   }
   else {
      num_kx = a_pp.countval("x_wave_number");
      m_x_wave_number.resize(num_kx);
      a_pp.queryarr("x_wave_number", m_x_wave_number, 0, num_kx);
   }

   if (!a_pp.contains("y_wave_number")) {
      LOKI_ABORT("Must supply y_wave_number");
   }
   else {
      num_ky = a_pp.countval("y_wave_number");
      m_y_wave_number.resize(num_ky);
      a_pp.queryarr("y_wave_number", m_y_wave_number, 0, num_ky);
   }

   if (!a_pp.contains("phase")) {
      LOKI_ABORT("Must supply phase");
   }
   else {
      num_phi = a_pp.countval("phase");
      m_phi.resize(num_phi);
      a_pp.queryarr("phase", m_phi, 0, num_phi);
   }

   if (num_amp != num_kx || num_amp != num_ky || num_amp != num_phi) {
      LOKI_ABORT("Number of amplitude, wave numbers, and phases no not match.");
   }
   m_num_waves = num_amp;
}


SimpleVELIC::~SimpleVELIC()
{
}


void
SimpleVELIC::set(
   ParallelArray& a_vz,
   const ProblemDomain& a_domain) const
{
   // Delegate evaluation to fortran.
   SET_SIMPLE_VELIC(*a_vz.getData(),
      BOX2D_TO_FORT(a_vz.dataBox()),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      m_num_waves,
      m_amp[0],
      m_x_wave_number[0],
      m_y_wave_number[0],
      m_phi[0]);
}


int
SimpleVELIC::numWaves()
{
   return m_num_waves;
}


void
SimpleVELIC::printParameters() const
{
   // Print all input parameters.
   Loki_Utilities::printF("  Using velocity field initialization:\n");
   Loki_Utilities::printF("    amplitude     = ");
   for (int i = 0; i < m_num_waves; ++i) {
      Loki_Utilities::printF("%e ", m_amp[i]);
   }
   Loki_Utilities::printF("\n");
   Loki_Utilities::printF("    x wave number = ");
   for (int i = 0; i < m_num_waves; ++i) {
      Loki_Utilities::printF("%e ", m_x_wave_number[i]);
   }
   Loki_Utilities::printF("\n");
   Loki_Utilities::printF("    y wave number = ");
   for (int i = 0; i < m_num_waves; ++i) {
      Loki_Utilities::printF("%e ", m_y_wave_number[i]);
   }
   Loki_Utilities::printF("\n");
   Loki_Utilities::printF("    phase         = ");
   for (int i = 0; i < m_num_waves; ++i) {
      Loki_Utilities::printF("%e ", m_phi[i]);
   }
   Loki_Utilities::printF("\n");
}

} // end namespace Loki
