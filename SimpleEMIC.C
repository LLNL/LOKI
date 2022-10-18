/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "SimpleEMIC.H"
#include "SimpleEMICF.H"
#include "Loki_Utilities.H"

namespace Loki {

const string SimpleEMIC::s_CLASS_NAME("SimpleEMIC");

bool
SimpleEMIC::isType(
   const string& a_name)
{
   if (a_name.compare(s_CLASS_NAME) == 0) {
      return true;
   }
   return false;
}

SimpleEMIC::SimpleEMIC(
   LokiInputParser& a_pp)
   : m_iparameters(NUM_IPARAMS)
{
   // Read which field we're initializing.
   string tmp;
   if (!a_pp.query("field", tmp)) {
      LOKI_ABORT("Must supply field");
   }
   else if (tmp.compare("E") == 0 || tmp.compare("e") == 0) {
      m_iparameters[FIELD] = E;
   }
   else if (tmp.compare("B") == 0 || tmp.compare("b") == 0) {
      m_iparameters[FIELD] = B;
   }
   else {
      LOKI_ABORT("Unknown field");
   }

   // All this input is required as there's no obvious defaults.
   int num_xamp, num_yamp, num_zamp, num_kx, num_ky, num_phi;
   if (!a_pp.contains("xamp")) {
      LOKI_ABORT("Must supply xamp");
   }
   else {
      num_xamp = a_pp.countval("xamp");
      m_x_amp.resize(num_xamp);
      a_pp.queryarr("xamp", m_x_amp, 0, num_xamp);
   }
   if (!a_pp.contains("yamp")) {
      LOKI_ABORT("Must supply yamp");
   }
   else {
      num_yamp = a_pp.countval("yamp");
      m_y_amp.resize(num_yamp);
      a_pp.queryarr("yamp", m_y_amp, 0, num_yamp);
   }
   if (!a_pp.contains("zamp")) {
      LOKI_ABORT("Must supply zamp");
   }
   else {
      num_zamp = a_pp.countval("zamp");
      m_z_amp.resize(num_zamp);
      a_pp.queryarr("zamp", m_z_amp, 0, num_zamp);
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

   if (num_xamp != num_yamp || num_xamp != num_zamp ||
       num_xamp != num_kx || num_xamp != num_ky || num_xamp != num_phi) {
      LOKI_ABORT("Number of amplitude, wave numbers, and phases no not match.");
   }
   m_num_waves = num_xamp;
}


SimpleEMIC::~SimpleEMIC()
{
}


void
SimpleEMIC::set(
   ParallelArray& a_u,
   const ProblemDomain& a_domain) const
{
   // Delegate evaluation to fortran.
   SET_SIMPLE_EMIC(BOX2D_TO_FORT(a_u.dataBox()),
      *a_u.getData(),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      m_num_waves,
      m_x_amp[0],
      m_y_amp[0],
      m_z_amp[0],
      m_x_wave_number[0],
      m_y_wave_number[0],
      m_phi[0],
      m_iparameters[0]);
}


bool
SimpleEMIC::initializesE()
{
   return m_iparameters[FIELD] == E;
}


bool
SimpleEMIC::initializesB()
{
   return m_iparameters[FIELD] == B;
}


int
SimpleEMIC::numWaves()
{
   return m_num_waves;
}


void
SimpleEMIC::printParameters() const
{
   // Print all input parameters.
   Loki_Utilities::printF("  Using EM field initialization:\n");
   Loki_Utilities::printF("    field         = %s\n",
      m_iparameters[FIELD] == E ? "E" : "B");
   Loki_Utilities::printF("    X amplitude   = ");
   for (int i = 0; i < m_num_waves; ++i ) {
      Loki_Utilities::printF("%e ", m_x_amp[i]);
   }
   Loki_Utilities::printF("\n");
   Loki_Utilities::printF("    Y amplitude   = ");
   for (int i = 0; i < m_num_waves; ++i ) {
      Loki_Utilities::printF("%e ", m_y_amp[i]);
   }
   Loki_Utilities::printF("\n");
   Loki_Utilities::printF("    Z amplitude   = ");
   for (int i = 0; i < m_num_waves; ++i ) {
      Loki_Utilities::printF("%e ", m_z_amp[i]);
   }
   Loki_Utilities::printF("\n");
   Loki_Utilities::printF("    x wave number = ");
   for (int i = 0; i < m_num_waves; ++i ) {
      Loki_Utilities::printF("%e ", m_x_wave_number[i]);
   }
   Loki_Utilities::printF("\n");
   Loki_Utilities::printF("    y wave number = ");
   for (int i = 0; i < m_num_waves; ++i ) {
      Loki_Utilities::printF("%e ", m_y_wave_number[i]);
   }
   Loki_Utilities::printF("\n");
   Loki_Utilities::printF("    phase         = ");
   for (int i = 0; i < m_num_waves; ++i ) {
      Loki_Utilities::printF("%e ", m_phi[i]);
   }
   Loki_Utilities::printF("\n");
}

} // end namespace Loki
