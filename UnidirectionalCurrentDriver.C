/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "UnidirectionalCurrentDriver.H"
#include "UnidirectionalCurrentDriverF.H"
#include "Loki_Defines.H"
#include "Loki_Utilities.H"
#include "Simulation.H"

namespace Loki {

const string UnidirectionalCurrentDriver::s_CLASS_NAME(
   "Unidirectional Current Driver");


bool
UnidirectionalCurrentDriver::isType(
   const string& a_name)
{
   if (a_name.compare(s_CLASS_NAME) == 0) {
      return true;
   }
   return false;
}


UnidirectionalCurrentDriver::UnidirectionalCurrentDriver(
   LokiInputParser& a_pp)
   : m_R(12),
     m_parameters(NUM_PARAMS)
{
   // The default is no shaping of the antenna.
   m_parameters[SHAPE] = 0.0;

   // Get what the user wants.
   parseParameters(a_pp);

   // Compute constant geometric terms from what the user has specified.
   m_parameters[J0_INPLANE] = m_parameters[J0]*cos(m_parameters[POLAR_ANGLE]);
   m_parameters[J0_OUTPLANE] = m_parameters[J0]*sin(m_parameters[POLAR_ANGLE]);
   m_parameters[L] = sqrt(pow(m_parameters[X1X]-m_parameters[X0X], 2)+
                          pow(m_parameters[X1Y]-m_parameters[X0Y], 2));
   m_parameters[PHI] = atan2(m_parameters[X1Y]-m_parameters[X0Y],
                             m_parameters[X1X]-m_parameters[X0X]);
   double pi = 4.0*atan(1.0);
   double theta = m_parameters[PHI] - pi/2.0;
   double sin_theta = sin(theta);
   double cos_theta = cos(theta);
   m_R[0] = -sin_theta*m_parameters[LIGHT_SPEED];
   m_R[1] = cos_theta*m_parameters[LIGHT_SPEED];
   m_R[2] = 0.0;
   m_R[3] = 0.0;
   m_R[4] = 0.0;
   m_R[5] = 1.0;
   m_R[6] = 0.0;
   m_R[7] = 0.0;
   m_R[8] = -m_parameters[LIGHT_SPEED];
   m_R[9] = -sin_theta;
   m_R[10] = cos_theta;
   m_R[11] = 0.0;
}


UnidirectionalCurrentDriver::~UnidirectionalCurrentDriver()
{
}


void
UnidirectionalCurrentDriver::evaluate(
   ParallelArray& a_antenna_source,
   const ParallelArray& a_omega_eff2,
   const ProblemDomain& a_domain,
   double a_time) const
{
   // If the driver is on at a_time, compute the antenna source by delegating
   // the work to fortran.
   if (m_parameters[T0] <= a_time &&
       a_time < m_parameters[T0] + m_parameters[TRAMPUP] +
                m_parameters[THOLD] + m_parameters[TRAMPDOWN]) {
      int status;
      evaluateUnidirectionalCurrentDriver(*a_antenna_source.getData(),
         status,
         *a_omega_eff2.getData(),
         BOX2D_TO_FORT(a_antenna_source.dataBox()),
         BOX2D_TO_FORT(a_antenna_source.interiorBox()),
         PROBLEMDOMAIN_TO_FORT(a_domain),
         a_time,
         m_parameters[0],
         m_R[0]);
      if (status != 0) {
         LOKI_ABORT("Antenna wave becomes evanescent.");
      }
   }
}


void
UnidirectionalCurrentDriver::printParameters() const
{
   // Write the driver parameters.
   Loki_Utilities::printF("  Using external current forcing:\n" );
   Loki_Utilities::printF("    omega       = %e\n", m_parameters[OMEGA]);
   Loki_Utilities::printF("    J0          = %e\n", m_parameters[J0]);
   Loki_Utilities::printF("    t0          = %e\n", m_parameters[T0]);
   Loki_Utilities::printF("    t_rampup    = %e\n", m_parameters[TRAMPUP]);
   Loki_Utilities::printF("    t_hold      = %e\n", m_parameters[THOLD]);
   Loki_Utilities::printF("    t_rampdown  = %e\n", m_parameters[TRAMPDOWN]);
   Loki_Utilities::printF("    polar_angle = %e\n", m_parameters[POLAR_ANGLE]);
   Loki_Utilities::printF("    from_pt     = (%e, %e)\n",
      m_parameters[X0X], m_parameters[X0Y]);
   Loki_Utilities::printF("    to_pt       = (%e, %e)\n",
      m_parameters[X1X], m_parameters[X1Y]);
   Loki_Utilities::printF("    fwhm        = %e\n",
      2.0*sqrt(2.0*log(2))/m_parameters[BETA]);
   if (m_parameters[SHAPE] != 0.0) {
      Loki_Utilities::printF("    shape       = %e\n", m_parameters[SHAPE]);
      Loki_Utilities::printF("    width       = %e\n", m_parameters[WIDTH]);
   }
   else {
      Loki_Utilities::printF("    no antenna shaping\n");
   }
}


bool
UnidirectionalCurrentDriver::needsDispersionRelation()
{
   return true;
}


void
UnidirectionalCurrentDriver::parseParameters(
   LokiInputParser& a_pp)
{
   // All parameters are required.
   if (a_pp.contains("omega")) {
      a_pp.get("omega", m_parameters[OMEGA]);
   }
   else {
      LOKI_ABORT("Must supply omega!");
   }
   if (a_pp.contains("J_0")) {
      a_pp.get("J_0", m_parameters[J0]);
   }
   else {
      LOKI_ABORT("Must supply J_0!");
   }
   if (a_pp.contains("t0")) {
      a_pp.get("t0", m_parameters[T0]);
   }
   else {
      LOKI_ABORT("Must supply t0!");
   }
   // Check for old t_ramp, t_off syntax.
   if (a_pp.contains("t_ramp")) {
      // Map old syntax into new.
      a_pp.get("t_ramp", m_parameters[TRAMPUP]);
      m_parameters[THOLD] = 0.0;
      if (a_pp.contains("t_off")) {
         a_pp.get("t_off", m_parameters[TRAMPDOWN]);
      }
      else {
         LOKI_ABORT("Must supply t_off!");
      }
   }
   else if (a_pp.contains("t_rampup")) {
      // Read new syntax.
      a_pp.get("t_rampup", m_parameters[TRAMPUP]);
      if (a_pp.contains("t_hold")) {
         a_pp.get("t_hold", m_parameters[THOLD]);
      }
      else {
         LOKI_ABORT("Must supply t_hold!");
      }
      if (a_pp.contains("t_rampdown")) {
         a_pp.get("t_rampdown", m_parameters[TRAMPDOWN]);
      }
      else {
         LOKI_ABORT("Must supply t_rampdown!");
      }
   }
   else {
      LOKI_ABORT("Must supply either old t_ramp or new t_rampup syntax!");
   }
   if (a_pp.contains("polar_angle")) {
      a_pp.get("polar_angle", m_parameters[POLAR_ANGLE]);
   }
   else {
      LOKI_ABORT("Must supply polar_angle!");
   }
   vector<double> tmp(2);
   if (!a_pp.queryarr("from_pt", tmp, 0, 2)) {
      LOKI_ABORT("Must supply from_pt!");
   }
   else if (a_pp.countval("from_pt") != 2) {
      LOKI_ABORT("from_pt must have exactly 2 entries.");
   }
   m_parameters[X0X] = tmp[0];
   m_parameters[X0Y] = tmp[1];
   if (!a_pp.queryarr("to_pt", tmp, 0, 2)) {
      LOKI_ABORT("Must supply to_pt!");
   }
   else if (a_pp.countval("to_pt") != 2) {
      LOKI_ABORT("to_pt must have exactly 2 entries.");
   }
   m_parameters[X1X] = tmp[0];
   m_parameters[X1Y] = tmp[1];
   if (a_pp.contains("fwhm")) {
      double fwhm;
      a_pp.get("fwhm", fwhm);
      m_parameters[BETA] = 2.0*sqrt(2.0*log(2))/fwhm;
   }
   else {
      LOKI_ABORT("Must supply fwhm");
   }

   // Cram the light speed into the array of parameters.
   m_parameters[LIGHT_SPEED] = Simulation::s_LIGHT_SPEED;

   // See if the user wants any shaping of the antenna.
   a_pp.query("shape", m_parameters[SHAPE]);

   // If there is shaping of the antenna then the user must specify the width
   // over which that shaping occurs.
   if (a_pp.contains("width")) {
      a_pp.get("width", m_parameters[WIDTH]);
   }
   else if (m_parameters[SHAPE] != 0.0) {
      LOKI_ABORT("Must supply width with antenna shaping!");
   }
}

} // end namespace Loki
