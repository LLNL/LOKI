/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "ShapedRampedCosineCurrentDriver.H"
#include "ShapedRampedCosineCurrentDriverF.H"
#include "Loki_Defines.H"
#include "Loki_Utilities.H"

namespace Loki {

const string ShapedRampedCosineCurrentDriver::s_CLASS_NAME(
   "Shaped Ramped Cosine Current Driver");


bool
ShapedRampedCosineCurrentDriver::isType(
   const string& a_name)
{
   if (a_name.compare(s_CLASS_NAME) == 0) {
      return true;
   }
   return false;
}


ShapedRampedCosineCurrentDriver::ShapedRampedCosineCurrentDriver(
   LokiInputParser& a_pp)
   : m_parameters(NUM_PARAMS)
{
   // Size the parameters and set some hopefully sane defaults then get what
   // the user actually wants from the input file.
   m_parameters[SHAPE]      = 1.0;
   m_parameters[OMEGA]      = 1.0;
   m_parameters[J0]         = 0.01;
   m_parameters[T0]         = 0.0;
   m_parameters[TRAMPUP]    = 10.0;
   m_parameters[THOLD]      = 10.0;
   m_parameters[TRAMPDOWN]  = 10.0;
   m_parameters[X0]         = 0.0;
   m_parameters[PLANE]      = 0.0;
   parseParameters(a_pp);
}


ShapedRampedCosineCurrentDriver::~ShapedRampedCosineCurrentDriver()
{
}


void
ShapedRampedCosineCurrentDriver::evaluate(
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
      evaluateShapedRampedCurrentDriver(*a_antenna_source.getData(),
         BOX2D_TO_FORT(a_antenna_source.dataBox()),
         BOX2D_TO_FORT(a_antenna_source.interiorBox()),
         PROBLEMDOMAIN_TO_FORT(a_domain),
         a_time,
         m_parameters[0]);
   }
}


void
ShapedRampedCosineCurrentDriver::printParameters() const
{
   // Write the driver parameters.
   Loki_Utilities::printF("  Using external current forcing:\n" );
   Loki_Utilities::printF("    width           = %e\n", m_parameters[WIDTH]);
   Loki_Utilities::printF("    apply direction = %i\n",
      static_cast<int>(m_parameters[APPLY_DIR]));
   Loki_Utilities::printF("    shape           = %e\n", m_parameters[SHAPE]);
   Loki_Utilities::printF("    omega           = %e\n", m_parameters[OMEGA]);
   Loki_Utilities::printF("    J0              = %e\n", m_parameters[J0]);
   Loki_Utilities::printF("    t0              = %e\n", m_parameters[T0]);
   Loki_Utilities::printF("    t_rampup        = %e\n", m_parameters[TRAMPUP]);
   Loki_Utilities::printF("    t_hold          = %e\n", m_parameters[THOLD]);
   Loki_Utilities::printF("    t_rampdown      = %e\n",
      m_parameters[TRAMPDOWN]);
   if (m_parameters[PLANE] == 0) {
      Loki_Utilities::printF("    x0              = %e\n", m_parameters[X0]);
   }
   else {
      Loki_Utilities::printF("    y0              = %e\n", m_parameters[X0]);
   }
   Loki_Utilities::printF("    plane           = %e\n", m_parameters[PLANE]);
}


bool
ShapedRampedCosineCurrentDriver::needsDispersionRelation()
{
   return false;
}


void
ShapedRampedCosineCurrentDriver::parseParameters(
   LokiInputParser& a_pp)
{
   // The width, apply_dir, and one of x0 or y0 are required.
   if (a_pp.contains("width")) {
      a_pp.get("width", m_parameters[WIDTH]);
   }
   else {
      LOKI_ABORT("Must supply width!");
   }
   bool contains_x0 = a_pp.contains("x0");
   bool contains_y0 = a_pp.contains("y0");
   if (contains_x0 || contains_y0) {
      if (contains_x0) {
         a_pp.get("x0", m_parameters[X0]);
         m_parameters[PLANE] = 0.0;
      }
      else {
         a_pp.get("y0", m_parameters[X0]);
         m_parameters[PLANE] = 1.0;
      }
   }
   else {
      LOKI_ABORT("Must supply one of x0 or y0!");
   }
   if (a_pp.contains("apply_dir")) {
      string tmp;
      a_pp.get("apply_dir", tmp);
      if (tmp.compare("x") == 0 || tmp.compare("X") == 0) {
         if (contains_x0) {
            LOKI_ABORT("Can not create x current in a y-z current sheet!");
         }
         else {
            m_parameters[APPLY_DIR] = 1.0;
         }
      }
      else if (tmp.compare("y") == 0 || tmp.compare("Y") == 0) {
         if (contains_y0) {
            LOKI_ABORT("Can not create y current in a x-z current sheet!");
         }
         else {
            m_parameters[APPLY_DIR] = 2.0;
         }
      }
      else if (tmp.compare("z") == 0 || tmp.compare("Z") == 0) {
         m_parameters[APPLY_DIR] = 3.0;
      }
      else {
         LOKI_ABORT("Unknown apply_dir!");
      }
   }
   else {
      LOKI_ABORT("Must supply apply_dir!");
   }

   // Everything else is not required which may be a bad idea as it's not clear
   // what valid defaults should be.
   // Check for old t_ramp, t_off syntax.
   if (a_pp.contains("t_ramp")) {
      // Read old synatax and map into new syntax.
      a_pp.query("t_ramp", m_parameters[TRAMPUP]);
      m_parameters[THOLD] = 0.0;
      a_pp.query("t_off", m_parameters[TRAMPDOWN]);
   }
   else {
      // Read new synatax.
      a_pp.query("t_rampup", m_parameters[TRAMPUP]);
      a_pp.query("t_hold", m_parameters[THOLD]);
      a_pp.query("t_rampdown", m_parameters[TRAMPDOWN]);
   }
   a_pp.query("shape",   m_parameters[SHAPE]);
   a_pp.query("omega",   m_parameters[OMEGA]);
   a_pp.query("J_0",     m_parameters[J0]);
   a_pp.query("t0",      m_parameters[T0]);
}

} // end namespace Loki
