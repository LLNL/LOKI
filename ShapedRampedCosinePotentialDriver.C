/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "ShapedRampedCosinePotentialDriver.H"
#include "ShapedRampedCosinePotentialDriverF.H"
#include "Loki_Defines.H"
#include "Loki_Utilities.H"

namespace Loki {

const string ShapedRampedCosinePotentialDriver::s_CLASS_NAME(
   "Shaped Ramped Cosine Potential Driver");


bool
ShapedRampedCosinePotentialDriver::isType(
   const string& a_name)
{
   if (a_name.compare(s_CLASS_NAME) == 0) {
      return true;
   }
   return false;
}


ShapedRampedCosinePotentialDriver::ShapedRampedCosinePotentialDriver(
   LokiInputParser& a_pp)
   : m_parameters(NUM_PARAMS)
{
   // Set some defaults which most likely don't make any sense.
   m_parameters[XWIDTH]    = 0.5;
   m_parameters[YWIDTH]    = 0.5;
   m_parameters[SHAPE]     = 1.0;
   m_parameters[OMEGA]     = 1.0;
   m_parameters[AMP]       = 0.01;
   m_parameters[T0]        = 0.0;
   m_parameters[TRAMPUP]   = 10.0;
   m_parameters[THOLD]     = 0.0;
   m_parameters[TRAMPDOWN] = 10.0;
   m_parameters[XSHAPE]    = 0.0;
   m_parameters[LWIDTH]    = 0.5;
   m_parameters[X0]        = 0.0;

   // Now get what the user really wants.
   parseParameters(a_pp);
}


ShapedRampedCosinePotentialDriver::~ShapedRampedCosinePotentialDriver()
{
}


void
ShapedRampedCosinePotentialDriver::evaluate(
   ParallelArray& a_phi,
   const ProblemDomain& a_domain,
   double a_time) const
{
   // Delegate evaluation to fortran.
   evaluateShapedRampedPotentialDriver(*a_phi.getData(),
      BOX2D_TO_FORT(a_phi.dataBox()),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      a_time,
      m_parameters[0]);
}


void
ShapedRampedCosinePotentialDriver::printParameters() const
{
   // Write all the input parameters.
   Loki_Utilities::printF("  Using external electric potential forcing:\n");
   Loki_Utilities::printF("    xwidth      = %e\n", m_parameters[XWIDTH]);
   Loki_Utilities::printF("    ywidth      = %e\n", m_parameters[YWIDTH]);
   Loki_Utilities::printF("    shape       = %e\n", m_parameters[SHAPE]);
   Loki_Utilities::printF("    omega       = %e\n", m_parameters[OMEGA]);
   Loki_Utilities::printF("    amp         = %e\n", m_parameters[AMP]);
   Loki_Utilities::printF("    t0          = %e\n", m_parameters[T0]);
   Loki_Utilities::printF("    t_rampup    = %e\n", m_parameters[TRAMPUP]);
   Loki_Utilities::printF("    t_hold      = %e\n", m_parameters[THOLD]);
   Loki_Utilities::printF("    t_rampdown  = %e\n", m_parameters[TRAMPDOWN]);
   Loki_Utilities::printF("    x_shape     = %e\n", m_parameters[XSHAPE]);
   Loki_Utilities::printF("    lwidth      = %e\n", m_parameters[LWIDTH]);
   Loki_Utilities::printF("    x0          = %e\n", m_parameters[X0]);
}


void
ShapedRampedCosinePotentialDriver::parseParameters(
   LokiInputParser& a_pp)
{
   // There's a bunch of backward compatibility here for older style syntax.
   // The parameters are all in terms of widths in configuration space
   // directions.  New syntax enters this directly.  Old syntax was in terms
   // of wave numbers and lengths which are converted to widths.

   // Read the information about x modulation.  None of this is required but it
   // probably should be as defaults are likely meaningless.
   if (a_pp.contains("kx") || a_pp.contains("Lx")) {
      if (a_pp.contains("xwidth")) {
         LOKI_ABORT("Mixed old and new syntax for kx/Lx");
      }
      else {
         double kx = 1.0;
         double lx = 1.0;
         a_pp.query("kx", kx);
         a_pp.query("Lx", lx);
         m_parameters[XWIDTH] = lx/(2.0*kx);
      }
   }
   else {
      a_pp.query("xwidth", m_parameters[XWIDTH]);
   }

   // Read the information about y modulation.  None of this is required but it
   // probably should be as defaults are likely meaningless.
   if (a_pp.contains("ky") || a_pp.contains("Ly")) {
      if (a_pp.contains("ywidth")) {
         LOKI_ABORT("Mixed old and new syntax for ky/Ly");
      }
      else {
         double ky = 1.0;
         double ly = 1.0;
         a_pp.query("ky", ky);
         a_pp.query("Ly", ly);
         m_parameters[YWIDTH] = ly/(2.0*ky);
      }
   }
   else {
      a_pp.query("ywidth", m_parameters[YWIDTH]);
   }

   // Read the information about x shape envelope.  None of this is required
   // but it probably should be as defaults are likely meaningless.
   if (a_pp.contains("kl") || a_pp.contains("Ll")) {
      if (a_pp.contains("lwidth")) {
         LOKI_ABORT("Mixed old and new syntax for kl/Ll");
      }
      else {
         double kl = 1.0;
         double ll = 1.0;
         a_pp.query("kl", kl);
         a_pp.query("Ll", ll);
         m_parameters[LWIDTH] = ll/(2.0*kl);
      }
   }
   else {
      a_pp.query("lwidth", m_parameters[LWIDTH]);
   }

   // These are not required but should be.
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
   a_pp.query("amp",     m_parameters[AMP]);
   a_pp.query("t0",      m_parameters[T0]);
   a_pp.query("x_shape", m_parameters[XSHAPE]);
   a_pp.query("x0",      m_parameters[X0]);
}

} // end namespace Loki
