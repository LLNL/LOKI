/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "GaussianCurrentDriver.H"
#include "GaussianCurrentDriverF.H"
#include "Loki_Defines.H"
#include "Loki_Utilities.H"
#include "Simulation.H"

namespace Loki {

const string GaussianCurrentDriver::s_CLASS_NAME(
   "Gaussian Current Driver");


bool
GaussianCurrentDriver::isType(
   const string& a_name)
{
   if (a_name.compare(s_CLASS_NAME) == 0) {
      return true;
   }
   return false;
}


GaussianCurrentDriver::GaussianCurrentDriver(
   LokiInputParser& a_pp,
   const ProblemDomain& a_domain)
   : m_parameters(NUM_PARAMS)
{
   // Read the input parameters.
   parseParameters(a_pp);

   // Compute the in-plane current shaping perpendicular to the current
   // direction.
   int shape_length;
   double x2Ds, ylo, dy, yCenter;
   if (m_parameters[PLANE] == 0) {
      if (m_parameters[X0] < a_domain.lower(0) ||
          m_parameters[X0] > a_domain.upper(0)) {
         LOKI_ABORT("Antenna is placed outside the domain.");
      }
      shape_length = a_domain.numberOfCells(1);
      m_shape_lo = a_domain.box().lower(1);
      m_shape_hi = a_domain.box().upper(1);
      m_x_apply = (m_parameters[X0]-a_domain.lower(0))/a_domain.dx(0);
      x2Ds = a_domain.lower(0) + a_domain.dx(0)*(m_x_apply+0.5) -
         m_parameters[FOCAL_PT];
      ylo = a_domain.lower(1);
      dy = a_domain.dx(1);
      yCenter = (a_domain.upper(1) + a_domain.lower(1)) / 2.0;
   }
   else {
      if (m_parameters[X0] < a_domain.lower(1) ||
          m_parameters[X0] > a_domain.upper(1)) {
         LOKI_ABORT("Antenna is placed outside the domain.");
      }
      shape_length = a_domain.numberOfCells(0);
      m_shape_lo = a_domain.box().lower(0);
      m_shape_hi = a_domain.box().upper(0);
      m_x_apply = (m_parameters[X0]-a_domain.lower(1))/a_domain.dx(1);
      x2Ds = a_domain.lower(1) + a_domain.dx(1)*(m_x_apply+0.5) -
         m_parameters[FOCAL_PT];
      ylo = a_domain.lower(0);
      dy = a_domain.dx(0);
      yCenter = (a_domain.upper(0) + a_domain.lower(0)) / 2.0;
   }
   m_shaping.resize(shape_length);
   m_phase.resize(shape_length);
   double x2Ds0 = x2Ds == 0 ? 1.0 : 0.0;
   double x2Dsn0 = x2Ds != 0 ? 1.0 : 0.0;
   double pi = 4.0*atan(1.0);
   double w0 = m_parameters[FWHM]/sqrt(2.0*log(2.0));
   double n = sqrt(1.0-m_parameters[NENC]);
   double xR = pi*pow(w0, 2)*n/m_parameters[LAMBDA];
   double k = 2.0*pi*n/m_parameters[LAMBDA];
   double wx = w0*sqrt(1.0+pow(x2Ds/xR, 2));
   double wx2 = pow(wx, 2);
   double Rx = x2Ds*(1.0+pow(xR/(x2Ds+x2Ds0), 2));
   double amp = sqrt(w0/wx);
   double phase_factor = 0.5*k*x2Dsn0/(Rx+x2Ds0);
   for (int i = 0; i < shape_length; ++i) {
      double y2D = ylo + dy*(i+0.5) - yCenter;
      double y2D2 = pow(y2D, 2);
      m_shaping[i] = amp*exp(-y2D2/wx2);
      m_phase[i] = y2D2*phase_factor;
   }
}


GaussianCurrentDriver::~GaussianCurrentDriver()
{
}


void
GaussianCurrentDriver::evaluate(
   ParallelArray& a_antenna_source,
   const ParallelArray& a_omega_eff2,
   const ProblemDomain& a_domain,
   double a_time) const
{
   // If the driver is on at a_time, compute the antenna source by delegating
   // the work to fortran.
   const ParallelArray::Box& ib = a_antenna_source.interiorBox();
   bool source_in_domain =
      (m_parameters[PLANE] == 0 &&
       (m_x_apply >= ib.lower(0) && m_x_apply <= ib.upper(0))) ||
      (m_parameters[PLANE] == 1 &&
       (m_x_apply >= ib.lower(1) && m_x_apply <= ib.upper(1)));
   bool time_for_antenna = m_parameters[T0] <= a_time &&
      a_time < m_parameters[T0] + m_parameters[TRAMPUP] +
               m_parameters[THOLD] + m_parameters[TRAMPDOWN];
   if (time_for_antenna && source_in_domain) {
      evaluateGaussianCurrentDriver(*a_antenna_source.getData(),
         BOX2D_TO_FORT(a_antenna_source.dataBox()),
         BOX2D_TO_FORT(ib),
         PROBLEMDOMAIN_TO_FORT(a_domain),
         a_time,
         m_shape_lo,
         m_shape_hi,
         m_parameters[0],
         m_shaping[0],
         m_phase[0],
         m_x_apply);
   }
}


void
GaussianCurrentDriver::printParameters() const
{
   // Write the driver parameters.
   Loki_Utilities::printF("  Using external current forcing:\n" );
   if (m_parameters[PLANE] == 0) {
      Loki_Utilities::printF("    x focal point               = %e\n",
         m_parameters[FOCAL_PT]);
      Loki_Utilities::printF("    y fwhm                      = %e\n",
         m_parameters[FWHM]);
   }
   else {
      Loki_Utilities::printF("    y focal point               = %e\n",
         m_parameters[FOCAL_PT]);
      Loki_Utilities::printF("    x fwhm                      = %e\n",
         m_parameters[FWHM]);
   }
   Loki_Utilities::printF("    frequency                   = %e\n",
      m_parameters[OMEGA]);
   Loki_Utilities::printF("    normalized electron density = %e\n",
      m_parameters[NENC]);
   Loki_Utilities::printF("    J0                          = %e\n",
      m_parameters[J0]);
   Loki_Utilities::printF("    t0                          = %e\n",
      m_parameters[T0]);
   Loki_Utilities::printF("    t_rampup                    = %e\n",
      m_parameters[TRAMPUP]);
   Loki_Utilities::printF("    t_hold                      = %e\n",
      m_parameters[THOLD]);
   Loki_Utilities::printF("    t_rampdown                  = %e\n",
      m_parameters[TRAMPDOWN]);
   if (m_parameters[PLANE] == 0) {
      Loki_Utilities::printF("    x0                          = %e\n",
         m_parameters[X0]);
   }
   else {
      Loki_Utilities::printF("    y0                          = %e\n",
         m_parameters[X0]);
   }
   Loki_Utilities::printF("    plane                          = %e\n",
      m_parameters[PLANE]);
}


bool
GaussianCurrentDriver::needsDispersionRelation()
{
   return false;
}


void
GaussianCurrentDriver::parseParameters(
   LokiInputParser& a_pp)
{
   // All parameters are required.
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

   if (a_pp.contains("focal_pt")) {
      a_pp.get("focal_pt", m_parameters[FOCAL_PT]);
   }
   else {
      LOKI_ABORT("Must supply focal_pt!");
   }

   if (a_pp.contains("fwhm")) {
      a_pp.get("fwhm", m_parameters[FWHM]);
   }
   else {
      LOKI_ABORT("Must supply fwhm!");
   }

   if (a_pp.contains("omega")) {
      a_pp.get("omega", m_parameters[OMEGA]);
   }
   else {
      LOKI_ABORT("Must supply omega!");
   }

   m_parameters[LAMBDA] =
      2.0*4.0*atan(1.0)*Simulation::s_LIGHT_SPEED/m_parameters[OMEGA];

   if (a_pp.contains("nenc")) {
      a_pp.get("nenc", m_parameters[NENC]);
   }
   else {
      LOKI_ABORT("Must supply nenc!");
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

   if (a_pp.contains("t_rampup")) {
      a_pp.get("t_rampup", m_parameters[TRAMPUP]);
   }
   else {
      LOKI_ABORT("Must supply t_rampup!");
   }

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

} // end namespace Loki
