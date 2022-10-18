/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "ShapedRampedCosineDriver.H"
#include "ShapedRampedCosineDriverF.H"
#include "Loki_Defines.H"
#include "Loki_Utilities.H"

#include <math.h>
#include <sstream>

namespace Loki {

const string ShapedRampedCosineDriver::s_CLASS_NAME(
   "Shaped Ramped Cosine Driver");


ShapedRampedCosineDriver::ShapedRampedCosineDriver(
   LokiInputParser& a_pp,
   int a_driver_num,
   int a_species_num)
   : m_parameters(NUM_PARAMS),
     m_do_random_phase(false),
     m_phase(0.0),
     m_phase_h(0.0),
     m_num_phase_evals(0),
     m_next(1),
     m_shape_type(0)
{
   // Construct a unique name for this driver which is the name of the
   // sub-database in the restart file containing its state.
   stringstream ss;
   ss << a_species_num << "_" << a_driver_num;
   string tmp(ss.str());
   m_name = s_CLASS_NAME + tmp;

   // Set some hopefully sane defaults.
   m_parameters[XWIDTH]                 = 0.5;
   m_parameters[YWIDTH]                 = 0.5;
   m_parameters[SHAPE]                  = 1.0;
   m_parameters[OMEGA]                  = 1.0;
   m_parameters[E0]                     = 0.01;
   m_parameters[T0]                     = 0.0;
   m_parameters[TRAMPUP]                = 10.0;
   m_parameters[THOLD]                  = 0.0;
   m_parameters[TRAMPDOWN]              = 10.0;
   m_parameters[XSHAPE]                 = 0.0;
   m_parameters[LWIDTH]                 = 0.5;
   m_parameters[X0]                     = 0.0;
   m_parameters[ALPHA]                  = 0.0;
   m_parameters[TRES]                   = 0.0;
   m_parameters[PHASE_DECAY_TIME_STEPS] = 0.0;
   m_parameters[FWHM]                   = 0.0;

   // Get what the user really wants.
   parseParameters(a_pp);
}


ShapedRampedCosineDriver::~ShapedRampedCosineDriver()
{
}


void
ShapedRampedCosineDriver::evaluate(
   ParallelArray& a_em_vars,
   ParallelArray& a_ext_efield,
   const ProblemDomain& a_domain,
   int a_sums_into,
   double a_time,
   double a_dt,
   bool a_first_rk_stage) const
{
   // If the driver is on at a_time then delegate its evaluation to fortran.
   if (m_parameters[T0] <= a_time &&
       a_time < m_parameters[T0] + m_parameters[TRAMPUP] +
                m_parameters[THOLD] + m_parameters[TRAMPDOWN]) {
      computePhase(a_dt, a_first_rk_stage);
      evaluateShapedRampedDriver(*a_em_vars.getData(),
         *a_ext_efield.getData(),
         BOX2D_TO_FORT(a_em_vars.dataBox()),
         a_em_vars.dataBox().numberOfCells(2),
         PROBLEMDOMAIN_TO_FORT(a_domain),
         a_sums_into,
         a_time,
         m_parameters[0],
         m_phase,
         m_shape_type);
   }
}


void
ShapedRampedCosineDriver::evaluateTimeEnvelope(
   double& a_envel,
   double a_time) const
{
   double t0 = m_parameters[T0];
   double t_rampup = m_parameters[TRAMPUP];
   double t_hold = m_parameters[THOLD];
   double t_rampdown = m_parameters[TRAMPDOWN];
   double E_0 = m_parameters[E0];
   if (t0 <= a_time && a_time < t0 + t_rampup + t_hold + t_rampdown) {
      evaluateShapedRampedDriverEnvelope(a_envel,
         a_time,
         t0,
         t_rampup,
         t_hold,
         t_rampdown,
         E_0);
   }
   else {
      a_envel = 0.0;
   }
}


void
ShapedRampedCosineDriver::printParameters() const
{
   // Write all the driver parameters.
   Loki_Utilities::printF("  Using external electric field forcing:\n" );
   Loki_Utilities::printF("    xwidth                 = %e\n",
      m_parameters[XWIDTH]);
   Loki_Utilities::printF("    ywidth                 = %e\n",
      m_parameters[YWIDTH]);
   Loki_Utilities::printF("    shape                  = %e\n",
      m_parameters[SHAPE]);
   Loki_Utilities::printF("    omega                  = %e\n",
      m_parameters[OMEGA]);
   Loki_Utilities::printF("    E0                     = %e\n",
      m_parameters[E0]);
   Loki_Utilities::printF("    t0                     = %e\n",
      m_parameters[T0]);
   Loki_Utilities::printF("    t_rampup               = %e\n",
      m_parameters[TRAMPUP]);
   Loki_Utilities::printF("    t_hold                 = %e\n",
      m_parameters[THOLD]);
   Loki_Utilities::printF("    t_rampdown             = %e\n",
      m_parameters[TRAMPDOWN]);
   Loki_Utilities::printF("    x_shape                = %e\n",
      m_parameters[XSHAPE]);
   if (m_shape_type == 0) {
      Loki_Utilities::printF("    shape type             = sin squared\n");
   }
   else {
      Loki_Utilities::printF("    shape type             = exp\n");
   }
   Loki_Utilities::printF("    lwidth                 = %e\n",
      m_parameters[LWIDTH]);
   Loki_Utilities::printF("    x0                     = %e\n",
      m_parameters[X0]);
   Loki_Utilities::printF("    alpha                  = %e\n",
      m_parameters[ALPHA]);
   Loki_Utilities::printF("    t_res                  = %e\n",
      m_parameters[TRES]);
   Loki_Utilities::printF("    phase decay time steps = %e\n",
      m_parameters[PHASE_DECAY_TIME_STEPS]);
   Loki_Utilities::printF("    fwhm                   = %e\n",
      m_parameters[FWHM]);
}


void
ShapedRampedCosineDriver::getFromDatabase(
   RestartReader& a_reader)
{
   // Get the subdirectory containing this object's state.
   a_reader.pushSubDir(m_name);

   // If this is a noisy driver, get the number of times computePhase was
   // called from the database.
   int times_called = 0;
   if (m_do_random_phase) {
      a_reader.readIntegerValue("num_phase_evals", times_called);

      // Now get the phase and half time phase.
      a_reader.readDoubleValue("phase", m_phase);
      a_reader.readDoubleValue("phase_h", m_phase_h);

      // Now run the random number generator times_called times so that it is
      // in the same state as when the restart file was generated.  Once this
      // is done the computation of the random phase will proceed as if there
      // had been no restart.
      for (int i = 0; i < times_called; ++i) {
         rand();
      }
      m_num_phase_evals = times_called;
   }
   a_reader.popSubDir();
}


void
ShapedRampedCosineDriver::putToDatabase(
   RestartWriter& a_writer,
   bool a_write_data) const
{
   // Make the subdirectory to contain this object's state.
   a_writer.pushSubDir(m_name);

   // Put m_num_phase_evals, m_phase, and, m_phase_h into the database.
   a_writer.writeIntegerValue("num_phase_evals",
      m_num_phase_evals,
      a_write_data);
   a_writer.writeDoubleValue("phase", m_phase, a_write_data);
   a_writer.writeDoubleValue("phase_h", m_phase_h, a_write_data);

   a_writer.popSubDir();
}


bool
ShapedRampedCosineDriver::isType(
   const string& a_name)
{
   if (a_name.compare(s_CLASS_NAME) == 0) {
      return true;
   }
   return false;
}


void
ShapedRampedCosineDriver::parseParameters(
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
   string test_str = "sin2";
   a_pp.query("shape_type", test_str);
   if (test_str.compare("sin2") == 0) {
      m_shape_type = 0;
   }
   else if (test_str.compare("exp") == 0) {
      m_shape_type = 1;
   }
   else {
      LOKI_ABORT("Unknown shape type");
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
   a_pp.query("E_0",     m_parameters[E0]);
   a_pp.query("t0",      m_parameters[T0]);
   a_pp.query("x_shape", m_parameters[XSHAPE]);
   a_pp.query("x0",      m_parameters[X0]);
   a_pp.query("alpha",   m_parameters[ALPHA]);
   a_pp.query("t_res",   m_parameters[TRES]);

   // These 2 parameters define the spectrum of the noisy phases.  So users
   // must specify either both or neither.
   bool contains_decay = false;
   bool contains_fwhm = false;
   if (a_pp.contains("phase_decay_time_steps")) {
      contains_decay = true;
   }
   if (a_pp.contains("fwhm")) {
      contains_fwhm = true;
   }
   if (contains_decay != contains_fwhm) {
      LOKI_ABORT("Both phase_decay_time_steps and fwhm must be specified for noisy drivers.");
   }
   // They want a noisy driver so compute the data members necessary to figure
   // out the phases.
   if (contains_decay) {
      a_pp.query("phase_decay_time_steps",
         m_parameters[PHASE_DECAY_TIME_STEPS]);
      a_pp.query("fwhm", m_parameters[FWHM]);
      m_factor_1 = 1.0 - 1.0/m_parameters[PHASE_DECAY_TIME_STEPS];
      m_factor_2 =
         sqrt(3.0/(m_parameters[PHASE_DECAY_TIME_STEPS]*log(2.0)))*m_parameters[FWHM];
      m_do_random_phase = true;
   }
}


void
ShapedRampedCosineDriver::computePhase(
   double a_dt,
   bool a_first_rk_stage) const
{
   // We only compute a new random phase for the 1st RK stage.  This phase
   // applies throughout the time step.
   // On initialization a_dt will be < 0.0.  In this case, if the simulation
   // has been restarted we want the last phase which was read from the restart
   // file, otherwise we want the initial phase which is 0.0.
   if (m_do_random_phase && a_first_rk_stage && a_dt > 0.0) {
      if (m_num_phase_evals > 0) {
         m_phase_h = m_factor_1*m_phase_h + m_factor_2*a_dt*(rand() - 0.5);
         m_phase += m_phase_h;
      }
      ++m_num_phase_evals;
   }
}

} // end namespace Loki
