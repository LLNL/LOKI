/*************************************************************************
 *
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * Written by Jeffrey Banks banksj3@rpi.edu (Rensselaer Polytechnic Institute,
 * Amos Eaton 301, 110 8th St., Troy, NY 12180); Jeffrey Hittinger
 * hittinger1@llnl.gov, William Arrighi arrighi2@llnl.gov, Richard Berger
 * berger5@llnl.gov, Thomas Chapman chapman29@llnl.gov (LLNL, P.O Box 808,
 * Livermore, CA 94551); Stephan Brunner stephan.brunner@epfl.ch (Ecole
 * Polytechnique Federale de Lausanne, EPFL SB SPC-TH, PPB 312, Station 13,
 * CH-1015 Lausanne, Switzerland).
 * CODE-744849
 *
 * All rights reserved.
 *
 * This file is part of Loki.  For details, see.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 ************************************************************************/
#include "ShapedRampedCosineDriver.H"
#include "ShapedRampedCosineDriverF.H"
#include "RestartManager.H"
#include "Loki_Utilities.H"

#include <math.h>
#include <sstream>

namespace Loki {

const aString ShapedRampedCosineDriver::s_CLASS_NAME(
   "Shaped Ramped Cosine Driver");


ShapedRampedCosineDriver::ShapedRampedCosineDriver(
   ParmParse& a_pp,
   int a_driver_num)
   : m_do_random_phase(false),
     m_phase(0.0),
     m_phase_h(0.0),
     m_num_phase_evals(0),
     m_next(1)
{
   // Construct a unique name for this driver which is the name of the
   // sub-database in the restart file containing its state.
   std::stringstream ss;
   ss << a_driver_num;
   aString tmp(ss.str());
   m_name = s_CLASS_NAME + tmp;

   // Size the parameters and set some hopefully sane defaults.
   m_parameters.resize(NUM_PARAMS);
   m_parameters(XWIDTH)                 = 0.5;
   m_parameters(YWIDTH)                 = 0.5;
   m_parameters(SHAPE)                  = 1.0;
   m_parameters(OMEGA)                  = 1.0;
   m_parameters(E0)                     = 0.01;
   m_parameters(T0)                     = 0.0;
   m_parameters(TRAMP)                  = 10.0;
   m_parameters(TOFF)                   = 10.0;
   m_parameters(XSHAPE)                 = 0.0;
   m_parameters(LWIDTH)                 = 0.5;
   m_parameters(X0)                     = 0.0;
   m_parameters(ALPHA)                  = 0.0;
   m_parameters(TRES)                   = 0.0;
   m_parameters(PHASE_DECAY_TIME_STEPS) = 0.0;
   m_parameters(FWHM)                   = 0.0;

   // Get what the user really wants.
   parseParameters(a_pp);

   // This object saves/writes state to restart so register it which the
   // RestartManager.
   RestartManager::getManager()->registerRestart(this);
}


ShapedRampedCosineDriver::~ShapedRampedCosineDriver()
{
}


void
ShapedRampedCosineDriver::evaluate(
   RealArray& a_Ex,
   RealArray& a_Ey,
   const tbox::Box& a_fill_box,
   const ProblemDomain& a_domain,
   real a_time,
   real a_dt,
   int a_stage) const
{
   // If the driver is on at a_time then delegate its evaluation to fortran.
   if (m_parameters(T0) <= a_time &&
       a_time < m_parameters(T0) + m_parameters(TRAMP) + m_parameters(TOFF)) {
      computePhase(a_dt, a_stage);
      evaluateShapedRampedDriver(*a_Ex.getDataPointer(),
         *a_Ey.getDataPointer(),
         BOX2D_TO_FORT(a_fill_box),
         PROBLEMDOMAIN_TO_FORT(a_domain),
         a_time,
         *m_parameters.getDataPointer(),
         m_phase);
   }
}


void
ShapedRampedCosineDriver::printParameters() const
{
   // Write all the driver parameters.
   printF("  Using external electric field forcing:\n" );
   printF("    xwidth                 = %e\n", m_parameters(XWIDTH));
   printF("    ywidth                 = %e\n", m_parameters(YWIDTH));
   printF("    shape                  = %e\n", m_parameters(SHAPE));
   printF("    omega                  = %e\n", m_parameters(OMEGA));
   printF("    E0                     = %e\n", m_parameters(E0));
   printF("    t0                     = %e\n", m_parameters(T0));
   printF("    t_ramp                 = %e\n", m_parameters(TRAMP));
   printF("    t_off                  = %e\n", m_parameters(TOFF));
   printF("    x_shape                = %e\n", m_parameters(XSHAPE));
   printF("    lwidth                 = %e\n", m_parameters(LWIDTH));
   printF("    x0                     = %e\n", m_parameters(X0));
   printF("    alpha                  = %e\n", m_parameters(ALPHA));
   printF("    t_res                  = %e\n", m_parameters(TRES));
   printF("    phase decay time steps = %e\n", m_parameters(PHASE_DECAY_TIME_STEPS));
   printF("    fwhm                   = %e\n", m_parameters(FWHM));
}


void
ShapedRampedCosineDriver::getFromRestart(
   const HDF_DataBase& a_db)
{
   // Get the subdatabase containing this object's state.
   HDF_DataBase sub_db;
   a_db.locate(sub_db, m_name);

   // If this is a noisy driver, get the number of times computePhase was
   // called from the database.
   int times_called = 0;
   if (sub_db.get(times_called, "num_phase_evals") == 0) {

      // Now get the phase and half time phase.
      sub_db.get(m_phase, "phase");
      sub_db.get(m_phase_h, "phase_h");

      // Now run the random number generator times_called times so that it is
      // in the same state as when the restart file was generated.  Once this
      // is done the computation of the random phase will proceed as if there
      // had been no restart.
      for (int i = 0; i < times_called; ++i) {
         rand();
      }
      m_num_phase_evals = times_called;
   }
}


void
ShapedRampedCosineDriver::putToRestart(
   HDF_DataBase& a_db,
   real a_time)
{
   NULL_USE(a_time);

   // Make the subdatabase to contain this object's state.
   HDF_DataBase sub_db;
   a_db.create(sub_db, m_name, "directory");

   // Put m_num_phase_evals, m_phase, and, m_phase_h into the database.
   sub_db.put(m_num_phase_evals, "num_phase_evals");
   sub_db.put(m_phase, "phase");
   sub_db.put(m_phase_h, "phase_h");
}


bool
ShapedRampedCosineDriver::isType(
   const aString& a_name)
{
   if (a_name.matches(s_CLASS_NAME)) {
      return true;
   }
   return false;
}


void
ShapedRampedCosineDriver::parseParameters(
   ParmParse& a_pp)
{
   // There's a bunch of backward compatibility here for older style syntax.
   // The parameters are all in terms of widths in configuration space
   // directions.  New syntax enters this directly.  Old syntax was in terms
   // of wave numbers and lengths which are converted to widths.

   // Read the information about x modulation.  None of this is required but it
   // probably should be as defaults are likely meaningless.
   if (a_pp.contains("kx") || a_pp.contains("Lx")) {
      if (a_pp.contains("xwidth")) {
         OV_ABORT("Mixed old and new syntax for kx/Lx");
      }
      else {
         real kx = 1.0;
         real lx = 1.0;
         a_pp.query("kx", kx);
         a_pp.query("Lx", lx);
         m_parameters(XWIDTH) = lx/(2.0*kx);
      }
   }
   else {
      a_pp.query("xwidth", m_parameters(XWIDTH));
   }

   // Read the information about y modulation.  None of this is required but it
   // probably should be as defaults are likely meaningless.
   if (a_pp.contains("ky") || a_pp.contains("Ly")) {
      if (a_pp.contains("ywidth")) {
         OV_ABORT("Mixed old and new syntax for ky/Ly");
      }
      else {
         real ky = 1.0;
         real ly = 1.0;
         a_pp.query("ky", ky);
         a_pp.query("Ly", ly);
         m_parameters(YWIDTH) = ly/(2.0*ky);
      }
   }
   else {
      a_pp.query("ywidth", m_parameters(YWIDTH));
   }

   // Read the information about x shape envelope.  None of this is required
   // but it probably should be as defaults are likely meaningless.
   if (a_pp.contains("kl") || a_pp.contains("Ll")) {
      if (a_pp.contains("lwidth")) {
         OV_ABORT("Mixed old and new syntax for kl/Ll");
      }
      else {
         real kl = 1.0;
         real ll = 1.0;
         a_pp.query("kl", kl);
         a_pp.query("Ll", ll);
         m_parameters(LWIDTH) = ll/(2.0*kl);
      }
   }
   else {
      a_pp.query("lwidth", m_parameters(LWIDTH));
   }

   // These are not required but should be.
   a_pp.query("shape",   m_parameters(SHAPE));
   a_pp.query("omega",   m_parameters(OMEGA));
   a_pp.query("E_0",     m_parameters(E0));
   a_pp.query("t0",      m_parameters(T0));
   a_pp.query("t_ramp",  m_parameters(TRAMP));
   a_pp.query("t_off",   m_parameters(TOFF));
   a_pp.query("x_shape", m_parameters(XSHAPE));
   a_pp.query("x0",      m_parameters(X0));
   a_pp.query("alpha",   m_parameters(ALPHA));
   a_pp.query("t_res",   m_parameters(TRES));

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
      OV_ABORT("Both phase_decay_time_steps and fwhm must be specified for noisy drivers.");
   }
   // They want a noisy driver so compute the data members necessary to figure
   // out the phases.
   if (contains_decay) {
      a_pp.query("phase_decay_time_steps",
         m_parameters(PHASE_DECAY_TIME_STEPS));
      a_pp.query("fwhm", m_parameters(FWHM));
      m_factor_1 = 1.0 - 1.0/m_parameters(PHASE_DECAY_TIME_STEPS);
      m_factor_2 =
         sqrt(3.0/(m_parameters(PHASE_DECAY_TIME_STEPS)*log(2.0)))*m_parameters(FWHM);
      m_do_random_phase = true;
   }
}


void
ShapedRampedCosineDriver::computePhase(
   real a_dt,
   int a_stage) const
{
   // We only compute a new random phase for the 1st RK stage.  This phase
   // applies throughout the time step.
   // On initialization a_dt will be < 0.0.  In this case, if the simulation
   // has been restarted we want the last phase which was read from the restart
   // file, otherwise we want the initial phase which is 0.0.
   if (m_do_random_phase && a_stage == 1 && a_dt > 0.0) {
      if (m_num_phase_evals > 0) {
         m_phase_h = m_factor_1*m_phase_h + m_factor_2*a_dt*(rand() - 0.5);
         m_phase += m_phase_h;
      }
      ++m_num_phase_evals;
   }
}

} // end namespace Loki
