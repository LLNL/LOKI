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
#include "ShapedRampedCosineCurrentDriver.H"
#include "ShapedRampedCosineCurrentDriverF.H"

namespace Loki {

const aString ShapedRampedCosineCurrentDriver::s_CLASS_NAME(
   "Shaped Ramped Cosine Current Driver");


bool
ShapedRampedCosineCurrentDriver::isType(
   const aString& a_name)
{
   if (a_name.matches(s_CLASS_NAME)) {
      return true;
   }
   return false;
}


ShapedRampedCosineCurrentDriver::ShapedRampedCosineCurrentDriver(
   ParmParse& a_pp)
{
   // Size the parameters and set some hopefully sane defaults then get what
   // the user actually wants from the input file.
   m_parameters.resize(NUM_PARAMS);
   m_parameters(SHAPE)  = 1.0;
   m_parameters(OMEGA)  = 1.0;
   m_parameters(J0)     = 0.01;
   m_parameters(T0)     = 0.0;
   m_parameters(TRAMP)  = 10.0;
   m_parameters(TOFF)   = 10.0;
   m_parameters(X0)     = 0.0;
   parseParameters(a_pp);
}


ShapedRampedCosineCurrentDriver::~ShapedRampedCosineCurrentDriver()
{
}


void
ShapedRampedCosineCurrentDriver::evaluate(
   RealArray& a_Jx,
   RealArray& a_Jy,
   RealArray& a_Jz,
   const tbox::Box& a_local_box,
   const ProblemDomain& a_domain,
   real a_time) const
{
   // If the driver is on at a_time, compute the current densities by delegating
   // the work to fortran.
   if (m_parameters(T0) <= a_time &&
       a_time < m_parameters(T0) + m_parameters(TRAMP) + m_parameters(TOFF)) {
      evaluateShapedRampedCurrentDriver(*a_Jx.getDataPointer(),
         *a_Jy.getDataPointer(),
         *a_Jz.getDataPointer(),
         BOX2D_TO_FORT(a_local_box),
         BOX2D_TO_FORT(a_domain.box()),
         PROBLEMDOMAIN_TO_FORT(a_domain),
         a_time,
         *m_parameters.getDataPointer());
   }
}


void
ShapedRampedCosineCurrentDriver::printParameters() const
{
   // Write the driver parameters.
   printF("  Using external current forcing:\n" );
   printF("    width           = %e\n", m_parameters(WIDTH));
   printF("    apply direction = %i\n", m_parameters(APPLY_DIR));
   printF("    shape           = %e\n", m_parameters(SHAPE));
   printF("    omega           = %e\n", m_parameters(OMEGA));
   printF("    J0              = %e\n", m_parameters(J0));
   printF("    t0              = %e\n", m_parameters(T0));
   printF("    t_ramp          = %e\n", m_parameters(TRAMP));
   printF("    t_off           = %e\n", m_parameters(TOFF));
   printF("    x0              = %e\n", m_parameters(X0));
}


void
ShapedRampedCosineCurrentDriver::parseParameters(
   ParmParse& a_pp)
{
   // The width and apply_dir are required.
   if (a_pp.contains("width")) {
      a_pp.get("width", m_parameters(WIDTH));
   }
   else {
      OV_ABORT("Must supply width!");
   }
   if (a_pp.contains("apply_dir")) {
      aString tmp;
      a_pp.get("apply_dir", tmp);
      if (tmp.matches("x") || tmp.matches("X")) {
         m_parameters(APPLY_DIR) = 1.0;
      }
      else if (tmp.matches("y") || tmp.matches("Y")) {
         m_parameters(APPLY_DIR) = 2.0;
      }
      else if (tmp.matches("z") || tmp.matches("Z")) {
         m_parameters(APPLY_DIR) = 3.0;
      }
      else {
         OV_ABORT("Unknown apply_dir!");
      }
   }
   else {
      OV_ABORT("Must supply apply_dir!");
   }

   // Everything else is not required which may be a bad idea as it's not clear
   // what valid defaults should be.
   a_pp.query("shape",   m_parameters(SHAPE));
   a_pp.query("omega",   m_parameters(OMEGA));
   a_pp.query("J_0",     m_parameters(J0));
   a_pp.query("t0",      m_parameters(T0));
   a_pp.query("t_ramp",  m_parameters(TRAMP));
   a_pp.query("t_off",   m_parameters(TOFF));
   a_pp.query("x0",      m_parameters(X0));
}

} // end namespace Loki
