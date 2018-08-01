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
#include "ShapedRampedCosinePotentialDriver.H"
#include "ShapedRampedCosinePotentialDriverF.H"

namespace Loki {

const aString ShapedRampedCosinePotentialDriver::s_CLASS_NAME(
   "Shaped Ramped Cosine Potential Driver");


bool
ShapedRampedCosinePotentialDriver::isType(
   const aString& a_name)
{
   if (a_name.matches(s_CLASS_NAME)) {
      return true;
   }
   return false;
}


ShapedRampedCosinePotentialDriver::ShapedRampedCosinePotentialDriver(
   ParmParse& a_pp)
{
   // Size the parameters and set some defaults which most likely don't
   // make any sense.
   m_parameters.resize(NUM_PARAMS);
   m_parameters(XWIDTH) = 0.5;
   m_parameters(YWIDTH) = 0.5;
   m_parameters(SHAPE)  = 1.0;
   m_parameters(OMEGA)  = 1.0;
   m_parameters(AMP)    = 0.01;
   m_parameters(T0)     = 0.0;
   m_parameters(TRAMP)  = 10.0;
   m_parameters(TOFF)   = 10.0;
   m_parameters(XSHAPE) = 0.0;
   m_parameters(LWIDTH) = 0.5;
   m_parameters(X0)     = 0.0;

   // Now get what the user really wants.
   parseParameters(a_pp);
}


ShapedRampedCosinePotentialDriver::~ShapedRampedCosinePotentialDriver()
{
}


void
ShapedRampedCosinePotentialDriver::evaluate(
   RealArray& a_phi,
   const tbox::Box& a_fill_box,
   const ProblemDomain& a_domain,
   real a_time) const
{
   // Delegate evaluation to fortran.
   evaluateShapedRampedPotentialDriver(*a_phi.getDataPointer(),
      BOX2D_TO_FORT(a_fill_box),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      a_time,
      *m_parameters.getDataPointer());
}


void
ShapedRampedCosinePotentialDriver::printParameters() const
{
   // Write all the input parameters.
   printF("  Using external electric potential forcing:\n");
   printF("    xwidth  = %e\n", m_parameters(XWIDTH));
   printF("    ywidth  = %e\n", m_parameters(YWIDTH));
   printF("    shape   = %e\n", m_parameters(SHAPE));
   printF("    omega   = %e\n", m_parameters(OMEGA));
   printF("    amp     = %e\n", m_parameters(AMP));
   printF("    t0      = %e\n", m_parameters(T0));
   printF("    t_ramp  = %e\n", m_parameters(TRAMP));
   printF("    t_off   = %e\n", m_parameters(TOFF));
   printF("    x_shape = %e\n", m_parameters(XSHAPE));
   printF("    lwidth  = %e\n", m_parameters(LWIDTH));
   printF("    x0      = %e\n", m_parameters(X0));
}


void
ShapedRampedCosinePotentialDriver::parseParameters(
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
   a_pp.query("amp",     m_parameters(AMP));
   a_pp.query("t0",      m_parameters(T0));
   a_pp.query("t_ramp",  m_parameters(TRAMP));
   a_pp.query("t_off",   m_parameters(TOFF));
   a_pp.query("x_shape", m_parameters(XSHAPE));
   a_pp.query("x0",      m_parameters(X0));
}

} // end namespace Loki
