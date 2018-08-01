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
#include "SimpleEMIC.H"
#include "SimpleEMICF.H"
#include "BoxOps.H"

namespace Loki {

const aString SimpleEMIC::s_CLASS_NAME("SimpleEMIC");

bool
SimpleEMIC::isType(
   const aString& a_name)
{
   if (a_name.matches(s_CLASS_NAME)) {
      return true;
   }
   return false;
}

SimpleEMIC::SimpleEMIC(
   ParmParse& a_pp)
{
   // Size the double and integer parameters.
   m_dparameters.resize(NUM_DPARAMS);
   m_iparameters.resize(NUM_IPARAMS);

   // Read which field we're initializing.
   aString tmp;
   if (!a_pp.query("field", tmp)) {
      OV_ABORT("Must supply field");
   }
   else if (tmp.matches("E") || tmp.matches("e")) {
      m_iparameters(FIELD) = E;
   }
   else if (tmp.matches("B") || tmp.matches("b")) {
      m_iparameters(FIELD) = B;
   }
   else {
      OV_ABORT("Unknown field");
   }

   // All this input is required as there's no obvious defaults.
   if (!a_pp.query("xamp", m_dparameters(XAMP))) {
      OV_ABORT("Must supply xamp");
   }
   if (!a_pp.query("yamp", m_dparameters(YAMP))) {
      OV_ABORT("Must supply yamp");
   }
   if (!a_pp.query("zamp", m_dparameters(ZAMP))) {
      OV_ABORT("Must supply zamp");
   }

   if (!a_pp.query("x_wave_number", m_dparameters(X_WAVE_NUMBER))) {
      OV_ABORT("Must supply x_wave_number");
   }

   if (!a_pp.query("y_wave_number", m_dparameters(Y_WAVE_NUMBER))) {
      OV_ABORT("Must supply y_wave_number");
   }

   if (!a_pp.query("phase", m_dparameters(PHI))) {
      OV_ABORT("Must supply phase");
   }
}


SimpleEMIC::~SimpleEMIC()
{
}


void
SimpleEMIC::set(
   RealArray& a_u,
   const ProblemDomain& a_domain) const
{
   // Delegate evaluation to fortran.
   tbox::Box local_box(BoxOps::getLocalBox(a_u));

   SET_SIMPLE_EMIC(BOX2D_TO_FORT(local_box),
      *a_u.getDataPointer(),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      *m_dparameters.getDataPointer(),
      *m_iparameters.getDataPointer());
}


bool
SimpleEMIC::initializesE()
{
   return m_iparameters(FIELD) == E;
}


bool
SimpleEMIC::initializesB()
{
   return m_iparameters(FIELD) == B;
}


void
SimpleEMIC::printParameters() const
{
   // Print all input parameters.
   printF("  Using EM field initialization:\n");
   printF("    field         = %s\n",
      m_iparameters(FIELD) == E ? "E" : "B");
   printF("    X amplitude   = %f\n", m_dparameters(XAMP));
   printF("    Y amplitude   = %f\n", m_dparameters(YAMP));
   printF("    Z amplitude   = %f\n", m_dparameters(ZAMP));
   printF("    x wave number = %f\n", m_dparameters(X_WAVE_NUMBER));
   printF("    y wave number = %f\n", m_dparameters(Y_WAVE_NUMBER));
   printF("    phase         = %f\n", m_dparameters(PHI));
}

} // end namespace Loki
