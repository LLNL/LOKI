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
#include "TrigTZSource.H"
#include "TZSourceF.H"
#include "BoxOps.H"

namespace Loki {

const aString TrigTZSource::s_CLASS_NAME("TrigTZSource");

bool
TrigTZSource::isType(
   const aString& a_name)
{
   if (a_name.matches(s_CLASS_NAME)) {
      return true;
   }
   return false;
}

TrigTZSource::TrigTZSource(
   ParmParse& a_pp)
{
   // Size the parameters array and get the single bit of user input for this
   // source which is required.
   m_dparameters.resize(NUM_DPARAMS);

   if (!a_pp.query("amp", m_dparameters(AMP))) {
      OV_ABORT("Must supply amp");
   }
}


TrigTZSource::~TrigTZSource()
{
}


void
TrigTZSource::set(
   RealArray& a_dist_func,
   const ProblemDomain& a_domain,
   real a_time) const
{
   // Delegate evaluation of the source to fortran.
   tbox::Box local_box(BoxOps::getLocalBox(a_dist_func));

   FORT_SET_TRIG_TZ_SOURCE(*a_dist_func.getDataPointer(),
      BOX4D_TO_FORT(local_box),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      a_time,
      *m_dparameters.getDataPointer());
}


void
TrigTZSource::computeError(
   RealArray& a_tz_error,
   const RealArray& a_dist_func,
   const ProblemDomain& a_domain,
   real a_time) const
{
   // Delegate evaluation of the error to fortran.
   tbox::Box local_box(BoxOps::getLocalBox(a_dist_func));

   FORT_COMPUTE_TRIG_TZ_SOURCE_ERROR(*a_tz_error.getDataPointer(),
      *a_dist_func.getDataPointer(),
      BOX4D_TO_FORT(local_box),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      a_time,
      *m_dparameters.getDataPointer());
}


void
TrigTZSource::printParameters() const
{
   // Print the only parameter.
   printF("  Using twilight zone source:\n");
   printF("    amplitude     = %f\n", m_dparameters(AMP));
}

} // end namespace Loki
