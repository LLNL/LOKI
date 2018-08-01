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
#include "PitchAngleCollisionOperator.H"
#include "PitchAngleCollisionOperatorF.H"

namespace Loki {

const aString PitchAngleCollisionOperator::s_CLASS_NAME(
   "Pitch Angle Collision Operator");


PitchAngleCollisionOperator::PitchAngleCollisionOperator(
   ParmParse& a_pp,
   int a_solution_order)
{
   // Size the double valued parameters.
   m_dparameters.resize(NUM_DPARAMS);

   // Size the integer valued parameters and set a default for the conservative
   // operator option.
   m_iparameters.resize(NUM_IPARAMS);
   m_iparameters(ICONS) = 1;
   m_iparameters(SOLN_ORDER) = a_solution_order;

   // Get what the user really wants.
   parseParameters(a_pp);

   if (m_iparameters(ICONS) != 1 && m_iparameters(SOLN_ORDER) == 6) {
      OV_ABORT("Non-conservative operator in 6th order not supported.");
   }
}


PitchAngleCollisionOperator::~PitchAngleCollisionOperator()
{
}


void
PitchAngleCollisionOperator::evaluate(
   RealArray& a_rhs,
   const RealArray& a_u,
   const tbox::Box& a_local_box,
   const tbox::Box& a_interior_box,
   const ProblemDomain& a_domain,
   real a_vflowx,
   real a_vflowy)
{
   // Set the flow velocities in the parameters.
   m_dparameters(VFLOWX) = a_vflowx;
   m_dparameters(VFLOWY) = a_vflowy;

   // Compute the thermal velocity based on the algorithm specified by the user.
   if (m_vthermal_method == LOCAL_VTHERMAL) {
      real vth = computeVthLocal(a_u,
         a_local_box,
         a_interior_box,
         a_domain,
         a_vflowx,
         a_vflowy);
      m_dparameters(VTHERMAL) = vth;
   }
   else if (m_vthermal_method == GLOBAL_VTHERMAL) {
      real vth = computeVthGlobal(a_u,
         a_local_box,
         a_interior_box,
         a_domain,
         a_vflowx,
         a_vflowy);
      m_dparameters(VTHERMAL) = vth;
   }

   // Delegate the operator to fortran.
   FORT_APPEND_PITCH_ANGLE_COLLISION(*a_rhs.getDataPointer(),
      *a_u.getDataPointer(),
      BOX4D_TO_FORT(a_local_box),
      BOX4D_TO_FORT(a_interior_box),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      *m_dparameters.getDataPointer(),
      *m_iparameters.getDataPointer());
}

real
PitchAngleCollisionOperator::computeRealLam(
   real a_dv) const
{
   real pi = 4.0*atan(1.0);
   return m_dparameters(NU)*pow(m_dparameters(VTHERMAL_DT), 3.0)*pi*pi/(a_dv*a_dv*std::max(a_dv, m_dparameters(VFLOOR)));
}


void
PitchAngleCollisionOperator::printParameters() const
{
   // Print the user selected input.
   printF("  Using pitch angle collision:\n" );
   printF("    vceil        = %e\n", m_dparameters(VCEIL));
   printF("    vfloor       = %e\n", m_dparameters(VFLOOR));
   if (m_vthermal_method == INPUT_VTHERMAL) {
      printF("    vthermal     = %e\n", m_dparameters(VTHERMAL));
   }
   else if (m_vthermal_method == LOCAL_VTHERMAL) {
      printF("    computing vthermal locally\n");
   }
   else if (m_vthermal_method == GLOBAL_VTHERMAL) {
      printF("    computing vthermal globally\n");
   }
   printF("    vthermal_dt  = %e\n", m_dparameters(VTHERMAL_DT));
   printF("    nuCoeff      = %e\n", m_dparameters(NU));
   printF("    conservative = %d\n", m_iparameters(ICONS));
}


bool
PitchAngleCollisionOperator::isType(
   const aString& a_name)
{
   if (a_name.matches(s_CLASS_NAME)) {
      return true;
   }
   return false;
}


void
PitchAngleCollisionOperator::parseParameters(
   ParmParse& a_pp)
{
   // vecil and vfloor are required
   if (!a_pp.query("collision_vceil", m_dparameters(VCEIL))) {
      OV_ABORT("Must supply collision_vceil");
   }
   if (!a_pp.query("collision_vfloor", m_dparameters(VFLOOR))) {
      OV_ABORT("Must supply collision_vfloor");
   }

   // This is kind of messy.  There are 2 basic ways that vthermal can be
   // "computed":
   // 1) Through one of 3 calculation methods:
   //    a) The "local vthermal" calculation
   //       In this case the user should not specify vthermal as that is
   //       contradictory
   //    b) The "global vthermal" calculation
   //       In this case the user should not specify vthermal as that is
   //       contradictory
   //    c) The "input vthermal"
   //       In this case the user must specify vthermal otherwise we know
   //       nothing.
   // 2) The user just gives us vthermal which is a shortcut the 1c above.
   bool has_vthermal_method = a_pp.contains("collision_vthermal_method");
   bool has_vthermal = a_pp.contains("collision_vthermal");
   if (!has_vthermal_method) {
      if (!has_vthermal) {
         OV_ABORT("If inputting collision_vthermal then it must be supplied");
         
      }
      else {
         a_pp.query("collision_vthermal", m_dparameters(VTHERMAL));
      }
      m_vthermal_method = INPUT_VTHERMAL;
   }
   else {
      aString vthermal_method("input vthermal");
      a_pp.query("collision_vthermal_method", vthermal_method);
      if (vthermal_method.matches("input vthermal")) {
         if (!has_vthermal) {
            OV_ABORT("If inputting collision_vthermal then it must be supplied");
         }
         else {
            a_pp.query("collision_vthermal", m_dparameters(VTHERMAL));
         }
         m_vthermal_method = INPUT_VTHERMAL;
      }
      else if (vthermal_method.matches("local vthermal")) {
         if (has_vthermal) {
            OV_ABORT("Either compute vthermal or specify vthermal but not both");
         }
         m_vthermal_method = LOCAL_VTHERMAL;
      }
      else if (vthermal_method.matches("global vthermal")) {
         if (has_vthermal) {
            OV_ABORT("Either compute vthermal or specify vthermal but not both");
         }
         m_vthermal_method = GLOBAL_VTHERMAL;
      }
      else {
         OV_ABORT("Unknown input for collision_vthermal_method");
      }
   }

   // These inputs are required.
   if (!a_pp.query("collision_vthermal_dt", m_dparameters(VTHERMAL_DT))) {
      OV_ABORT("Must supply collision_vthermal_dt");
   }
   if (!a_pp.query("collision_nuCoeff", m_dparameters(NU))) {
      OV_ABORT("Must supply collision_nuCoeff");
   }

   // The default is the conservative operator.
   a_pp.query("collision_conservative", m_iparameters(ICONS));
}

} // end namespace Loki
