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
#include "InterspeciesCollisionOperator.H"
#include "InterspeciesCollisionOperatorF.H"
#include "IntegralOp.H"

namespace Loki {

const aString InterspeciesCollisionOperator::s_CLASS_NAME(
   "Interspecies Collision Operator");


InterspeciesCollisionOperator::InterspeciesCollisionOperator(
   ParmParse& a_pp,
   int a_solution_order)
   : m_back_reaction(false),
     m_compute_moments(true),
     m_compute_constants(true)
{
   // Size the double and integer operator parameters.  Set what we know and
   // read from user input what we do not know.
   m_dparameters.resize(NUM_DPARAMS);
   m_iparameters.resize(NUM_IPARAMS);
   m_iparameters(SOLN_ORDER) = a_solution_order;
   m_iparameters(NUM_GHOSTS) = (a_solution_order == 4) ? 2 : 3;
   parseParameters(a_pp);
}


InterspeciesCollisionOperator::~InterspeciesCollisionOperator()
{
}


void
InterspeciesCollisionOperator::evaluate(
   RealArray& a_rhs,
   const RealArray& a_u,
   const tbox::Box& a_local_box,
   const tbox::Box& a_interior_box,
   const ProblemDomain& a_domain,
   real a_vflowx,
   real a_vflowy)
{
   // Compute vthermal and alpha if the user specified a computation method.
   m_dparameters(VFLOWX) = a_vflowx;
   m_dparameters(VFLOWY) = a_vflowy;
   if (m_vthermal_method == LOCAL_VTHERMAL) {
      real vth = computeVthLocal(a_u,
         a_local_box,
         a_interior_box,
         a_domain,
         a_vflowx,
         a_vflowy);
      m_dparameters(VTHERMAL) = vth;
      m_dparameters(ALPHA) = 1.0/(vth*vth);
   }
   else if (m_vthermal_method == GLOBAL_VTHERMAL) {
      real vth = computeVthGlobal(a_u,
         a_local_box,
         a_interior_box,
         a_domain,
         a_vflowx,
         a_vflowy);
      m_dparameters(VTHERMAL) = vth;
      m_dparameters(ALPHA) = 1.0/(vth*vth);
   }

   // We need to compute the constants once.  At the time of construction,
   // the interior box is not known so we do it at the first evaluation of
   // the operator.
   if (m_compute_constants) {
      computeConstants(a_interior_box, a_domain);
      m_compute_constants = false;
   }

   // If the back reaction is computed then there's a bunch of other things that
   // need to be computed.
   if (m_back_reaction) {
      // Again, we need to compute the various moments once but at the time of
      // construction we don't know the local or interior boxes so we do it at
      // the first evaluation of the operator.
      if (m_compute_moments) {
         constructBackReactionArrays(a_local_box);
         constructBackReactionReductionSchedules(a_interior_box, a_domain);
         computeBackReactionMoments(a_local_box, a_interior_box, a_domain);
         m_compute_moments = false;
      }

      // Compute back reaction terms.
      m_rMomx = 0.0;
      m_rMomy = 0.0;
      m_rKE = 0.0;
      FORT_COMPUTE_INTERSPECIES_NUMERATORS(*m_rMomx.getDataPointer(),
         *m_rMomy.getDataPointer(),
         *m_rKE.getDataPointer(),
         *m_cMomx.getDataPointer(),
         *m_cMomy.getDataPointer(),
         *m_cKE.getDataPointer(),
         *m_fM.getDataPointer(),
         *a_u.getDataPointer(),
         BOX4D_TO_FORT(a_local_box),
         BOX4D_TO_FORT(a_interior_box));
      m_IMomxN = 0.0;
      m_IMomyN = 0.0;
      m_IKEN = 0.0;
      m_x_mom_reduction->execute(m_IMomxN);
      m_y_mom_reduction->execute(m_IMomyN);
      m_ke_reduction->execute(m_IKEN);
   }

   // Compute collision operator
   FORT_APPEND_INTERSPECIES_COLLISION(*a_rhs.getDataPointer(),
      *a_u.getDataPointer(),
      *m_vx.getDataPointer(),
      *m_vy.getDataPointer(),
      *m_khat.getDataPointer(),
      *m_hhat.getDataPointer(),
      BOX4D_TO_FORT(a_local_box),
      BOX4D_TO_FORT(a_interior_box),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      *m_dparameters.getDataPointer(),
      *m_iparameters.getDataPointer());

   if (m_back_reaction) {
      // Fold back reaction terms into collision operator.
      FORT_APPEND_INTERSPECIES_BACK_REACTION(*a_rhs.getDataPointer(),
         *m_cMomx.getDataPointer(),
         *m_cMomy.getDataPointer(),
         *m_cKE.getDataPointer(),
         *m_IMomxN.getDataPointer(),
         *m_IMomyN.getDataPointer(),
         *m_IKEN.getDataPointer(),
         *m_IMomxD.getDataPointer(),
         *m_IMomyD.getDataPointer(),
         *m_IKED.getDataPointer(),
         BOX4D_TO_FORT(a_local_box),
         BOX4D_TO_FORT(a_interior_box),
         *m_dparameters.getDataPointer());
   }
}

real
InterspeciesCollisionOperator::computeRealLam(
   real a_dv) const
{
   real pi = 4.0*atan(1.0);
   return m_dparameters(NU)*pow(m_dparameters(VTHERMAL_DT), 3.0)*pi*pi/(a_dv*a_dv*std::max(a_dv, m_dparameters(VFLOOR)));
}


void
InterspeciesCollisionOperator::printParameters() const
{
   // Print all operator parameters in play.
   printF("  Using Interspecies collision:\n" );
   printF("    vceil       = %e\n", m_dparameters(VCEIL));
   printF("    vfloor      = %e\n", m_dparameters(VFLOOR));
   if (m_vthermal_method == INPUT_VTHERMAL) {
      printF("    vthermal     = %e\n", m_dparameters(VTHERMAL));
   }
   else if (m_vthermal_method == LOCAL_VTHERMAL) {
      printF("    computing vthermal locally\n");
   }
   else if (m_vthermal_method == GLOBAL_VTHERMAL) {
      printF("    computing vthermal globally\n");
   }
   printF("    vthermal_dt = %e\n", m_dparameters(VTHERMAL_DT));
   printF("    nuCoeff     = %e\n", m_dparameters(NU));
   if (m_vthermal_method == INPUT_VTHERMAL) {
     printF("    alpha       = %e\n", m_dparameters(ALPHA));
   }
   else if (m_vthermal_method == LOCAL_VTHERMAL) {
      printF("    computing alpha locally\n");
   }
   else if (m_vthermal_method == GLOBAL_VTHERMAL) {
      printF("    computing alpha globally\n");
   }
   if (m_back_reaction) {
      printF("    massR       = %e\n", m_dparameters(MASSR));
      printF("    including back reaction\n");
   }
   else {
      printF("    not including back reaction\n");
   }
}


bool
InterspeciesCollisionOperator::isType(
   const aString& a_name)
{
   if (a_name.matches(s_CLASS_NAME)) {
      return true;
   }
   return false;
}


void
InterspeciesCollisionOperator::parseParameters(
   ParmParse& a_pp)
{
   // The collision_vceil and collision_vfloor are required.
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
   // There is also the alpha parameter which should be 1/vth**2.  There's
   // a bunch of logic here to check for user double talk about what they
   // want to do WRT alpha.  Truthfully, the ALPHA parameter should be
   // eliminated as it is really a derived quantity.
   bool has_vthermal_method = a_pp.contains("collision_vthermal_method");
   bool has_vthermal = a_pp.contains("collision_vthermal");
   bool has_alpha = a_pp.contains("collision_alpha");
   if (!has_vthermal_method) {
      if (!has_vthermal) {
         OV_ABORT("If inputting collision_vthermal then it must be supplied");
         
      }
      else if (!has_alpha) {
         OV_ABORT("If inputting collision_vthermal then collision_alpha must be supplied");
         
      }
      else {
         a_pp.query("collision_vthermal", m_dparameters(VTHERMAL));
         a_pp.query("collision_alpha", m_dparameters(ALPHA));
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
         else if (!has_alpha) {
            OV_ABORT("If inputting collision_vthermal then collision_alpha must be supplied");
         
         }
         else {
            a_pp.query("collision_vthermal", m_dparameters(VTHERMAL));
            a_pp.query("collision_alpha", m_dparameters(ALPHA));
         }
         m_vthermal_method = INPUT_VTHERMAL;
      }
      else if (vthermal_method.matches("local vthermal")) {
         if (has_vthermal) {
            OV_ABORT("Either compute vthermal or specify vthermal but not both");
         }
         else if (has_alpha) {
            OV_ABORT("Either compute alpha or specify alpha but not both");
         }
         m_vthermal_method = LOCAL_VTHERMAL;
      }
      else if (vthermal_method.matches("global vthermal")) {
         if (has_vthermal) {
            OV_ABORT("Either compute vthermal or specify vthermal but not both");
         }
         else if (has_alpha) {
            OV_ABORT("Either compute alpha or specify alpha but not both");
         }
         m_vthermal_method = GLOBAL_VTHERMAL;
      }
      else {
         OV_ABORT("Unknown input for collision_vthermal_method");
      }
   }

   // The collision_vthermal_dt and collision_nuCoeff are required.
   if (!a_pp.query("collision_vthermal_dt", m_dparameters(VTHERMAL_DT))) {
      OV_ABORT("Must supply collision_vthermal_dt");
   }
   if (!a_pp.query("collision_nuCoeff", m_dparameters(NU))) {
      OV_ABORT("Must supply collision_nuCoeff");
   }

   // Read if the back reaction should be computed.  If it is, then the mass
   // ratio must be supplied.
   aString tmp("false");
   a_pp.query("collision_back_reaction", tmp);
   m_back_reaction = tmp.matches("false") ? false : true;
   if (m_back_reaction &&
       !a_pp.query("collision_massR", m_dparameters(MASSR))) {
      OV_ABORT("Must supply collision_massR");
   }
}


void
InterspeciesCollisionOperator::computeConstants(
   const tbox::Box& a_interior_box,
   const ProblemDomain& a_domain)
{
   // Dimension constants then delegate the computation of them to fortran.
   int n3a = a_interior_box.lower(2);
   int n3b = a_interior_box.upper(2);
   int n4a = a_interior_box.lower(3);
   int n4b = a_interior_box.upper(3);
   int nGhosts = m_iparameters(NUM_GHOSTS);
   Range R1hat(n3a-nGhosts, n3b+nGhosts);
   Range R2hat(n4a-nGhosts, n4b+nGhosts);
   m_khat.redim(R1hat, R2hat);
   m_hhat.redim(R1hat, R2hat);
   m_vx.redim(R1hat, R2hat);
   m_vy.redim(R1hat, R2hat);
   FORT_COMPUTE_INTERSPECIES_CONSTANTS(*m_khat.getDataPointer(),
      *m_hhat.getDataPointer(),
      *m_vx.getDataPointer(),
      *m_vy.getDataPointer(),
      BOX4D_TO_FORT(a_interior_box),
      nGhosts,
      PROBLEMDOMAIN_TO_FORT(a_domain),
      *m_dparameters.getDataPointer(),
      *m_iparameters.getDataPointer());
}


void
InterspeciesCollisionOperator::constructBackReactionArrays(
   const tbox::Box& a_local_box)
{
   // Local box limits.
   int nd1a = a_local_box.lower(0);
   int nd1b = a_local_box.upper(0);
   int nd2a = a_local_box.lower(1);
   int nd2b = a_local_box.upper(1);
   int nd3a = a_local_box.lower(2);
   int nd3b = a_local_box.upper(2);
   int nd4a = a_local_box.lower(3);
   int nd4b = a_local_box.upper(3);

   // Ranges for local arrays.
   Range R1(nd1a, nd1b);
   Range R2(nd2a, nd2b);
   Range R3(nd3a, nd3b);
   Range R4(nd4a, nd4b);

   // Dimension persistent 4D momenta, kinetic energy, and maxwellian
   // distribution.
   m_cMomx.redim(R1, R2, R3, R4);
   m_cMomx = 0.0;
   m_cMomy.redim(R1, R2, R3, R4);
   m_cMomy = 0.0;
   m_cKE.redim(R1, R2, R3, R4);
   m_cKE = 0.0;
   m_fM.redim(R1, R2, R3, R4);

   // Dimension intermediate 2D locally reduced momenta and kinetic energy.
   m_rMomx.redim(R1, R2);
   m_rMomx = 0.0;
   m_rMomy.redim(R1, R2);
   m_rMomy = 0.0;
   m_rKE.redim(R1, R2);
   m_rKE = 0.0;

   // Dimension persistent 2D globally reduced momenta and kinetic energy
   // denominator terms.
   m_IMomxD.redim(R1, R2);
   m_IMomxD = 0.0;
   m_IMomyD.redim(R1, R2);
   m_IMomyD = 0.0;
   m_IKED.redim(R1, R2);
   m_IKED = 0.0;

   // Dimension persistent 2D globally reduced momenta and kinetic energy
   // numerator terms.
   m_IMomxN.redim(R1, R2);
   m_IMomyN.redim(R1, R2);
   m_IKEN.redim(R1, R2);
}


void
InterspeciesCollisionOperator::constructBackReactionReductionSchedules(
   const tbox::Box& a_interior_box,
   const ProblemDomain& a_domain)
{
   // Build reduction schedules for momenta and kinetic energy denominators.
   Array<bool> collapse_dir(a_domain.dim(), false);
   collapse_dir[V1] = true;
   collapse_dir[V2] = true;
   double measure(a_domain.dx(2) * a_domain.dx(3));

   // x-momentum reduction schedule
   tbox::Pointer<ReductionOp> int_op_x_mom(new IntegralOp(measure));
   m_x_mom_reduction = new ReductionSchedule4D(m_rMomx,
      a_interior_box,
      a_domain.box(),
      collapse_dir,
      int_op_x_mom,
      m_processor_range,
      m_comm);

   // y-momentum reduction schedule
   tbox::Pointer<ReductionOp> int_op_y_mom(new IntegralOp(measure));
   m_y_mom_reduction = new ReductionSchedule4D(m_rMomy,
      a_interior_box,
      a_domain.box(),
      collapse_dir,
      int_op_y_mom,
      m_processor_range,
      m_comm);

   // KE reduction schedule
   tbox::Pointer<ReductionOp> int_op_ke(new IntegralOp(measure));
   m_ke_reduction = new ReductionSchedule4D(m_rKE,
      a_interior_box,
      a_domain.box(),
      collapse_dir,
      int_op_ke,
      m_processor_range,
      m_comm);
}


void
InterspeciesCollisionOperator::computeBackReactionMoments(
   const tbox::Box& a_local_box,
   const tbox::Box& a_interior_box,
   const ProblemDomain& a_domain)
{
   // Construct 4D momenta and kinetic energy temps.
   // Local box limits.
   int nd1a = a_local_box.lower(0);
   int nd1b = a_local_box.upper(0);
   int nd2a = a_local_box.lower(1);
   int nd2b = a_local_box.upper(1);
   int nd3a = a_local_box.lower(2);
   int nd3b = a_local_box.upper(2);
   int nd4a = a_local_box.lower(3);
   int nd4b = a_local_box.upper(3);

   // Ranges for local arrays.
   Range R1(nd1a, nd1b);
   Range R2(nd2a, nd2b);
   Range R3(nd3a, nd3b);
   Range R4(nd4a, nd4b);

   RealArray uMomx(R1, R2, R3, R4);
   RealArray uMomy(R1, R2, R3, R4);
   RealArray uKE(R1, R2, R3, R4);

   // Compute temporary 4D momenta and kinetic energy and persistent 4D
   // maxwellian distribution.
   FORT_COMPUTE_INTERSPECIES_TEMPS(*m_fM.getDataPointer(),
      *uMomx.getDataPointer(),
      *uMomy.getDataPointer(),
      *uKE.getDataPointer(),
      BOX4D_TO_FORT(a_local_box),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      *m_dparameters.getDataPointer());

   // Compute persistent 4D momenta and kinetic energy.
   FORT_APPEND_INTERSPECIES_COLLISION(*m_cMomx.getDataPointer(),
      *uMomx.getDataPointer(),
      *m_vx.getDataPointer(),
      *m_vy.getDataPointer(),
      *m_khat.getDataPointer(),
      *m_hhat.getDataPointer(),
      BOX4D_TO_FORT(a_local_box),
      BOX4D_TO_FORT(a_interior_box),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      *m_dparameters.getDataPointer(),
      *m_iparameters.getDataPointer());
   FORT_APPEND_INTERSPECIES_COLLISION(*m_cMomy.getDataPointer(),
      *uMomy.getDataPointer(),
      *m_vx.getDataPointer(),
      *m_vy.getDataPointer(),
      *m_khat.getDataPointer(),
      *m_hhat.getDataPointer(),
      BOX4D_TO_FORT(a_local_box),
      BOX4D_TO_FORT(a_interior_box),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      *m_dparameters.getDataPointer(),
      *m_iparameters.getDataPointer());
   FORT_APPEND_INTERSPECIES_COLLISION(*m_cKE.getDataPointer(),
      *uKE.getDataPointer(),
      *m_vx.getDataPointer(),
      *m_vy.getDataPointer(),
      *m_khat.getDataPointer(),
      *m_hhat.getDataPointer(),
      BOX4D_TO_FORT(a_local_box),
      BOX4D_TO_FORT(a_interior_box),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      *m_dparameters.getDataPointer(),
      *m_iparameters.getDataPointer());

   // Fill m_rMomx, m_rMomy, and m_rKE with quantities to be reduced to form
   // the back reaction denominator terms.
   FORT_COMPUTE_INTERSPECIES_DENOMS(*m_rMomx.getDataPointer(),
      *m_rMomy.getDataPointer(),
      *m_rKE.getDataPointer(),
      *m_cMomx.getDataPointer(),
      *m_cMomy.getDataPointer(),
      *m_cKE.getDataPointer(),
      BOX4D_TO_FORT(a_local_box),
      BOX4D_TO_FORT(a_interior_box),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      *m_dparameters.getDataPointer());

   // Reduce m_rMomx. m_rMomy, and m_rKE to get the back reaction denominator
   // terms.
   m_x_mom_reduction->execute(m_IMomxD);
   m_y_mom_reduction->execute(m_IMomyD);
   m_ke_reduction->execute(m_IKED);
}

} // end namespace Loki
