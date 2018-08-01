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
#include "PerturbedMaxwellianIC.H"
#include "Array.H"
#include "BoxOps.H"
#include "PerturbedMaxwellianICF.H"
#include "Loki_Utilities.H"

namespace Loki {

bool
PerturbedMaxwellianIC::isType(
   const aString& a_name)
{
   // Kind of funky but the concept of name is a bit overloaded here.  There
   // are essentially 3 related initial conditions covered by this class.
   if (a_name.matches("Perturbed Maxwellian") ||
       a_name.matches("Landau damping") ||
       a_name.matches("Maxwellian with noise")) {
      return true;
   }
   return false;
}


PerturbedMaxwellianIC::PerturbedMaxwellianIC(
   ParmParse& a_pp,
   real a_vflowinitx,
   real a_vflowinity)
   : m_ic_option(1)
{
   // Set the default parameters to some hopefully sane values.
   // The initial flow velocities are passed in from the KineticSpecies as
   // these are needed by both the species and its initializer.
   m_parameters.resize(NUM_PARAMS);
   m_parameters(ALPHA)         = 1.0;
   m_parameters(BETA)          = 1.0;
   m_parameters(VX0)           = 0.0;
   m_parameters(VY0)           = 0.0;
   m_parameters(FRAC)          = 1.0;
   m_parameters(A)             = 0.0;
   m_parameters(KX1)           = 0.5;
   m_parameters(KY1)           = 0.5;
   m_parameters(B)             = 0.0;
   m_parameters(KX2)           = 0.5;
   m_parameters(C)             = 0.0;
   m_parameters(KY2)           = 0.5; 
   m_parameters(X_WAVE_NUMBER) = 0.0;
   m_parameters(Y_WAVE_NUMBER) = 0.0;
   m_parameters(FLOW_VEL_PHI)  = 0.0;
   m_parameters(SPATIAL_PHI)   = 0.0;
   m_parameters(VFLOWINITX)    = a_vflowinitx;
   m_parameters(VFLOWINITY)    = a_vflowinity;

   m_noise_amp.resize(0);
   m_noise_phase.resize(0);

   // Get the parameters that the user really wants.
   parseParameters(a_pp);
}


PerturbedMaxwellianIC::~PerturbedMaxwellianIC()
{
}


void
PerturbedMaxwellianIC::set(
   RealArray& a_u,
   const ProblemDomain& a_domain,
   const tbox::Box& a_grown_global_box,
   real a_time) const
{
   NULL_USE(a_time);
   NULL_USE(a_grown_global_box);

   // This is just delegated to fortran.
   tbox::Box local_box(BoxOps::getLocalBox(a_u));
   setIC_phase_4D(BOX4D_TO_FORT(local_box),
      *a_u.getDataPointer(),
      m_ic_option,
      *m_parameters.getDataPointer(),
      PROBLEMDOMAIN_TO_FORT(a_domain),
      INTVECTOR4D_TO_FORT(a_domain.numberOfCells()),
      *m_noise_amp.getDataPointer(),
      *m_noise_phase.getDataPointer(),
      m_noise_amp.elementCount());
}


void
PerturbedMaxwellianIC::printParameters() const
{
   // Write out the input parameters.  Note that not all parameters are used
   // by all 3 of the initial conditions embedded in this class.
   if (m_ic_option == 1) {
      printF("\n  Using built in initial conditions:\n");
   }
   else if (m_ic_option == 2) {
      printF("\n  Using Landau damping initial conditions:\n");
   }
   else if (m_ic_option == 3) {
      printF("\n  Using Maxwellian initial conditions with noise:\n");
   }

   printF("    alpha               = %e\n", m_parameters(ALPHA));
   printF("    beta                = %e\n", m_parameters(BETA));
   printF("    vx0                 = %e\n", m_parameters(VX0));
   printF("    vy0                 = %e\n", m_parameters(VY0));
   printF("    frac                = %e\n", m_parameters(FRAC));

   if (m_ic_option == 1 || m_ic_option == 2) {
      printF("    A                   = %e\n", m_parameters(A));
      printF("    kx1                 = %e\n", m_parameters(KX1));
      printF("    ky1                 = %e\n", m_parameters(KY1));
      if (m_ic_option == 1) {
         printF("    B                   = %e\n", m_parameters(B));
         printF("    kx2                 = %e\n", m_parameters(KX2));
         printF("    C                   = %e\n", m_parameters(C));
         printF("    ky2                 = %e\n", m_parameters(KY2));
      }
   }
   printF("    x wave number       = %e\n", m_parameters(X_WAVE_NUMBER));
   printF("    y wave number       = %e\n", m_parameters(Y_WAVE_NUMBER));
   printF("    flow velocity phase = %e\n", m_parameters(FLOW_VEL_PHI));
   printF("    spatial phase       = %e\n", m_parameters(SPATIAL_PHI));

   if (m_ic_option == 3) {
      printF("    perturbing %i modes\n", m_noise_amp.elementCount());
      printF("    mode #: amplitude, phase\n");
      for (int k(0); k < m_noise_amp.elementCount(); ++k) {
         printF("    mode %i: %e, %e\n", k, m_noise_amp(k), m_noise_phase(k));
      }
   }
}


void
arrayToRealArray(
   RealArray& a_real_array,
   const Array<double>& a_array)
{
   // The ParmParse interface supports reading Arrays but our data members are
   // RealArrays.  They could be Arrays as their API does everything that we
   // need a RealArray to do but that's the way this class has been for ages.
   a_real_array.resize(static_cast<int>(a_array.length()));
   for (int i(0); i < a_array.length(); ++i) {
      a_real_array(i) = a_array[i];
   }
}


void
PerturbedMaxwellianIC::parseParameters(
   ParmParse& a_pp)
{
   // Figure out which variant of this initializer the user wants.  This is
   // required.
   aString ic_name("Perturbed Maxwellian");
   a_pp.query("name", ic_name);
   if (ic_name.matches("Perturbed Maxwellian")) {
      m_ic_option = 1;
   }
   else if (ic_name.matches("Landau damping")) {
      m_ic_option = 2;
   }
   else if (ic_name.matches("Maxwellian with noise")) {
      m_ic_option = 3;
   }
   else {
      OV_ABORT("Unknown name input");
   }

   // None of the rest of the input is required so if the user doesn't specify
   // something they get the default which may not be a sane value.  It might
   // be wise to make some or all of these required.  Different subsets of
   // these inputs are relevant for a given value of m_ic_option.
   int num_noise(0);
   if (a_pp.query("number_of_noisy_modes", num_noise)) {
      Array<double> noise_amp(num_noise);
      a_pp.getarr("noise_amplitudes", noise_amp, 0, num_noise);
      arrayToRealArray(m_noise_amp, noise_amp);

      Array<double> noise_phase(num_noise);
      a_pp.getarr("noise_phases", noise_phase, 0, num_noise);
      arrayToRealArray(m_noise_phase, noise_phase);
   }

   a_pp.query("alpha",         m_parameters(ALPHA));
   a_pp.query("beta",          m_parameters(BETA));
   a_pp.query("vx0",           m_parameters(VX0));
   a_pp.query("vy0",           m_parameters(VY0));
   a_pp.query("frac",          m_parameters(FRAC));
   a_pp.query("A",             m_parameters(A));
   a_pp.query("B",             m_parameters(B));
   a_pp.query("C",             m_parameters(C));
   a_pp.query("kx1",           m_parameters(KX1));
   a_pp.query("kx2",           m_parameters(KX2));
   a_pp.query("ky1",           m_parameters(KY1));
   a_pp.query("ky2",           m_parameters(KY2));
   a_pp.query("x_wave_number", m_parameters(X_WAVE_NUMBER));
   a_pp.query("y_wave_number", m_parameters(Y_WAVE_NUMBER));
   a_pp.query("phase",         m_parameters(FLOW_VEL_PHI));
   a_pp.query("spatial_phase", m_parameters(SPATIAL_PHI));
}

} // end namespace Loki
