/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "PerturbedMaxwellianIC.H"
#include "JuttnerThermal.H"
#include "MaxwellianThermal.H"
#include "Directions.H"
#include "getMomentumF.H"
#include "Simulation.H"
#include "Loki_Defines.H"
#include "Loki_Utilities.H"

namespace Loki {

bool
PerturbedMaxwellianIC::isType(
   const string& a_name)
{
   // Kind of funky but the concept of name is a bit overloaded here.  There
   // are essentially 3 related initial conditions covered by this class.
   if (a_name.compare("Perturbed Maxwellian") == 0 ||
       a_name.compare("Landau damping") == 0 ||
       a_name.compare("Maxwellian with noise") == 0) {
      return true;
   }
   return false;
}


PerturbedMaxwellianIC::PerturbedMaxwellianIC(
   LokiInputParser& a_pp,
   const ProblemDomain& a_domain,
   double a_mass)
  : ICInterface(a_domain),
    m_ic_option(1)
{
   // Set the default parameters to some hopefully sane values.
   // The initial flow velocities are passed in from the KineticSpecies as
   // these are needed by both the species and its initializer.
   m_tx = 1.0;
   m_ty = 1.0;
   m_vx0 = 0.0;
   m_vy0 = 0.0;
   m_x_wave_number = 0.0;
   m_y_wave_number = 0.0;
   m_flow_vel_phi = 0.0;
   m_frac = 1.0;
   m_a = 0.0;
   m_kx1 = 0.5;
   m_ky1 = 0.5;
   m_b = 0.0;
   m_kx2 = 0.5;
   m_c = 0.0;
   m_ky2 = 0.5;
   m_spatial_phi = 0.0;

   // Get the parameters that the user really wants.
   parseParameters(a_pp);

   m_Lx = m_domain.dx(X1)*m_domain.numberOfCells(X1);

   if (m_maxwellian_thermal) {
      m_thermal = new MaxwellianThermal(m_vflowinitx,
         m_vflowinity,
         m_vx0,
         m_vy0,
         a_mass,
         m_tx,
         m_ty);
   }
   else {
      m_thermal = new JuttnerThermal(m_vflowinitx,
         m_vflowinity,
         m_vx0,
         m_vy0,
         a_mass,
         m_tx,
         m_ty);
   }
}


PerturbedMaxwellianIC::~PerturbedMaxwellianIC()
{
}


void
PerturbedMaxwellianIC::cache(
   const ParallelArray& a_u)
{
   const ParallelArray::Box& db = a_u.dataBox();
   const ParallelArray::Box& ib = a_u.interiorBox();
   double pi = 4.0*atan(1.0);

   if (m_factorable) {
      // Build m_fx and m_fv.
      ParallelArray::Box factored_box(CDIM);
      vector<int> num_global_cells(CDIM);
      for (int i = 0; i < CDIM; ++i) {
         factored_box.lower(i) = ib.lower(i);
         factored_box.upper(i) = ib.upper(i);
         num_global_cells[i] = m_domain.numberOfCells(i);
      }
      m_fx.partition(factored_box, CDIM, a_u.numGhosts(), num_global_cells);
      for (int i = CDIM; i < PDIM; ++i) {
         factored_box.lower(i-CDIM) = ib.lower(i);
         factored_box.upper(i-CDIM) = ib.upper(i);
         num_global_cells[i-CDIM] = m_domain.numberOfCells(i);
      }
      m_fv.partition(factored_box, CDIM, a_u.numGhosts(), num_global_cells);

      // Compute m_fx.
      for (int i2 = db.lower(X2); i2 <= db.upper(X2); ++i2) {
         double x2 = m_domain.lower(X2)+(i2+0.5)*m_domain.dx(X2);
         // Enforce periodicity.
         if (m_domain.isPeriodic(X2)) {
            if (x2 < m_domain.lower(X2)) {
               x2 += m_domain.numberOfCells(X2)*m_domain.dx(X2);
            }
            else if (x2 > m_domain.upper(X2)) {
               x2 -= m_domain.numberOfCells(X2)*m_domain.dx(X2);
            }
         }
         for (int i1 = db.lower(X1); i1 <= db.upper(X1); ++i1) {
            double fx;
            double x1 = m_domain.lower(X1)+(i1+0.5)*m_domain.dx(X1);
            // Enforce periodicity.
            if (m_domain.isPeriodic(X1)) {
               if (x1 < m_domain.lower(X1)) {
                  x1 += m_domain.numberOfCells(X1)*m_domain.dx(X1);
               }
               else if (x1 > m_domain.upper(X1)) {
                  x1 -= m_domain.numberOfCells(X1)*m_domain.dx(X1);
               }
            }

            // Compute m_fx at this point in configuration space.
            if (m_ic_option == 1) {
               // Perturbed Maxwellian
               fx = 1.0+m_a*cos(m_kx1*x1+m_spatial_phi)*cos(m_ky1*x2+m_spatial_phi)+
                        m_b*cos(m_kx2*x1+m_spatial_phi)+
                        m_c*cos(m_ky2*x2+m_spatial_phi);
            }
            else if (m_ic_option == 2) {
               // Landau damping
               fx = 1.0+m_a*cos(m_kx1*x1+m_ky1*x2+m_spatial_phi);
            }
            else if (m_ic_option == 3) {
               // Maxwellian with noise
               fx = 1.0;
               for (int k = 1; k <= static_cast<int>(m_noise_amp.size()); ++k) {
                  fx += m_noise_amp[k-1]*
                     cos(2.0*pi*k*(x1+m_noise_phase[k-1])/m_Lx+m_spatial_phi);
               }
            }
            m_fx(i1, i2) = fx;
         }
      }

      // Compute m_fv.
      for (int i4 = db.lower(V2); i4 <= db.upper(V2); ++i4) {
         double x4 = m_domain.lower(V2)+(i4+0.5)*m_domain.dx(V2);
         for (int i3 = db.lower(V1); i3 <= db.upper(V1); ++i3) {
            double x3 = m_domain.lower(V1)+(i3+0.5)*m_domain.dx(V1);
            // Compute m_fv at this point in velocity space.
            m_fv(i3, i4) = m_thermal->thermalFactor(0.0, x3, x4);
         }
      }
   }
   else {
      // Build m_f.
      vector<int> num_global_cells(PDIM);
      for (int i = 0; i < PDIM; ++i) {
         num_global_cells[i] = m_domain.numberOfCells(i);
      }
      m_f.partition(ib, PDIM, a_u.numGhosts(), num_global_cells);

      // Compute m_f.
      double fnorm = m_thermal->fnorm();
      for (int i4 = db.lower(V2); i4 <= db.upper(V2); ++i4) {
         double x4 = m_domain.lower(V2)+(i4+0.5)*m_domain.dx(V2);
         for (int i3 = db.lower(V1); i3 <= db.upper(V1); ++i3) {
            double x3 = m_domain.lower(V1)+(i3+0.5)*m_domain.dx(V1);
            for (int i2 = db.lower(X2); i2 <= db.upper(X2); ++i2) {
               double x2 = m_domain.lower(X2)+(i2+0.5)*m_domain.dx(X2);
               // Enforce periodicity.
               if (m_domain.isPeriodic(X2)) {
                  if (x2 < m_domain.lower(X2)) {
                     x2 += m_domain.numberOfCells(X2)*m_domain.dx(X2);
                  }
                  else if (x2 > m_domain.upper(X2)) {
                     x2 -= m_domain.numberOfCells(X2)*m_domain.dx(X2);
                  }
               }
               for (int i1 = db.lower(X1); i1 <= db.upper(X1); ++i1) {
                  double x1 = m_domain.lower(X1)+(i1+0.5)*m_domain.dx(X1);
                  // Enforce periodicity.
                  if (m_domain.isPeriodic(X1)) {
                     if (x1 < m_domain.lower(X1)) {
                        x1 += m_domain.numberOfCells(X1)*m_domain.dx(X1);
                     }
                     else if (x1 > m_domain.upper(X1)) {
                        x1 -= m_domain.numberOfCells(X1)*m_domain.dx(X1);
                     }
                  }

                  // Compute m_f at this point in phase space.
                  double spatial_factor =
                     cos(m_x_wave_number*x1+m_y_wave_number*x2+m_flow_vel_phi);
                  double fv = m_thermal->thermalFactor(spatial_factor, x3, x4);
                  if (m_ic_option == 1) {
                     // Perturbed Maxwellian
                     double fx =
                        1.0+m_a*cos(m_kx1*x1+m_spatial_phi)*cos(m_ky1*x2+m_spatial_phi)+
                        m_b*cos(m_kx2*x1+m_spatial_phi)+
                        m_c*cos(m_ky2*x2+m_spatial_phi);
                     m_f(i1, i2, i3, i4) = fnorm*fv*fx*m_frac;
                  }
                  else if (m_ic_option == 2) {
                     // Landau damping
                     double fx = 1.0+m_a*cos(m_kx1*x1+m_ky1*x2+m_spatial_phi);
                     m_f(i1, i2, i3, i4) = fnorm*fv*fx*m_frac;
                  }
                  else if (m_ic_option == 3) {
                     // Maxwellian with noise
                     double fx = 1.0;
                     for (int k = 1;
                          k <= static_cast<int>(m_noise_amp.size()); ++k) {
                        fx += m_noise_amp[k-1]*
                           cos(2.0*pi*k*(x1+m_noise_phase[k-1])/m_Lx+m_spatial_phi);
                     }
                     m_f(i1, i2, i3, i4) = fx*fnorm*fv*m_frac;
                  }
               }
            }
         }
      }
   }
}


void
PerturbedMaxwellianIC::set(
   ParallelArray& a_u) const
{
   const ParallelArray::Box& data_box = a_u.dataBox();
   for (int i4 = data_box.lower(V2); i4 <= data_box.upper(V2); ++i4) {
      for (int i3 = data_box.lower(V1); i3 <= data_box.upper(V1); ++i3) {
         for (int i2 = data_box.lower(X2); i2 <= data_box.upper(X2); ++i2) {
            for (int i1 = data_box.lower(X1); i1 <= data_box.upper(X1); ++i1) {
               a_u(i1, i2, i3, i4) = getIC_At_Pt(i1, i2, i3, i4);
            }
         }
      }
   }
}


double
PerturbedMaxwellianIC::getIC_At_Pt(
   int a_i1,
   int a_i2,
   int a_i3,
   int a_i4) const
{
   // To avoid order of operation differences with earlier versions of this
   // initial condition that did not cache their result the cached factored
   // results are multiplied in this specific order.  Do not change.
   if (m_factorable) {
      if (m_ic_option == 3) {
         return m_fx(a_i1, a_i2)*m_thermal->fnorm()*m_fv(a_i3, a_i4)*m_frac;
      }
      else {
         return m_thermal->fnorm()*m_fv(a_i3, a_i4)*m_fx(a_i1, a_i2)*m_frac;
      }
   }
   else {
      return m_f(a_i1, a_i2, a_i3, a_i4);
   }
}


void
PerturbedMaxwellianIC::printParameters() const
{
   // Write out the input parameters.  Note that not all parameters are used
   // by all 3 of the initial conditions embedded in this class.
   if (m_ic_option == 1) {
      Loki_Utilities::printF("\n  Using built in initial conditions:\n");
   }
   else if (m_ic_option == 2) {
      Loki_Utilities::printF("\n  Using Landau damping initial conditions:\n");
   }
   else if (m_ic_option == 3) {
      Loki_Utilities::printF("\n  Using Maxwellian initial conditions with noise:\n");
   }

   if (m_maxwellian_thermal) {
      Loki_Utilities::printF("    Maxwellian thermal\n");
   }
   else {
      Loki_Utilities::printF("    Juttner thermal\n");
   }
   Loki_Utilities::printF("    tx                  = %e\n", m_tx);
   Loki_Utilities::printF("    ty                  = %e\n", m_ty);
   Loki_Utilities::printF("    vx0                 = %e\n", m_vx0);
   Loki_Utilities::printF("    vy0                 = %e\n", m_vy0);
   Loki_Utilities::printF("    vflowinitx          = %e\n", m_vflowinitx);
   Loki_Utilities::printF("    vflowinity          = %e\n", m_vflowinity);
   Loki_Utilities::printF("    frac                = %e\n", m_frac);

   if (m_ic_option == 1 || m_ic_option == 2) {
      Loki_Utilities::printF("    A                   = %e\n", m_a);
      Loki_Utilities::printF("    kx1                 = %e\n", m_kx1);
      Loki_Utilities::printF("    ky1                 = %e\n", m_ky1);
      if (m_ic_option == 1) {
         Loki_Utilities::printF("    B                   = %e\n", m_b);
         Loki_Utilities::printF("    kx2                 = %e\n", m_kx2);
         Loki_Utilities::printF("    C                   = %e\n", m_c);
         Loki_Utilities::printF("    ky2                 = %e\n", m_ky2);
      }
   }
   Loki_Utilities::printF("    x wave number       = %e\n", m_x_wave_number);
   Loki_Utilities::printF("    y wave number       = %e\n", m_y_wave_number);
   Loki_Utilities::printF("    flow velocity phase = %e\n", m_flow_vel_phi);
   Loki_Utilities::printF("    spatial phase       = %e\n", m_spatial_phi);

   if (m_ic_option == 3) {
      int num_modes = static_cast<int>(m_noise_amp.size());
      Loki_Utilities::printF("    perturbing %i modes\n", num_modes);
      Loki_Utilities::printF("    mode #: amplitude, phase\n");
      for (int k(0); k < num_modes; ++k) {
         Loki_Utilities::printF("    mode %i: %e, %e\n",
            k,
            m_noise_amp[k],
            m_noise_phase[k]);
      }
   }
}


void
PerturbedMaxwellianIC::parseParameters(
   LokiInputParser& a_pp)
{
   // Figure out which variant of this initializer the user wants.  This is
   // required.
   string ic_name("Perturbed Maxwellian");
   a_pp.query("name", ic_name);
   if (ic_name.compare("Perturbed Maxwellian") == 0) {
      m_ic_option = 1;
   }
   else if (ic_name.compare("Landau damping") == 0) {
      m_ic_option = 2;
   }
   else if (ic_name.compare("Maxwellian with noise") == 0) {
      m_ic_option = 3;
   }
   else {
      LOKI_ABORT("Unknown name input");
   }

   // Determine if the user wants Maxwellian or Juttner thermal behavior.
   string test_str = "true";
   a_pp.query("maxwellian_thermal", test_str);
   m_maxwellian_thermal = test_str.compare("true") == 0 ? true : false;

   // None of the rest of the input is required so if the user doesn't specify
   // something they get the default which may not be a sane value.  It might
   // be wise to make some or all of these required.  Different subsets of
   // these inputs are relevant for a given value of m_ic_option.
   int num_noise(0);
   if (a_pp.query("number_of_noisy_modes", num_noise)) {
      m_noise_amp.resize(num_noise);
      a_pp.getarr("noise_amplitudes", m_noise_amp, 0, num_noise);

      m_noise_phase.resize(num_noise);
      a_pp.getarr("noise_phases", m_noise_phase, 0, num_noise);
   }

   if (a_pp.contains("alpha") || a_pp.contains("beta")) {
      LOKI_ABORT("Parameters alpha and beta are deprecated in favor of tx and ty.");
   }
   a_pp.query("tx",            m_tx);
   a_pp.query("ty",            m_ty);
   a_pp.query("vx0",           m_vx0);
   a_pp.query("vy0",           m_vy0);
   if (m_vx0 == 0.0 && m_vy0 == 0) {
      m_factorable = true;
   }
   a_pp.query("frac",          m_frac);
   a_pp.query("vflowinitx",    m_vflowinitx);
   a_pp.query("vflowinity",    m_vflowinity);
   a_pp.query("A",             m_a);
   a_pp.query("B",             m_b);
   a_pp.query("C",             m_c);
   a_pp.query("kx1",           m_kx1);
   a_pp.query("kx2",           m_kx2);
   a_pp.query("ky1",           m_ky1);
   a_pp.query("ky2",           m_ky2);
   a_pp.query("x_wave_number", m_x_wave_number);
   a_pp.query("y_wave_number", m_y_wave_number);
   a_pp.query("phase",         m_flow_vel_phi);
   a_pp.query("spatial_phase", m_spatial_phi);
}

} // end namespace Loki
