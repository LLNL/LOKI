/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "InterpenetratingStreamIC.H"
#include "JuttnerThermal.H"
#include "MaxwellianThermal.H"
#include "Directions.H"
#include "getMomentumF.H"
#include "Simulation.H"
#include "Loki_Defines.H"
#include "Loki_Utilities.H"

namespace Loki {

bool
InterpenetratingStreamIC::isType(
   const string& a_name)
{
   if (a_name.compare("Interpenetrating Stream") == 0) {
      return true;
   }
   return false;
}


InterpenetratingStreamIC::InterpenetratingStreamIC(
   LokiInputParser& a_pp,
   const ProblemDomain& a_domain,
   double a_mass,
   int a_n_ghosts)
  : ICInterface(a_domain),
    m_box_syntax(false),
    m_floor(0.0)
{
   // This initial condition is factorable by construction.
   m_factorable = true;

   // Get user input.
   parseParameters(a_pp);

   // Precompute some constants derived from the input parameters.
   if (m_box_syntax) {
      if (m_maxwellian_thermal) {
         m_thermal = new MaxwellianThermal(m_vflowinitx,
            m_vflowinity,
            0.0,
            0.0,
            a_mass,
            m_tx,
            m_ty);
      }
      else {
         m_thermal = new JuttnerThermal(m_vflowinitx,
            m_vflowinity,
            0.0,
            0.0,
            a_mass,
            m_tx,
            m_ty);
      }
   }
   else {
      if (m_maxwellian_thermal) {
         m_thermal = new MaxwellianThermal(m_vl0,
            m_vt0,
            0.0,
            0.0,
            a_mass,
            m_tl,
            m_tt);
         if (m_two_sided) {
            m_thermal2 = new MaxwellianThermal(-m_vl0,
               m_vt0,
               0.0,
               0.0,
               a_mass,
               m_tl,
               m_tt);
         }
      }
      else {
         m_thermal = new JuttnerThermal(m_vl0,
            m_vt0,
            0.0,
            0.0,
            a_mass,
            m_tl,
            m_tt);
         if (m_two_sided) {
            m_thermal2 = new JuttnerThermal(-m_vl0,
               m_vt0,
               0.0,
               0.0,
               a_mass,
               m_tl,
               m_tt);
         }
      }
   }
}


InterpenetratingStreamIC::~InterpenetratingStreamIC()
{
}


void
InterpenetratingStreamIC::cache(
   const ParallelArray& a_u)
{
   const ParallelArray::Box& db = a_u.dataBox();
   const ParallelArray::Box& ib = a_u.interiorBox();

   // Build m_fx and m_fv.
   ParallelArray::Box factored_box(CDIM);
   vector<int> num_global_cells(CDIM);
   for (int i = 0; i < CDIM; ++i) {
      factored_box.lower(i) = ib.lower(i);
      factored_box.upper(i) = ib.upper(i);
      num_global_cells[i] = m_domain.numberOfCells(i);
   }
   m_fx.partition(factored_box, CDIM, a_u.numGhosts(), num_global_cells);
   // If necessary, build m_fx2.
   if (!m_box_syntax && (m_two_sided || m_centered)) {
      m_fx2.partition(factored_box, CDIM, a_u.numGhosts(), num_global_cells);
   }
   for (int i = CDIM; i < PDIM; ++i) {
      factored_box.lower(i-CDIM) = ib.lower(i);
      factored_box.upper(i-CDIM) = ib.upper(i);
      num_global_cells[i-CDIM] = m_domain.numberOfCells(i);
   }
   m_fv.partition(factored_box, CDIM, a_u.numGhosts(), num_global_cells);
   int nd3a = m_fv.dataBox().lower(V1-CDIM);
   int nd3b = m_fv.dataBox().upper(V1-CDIM);
   int nd4a = m_fv.dataBox().lower(V2-CDIM);
   int nd4b = m_fv.dataBox().upper(V2-CDIM);
   // If necessary, build m_fv2.
   if (!m_box_syntax && m_two_sided) {
      m_fv2.partition(factored_box, CDIM, a_u.numGhosts(), num_global_cells);
   }

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
         if (m_box_syntax) {
            // When boxes are specified, a spatial weighting due to the the
            // Heaviside for each box is needed at each point in configuration
            // space.
            double spatial_wt = 0.0;
            for (int box = 0; box < m_num_boxes; ++box) {
               double h = 1.0;
               double edge = m_box_x_hi[box];
               if (edge != m_domain.upper(X1)) {
                  h = 0.5*(1.0+erf(-m_beta*(x1-edge)));
               }
               edge = m_box_x_lo[box];
               if (edge != m_domain.lower(X1)) {
                  h *= 0.5*(1.0+erf(m_beta*(x1-edge)));
               }
               edge = m_box_y_hi[box];
               if (edge != m_domain.upper(X2)) {
                  h *= 0.5*(1.0+erf(-m_beta*(x2-edge)));
               }
               edge = m_box_y_lo[box];
               if (edge != m_domain.lower(X2)) {
                  h *= 0.5*(1.0+erf(m_beta*(x2-edge)));
               }
               spatial_wt += m_box_frac[box]*h;
            }
            m_fx(i1, i2) = spatial_wt;
         }
         else {
            double xi0 = -m_d;
            double xi = x1*cos(m_theta)+x2*sin(m_theta);
            if (m_two_sided) {
               m_fx(i1, i2) = m_floor/2.0 +
                  (m_frac-m_floor)*0.5*(1.0+erf(-m_beta*(xi-xi0)));
               m_fx2(i1, i2) = m_floor/2.0 +
                  (m_frac2-m_floor)*0.5*(1.0+erf(m_beta*(xi+xi0)));
            }
            else if (m_centered) {
               m_fx(i1, i2) = m_floor +
                  (sqrt(m_frac)-m_floor)*0.5*(1.0+erf(-m_beta*(xi-xi0)));
               m_fx2(i1, i2) = m_floor +
                  (sqrt(m_frac)-m_floor)*0.5*(1.0+erf(m_beta*(xi+xi0)));
            }
            else {
               m_fx(i1, i2) = m_floor +
                  (m_frac-m_floor)*0.5*(1.0+erf(-m_beta*(xi-xi0)));
            }
         }
      }
   }

   // Compute m_fv.
   for (int i4 = db.lower(V2); i4 <= db.upper(V2); ++i4) {
      double x4 = m_domain.lower(V2)+(i4+0.5)*m_domain.dx(V2);
      for (int i3 = db.lower(V1); i3 <= db.upper(V1); ++i3) {
         double x3 = m_domain.lower(V1)+(i3+0.5)*m_domain.dx(V1);
         if (m_box_syntax) {
            m_fv(i3, i4) =
               m_thermal->fnorm()*m_thermal->thermalFactor(0.0, x3, x4);
         }
         else {
            double vl = x3*cos(m_theta) + x4*sin(m_theta);
            double vt = -x3*sin(m_theta) + x4*cos(m_theta);
            m_fv(i3, i4) =
               m_thermal->fnorm()*m_thermal->thermalFactor(0.0, vl, vt);
            if (m_two_sided) {
               m_fv2(i3, i4) =
                  m_thermal2->fnorm()*m_thermal2->thermalFactor(0.0, vl, vt);
            }
         }
      }
   }
}


void
InterpenetratingStreamIC::set(
   ParallelArray& a_u) const
{
   // For each zone, compute the initial condition.
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
InterpenetratingStreamIC::getIC_At_Pt(
   int a_i1,
   int a_i2,
   int a_i3,
   int a_i4) const
{
   // To avoid order of operation differences with earlier versions of this
   // initial condition that did not cache their result the cached factored
   // results are multiplied in this specific order.  Do not change.
   if (!m_box_syntax && m_two_sided) {
      return (m_fx(a_i1, a_i2)*m_fv(a_i3, a_i4)+
              m_fx2(a_i1, a_i2)*m_fv2(a_i3, a_i4));
   }
   else if (!m_box_syntax && m_centered) {
      return m_fv(a_i3, a_i4)*m_fx(a_i1, a_i2)*m_fx2(a_i1, a_i2);
   }
   else {
      return m_fv(a_i3, a_i4)*m_fx(a_i1, a_i2);
   }
}


void
InterpenetratingStreamIC::printParameters() const
{
   if (m_box_syntax) {
      printParametersBox();
   }
   else {
      printParametersHalfPlane();
   }
}


void
InterpenetratingStreamIC::printParametersBox() const
{
   Loki_Utilities::printF("\n  Using box syntax slab initial conditions:\n");
   if (m_maxwellian_thermal) {
      Loki_Utilities::printF("    Maxwellian thermal\n");
   }
   else {
      Loki_Utilities::printF("    Juttner thermal\n");
   }
   Loki_Utilities::printF("    X temperature                  = %e\n", m_tx);
   Loki_Utilities::printF("    Y temperature                  = %e\n", m_ty);
   Loki_Utilities::printF("    X drift speed                  = %e\n",
      m_vflowinitx);
   Loki_Utilities::printF("    Y drift speed                  = %e\n",
      m_vflowinity);
   Loki_Utilities::printF("    Transition sharpness parameter = %e\n",m_beta);
   Loki_Utilities::printF("    Number of boxes                = %d\n",
      m_num_boxes);
   for (int i = 0; i < m_num_boxes; ++i) {
      Loki_Utilities::printF("        Box %d x_lo                = %e\n",
         m_box_x_lo[i]);
      Loki_Utilities::printF("        Box %d x_hi                = %e\n",
         m_box_x_hi[i]);
      Loki_Utilities::printF("        Box %d y_lo                = %e\n",
         m_box_y_lo[i]);
      Loki_Utilities::printF("        Box %d y_hi                = %e\n",
         m_box_y_hi[i]);
      Loki_Utilities::printF("        Box %d frac                = %e\n",
         m_box_frac[i]);
   }
}


void
InterpenetratingStreamIC::printParametersHalfPlane() const
{
   Loki_Utilities::printF("\n  Using half-plane syntax slab initial conditions:\n");
   if (m_maxwellian_thermal) {
      Loki_Utilities::printF("    Maxwellian thermal\n");
   }
   else {
      Loki_Utilities::printF("    Juttner thermal\n");
   }
   Loki_Utilities::printF("    Longitudinal temperature                  = %e\n",
      m_tl);
   Loki_Utilities::printF("    Transverse temperature                    = %e\n",
      m_tt);
   Loki_Utilities::printF("    Drift direction                           = %e\n",
      m_theta);
   Loki_Utilities::printF("    Longitudinal drift speed                  = %e\n",
      m_vl0);
   Loki_Utilities::printF("    Transverse drift speed                    = %e\n",
      m_vt0);
   Loki_Utilities::printF("    Distance from origin                      = %e\n",
      m_d);
   Loki_Utilities::printF("    Transition sharpness parameter            = %e\n",
      m_beta);
   Loki_Utilities::printF("    Species minimum relative weight           = %e\n",
      m_floor);
   Loki_Utilities::printF("    Species maximum relative weight           = %e\n",
      m_frac);
   if (m_two_sided) {
      Loki_Utilities::printF("    Two sided species maximum relative weight = %e\n",
         m_frac2);
   }
   if (m_centered) {
      Loki_Utilities::printF("    Centered slab\n");
   }
}


void
InterpenetratingStreamIC::parseParameters(
   LokiInputParser& a_pp)
{
   // See if the syntax type is specified.  The half plane syntax is the
   // default.
   if (a_pp.contains("syntax")) {
      string input_syntax("half plane");
      a_pp.query("syntax", input_syntax);
      if (input_syntax.compare("half plane") == 0) {
         m_box_syntax = false;
      }
      else if (input_syntax.compare("box") == 0) {
         m_box_syntax = true;
      }
      else {
         LOKI_ABORT("Unknown input syntax.");
      }
   }

   // Determine if the user wants Maxwellian or Juttner thermal behavior.
   string test_str = "true";
   a_pp.query("maxwellian_thermal", test_str);
   m_maxwellian_thermal = test_str.compare("true") == 0 ? true : false;

   if (m_box_syntax) {
      parseParametersBox(a_pp);
   }
   else {
      if (Simulation::s_DO_RELATIVITY) {
         LOKI_ABORT("Currently can not use half plane syntax with relativistic problems.");
      }
      parseParametersHalfPlane(a_pp);
   }
}


void
InterpenetratingStreamIC::parseParametersBox(
   LokiInputParser& a_pp)
{
   // Read the input.  The box bounds are the domain bounds unless the user
   // specifies.  This makes it easy to place boxes at the physical domain.
   if (!a_pp.query("num_boxes", m_num_boxes)) {
      LOKI_ABORT("Number of boxes is required.");
   }
   m_box_x_hi.resize(m_num_boxes);
   m_box_x_lo.resize(m_num_boxes);
   m_box_y_hi.resize(m_num_boxes);
   m_box_y_lo.resize(m_num_boxes);
   m_box_frac.resize(m_num_boxes);
   for (int i = 0; i < m_num_boxes; ++i) {
      ostringstream input_string;
      input_string << "box_" << i+1 << "_x_hi";
      if (a_pp.contains(input_string.str().c_str())) {
         a_pp.query(input_string.str().c_str(), m_box_x_hi[i]);
      }
      else {
         m_box_x_hi[i] = m_domain.upper(X1);
      }
      input_string.str("");
      input_string << "box_" << i+1 << "_x_lo";
      if (a_pp.contains(input_string.str().c_str())) {
         a_pp.query(input_string.str().c_str(), m_box_x_lo[i]);
      }
      else {
         m_box_x_lo[i] = m_domain.lower(X1);
      }
      input_string.str("");
      input_string << "box_" << i+1 << "_y_hi";
      if (a_pp.contains(input_string.str().c_str())) {
         a_pp.query(input_string.str().c_str(), m_box_y_hi[i]);
      }
      else {
         m_box_y_hi[i] = m_domain.upper(X2);
      }
      input_string.str("");
      input_string << "box_" << i+1 << "_y_lo";
      if (a_pp.contains(input_string.str().c_str())) {
         a_pp.query(input_string.str().c_str(), m_box_y_lo[i]);
      }
      else {
         m_box_y_lo[i] = m_domain.lower(X2);
      }
      input_string.str("");
      input_string << "box_" << i+1 << "_frac";
      if (a_pp.contains(input_string.str().c_str())) {
         a_pp.query(input_string.str().c_str(), m_box_frac[i]);
      }
      else {
         LOKI_ABORT("Box fraction is required.");
      }
      input_string.str("");
   }
   // The sign of the transition sharpness parameter varies depending on which
   // box edge is being formed so make it postive when read in.
   if (!a_pp.query("beta", m_beta)) {
      LOKI_ABORT("Transition sharpness parameter is required.");
   }
   else {
      m_beta = abs(m_beta);
   }
   if (!a_pp.query("tx", m_tx)) {
      LOKI_ABORT("X temperature is required.");
   }
   if (!a_pp.query("ty", m_ty)) {
      LOKI_ABORT("Y temperature is required.");
   }
   a_pp.query("vflowinitx", m_vflowinitx);
   a_pp.query("vflowinity", m_vflowinity);
   if (a_pp.contains("vl0")) {
      LOKI_ABORT("InterpenetratingStream vl0 input no longer used in favor of vflowinitx.");
   }
   if (a_pp.contains("vt0")) {
      LOKI_ABORT("InterpenetratingStream vt0 input no longer used in favor of vflowinity.");
   }
}


void
InterpenetratingStreamIC::parseParametersHalfPlane(
   LokiInputParser& a_pp)
{
   // Read the input.  All parameters are required as there are no sane
   // defaults.
   if (!a_pp.query("tl", m_tl)) {
      LOKI_ABORT("Longitudinal temperature is required.");
   }
   if (!a_pp.query("tt", m_tt)) {
      LOKI_ABORT("Transverse temperature is required.");
   }
   if (!a_pp.query("theta", m_theta)) {
      LOKI_ABORT("Drift direction is required.");
   }
   a_pp.query("vflowinitx", m_vflowinitx);
   a_pp.query("vflowinity", m_vflowinity);
   if (a_pp.contains("vl0")) {
      LOKI_ABORT("InterpenetratingStream vl0 input no longer used in favor of vflowinitx.");
   }
   if (a_pp.contains("vt0")) {
      LOKI_ABORT("InterpenetratingStream vt0 input no longer used in favor of vflowinity.");
   }
   m_vl0 = m_vflowinitx*cos(m_theta) + m_vflowinity*sin(m_theta);
   m_vt0 = -m_vflowinitx*sin(m_theta) + m_vflowinity*cos(m_theta);
   if (!a_pp.query("d", m_d)) {
      LOKI_ABORT("Distance from origin is required.");
   }
   if (!a_pp.query("beta", m_beta)) {
      LOKI_ABORT("Transition sharpness parameter is required.");
   }
   if (!a_pp.query("frac", m_frac)) {
      LOKI_ABORT("Species relative weight parameter is required.");
   }
   a_pp.query("floor", m_floor);
   string test_str = "false";
   a_pp.query("two_sided", test_str);
   m_two_sided = test_str.compare("true") == 0 ? true : false;
   string test_str2 = "false";
   a_pp.query("centered", test_str2);
   m_centered = test_str2.compare("true") == 0 ? true : false;
   if (m_two_sided && m_centered) {
      LOKI_ABORT("Only one of a two sided or centered slab may be specified.");
   }
   if (m_two_sided && !a_pp.query("frac2", m_frac2)) {
      LOKI_ABORT("Two sided species relative weight parameter is required.");
   }
}

} // end namespace Loki
