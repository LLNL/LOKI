/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "ExternalDistKrook.H"
#include "Directions.H"
#include "Loki_Defines.H"
#include "Loki_Utilities.H"

namespace Loki {

ExternalDistKrookLayer::ExternalDistKrookLayer(
   const tbox::Dimension& a_dim,
   LokiInputParser& a_pp,
   const ProblemDomain& a_domain)
   : m_dim(a_dim),
     m_x_lo_lo(a_dim),
     m_x_hi_lo(a_dim),
     m_x_lo_hi(a_dim),
     m_x_hi_hi(a_dim),
     m_has_layer_lo(a_dim, false),
     m_has_layer_hi(a_dim, false),
     m_has_layer(false),
     m_overlaps(false),
     m_coefficient(1.0),
     m_domain(a_domain),
     m_external_distribution_ic(0)
{
   // Initialize upper and lower ends of ExternalDistKrookLayer to the full
   // domain and that there is no layer in any direction.
   for (int d = 0; d < m_dim; ++d) {
      m_x_lo_lo[d] = m_domain.lower(d);
      m_x_hi_lo[d] = m_domain.lower(d);
      m_x_lo_hi[d] = m_domain.upper(d);
      m_x_hi_hi[d] = m_domain.upper(d);
   }

   // Read the ExternalDistKrookLayer specified by the user and make sure the
   // user is sane.
   parseParameters(a_pp);
   checkParameters();
}


ExternalDistKrookLayer::~ExternalDistKrookLayer()
{
   if (m_external_distribution_ic != 0) {
      delete m_external_distribution_ic;
   }
}


void
ExternalDistKrookLayer::initialize(
   const ParallelArray& a_dist_func,
   int a_spatial_solution_order)
{
   // Determine if this domain overlaps the layer.  If it does, then create the
   // external distribution initial condition object and have it cache the
   // initial condition.  Then compute the constant layer coefficients, m_nu.
   double xmin = m_domain.lower(X1);
   double ymin = m_domain.lower(X2);
   double dx = m_domain.dx(X1);
   double dy = m_domain.dx(X2);
   double xlolo = lolo(X1);
   double xhilo = hilo(X1);
   double xlohi = lohi(X1);
   double xhihi = hihi(X1);
   double ylolo = lolo(X2);
   double yhilo = hilo(X2);
   double ylohi = lohi(X2);
   double yhihi = hihi(X2);
   const ParallelArray::Box& ib = a_dist_func.interiorBox();
   int n1a = ib.lower(X1);
   int n1b = ib.upper(X1);
   int n2a = ib.lower(X2);
   int n2b = ib.upper(X2);
   m_overlaps =
      (xmin+(0.5+n1b)*dx > lolo(X1) && xmin+(0.5+n1a)*dx < hihi(X1)) ||
      (ymin+(0.5+n2b)*dy > lolo(X2) && ymin+(0.5+n2a)*dy < hihi(X2));
   if (m_has_layer && m_overlaps) {
      // Create the initial condition object and cache the initial condition
      // values.
      m_external_distribution_ic =
         new External2VIC(m_ext_2v_file, m_domain, a_dist_func.numGhosts());
      m_external_distribution_ic->cache(a_dist_func);

      // Allocate and compute m_nu.
      ParallelArray::Box base_space(CDIM);
      vector<int> num_cells(CDIM);
      for (int d = 0; d < CDIM; ++d) {
         base_space.lower(d) = ib.lower(d);
         base_space.upper(d) = ib.upper(d);
         num_cells[d] = m_domain.numberOfCells(X1);
      }
      m_nu.partition(base_space, CDIM, a_dist_func.numGhosts(), num_cells);
      for (int i2 = n2a; i2 <= n2b; ++i2) {
         double y = ymin+(0.5+i2)*dy;
         double ya, yb, eta, nuy;
         if (y <= ylolo) {
            eta = 0.0;
         }
         else if (y < yhilo) {
            ya = ylolo;
            yb = yhilo;
            eta = (y-ya)/(yb-ya);
         }
         else if (y <= ylohi) {
            eta = 1.0;
         }
         else if (y < yhihi) {
            ya = yhihi;
            yb = ylohi;
            eta = (y-ya)/(yb-ya);
         }
         else {
            eta = 0.0;
         }
         if (a_spatial_solution_order == 4) {
            nuy = -pow(eta, 4)*(
               +20.0*pow(eta, 3)
               -70.0*pow(eta, 2)
               +84.0*eta
               -35.0);
         }
         else {
            nuy = -pow(eta, 6)*(
               +252.0*pow(eta, 5)
               -1386.0*pow(eta, 4)
               +3080.0*pow(eta, 3)
               -3465.0*pow(eta, 2)
               +1980.0*eta
               -462.0);
         }
         for (int i1 = n1a; i1 <= n1b; ++i1) {
            double x = xmin+(0.5+i1)*dx;
            double xa, xb, xi, nux;
            if (x <= xlolo) {
               xi = 0.0;
            }
            else if (x < xhilo) {
               xa = xlolo;
               xb = xhilo;
               xi = (x-xa)/(xb-xa);
            }
            else if (x <= xlohi) {
               xi = 1.0;
            }
            else if (x < xhihi) {
               xa = xhihi;
               xb = xlohi;
               xi = (x-xa)/(xb-xa);
            }
            else {
               xi = 0.0;
            }
            if (a_spatial_solution_order == 4) {
               nux = -pow(xi, 4)*(
                  +20.0*pow(xi, 3)
                  -70.0*pow(xi, 2)
                  +84.0*xi
                  -35.0);
            }
            else {
               nux = -pow(xi, 6)*(
                  +252.0*pow(xi, 5)
                  -1386.0*pow(xi, 4)
                  +3080.0*pow(xi, 3)
                  -3465.0*pow(xi, 2)
                  +1980.0*xi
                  -462.0);
            }
            m_nu(i1, i2) = m_coefficient*nux*nuy;
         }
      }
   }
}


void
ExternalDistKrookLayer::parseParameters(
   LokiInputParser& a_pp)
{
   // Read the overall scaling factor for the layer.
   a_pp.query("coefficient", m_coefficient);

   // For each direction read upper and lower ends of the layer.
   for (int d(0); d < m_dim; ++d) {
      ostringstream name;

      name << "outer_x" << d+1 << "a";
      if (a_pp.query(name.str().c_str(), m_x_lo_lo[d])) {
         m_has_layer_lo[d] = true;
         m_has_layer = true;
      }
      name.str("");

      name << "outer_x" << d+1 << "b";
      if (a_pp.query(name.str().c_str(), m_x_hi_hi[d])) {
         m_has_layer_hi[d] = true;
         m_has_layer = true;
      }
      name.str("");

      name << "inner_x" << d+1 << "a";
      if (a_pp.query(name.str().c_str(), m_x_hi_lo[d])) {
         m_has_layer_lo[d] = true;
         m_has_layer = true;
      }
      name.str("");

      name << "inner_x" << d+1 << "b";
      if (a_pp.query(name.str().c_str(), m_x_lo_hi[d])) {
         m_has_layer_hi[d] = true;
         m_has_layer = true;
      }
      name.str("");
   }

   // Read the name of the file containing the external velocity distribution.
   if (m_has_layer &&
       !a_pp.query("external_distribution_file", m_ext_2v_file)) {
      LOKI_ABORT("Must supply name of external distribution file.");
   }
}


void
ExternalDistKrookLayer::checkParameters()
{
   if (m_coefficient<0.0) {
      LOKI_ABORT("ExternalDistKrookLayer must have non-negative coefficient");
   }

   for (int d(0); d < m_dim; ++d) {
      if (m_domain.lower(d) > m_x_lo_lo[d]) {
         ostringstream msg;
         msg << "ExternalDistKrookLayer must be inside physical domain in dimension" << d;
         LOKI_ABORT(msg.str());
      }

      if (m_x_hi_lo[d] < m_x_lo_lo[d]) {
         ostringstream msg;
         msg << "ExternalDistKrookLayer inner box must be inside outer box in dimension" << d;
         LOKI_ABORT(msg.str());
      }

      if (m_x_hi_lo[d] > m_x_lo_hi[d]) {
         ostringstream msg;
         msg << "ExternalDistKrookLayer inner box lower bound must be < inner box upper bound in dimension" << d;
         LOKI_ABORT(msg.str());
      }

      if (m_x_lo_hi[d] > m_x_hi_hi[d]) {
         ostringstream msg;
         msg << "ExternalDistKrookLayer inner box must be inside outer box in dimension" << d;
         LOKI_ABORT(msg.str());
      }

      if (m_domain.upper(d) < m_x_hi_hi[d]) {
         ostringstream msg;
         msg << "ExternalDistKrookLayer must be inside physical domain in dimension" << d;
         LOKI_ABORT(msg.str());
      }
   }
}


void ExternalDistKrookLayer::printParameters() const
{
   // Print the layer upper and lower ends along with the exponent and scaling.
   if (m_has_layer) {
      Loki_Utilities::printF("\n  ExternalDistKrookLayer parameters:\n");

      for (int d(0); d < m_dim; ++d) {
         if (m_has_layer_lo[d]) {
            Loki_Utilities::printF("    outer_x%da                 = %e\n",
               d+1,
               m_x_lo_lo[d]);
            Loki_Utilities::printF("    inner_x%da                 = %e\n",
               d+1,
               m_x_hi_lo[d]);
         }
         if (m_has_layer_hi[d]) {
            Loki_Utilities::printF("    inner_x%db                 = %e\n",
               d+1,
               m_x_lo_hi[d]);
            Loki_Utilities::printF("    outer_x%db                 = %e\n",
               d+1,
               m_x_hi_hi[d]);
         }
      }

      Loki_Utilities::printF("    coefficient                = %e\n",
         m_coefficient);
      Loki_Utilities::printF("    external distribution file = %s\n",
         m_ext_2v_file.c_str());
   }
}


void
ExternalDistKrookLayer::putToDatabase(
   RestartWriter& a_writer,
   bool a_write_data) const
{
   // Write all the info about the layer.
   a_writer.writeDoubleArray("x_lo_lo_external_dist_krook",
      &m_x_lo_lo[0],
      m_dim,
      a_write_data);
   a_writer.writeDoubleArray("x_hi_lo_external_dist_krook",
      &m_x_hi_lo[0],
      m_dim,
      a_write_data);
   a_writer.writeDoubleArray("x_lo_hi_external_dist_krook",
      &m_x_lo_hi[0],
      m_dim,
      a_write_data);
   a_writer.writeDoubleArray("x_hi_hi_external_dist_krook",
      &m_x_hi_hi[0],
      m_dim,
      a_write_data);
   a_writer.writeDoubleValue("externalDistKrookCoeff",
      m_coefficient,
      a_write_data);

   // These could be intuited in getFromDatabase based on the other stuff but
   // we're writing it all out just the same.
   a_writer.writeIntegerValue("externalDistKrookHasLayer",
      static_cast<int>(m_has_layer),
      a_write_data);
   a_writer.writeIntegerValue("externalDistKrookHasLayer_lo_0",
      static_cast<int>(m_has_layer_lo[X1]),
      a_write_data);
   a_writer.writeIntegerValue("externalDistKrookHasLayer_lo_1",
      static_cast<int>(m_has_layer_lo[X2]),
      a_write_data);
   a_writer.writeIntegerValue("externalDistKrookHasLayer_hi_0",
      static_cast<int>(m_has_layer_hi[X1]),
      a_write_data);
   a_writer.writeIntegerValue("externalDistKrookHasLayer_hi_1",
      static_cast<int>(m_has_layer_hi[X2]),
      a_write_data);
}


void
ExternalDistKrookLayer::getFromDatabase(
   RestartReader& a_reader)
{
   // Read all the info about the layer.
   a_reader.readDoubleArray("x_lo_lo_external_dist_krook",
      &m_x_lo_lo[0],
      m_dim);
   a_reader.readDoubleArray("x_hi_lo_external_dist_krook",
      &m_x_hi_lo[0],
      m_dim);
   a_reader.readDoubleArray("x_lo_hi_external_dist_krook",
      &m_x_lo_hi[0],
      m_dim);
   a_reader.readDoubleArray("x_hi_hi_external_dist_krook",
      &m_x_hi_hi[0],
      m_dim);
   a_reader.readDoubleValue("externalDistKrookCoeff", m_coefficient);

   // We actually have enough info at this point to figure the rest out but its
   // in the file so we'll just read it.
   int input_val;
   a_reader.readIntegerValue("externalDistKrookHasLayer", input_val);
   m_has_layer = static_cast<bool>(input_val);
   a_reader.readIntegerValue("externalDistKrookHasLayer_lo_0", input_val);
   m_has_layer_lo[X1] = static_cast<bool>(input_val);
   a_reader.readIntegerValue("externalDistKrookHasLayer_lo_1", input_val);
   m_has_layer_lo[X2] = static_cast<bool>(input_val);
   a_reader.readIntegerValue("externalDistKrookHasLayer_hi_0", input_val);
   m_has_layer_hi[X1] = static_cast<bool>(input_val);
   a_reader.readIntegerValue("externalDistKrookHasLayer_hi_1", input_val);
   m_has_layer_hi[X2] = static_cast<bool>(input_val);
}

} // end namespace Loki
