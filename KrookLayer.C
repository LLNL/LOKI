/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "KrookLayer.H"
#include "Directions.H"
#include "Loki_Defines.H"
#include "Loki_Utilities.H"

namespace Loki {

KrookLayer::KrookLayer(
   const tbox::Dimension& a_dim,
   LokiInputParser& a_pp,
   const ProblemDomain& a_domain)
   : m_dim(a_dim),
     m_has_layer(false),
     m_overlaps(false),
     m_power(3),
     m_coefficient(1.0),
     m_domain(a_domain)
{
   // Initialize upper and lower ends of Krook layer to the full domain and
   // that there is no layer in any direction.
   m_x_lo = m_domain.lower();
   m_x_hi = m_domain.upper();

   m_has_layer_lo.resize(m_dim);
   m_has_layer_hi.resize(m_dim);
   for (int d(0); d < m_dim; ++d) {
      m_has_layer_lo[d] = false;
      m_has_layer_hi[d] = false;
   }

   // Read the Krook layer specified by the user and make sure the user is sane.
   parseParameters(a_pp);
   checkParameters();
}


KrookLayer::~KrookLayer()
{
}


void
KrookLayer::initialize(
   const ParallelArray& a_dist_func,
   int a_spatial_solution_order)
{
   // Determine if this domain overlaps the layer.  If it does, then compute the
   // constant layer coefficients, m_nu.
   double xmin = m_domain.lower(X1);
   double xmax = m_domain.upper(X1);
   double ymin = m_domain.lower(X2);
   double ymax = m_domain.upper(X2);
   double dx = m_domain.dx(X1);
   double dy = m_domain.dx(X2);
   double xlo = lower(X1);
   double xhi = upper(X1);
   double ylo = lower(X2);
   double yhi = upper(X2);
   const ParallelArray::Box& ib = a_dist_func.interiorBox();
   int n1a = ib.lower(X1);
   int n1b = ib.upper(X1);
   int n2a = ib.lower(X2);
   int n2b = ib.upper(X2);
   m_overlaps = xmin+(0.5+n1a)*dx < xlo ||
      xmin+(0.5+n1b)*dx > xhi ||
      ymin+(0.5+n2a)*dy < ylo ||
      ymin+(0.5+n2b)*dy > yhi;
   if (m_has_layer && m_overlaps) {
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
         if (y < ylo) {
            ya = ylo;
            yb = ymin;
            eta = (y-ya)/(yb-ya);
//            nuy = pow(eta, m_power);
         }
         else if (y > yhi) {
            ya = yhi;
            yb = ymax;
            eta = (y-ya)/(yb-ya);
//            nuy = pow(eta, m_power);
         }
         else {
//            nuy = 0.0;
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
            if (x < xlo) {
               xa = xlo;
               xb = xmin;
               xi = (x-xa)/(xb-xa);
//               nux = pow(xi, m_power);
            }
            else if (x > xhi) {
               xa = xhi;
               xb = xmax;
               xi = (x-xa)/(xb-xa);
//               nuy = pow(xi, m_power);
            }
            else {
//               nuy = 0.0;
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
            m_nu(i1, i2) = m_coefficient*((1.0-nuy)*nux+nuy);
         }
      }
   }
}


void
KrookLayer::parseParameters(
   LokiInputParser& a_pp)
{
   // Read the exponent and overall scaling factor for the layer.
   a_pp.query("power", m_power);
   a_pp.query("coefficient", m_coefficient);

   // For each direction read upper and lower ends of the layer.
   for (int d(0); d < m_dim; ++d) {
      ostringstream name;

      name << "x" << d+1 << "a";
      if (a_pp.query(name.str().c_str(), m_x_lo[d])) {
         m_has_layer_lo[d] = true;
         m_has_layer = true;
      }
      name.str("");

      name << "x" << d+1 << "b";
      if (a_pp.query(name.str().c_str(), m_x_hi[d])) {
         m_has_layer_hi[d] = true;
         m_has_layer = true;
      }
      name.str("");
   }
}


void
KrookLayer::checkParameters()
{
   if (m_power<=0.0) {
      LOKI_ABORT("Krook layer must have positive power");
   }

   if (m_coefficient<0.0) {
      LOKI_ABORT("Krook layer must have non-negative coefficient");
   }

   for (int d(0); d < m_dim; ++d) {
      if (m_domain.lower(d) > m_x_lo[d]) {
         ostringstream msg;
         msg << "Krook layer must be inside physical domain (x" << d+1 << "a)";
         LOKI_ABORT(msg.str());
      }

      if (m_domain.upper(d) < m_x_hi[d]) {
         ostringstream msg;
         msg << "Krook layer must be inside physical domain (x" << d+1 << "b)";
         LOKI_ABORT(msg.str());
      }

      if (m_x_lo[d] >= m_x_hi[d]) {
         ostringstream msg;
         msg << "Krook layer lower bound (x" << d+1
             << "a must be < upper bound (x" << d+1 << "b)";
         LOKI_ABORT(msg.str());
      }
   }
}


void KrookLayer::printParameters() const
{
   // Print the layer upper and lower ends along with the exponent and scaling.
   if (m_has_layer) {
      Loki_Utilities::printF("\n  Krook layer parameters:\n");

      for (int d(0); d < m_dim; ++d) {
         if (m_has_layer_lo[d]) {
            Loki_Utilities::printF("    x%da krook  = %e\n", d+1, m_x_lo[d]);
         }
         if (m_has_layer_hi[d]) {
            Loki_Utilities::printF("    x%db krook  = %e\n", d+1, m_x_hi[d]);
         }
      }

      Loki_Utilities::printF("    krook power = %e, krook coefficient = %e\n",
         m_power,
         m_coefficient);
   }
}


void
KrookLayer::putToDatabase(
   RestartWriter& a_writer,
   bool a_write_data) const
{
   // Write all the info about the layer.
   a_writer.writeDoubleArray("x_lo_krook",
      &m_x_lo[0],
      static_cast<int>(m_x_lo.size()),
      a_write_data);
   a_writer.writeDoubleArray("x_hi_krook",
      &m_x_hi[0],
      static_cast<int>(m_x_hi.size()),
      a_write_data);
   a_writer.writeDoubleValue("krookPower", m_power, a_write_data);
   a_writer.writeDoubleValue("krookCoeff", m_coefficient, a_write_data);

   // These could be intuited in getFromDatabase based on the other stuff but
   // we're writing it all out just the same.
   a_writer.writeIntegerValue("krookHasLayer",
      static_cast<int>(m_has_layer),
      a_write_data);
   a_writer.writeIntegerValue("krookHasLayer_lo_0",
      static_cast<int>(m_has_layer_lo[X1]),
      a_write_data);
   a_writer.writeIntegerValue("krookHasLayer_lo_1",
      static_cast<int>(m_has_layer_lo[X2]),
      a_write_data);
   a_writer.writeIntegerValue("krookHasLayer_hi_0",
      static_cast<int>(m_has_layer_hi[X1]),
      a_write_data);
   a_writer.writeIntegerValue("krookHasLayer_hi_1",
      static_cast<int>(m_has_layer_hi[X2]),
      a_write_data);
}


void
KrookLayer::getFromDatabase(
   RestartReader& a_reader)
{
   // Read all the info about the layer.
   a_reader.readDoubleArray("x_lo_krook",
      &m_x_lo[0],
      static_cast<int>(m_x_lo.size()));
   a_reader.readDoubleArray("x_hi_krook",
      &m_x_hi[0],
      static_cast<int>(m_x_hi.size()));
   a_reader.readDoubleValue("krookPower", m_power);
   a_reader.readDoubleValue("krookCoeff", m_coefficient);

   // We actually have enough info at this point to figure the rest out but its
   // in the file so we'll just read it.
   int input_val;
   a_reader.readIntegerValue("krookHasLayer", input_val);
   m_has_layer = static_cast<bool>(input_val);
   a_reader.readIntegerValue("krookHasLayer_lo_0", input_val);
   m_has_layer_lo[X1] = static_cast<bool>(input_val);
   a_reader.readIntegerValue("krookHasLayer_lo_1", input_val);
   m_has_layer_lo[X2] = static_cast<bool>(input_val);
   a_reader.readIntegerValue("krookHasLayer_hi_0", input_val);
   m_has_layer_hi[X1] = static_cast<bool>(input_val);
   a_reader.readIntegerValue("krookHasLayer_hi_1", input_val);
   m_has_layer_hi[X2] = static_cast<bool>(input_val);
}

} // end namespace Loki
