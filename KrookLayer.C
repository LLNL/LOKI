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
#include "KrookLayer.H"

#include "Directions.H"

namespace Loki {

KrookLayer::KrookLayer(
   const tbox::Dimension& a_dim,
   ParmParse& a_pp,
   const ProblemDomain& a_domain)
   : m_dim(a_dim),
     m_has_layer(false),
     m_power(3),
     m_coefficient(1.0)
{
   // Initialize upper and lower ends of Krook layer to the full domain and
   // that there is no layer in any direction.
   m_x_lo = a_domain.lower();
   m_x_hi = a_domain.upper();

   m_has_layer_lo.resize(a_dim);
   m_has_layer_hi.resize(a_dim);
   for (int d(0); d < m_dim; ++d) {
      m_has_layer_lo[d] = false;
      m_has_layer_hi[d] = false;
   }

   // Read the Krook layer specified by the user and make sure the user is sane.
   parseParameters(a_pp);
   checkParameters(a_domain);
}


KrookLayer::~KrookLayer()
{
}


void
KrookLayer::parseParameters(
   ParmParse& a_pp)
{
   // Read the exponent and overall scaling factor for the layer.
   a_pp.query("power", m_power);
   a_pp.query("coefficient", m_coefficient);

   // For each direction read upper and lower ends of the layer.
   for (int d(0); d < m_dim; ++d) {
      char name[5];

      sprintf(name, "x%da", d+1);
      if (a_pp.query(name, m_x_lo[d])) {
         m_has_layer_lo[d] = true;
         m_has_layer = true;
      }

      sprintf(name, "x%db", d+1);
      if (a_pp.query(name, m_x_hi[d])) {
         m_has_layer_hi[d] = true;
         m_has_layer = true;
      }
   }
}


void
KrookLayer::checkParameters(
   const ProblemDomain& a_domain)
{
   if (m_power<=0.0) {
      OV_ABORT("Krook layer must have positive power");
   }

   if (m_coefficient<=0.0) {
      OV_ABORT("Krook layer must have positive coefficient");
   }

   for (int d(0); d < m_dim; ++d) {
      if (a_domain.lower(d) > m_x_lo[d]) {
         char msg[80];
         sprintf(msg,
            "Krook layer must be inside physical domain (x%da)",
            d+1);
         OV_ABORT(msg);
      }

      if (a_domain.upper(d) < m_x_hi[d]) {
         char msg[80];
         sprintf(msg,
            "Krook layer must be inside physical domain (x%db)",
            d+1);
         OV_ABORT(msg);
      }

      if ((m_has_layer_lo[d] || m_has_layer_hi[d]) && a_domain.isPeriodic(d)) {
         char msg[80];
         sprintf(msg,
            "Cannot specify Krook layer in periodic direction (d=%d)",
            d);
         OV_ABORT(msg);
      }
   }
}


void KrookLayer::printParameters() const
{
   // Print the layer upper and lower ends along with the exponent and scaling.
   if (m_has_layer) {
      printF("\n  Krook layer parameters:\n");

      for (int d(0); d < m_dim; ++d) {
         if (m_has_layer_lo[d]) {
            printF("    x%da krook = %e\n", d+1, m_x_lo[d]);
         }
         if (m_has_layer_hi[d]) {
            printF("    x%db krook = %e\n", d+1, m_x_hi[d]);
         }
      }

      printF("    krook power = %e, krook coefficient = %e\n",
         m_power,
         m_coefficient);
   }
}


void
KrookLayer::putToDatabase(
   HDF_DataBase& a_db) const
{
   // Write all the info about the layer.
   a_db.put(m_x_lo.dataPtr(), "x_lo_krook", static_cast<int>(m_x_lo.length()));
   a_db.put(m_x_hi.dataPtr(), "x_hi_krook", static_cast<int>(m_x_hi.length()));
   a_db.put(m_power, "krookPower");
   a_db.put(m_coefficient, "krookCoeff");

   // These could be intuited in getFromDatabase based on the other stuff but
   // we're writting it all out just the same.
   a_db.put(static_cast<int>(m_has_layer), "krookHasLayer");
   a_db.put(static_cast<int>(m_has_layer_lo[X1]), "krookHasLayer_lo_0");
   a_db.put(static_cast<int>(m_has_layer_lo[X2]), "krookHasLayer_lo_1");
   a_db.put(static_cast<int>(m_has_layer_hi[X1]), "krookHasLayer_hi_0");
   a_db.put(static_cast<int>(m_has_layer_hi[X2]), "krookHasLayer_hi_1");
}


void
KrookLayer::getFromDatabase(
   const HDF_DataBase& a_db)
{
   // Read all the info about the layer.
   a_db.get(m_x_lo.dataPtr(), "x_lo_krook", static_cast<int>(m_x_lo.length()));
   a_db.get(m_x_hi.dataPtr(), "x_hi_krook", static_cast<int>(m_x_hi.length()));
   a_db.get(m_power, "krookPower");
   a_db.get(m_coefficient, "krookCoeff");

   // We actually have enough info at this point to figure the rest out but its
   // in the file so we'll just read it.
   int input_val;
   a_db.get(input_val, "krookHasLayer");
   m_has_layer = static_cast<bool>(input_val);
   a_db.get(input_val, "krookHasLayer_lo_0");
   m_has_layer_lo[X1] = static_cast<bool>(input_val);
   a_db.get(input_val, "krookHasLayer_lo_1");
   m_has_layer_lo[X2] = static_cast<bool>(input_val);
   a_db.get(input_val, "krookHasLayer_hi_0");
   m_has_layer_hi[X1] = static_cast<bool>(input_val);
   a_db.get(input_val, "krookHasLayer_hi_1");
   m_has_layer_hi[X2] = static_cast<bool>(input_val);
}

} // end namespace Loki
