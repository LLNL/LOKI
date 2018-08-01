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
#include "ProblemDomain.H"

#include "Directions.H"

namespace Loki {

// Helper function to compute m_dx and m_box once we know m_dim, m_x_hi, m_x_lo,
// and m_n_cells.
void
ProblemDomain::define()
{
   for (int d(0); d < m_dim; ++d) {
      m_dx[d] = (m_x_hi[d] - m_x_lo[d]) / m_n_cells[d];
   }

   m_box = tbox::Box(tbox::IntVector::Zero(m_dim), m_n_cells - 1);
}


ProblemDomain::ProblemDomain(
   const tbox::Dimension& a_dim,
   HDF_DataBase& a_db)
   : m_dim( a_dim ),
     m_x_lo( a_dim ),
     m_x_hi( a_dim ),
     m_dx( a_dim ),
     m_n_cells( a_dim ),
     m_box( a_dim ),
     m_is_periodic( a_dim )
{
   // Read from the restart file.
   getFromDatabase( a_db );
}


ProblemDomain::ProblemDomain(
   const tbox::Dimension& a_dim,
   int a_spatial_solution_order,
   ParmParse& a_pp)
   : m_dim(a_dim),
     m_x_lo(a_dim),
     m_x_hi(a_dim),
     m_dx(a_dim),
     m_n_cells(a_dim),
     m_box(a_dim),
     m_is_periodic(a_dim)
{
   // Read the fundamental quantities from user input and compute derived
   // quantities.
   parseParameters(a_spatial_solution_order, a_pp);
   define();
}


ProblemDomain::ProblemDomain(
   const tbox::Dimension& a_dim,
   const tbox::IntVector& a_n_cells,
   const Array<double>&   a_x_lo,
   const Array<double>&   a_x_hi,
   const Array<bool>&     a_is_periodic)
   : m_dim(a_dim),
     m_x_lo(a_x_lo),
     m_x_hi(a_x_hi),
     m_dx(a_dim),
     m_n_cells(a_n_cells),
     m_box(a_dim),
     m_is_periodic(a_is_periodic)
{
   // Fundamental quantities supplied just need to compute derived quantities.
   define();
}


ProblemDomain::ProblemDomain(
   const ProblemDomain& a_other)
   : m_dim(a_other.m_dim),
     m_x_lo(a_other.m_x_lo),
     m_x_hi(a_other.m_x_hi),
     m_dx(a_other.m_dx),
     m_n_cells(a_other.m_n_cells),
     m_box(a_other.m_box),
     m_is_periodic(a_other.m_is_periodic)
{
}


ProblemDomain::~ProblemDomain()
{
}


void
ProblemDomain::parseParameters(
   int a_spatial_solution_order,
   ParmParse& a_pp)
{
   // Size the limits and put some bogus values into it.  We should probably
   // require the domain_limits input.
   Array<double> limits(2 * m_dim);
   for (int d(0); d < m_dim; ++d) {
      limits[2*d] = 0.0;
      limits[2*d+1] = 1.0;
   }
   a_pp.queryarr("domain_limits",
      limits,
      0,
      static_cast<int>(limits.length()));
   for (int d(0); d < m_dim; ++d) {
      m_x_lo[d] = limits[2*d];
      m_x_hi[d] = limits[2*d+1];
   }

   // Read the periodic_dir strings and convert to booleans.  This is optional.
   Array<aString> periodic_str(m_dim);
   for (int d(0); d < m_dim; ++d) {
      periodic_str[d] = aString("false");
   }
   a_pp.queryarr("periodic_dir",
      periodic_str,
      0,
      static_cast<int>(periodic_str.length()));
   for (int d(0); d < m_dim; ++d) {
      m_is_periodic[d] = periodic_str[d].matches("false") ? false : true;
   }

   // Read the number of cells in each direction.  Again, this should probably
   // be required.
   Array<int> tmp(m_dim);
   for (int d(0); d < m_dim; ++d) {
      tmp[d] = m_n_cells[d];
   }
   a_pp.queryarr("N", tmp, 0, static_cast<int>(tmp.length()));
   for (int d(0); d < m_dim; ++d) {
      m_n_cells[d] = tmp[d];
      // This is supposed to check that the configuration space dimensions
      // are large enough for Poisson and Maxwell.  It assumes that there is
      // a single Poisson or Maxwell processor.
      if (m_n_cells[d] < a_spatial_solution_order+1) {
         OV_ABORT("Number of cells must be at least stencil width");
      }
   }
}


void
ProblemDomain::printParameters() const
{
   // Write the basic info.
   printF("\n  ProblemDomain: Using the following parameters:\n");
   for (int d(0); d < m_dim; ++d) {
      int i(d+1);
      printF("    x%da = %e, x%db = %e, NX%d = %i\n",
         i, m_x_lo[d],
         i, m_x_hi[d],
         i, m_n_cells[d]);
   }
}


void
ProblemDomain::putToDatabase(
   HDF_DataBase& a_db) const
{
   // Write everything to the restart file including dx which is derivable.
   m_n_cells.putToDatabase(a_db, "N");
   a_db.put(m_x_lo.dataPtr(), "x_lo", static_cast<int>(m_x_lo.length()));
   a_db.put(m_x_hi.dataPtr(), "x_hi", static_cast<int>(m_x_hi.length()));
   a_db.put(m_dx.dataPtr(), "dx", static_cast<int>(m_dx.length()));
   a_db.put(static_cast<int>(m_is_periodic[X1]), "isPeriodic_0");
   a_db.put(static_cast<int>(m_is_periodic[X2]), "isPeriodic_1");
}


void
ProblemDomain::getFromDatabase(
   const HDF_DataBase& a_db)
{
   // Read everything to the restart file including dx which is derivable.
   m_n_cells.getFromDatabase(a_db, "N");
   a_db.get(m_x_lo.dataPtr(), "x_lo", static_cast<int>(m_x_lo.length()));
   a_db.get(m_x_hi.dataPtr(), "x_hi", static_cast<int>(m_x_hi.length()));
   a_db.get(m_dx.dataPtr(), "dx", static_cast<int>(m_dx.length()));
   int input_val;
   a_db.get(input_val, "isPeriodic_0");
   m_is_periodic[X1] = static_cast<bool>(input_val);
   a_db.get(input_val, "isPeriodic_1");
   m_is_periodic[X2] = static_cast<bool>(input_val);
}

} // end namespace Loki
