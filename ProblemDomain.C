/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "ProblemDomain.H"
#include "Directions.H"
#include "Loki_Defines.H"
#include "Loki_Utilities.H"

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
   RestartReader& a_reader)
   : m_dim(a_dim),
     m_x_lo(a_dim),
     m_x_hi(a_dim),
     m_dx(a_dim),
     m_n_cells(a_dim),
     m_box(a_dim),
     m_is_periodic(a_dim)
{
   // Read from the restart file.
   getFromDatabase(a_reader);
}


ProblemDomain::ProblemDomain(
   const tbox::Dimension& a_dim,
   int a_spatial_solution_order,
   LokiInputParser& a_pp)
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
   const vector<double>&  a_x_lo,
   const vector<double>&  a_x_hi,
   const deque<bool>&     a_is_periodic)
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
   LokiInputParser& a_pp)
{
   // Size the limits and put some bogus values into it.  We should probably
   // require the domain_limits input.
   vector<double> limits(2 * m_dim);
   for (int d(0); d < m_dim; ++d) {
      limits[2*d] = 0.0;
      limits[2*d+1] = 1.0;
   }
   a_pp.queryarr("domain_limits",
      limits,
      0,
      static_cast<int>(limits.size()));
   for (int d(0); d < m_dim; ++d) {
      if (limits[2*d] >= limits[2*d+1]) {
         LOKI_ABORT("Configuration space lower bound >= upper bound.");
      }
      m_x_lo[d] = limits[2*d];
      m_x_hi[d] = limits[2*d+1];
   }

   // Read the periodic_dir strings and convert to booleans.  This is optional.
   vector<string> periodic_str(m_dim);
   for (int d(0); d < m_dim; ++d) {
      periodic_str[d] = string("false");
   }
   a_pp.queryarr("periodic_dir",
      periodic_str,
      0,
      static_cast<int>(periodic_str.size()));
   for (int d(0); d < m_dim; ++d) {
      m_is_periodic[d] = periodic_str[d].compare("false") == 0 ? false : true;
   }

   // Read the number of cells in each direction.  Again, this should probably
   // be required.
   vector<int> tmp(m_dim);
   for (int d(0); d < m_dim; ++d) {
      tmp[d] = m_n_cells[d];
   }
   a_pp.queryarr("N", tmp, 0, static_cast<int>(tmp.size()));
   for (int d(0); d < m_dim; ++d) {
      m_n_cells[d] = tmp[d];
      // This is supposed to check that the configuration space dimensions
      // are large enough for Poisson and Maxwell.  It assumes that there is
      // a single Poisson or Maxwell processor.
      if (m_n_cells[d] < a_spatial_solution_order+1) {
         LOKI_ABORT("Number of cells must be at least stencil width");
      }
   }
}


void
ProblemDomain::printParameters() const
{
   // Write the basic info.
   Loki_Utilities::printF("\n  ProblemDomain: Using the following parameters:\n");
   for (int d(0); d < m_dim; ++d) {
      int i(d+1);
      Loki_Utilities::printF("    x%da = %e, x%db = %e, NX%d = %i\n",
         i, m_x_lo[d],
         i, m_x_hi[d],
         i, m_n_cells[d]);
   }
}


void
ProblemDomain::putToDatabase(
   RestartWriter& a_writer,
   bool a_write_data) const
{
   // Write everything to the restart file including dx which is derivable.
   m_n_cells.putToDatabase(a_writer, "N", a_write_data);
   a_writer.writeDoubleArray("x_lo",
      &m_x_lo[0],
      static_cast<int>(m_x_lo.size()),
      a_write_data);
   a_writer.writeDoubleArray("x_hi",
      &m_x_hi[0],
      static_cast<int>(m_x_hi.size()),
      a_write_data);
   a_writer.writeDoubleArray("dx",
      &m_dx[0],
      static_cast<int>(m_dx.size()),
      a_write_data);
   a_writer.writeIntegerValue("isPeriodic_0",
      static_cast<int>(m_is_periodic[X1]),
      a_write_data);
   a_writer.writeIntegerValue("isPeriodic_1",
      static_cast<int>(m_is_periodic[X2]),
      a_write_data);
}


void
ProblemDomain::getFromDatabase(
   RestartReader& a_reader)
{
   // Read everything from the restart file including dx which is derivable.
   m_n_cells.getFromDatabase(a_reader, "N");
   a_reader.readDoubleArray("x_lo",
      &m_x_lo[0],
      static_cast<int>(m_x_lo.size()));
   a_reader.readDoubleArray("x_hi",
      &m_x_hi[0],
      static_cast<int>(m_x_hi.size()));
   a_reader.readDoubleArray("dx",
      &m_dx[0],
      static_cast<int>(m_dx.size()));
   int input_val;
   a_reader.readIntegerValue("isPeriodic_0", input_val);
   m_is_periodic[X1] = static_cast<bool>(input_val);
   a_reader.readIntegerValue("isPeriodic_1", input_val);
   m_is_periodic[X2] = static_cast<bool>(input_val);
}

} // end namespace Loki
