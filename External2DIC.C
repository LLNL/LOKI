/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "External2DIC.H"
#include "JuttnerThermal.H"
#include "MaxwellianThermal.H"
#include "Directions.H"
#include "getMomentumF.H"
#include "Simulation.H"
#include "Loki_Utilities.H"
#include "hdf5.h"

namespace Loki {

bool
External2DIC::isType(
   const string& a_name)
{
   if (a_name.compare("External 2D") == 0) {
      return true;
   }
   return false;
}


External2DIC::External2DIC(
   LokiInputParser& a_pp,
   const ProblemDomain& a_domain,
   double a_mass,
   int a_num_ghosts,
   const tbox::Box& a_global_box)
  : ICInterface(a_domain),
    m_num_ghosts(a_num_ghosts)
{
   // Set some hopefully sane input parameter defaults.
   m_tx = 1.0;
   m_ty = 1.0;
   m_vx0 = 0.0;
   m_vy0 = 0.0;
   m_x_wave_number = 0.0;
   m_y_wave_number = 0.0;
   m_flow_vel_phi = 0.0;
   m_frac = 1.0;

   // Now read what the user really wants.
   parseParameters(a_pp);

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

   // Open 2D initial condition distribution file.
   hid_t file_id = H5Fopen(m_ext_2d_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   if (file_id < 0) {
      LOKI_ABORT("Unable to open external 2D distribution file.");
   }

   // Open the dataset containing the 2D distribution.
   hid_t dataset_id;
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
   dataset_id = H5Dopen(file_id, "2D_dist");
#else
   dataset_id = H5Dopen(file_id, "2D dist", H5P_DEFAULT);
#endif
#else
   dataset_id = H5Dopen(file_id, "2D_dist");
#endif
   if (dataset_id > 0) {
      // Check that the dataset is the right size.
      m_nx = a_global_box.numberCells(0);
      int data_size = m_nx*a_global_box.numberCells(1);
      hid_t dspace = H5Dget_space(dataset_id);
      if (dspace < 0) {
         LOKI_ABORT("Can not get dataspace for dataset \"2D dist\".");
      }
      hsize_t nsel = H5Sget_select_npoints(dspace);
      if (static_cast<int>(nsel) != data_size) {
         LOKI_ABORT("Distribution size does not match configuration space.");
      }

      // Now that everything checks out, read the dataset.
      m_data.resize(data_size);
      herr_t errf = H5Dread(dataset_id,
                            H5T_NATIVE_DOUBLE,
                            H5S_ALL,
                            H5S_ALL,
                            H5P_DEFAULT,
                            &m_data[0]);
      if (errf < 0) {
         LOKI_ABORT("Can not read dataset \"2D dist\".");
      }

      // Close the dataspace.
      errf = H5Sclose(dspace);
      if (errf < 0) {
         LOKI_ABORT("Can not close dataspace for dataset \"2D dist\".");
      }

      // Close the dataset containing the 2D distribution.
      errf = H5Dclose(dataset_id);
      if (errf < 0) {
         LOKI_ABORT("Can not close dataset \"2D dist\".");
      }
   }
   else {
      LOKI_ABORT("Can not open dataset \"2D dist\".");
   }
}


External2DIC::~External2DIC()
{
}


void
External2DIC::cache(
   const ParallelArray& a_u)
{
   const ParallelArray::Box& db = a_u.dataBox();
   const ParallelArray::Box& ib = a_u.interiorBox();

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
         int idx2 = i2;
         // Enforce periodicity.
         if (m_domain.isPeriodic(X2)) {
            if (i2 < m_domain.box().lower(X2)) {
               idx2 += m_domain.numberOfCells(X2);
            }
            else if (i2 > m_domain.box().upper(X2)) {
               idx2 -= m_domain.numberOfCells(X2);
            }
         }
         for (int i1 = db.lower(X1); i1 <= db.upper(X1); ++i1) {
            int idx1 = i1;
            // Enforce periodicity.
            if (m_domain.isPeriodic(X1)) {
               if (i1 < m_domain.box().lower(X1)) {
                  idx1 += m_domain.numberOfCells(X1);
               }
               else if (i1 > m_domain.box().upper(X1)) {
                  idx1 -= m_domain.numberOfCells(X1);
               }
            }
            // Compute m_fx at this point in configuration space.
            m_fx(i1, i2) = m_data[(idx1+m_num_ghosts)+m_nx*(idx2+m_num_ghosts)];
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
               int idx2 = i2;
               // Enforce periodicity.
               if (m_domain.isPeriodic(X2)) {
                  if (i2 < m_domain.box().lower(X2)) {
                     idx2 += m_domain.numberOfCells(X2);
                  }
                  else if (i2 > m_domain.box().upper(X2)) {
                     idx2 -= m_domain.numberOfCells(X2);
                  }
               }
               for (int i1 = db.lower(X1); i1 <= db.upper(X1); ++i1) {
                  int idx1 = i1;
                  // Enforce periodicity.
                  if (m_domain.isPeriodic(X1)) {
                     if (i1 < m_domain.box().lower(X1)) {
                        idx1 += m_domain.numberOfCells(X1);
                     }
                     else if (i1 > m_domain.box().upper(X1)) {
                        idx1 -= m_domain.numberOfCells(X1);
                     }
                  }

                  // Compute m_f at this point in phase space.
                  double fx =
                     m_data[(idx1+m_num_ghosts)+m_nx*(idx2+m_num_ghosts)];
                  double spatial_factor =
                     cos(m_x_wave_number*idx1+m_y_wave_number*idx2+m_flow_vel_phi);
                  double fv = m_thermal->thermalFactor(spatial_factor, x3, x4);
                  m_f(i1, i2, i3, i4) = fv*fx;
	       }
            }
         }
      }
   }
   m_data.clear();
}


void
External2DIC::set(
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
External2DIC::getIC_At_Pt(
   int a_i1,
   int a_i2,
   int a_i3,
   int a_i4) const
{
   // To avoid order of operation differences with earlier versions of this
   // initial condition that did not cache their result the cached factored
   // results are multiplied in this specific order.  Do not change.
   if (m_factorable) {
      return m_frac*m_thermal->fnorm()*m_fv(a_i3, a_i4)*m_fx(a_i1, a_i2);
   }
   else {
      return m_frac*m_thermal->fnorm()*m_f(a_i1, a_i2, a_i3, a_i4);
   }
}


void
External2DIC::printParameters() const
{
   // Print all the user input parameters.
   Loki_Utilities::printF("\n  Using external 2D distribution:\n");
   Loki_Utilities::printF("    external 2D file        = %s\n",
      m_ext_2d_file.c_str());
   if (m_maxwellian_thermal) {
      Loki_Utilities::printF("    Maxwellian thermal\n");
   }
   else {
      Loki_Utilities::printF("    Juttner thermal\n");
   }
   Loki_Utilities::printF("    tx                      = %e\n", m_tx);
   Loki_Utilities::printF("    ty                      = %e\n", m_ty);
   Loki_Utilities::printF("    vx0                     = %e\n", m_vx0);
   Loki_Utilities::printF("    vy0                     = %e\n", m_vy0);
   Loki_Utilities::printF("    vflowinitx              = %e\n", m_vflowinitx);
   Loki_Utilities::printF("    vflowinity              = %e\n", m_vflowinity);
   Loki_Utilities::printF("    x wave number           = %e\n",
      m_x_wave_number);
   Loki_Utilities::printF("    y wave number           = %e\n",
      m_y_wave_number);
   Loki_Utilities::printF("    flow velocity phase     = %e\n", m_flow_vel_phi);
   Loki_Utilities::printF("    species relative weight = %e\n", m_frac);
}


void
External2DIC::parseParameters(
   LokiInputParser& a_pp)
{
   // The file name is required.
   string tmp;
   if (!a_pp.query("file_name", tmp)) {
      LOKI_ABORT("Must supply name of external 2D distribution file.");
   }
   m_ext_2d_file = tmp;

   // Determine if the user wants Maxwellian or Juttner thermal behavior.
   string test_str = "true";
   a_pp.query("maxwellian_thermal", test_str);
   m_maxwellian_thermal = test_str.compare("true") == 0 ? true : false;

   // Currently none of these inputs is required so the defaults better be
   // sane.
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
   a_pp.query("vflowinitx",    m_vflowinitx);
   a_pp.query("vflowinity",    m_vflowinity);
   a_pp.query("x_wave_number", m_x_wave_number);
   a_pp.query("y_wave_number", m_y_wave_number);
   a_pp.query("phase",         m_flow_vel_phi);
   a_pp.query("frac",          m_frac);
}
} // end namespace Loki
