/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "External2VIC.H"
#include "Directions.H"
#include "getMomentumF.H"
#include "Simulation.H"
#include "Loki_Utilities.H"
#include "hdf5.h"

namespace Loki {

bool
External2VIC::isType(
   const string& a_name)
{
   if (a_name.compare("External 2V") == 0) {
      return true;
   }
   return false;
}


External2VIC::External2VIC(
   const string& a_file,
   const ProblemDomain& a_domain,
   int a_num_ghosts)
  : ICInterface(a_domain),
    m_nvx(a_domain.numberOfCells(2)+2*a_num_ghosts),
    m_num_ghosts(a_num_ghosts)
{
   // This initial condition is factorable by construction.
   m_factorable = true;

   // Set m_frac.  For this initial condition, m_frac is 1.
   m_frac = 1.0;

   // Open 2V initial condition distribution file.
   hid_t file_id = H5Fopen(a_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   if (file_id < 0) {
      LOKI_ABORT("Unable to open external 2V distribution file.");
   }

   // Read this velocity space.
   hid_t dataset_id;
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
   dataset_id = H5Dopen(file_id, "vspace");
#else
   dataset_id = H5Dopen(file_id, "vspace", H5P_DEFAULT);
#endif
#else
   dataset_id = H5Dopen(file_id, "vspace");
#endif
   if (dataset_id > 0) {
      // Check that the dataset is the right size.
      int data_size = m_nvx*(m_domain.numberOfCells(3) + 2*m_num_ghosts);
      hid_t dspace = H5Dget_space(dataset_id);
      if (dspace < 0) {
         LOKI_ABORT("Can not get dataspace for dataset of a vspace.");
      }
      hsize_t nsel = H5Sget_select_npoints(dspace);
      if (static_cast<int>(nsel) != data_size) {
         LOKI_ABORT("Distribution size does not match velocity space.");
      }

      // Now that everything checks out, read the dataset.
      m_vspace_data.resize(data_size);
      herr_t errf = H5Dread(dataset_id,
         H5T_NATIVE_DOUBLE,
         H5S_ALL,
         H5S_ALL,
         H5P_DEFAULT,
         &m_vspace_data[0]);
      if (errf < 0) {
         LOKI_ABORT("Can not read dataset for a vspace.");
      }

      // Close the dataspace.
      errf = H5Sclose(dspace);
      if (errf < 0) {
         LOKI_ABORT("Can not close dataspace for dataset of a vspace.");
      }

      // Close the dataset containing the velocity space.
      errf = H5Dclose(dataset_id);
      if (errf < 0) {
         LOKI_ABORT("Can not close dataset for a vspace.");
      }
   }
   else {
      LOKI_ABORT("Can not open dataset for a vspace.");
   }

   // Close the 2V initial condition distribution file.
   H5Fclose(file_id);
}


External2VIC::~External2VIC()
{
}


void
External2VIC::cache(
   const ParallelArray& a_u)
{
   // Build m_fx and m_fv.
   const ParallelArray::Box& db = a_u.dataBox();
   const ParallelArray::Box& ib = a_u.interiorBox();
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
   m_fx = 1.0;

   // Compute m_fv.
   for (int i4 = db.lower(V2); i4 <= db.upper(V2); ++i4) {
      for (int i3 = db.lower(V1); i3 <= db.upper(V1); ++i3) {
         int vspace_idx = (i3+m_num_ghosts)+m_nvx*(i4+m_num_ghosts);
         m_fv(i3, i4) = m_vspace_data[vspace_idx];
      }
   }
   m_vspace_data.clear();
}


void
External2VIC::set(
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
External2VIC::getIC_At_Pt(
   int a_i1,
   int a_i2,
   int a_i3,
   int a_i4) const
{
   return m_fv(a_i3, a_i4);
}


void
External2VIC::printParameters() const
{
}
} // end namespace Loki
