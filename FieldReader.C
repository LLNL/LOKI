/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "FieldReader.H"
#include "Loki_Defines.H"

namespace Loki {

FieldReader::FieldReader(
   const string& a_name)
   : ReaderWriterBase()
{
   // Open this field file and the root group.
   openFileAndRoot(a_name, m_file, m_root);
}


FieldReader::~FieldReader()
{
   // Close the file and root group.
   closeFileAndRoot(m_file, m_root);
}


void
FieldReader::readTime(
   const string& a_name,
   double& a_time)
{
   // Open the dataset for the time of this time slice.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
   hid_t dsid = H5Dopen(m_root, a_name.c_str());
#else
   hid_t dsid = H5Dopen(m_root, a_name.c_str(), H5P_DEFAULT);
#endif
#else
   hid_t dsid = H5Dopen(m_root, a_name.c_str());
#endif
   if (dsid < 0) {
      LOKI_ABORT("Can not open the dataset for the time of a time slice.");
   }

   // Read the time.
   herr_t errf = H5Dread(dsid,
      H5T_NATIVE_DOUBLE,
      H5S_ALL,
      H5S_ALL,
      H5P_DEFAULT,
      &a_time);
   if (errf < 0) {
      LOKI_ABORT("Can not read the time of a field time slice.");
   }

   // Close the dataset for the time.
   errf = H5Dclose(dsid);
   if (errf < 0) {
      LOKI_ABORT("Can not close dataset for the time of a field time slice.");
   }
}


void
FieldReader::readCoords(
   vector<double>& a_x_coords,
   vector<double>& a_y_coords)
{
   readDoubleArray("x", m_root, a_x_coords);
   readDoubleArray("y", m_root, a_y_coords);
}


void
FieldReader::readField(
   const string& a_field_name,
   vector<double>& a_field)
{
   readDoubleArray(a_field_name, m_root, a_field);
}


void
FieldReader::readTotalNumTimeSlices(
   int& a_total_num_time_slices)
{
   readIntegerValue("total_num_time_slices", m_root, a_total_num_time_slices);
}


void
FieldReader::readNumTimeSlicesInFile(
   int& a_num_time_slices_in_file)
{
   readIntegerValue("num_time_slices_in_this_file",
      m_root,
      a_num_time_slices_in_file);
}
} // end namespace Loki
