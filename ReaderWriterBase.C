/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "ReaderWriterBase.H"
#include "Loki_Defines.H"
#include <sstream>

namespace Loki {

ReaderWriterBase::ReaderWriterBase()
{
}


ReaderWriterBase::~ReaderWriterBase()
{
}


void
ReaderWriterBase::createFileAndRoot(
   const string& a_file_name,
   hid_t& a_file,
   hid_t& a_root,
   bool a_parallel_io)
{
   // Create the file.
   hid_t plistID;
   if (a_parallel_io) {
      plistID = H5Pcreate(H5P_FILE_ACCESS);
      herr_t errf = H5Pset_fapl_mpio(plistID, MPI_COMM_WORLD, MPI_INFO_NULL);
      if (errf < 0) {
         LOKI_ABORT("Unable to set file access property list.");
      }
   }
   else {
      plistID = H5P_DEFAULT;
   }
   a_file = H5Fcreate(a_file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plistID);
   H5Pclose(plistID);
   if (a_file < 0) {
      ostringstream os;
      os << "Unable to create file " << a_file_name << ".";
      LOKI_ABORT(os.str().c_str());
   }

   // Create the root group.
   createGroup("root", a_root, a_file);
}


void
ReaderWriterBase::openFileAndRoot(
   const string& a_file_name,
   hid_t& a_file,
   hid_t& a_root)
{
   // Open the file.
   a_file = H5Fopen(a_file_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
   if (a_file < 0) {
      ostringstream os;
      os << "Unable to open file " << a_file_name << ".";
      LOKI_ABORT(os.str().c_str());
   }

   // Open the root group.
   openGroup("root", a_root, a_file);
}


void
ReaderWriterBase::closeFileAndRoot(
   hid_t& a_file,
   hid_t& a_root)
{
   // Close the group.
   closeGroup(a_root);

   // Close the file.
   herr_t errf = H5Fclose(a_file);
   if (errf < 0) {
      LOKI_ABORT("Unable to close file.");
   }
}


void
ReaderWriterBase::createGroup(
   const string& a_group_name,
   hid_t& a_group,
   hid_t& a_home)
{
   // Create the group.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Gcreate_vers) && H5Gcreate_vers == 1
   a_group = H5Gcreate(a_home, a_group_name.c_str(), 0);
#else
   a_group = H5Gcreate(a_home,
      a_group_name.c_str(),
      H5P_DEFAULT,
      H5P_DEFAULT,
      H5P_DEFAULT);
#endif
#else
   a_group = H5Gcreate(a_home, a_group_name.c_str(), 0);
#endif
   if (a_group < 0) {
      LOKI_ABORT("Can not create group.");
   }

   // Give the group the "directory" attribute.
   hsize_t dims[1]={10};
   hsize_t rank=1;
   hid_t dataspaceID = H5Screate_simple(rank, dims, NULL);
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Acreate_vers) && H5Acreate_vers == 1
   hid_t classNameAttribID = H5Acreate(a_group,
      "className",
      H5T_NATIVE_UCHAR,
      dataspaceID,
      H5P_DEFAULT);
#else
   hid_t classNameAttribID = H5Acreate(a_group,
      "className",
      H5T_NATIVE_UCHAR,
      dataspaceID,
      H5P_DEFAULT,
      H5P_DEFAULT);
#endif
#else
   hid_t classNameAttribID = H5Acreate(a_group,
      "className",
      H5T_NATIVE_UCHAR,
      dataspaceID,
      H5P_DEFAULT);
#endif
   H5Sclose(dataspaceID);
   int istat = H5Awrite(classNameAttribID, H5T_NATIVE_UCHAR, "directory");
   if( istat<0 ) {
      LOKI_ABORT("Unable to write group attribute");
   }
   H5Aclose(classNameAttribID);
}


void
ReaderWriterBase::openGroup(
   const string& a_group_name,
   hid_t& a_group,
   hid_t& a_home)
{
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Gopen_vers) && H5Gopen_vers == 1
   a_group = H5Gopen(a_home, a_group_name.c_str());
#else
   a_group = H5Gopen(a_home, a_group_name.c_str(), H5P_DEFAULT);
#endif
#else
   a_group = H5Gopen(a_home, a_group_name.c_str());
#endif
   if (a_group < 0) {
      LOKI_ABORT("Can not open group.");
   }
}


void
ReaderWriterBase::closeGroup(
   hid_t& a_group)
{
   herr_t errf = H5Gclose(a_group);
   if (errf < 0) {
      LOKI_ABORT("Unable to close group.");
   }
}


void
ReaderWriterBase::writeDoubleArray(
   const string& a_name,
   hid_t a_home,
   hid_t a_dspace,
   const double* a_vals)
{
   // Create the dataset for the array of doubles.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dcreate_vers) && H5Dcreate_vers == 1
   hid_t dset = H5Dcreate(a_home,
      a_name.c_str(),
      H5T_NATIVE_DOUBLE,
      a_dspace,
      H5P_DEFAULT);
#else
   hid_t dset = H5Dcreate(a_home,
      a_name.c_str(),
      H5T_NATIVE_DOUBLE,
      a_dspace,
      H5P_DEFAULT,
      H5P_DEFAULT,
      H5P_DEFAULT);
#endif
#else
   hid_t dset = H5Dcreate(a_home,
      a_name.c_str(),
      H5T_NATIVE_DOUBLE,
      a_dspace,
      H5P_DEFAULT);
#endif
   if (dset < 0) {
      LOKI_ABORT("Unable to create dataset for array of doubles.");
   }

   // Write the array of doubles.
   herr_t errf = H5Dwrite(dset,
      H5T_NATIVE_DOUBLE,
      H5S_ALL,
      H5S_ALL,
      H5P_DEFAULT,
      a_vals);
   if (errf < 0) {
      LOKI_ABORT("Unable to write to dataset for array of doubles.");
   }

   // Close the dataset for array of doubles.
   errf = H5Dclose(dset);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataset for array of doubles.");
   }
}


void
ReaderWriterBase::writeDoubleValue(
   const string& a_name,
   hid_t a_home,
   double a_val)
{
   // Create the dataspace for a single value.
   hsize_t dims = 1;
   hid_t dspace = H5Screate_simple(1, &dims, NULL);
   if (dspace < 0) {
      LOKI_ABORT("Unable to create dataspace for a single double.");
   }

   // Create the dataset for the double value.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dcreate_vers) && H5Dcreate_vers == 1
   hid_t dset = H5Dcreate(a_home,
      a_name.c_str(),
      H5T_NATIVE_DOUBLE,
      dspace,
      H5P_DEFAULT);
#else
   hid_t dset = H5Dcreate(a_home,
      a_name.c_str(),
      H5T_NATIVE_DOUBLE,
      dspace,
      H5P_DEFAULT,
      H5P_DEFAULT,
      H5P_DEFAULT);
#endif
#else
   hid_t dset =
      H5Dcreate(a_home, a_name.c_str(), H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT);
#endif
   if (dset < 0) {
      LOKI_ABORT("Unable to create dataset for a single double.");
   }

   // Write the double value.
   herr_t errf =
      H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &a_val);
   if (errf < 0) {
      LOKI_ABORT("Unable to write to dataset for a single double.");
   }

   // Close the dataset and dataspace for a single double.
   errf = H5Dclose(dset);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataset for a single double.");
   }
   errf = H5Sclose(dspace);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataspace for a single double.");
   }
}


void
ReaderWriterBase::writeIntegerArray(
   const string& a_name,
   hid_t a_home,
   hid_t a_dspace,
   const int* a_vals)
{
   // Create the dataset for the array of integers.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dcreate_vers) && H5Dcreate_vers == 1
   hid_t dset = H5Dcreate(a_home,
      a_name.c_str(),
      H5T_STD_I32BE,
      a_dspace,
      H5P_DEFAULT);
#else
   hid_t dset = H5Dcreate(a_home,
      a_name.c_str(),
      H5T_STD_I32BE,
      a_dspace,
      H5P_DEFAULT,
      H5P_DEFAULT,
      H5P_DEFAULT);
#endif
#else
   hid_t dset = H5Dcreate(a_home,
      a_name.c_str(),
      H5T_STD_I32BE,
      a_dspace,
      H5P_DEFAULT);
#endif
   if (dset < 0) {
      LOKI_ABORT("Unable to create dataset for array of integers.");
   }

   // Write the array of integers.
   herr_t errf =
      H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, a_vals);
   if (errf < 0) {
      LOKI_ABORT("Unable to write to dataset for array of integers.");
   }

   // Close the dataset for array of integers.
   errf = H5Dclose(dset);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataset for array of integers.");
   }
}


void
ReaderWriterBase::writeIntegerValue(
   const string& a_name,
   hid_t a_home,
   int a_val)
{
   // Create the dataspace for a single value.
   hsize_t dims = 1;
   hid_t dspace = H5Screate_simple(1, &dims, NULL);
   if (dspace < 0) {
      LOKI_ABORT("Unable to create dataspace for a single integer.");
   }

   // Create the dataset for the integer value.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dcreate_vers) && H5Dcreate_vers == 1
   hid_t dset = H5Dcreate(a_home,
      a_name.c_str(),
      H5T_STD_I32BE,
      dspace,
      H5P_DEFAULT);
#else
   hid_t dset = H5Dcreate(a_home,
      a_name.c_str(),
      H5T_STD_I32BE,
      dspace,
      H5P_DEFAULT,
      H5P_DEFAULT,
      H5P_DEFAULT);
#endif
#else
   hid_t dset = H5Dcreate(a_home,
      a_name.c_str(),
      H5T_STD_I32BE,
      dspace,
      H5P_DEFAULT);
#endif
   if (dset < 0) {
      LOKI_ABORT("Unable to create dataset for a single integer.");
   }

   // Write the integer value.
   herr_t errf =
      H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &a_val);
   if (errf < 0) {
      LOKI_ABORT("Unable to write to dataset for a single integer.");
   }

   // Close the dataset and dataspace for a single integer.
   errf = H5Dclose(dset);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataset for a single integer.");
   }
   errf = H5Sclose(dspace);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataspace for a single integer.");
   }
}


void
ReaderWriterBase::readDoubleArray(
   const string& a_name,
   hid_t a_home,
   vector<double>& a_vals)
{
   // Open the dataset for the double array.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
   hid_t dsid = H5Dopen(a_home, a_name.c_str());
#else
   hid_t dsid = H5Dopen(a_home, a_name.c_str(), H5P_DEFAULT);
#endif
#else
   hid_t dsid = H5Dopen(a_home, a_name.c_str());
#endif
   if (dsid < 0) {
      LOKI_ABORT("Unable to open dataset for double array.");
   }

   // Get the size of the dataset.
   hid_t dspid = H5Dget_space(dsid);
   if (dspid < 0) {
      LOKI_ABORT("Can not get dataspace size for double array.");
   }
   hsize_t nsel = H5Sget_select_npoints(dspid);

   // Allocate space for the double array.
   a_vals.resize(nsel);

   // Read the double array.
   herr_t errf = H5Dread(dsid,
      H5T_NATIVE_DOUBLE,
      H5S_ALL,
      H5S_ALL,
      H5P_DEFAULT,
      &a_vals[0]);
   if (errf < 0) {
      LOKI_ABORT("Unable to read from dataset for double array.");
   }

   // Close the dataset for the double array.
   errf = H5Dclose(dsid);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataset for double array.");
   }
}


void
ReaderWriterBase::readDoubleValue(
   const string& a_name,
   hid_t a_home,
   double& a_val)
{
   // Open the dataset for the double value.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
   hid_t dset = H5Dopen(a_home, a_name.c_str());
#else
   hid_t dset = H5Dopen(a_home, a_name.c_str(), H5P_DEFAULT);
#endif
#else
   hid_t dset = H5Dopen(a_home, a_name.c_str());
#endif
   if (dset < 0) {
      LOKI_ABORT("Unable to open dataset for a single double.");
   }

   // Read the double value.
   herr_t errf =
      H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &a_val);
   if (errf < 0) {
      LOKI_ABORT("Unable to read from dataset for a single double.");
   }

   // Close the dataset for a single double.
   errf = H5Dclose(dset);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataset for a single double.");
   }
}


void
ReaderWriterBase::readIntegerValue(
   const string& a_name,
   hid_t a_home,
   int& a_val)
{
   // Open the dataset for the integer value.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
   hid_t dset = H5Dopen(a_home, a_name.c_str());
#else
   hid_t dset = H5Dopen(a_home, a_name.c_str(), H5P_DEFAULT);
#endif
#else
   hid_t dset = H5Dopen(a_home, a_name.c_str());
#endif
   if (dset < 0) {
      LOKI_ABORT("Unable to open dataset for a single integer.");
   }

   // Read the integer value.
   herr_t errf =
      H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &a_val);
   if (errf < 0) {
      LOKI_ABORT("Unable to read from dataset for a single integer.");
   }

   // Close the dataset for a single integer.
   errf = H5Dclose(dset);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataset for a single integer.");
   }
}
} // end namespace Loki
