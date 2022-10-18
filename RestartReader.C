/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "RestartReader.H"
#include "Loki_Defines.H"
#include "Loki_Utilities.H"

#include <sstream>

namespace Loki {

const int RestartReader::s_TAG_BATON = 4379;

RestartReader::RestartReader(
   const string& a_base_name,
   int a_max_num_files,
   bool a_post_processing)
   : RestartReaderWriterBase(a_base_name, a_max_num_files, a_post_processing)
{
   // Open the metadata file and its root group.
   m_metadata_groups.push(0);
   m_metadata_group_names.push("root");
   openFileAndRoot(a_base_name, m_metadata_file, m_metadata_groups.top());
}


RestartReader::~RestartReader()
{
}


void
RestartReader::readIntegerValue(
   const string& a_name,
   int& a_val)
{
   // Open the dataset for the integer value.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
   hid_t dset = H5Dopen(m_metadata_groups.top(), a_name.c_str());
#else
   hid_t dset = H5Dopen(m_metadata_groups.top(), a_name.c_str(), H5P_DEFAULT);
#endif
#else
   hid_t dset = H5Dopen(m_metadata_groups.top(), a_name.c_str());
#endif
   if (dset < 0) {
      LOKI_ABORT("Unable to open dataset for a single integer.");
   }

   // Get the dataspace for a single integer value.
   hid_t dspace = H5Dget_space(dset);

   // Set the data transfer mode.
   hid_t plistID = H5Pcreate(H5P_DATASET_XFER);
   if (plistID < 0) {
      LOKI_ABORT("Could not create xfer property list for a single integer.");
   }
   H5Pset_dxpl_mpio(plistID, H5FD_MPIO_INDEPENDENT);

   // Read the integer value.
   herr_t errf =
      H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &a_val);
   if (errf < 0) {
      LOKI_ABORT("Unable to read from dataset for a single integer.");
   }

   H5Pclose(plistID);
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
RestartReader::readIntegerArray(
   const string& a_name,
   int* a_vals,
   int a_num_vals)
{
   // Open the dataset for the integer array.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
   hid_t dset = H5Dopen(m_metadata_groups.top(), a_name.c_str());
#else
   hid_t dset = H5Dopen(m_metadata_groups.top(), a_name.c_str(), H5P_DEFAULT);
#endif
#else
   hid_t dset = H5Dopen(m_metadata_groups.top(), a_name.c_str());
#endif
   if (dset < 0) {
      LOKI_ABORT("Unable to open dataset for an integer array");
   }

   // Get the dataspace for the integer array.
   hsize_t dims[1];
   hid_t dspace = H5Dget_space(dset);
   H5Sget_simple_extent_dims(dspace, dims, NULL);
   if (dims[0] != a_num_vals) {
      LOKI_ABORT("Requested read of integer array does not match size in file.");
   }

   // Describe the array in memory.
   hid_t memspace = H5Screate_simple(1, dims, NULL);

   // Set the data transfer mode.
   hid_t plistID = H5Pcreate(H5P_DATASET_XFER);
   if (plistID < 0) {
      LOKI_ABORT("Could not create xfer property list for an integer array.");
   }
   H5Pset_dxpl_mpio(plistID, H5FD_MPIO_INDEPENDENT);

   // Read the integer array.
   herr_t errf =
      H5Dread(dset, H5T_NATIVE_INT, memspace, dspace, plistID, a_vals);
   if (errf < 0) {
      LOKI_ABORT("Unable to read from dataset for a integer array.");
   }

   H5Pclose(plistID);
   H5Sclose(memspace);
   // Close the dataset and dataspace for the integer array.
   errf = H5Dclose(dset);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataset for a integer array.");
   }
   errf = H5Sclose(dspace);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataspace for a integer array.");
   }
}


void
RestartReader::readDoubleValue(
   const string& a_name,
   double& a_val)
{
   // Open the dataset for the double value.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
   hid_t dset = H5Dopen(m_metadata_groups.top(), a_name.c_str());
#else
   hid_t dset = H5Dopen(m_metadata_groups.top(), a_name.c_str(), H5P_DEFAULT);
#endif
#else
   hid_t dset = H5Dopen(m_metadata_groups.top(), a_name.c_str());
#endif
   if (dset < 0) {
      LOKI_ABORT("Unable to open dataset for a single double.");
   }

   // Get the dataspace for a single double value.
   hid_t dspace = H5Dget_space(dset);

   // Set the data transfer mode.
   hid_t plistID = H5Pcreate(H5P_DATASET_XFER);
   if (plistID < 0) {
      LOKI_ABORT("Could not create xfer property list for a single double.");
   }
   H5Pset_dxpl_mpio(plistID, H5FD_MPIO_INDEPENDENT);

   // Read the double value.
   herr_t errf =
      H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &a_val);
   if (errf < 0) {
      LOKI_ABORT("Unable to read from dataset for a single double.");
   }

   H5Pclose(plistID);
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
RestartReader::readBulkDoubleValue(
   const string& a_name,
   double& a_val)
{
   // Read the value from the bulkdata file.

   // Get this processor's bulkdata file name.
   int file_num =
      Loki_Utilities::s_my_id*m_max_num_files/Loki_Utilities::s_num_procs;
   ostringstream bulkdata_file_name;
   bulkdata_file_name << m_base_name << ".g" << file_num;

   // If this processor is the first to read from its bulkdata file then it can
   // just open it.  If this processor is not the first to read from its
   // bulkdata file then it must wait for the previous reader to complete
   // reading from the file and then open the bulkdata file.
   if (!m_first_reader_writer) {
      MPI_Status status;
      int info;
      int prev_reader = Loki_Utilities::s_my_id-1;
      int err = MPI_Recv(&info,
         1,
         MPI_INT,
         prev_reader,
         s_TAG_BATON,
         MPI_COMM_WORLD,
         &status);
      if (err != MPI_SUCCESS) {
         LOKI_ABORT("Could not get OK from previous reader.");
      }
   }
   openFileAndRoot(bulkdata_file_name.str(),
      m_bulkdata_file,
      m_bulkdata_root);

   // Read the value from the bulkdata file.
   ostringstream array_name;
   array_name << m_metadata_group_names.top() << "\\" << a_name
              << ".p" << Loki_Utilities::s_my_id;
   ReaderWriterBase::readDoubleValue(array_name.str(),
      m_bulkdata_root,
      a_val);

   // Close the bulkdata file and its root group.
   closeFileAndRoot(m_bulkdata_file, m_bulkdata_root);

   // If this processor is not the last to read from its bulkdata file then send
   // a message to the next processor so that it may proceed with opening the
   // bulkdata file and writing to it.
   if (!m_last_reader_writer) {
      int info = 0;
      int next_reader = Loki_Utilities::s_my_id+1;
      MPI_Send(&info, 1, MPI_INT, next_reader, s_TAG_BATON, MPI_COMM_WORLD);
   }
}


void
RestartReader::readDoubleArray(
   const string& a_name,
   double* a_vals,
   int a_num_vals)
{
   // Open the dataset for the double array.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
   hid_t dset = H5Dopen(m_metadata_groups.top(), a_name.c_str());
#else
   hid_t dset = H5Dopen(m_metadata_groups.top(), a_name.c_str(), H5P_DEFAULT);
#endif
#else
   hid_t dset = H5Dopen(m_metadata_groups.top(), a_name.c_str());
#endif
   if (dset < 0) {
      LOKI_ABORT("Unable to open dataset for a double array");
   }

   // Get the dataspace for the double array.
   hsize_t dims[1];
   hid_t dspace = H5Dget_space(dset);
   H5Sget_simple_extent_dims(dspace, dims, NULL);
   if (dims[0] != a_num_vals) {
      LOKI_ABORT("Requested read of double array does not match size in file.");
   }

   // Describe the array in memory.
   hid_t memspace = H5Screate_simple(1, dims, NULL);

   // Set the data transfer mode.
   hid_t plistID = H5Pcreate(H5P_DATASET_XFER);
   if (plistID < 0) {
      LOKI_ABORT("Could not create xfer property list for a double array.");
   }
   H5Pset_dxpl_mpio(plistID, H5FD_MPIO_INDEPENDENT);

   // Read the double array.
   herr_t errf =
      H5Dread(dset, H5T_NATIVE_DOUBLE, memspace, dspace, plistID, a_vals);
   if (errf < 0) {
      LOKI_ABORT("Unable to read to dataset for a double array.");
   }

   H5Pclose(plistID);
   H5Sclose(memspace);
   // Close the dataset and dataspace for the double array.
   errf = H5Dclose(dset);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataset for a double array.");
   }
   errf = H5Sclose(dspace);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataspace for a double array.");
   }
}


void
RestartReader::readString(
   const string& a_name,
   string& a_val)
{
   // Open the dataset for the string.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
   hid_t dset = H5Dopen(m_metadata_groups.top(), a_name.c_str());
#else
   hid_t dset = H5Dopen(m_metadata_groups.top(), a_name.c_str(), H5P_DEFAULT);
#endif
#else
   hid_t dset = H5Dopen(m_metadata_groups.top(), a_name.c_str());
#endif
   if (dset < 0) {
      LOKI_ABORT("Unable to open dataset for a string");
   }

   // Get the dataspace for the string.
   hsize_t dims[1];
   hid_t dspace = H5Dget_space(dset);
   H5Sget_simple_extent_dims(dspace, dims, NULL);

   // Read the string.
   char* temp = new char[dims[0]];
   herr_t errf =
      H5Dread(dset, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
   if (errf < 0) {
      LOKI_ABORT("Unable to read from dataset for a string.");
   }
   a_val = temp;
   delete [] temp;

   // Close the dataset and dataspace for the string.
   errf = H5Dclose(dset);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataset for a string.");
   }
   errf = H5Sclose(dspace);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataspace for a string.");
   }
}


void
RestartReader::readParallelArray(
   const string& a_name,
   ParallelArray& a_array)
{
   // This is only for the use of the main code.
   if (m_post_processing) {
      LOKI_ABORT("Post processor using improper method to read ParallelArray.");
   }

   // We should only be reading 2D, 3D (2 spatial dims + components) or 4D data.
   int dim = a_array.dim();
   if (dim != 2 && dim != 3 && dim != 4) {
      LOKI_ABORT("Trying to read ParallelArray that is not 2D or 4D.");
   }

   // Get this processor's bulkdata file name.
   int file_num =
      Loki_Utilities::s_my_id*m_max_num_files/Loki_Utilities::s_num_procs;
   ostringstream bulkdata_file_name;
   bulkdata_file_name << m_base_name << ".g" << file_num;

   // Read the ParallelArray from the bulkdata file.

   // If this processor is the first to read from its bulkdata file then it can
   // just open it.  If this processor is not the first to read from its
   // bulkdata file then it must wait for the previous reader to complete
   // reading from the file and then open the bulkdata file.
   if (!m_first_reader_writer) {
      MPI_Status status;
      int info;
      int prev_reader = Loki_Utilities::s_my_id-1;
      int err = MPI_Recv(&info,
         1,
         MPI_INT,
         prev_reader,
         s_TAG_BATON,
         MPI_COMM_WORLD,
         &status);
      if (err != MPI_SUCCESS) {
         LOKI_ABORT("Could not get OK from previous reader.");
      }
   }
   openFileAndRoot(bulkdata_file_name.str(),
      m_bulkdata_file,
      m_bulkdata_root);

   // Check to see if a_array actually contains any data.  If it does try to
   // read its dataset.  Recall from the writer that early versions of HDF5 do
   // not allow writing empty datasets.  So we write a single scalar in place of
   // empty datasets.  So we can't use just read the empty dataset which would
   // be a no-op since it's not actually empty.
   if (!a_array.dataBox().empty()) {
      // Open the dataset for this rank's ParallelArray.
      ostringstream array_name;
      array_name << m_metadata_group_names.top() << "\\" << a_name
                 << ".p" << Loki_Utilities::s_my_id;
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
      hid_t dset = H5Dopen(m_bulkdata_root, array_name.str().c_str());
#else
      hid_t dset = H5Dopen(m_bulkdata_root, array_name.str().c_str(), H5P_DEFAULT);
#endif
#else
      hid_t dset = H5Dopen(m_bulkdata_root, array_name.str().c_str());
#endif
      if (dset < 0) {
         LOKI_ABORT("Unable to open dataset for a ParallelArray.");
      }

      // Get the dataspace for the ParallelArray.
      hsize_t dims[dim];
      hid_t dspace = H5Dget_space(dset);
      H5Sget_simple_extent_dims(dspace, dims, NULL);

      // Describe the array in memory.
      hid_t memspace = H5Screate_simple(dim, dims, NULL);

      // Read the ParallelArray.
      herr_t errf = H5Dread(dset,
         H5T_NATIVE_DOUBLE,
         memspace,
         dspace,
         H5P_DEFAULT,
         a_array.getData());
      if (errf < 0) {
         LOKI_ABORT("Unable to read to dataset for a double array.");
      }

      H5Sclose(memspace);
      // Close the dataset and dataspace for the ParallelArray.
      errf = H5Dclose(dset);
      if (errf < 0) {
         LOKI_ABORT("Unable to close dataset for a ParallelArray.");
      }
      errf = H5Sclose(dspace);
      if (errf < 0) {
         LOKI_ABORT("Unable to close dataspace for a ParallelArray.");
      }
   }

   // Close the bulkdata file and its root group.
   closeFileAndRoot(m_bulkdata_file, m_bulkdata_root);

   // If this processor is not the last to read from its bulkdata file then send
   // a message to the next processor so that it may proceed with opening the
   // bulkdata file and writing to it.
   if (!m_last_reader_writer) {
      int info = 0;
      int next_reader = Loki_Utilities::s_my_id+1;
      MPI_Send(&info, 1, MPI_INT, next_reader, s_TAG_BATON, MPI_COMM_WORLD);
   }
}


void
RestartReader::readParallelArray(
   const string& a_name,
   ParallelArray& a_array,
   int a_num_generating_processes)
{
   // This is only for the use of the post processor.
   if (!m_post_processing) {
      LOKI_ABORT("Main code using improper method to read ParallelArray.");
   }

   // We should only be reading 4D data for now.
   if (a_array.dim() != 4) {
      LOKI_ABORT("Trying to read ParallelArray that is not 4D.");
   }

   // Get the global serial ParallelArray's extent and data buffer.
   const ParallelArray::Box& serial_box = a_array.dataBox();
   int ns0 = serial_box.numberOfCells(0);
   int ns1 = serial_box.numberOfCells(1);
   int ns2 = serial_box.numberOfCells(2);
   int ns3 = serial_box.numberOfCells(3);
   double* serialBuff = a_array.getData();

   // Create a box which will hold the bounds of each local distributed
   // ParallelArray to be read.
   tbox::Dimension dim(4);
   tbox::Box data_box(dim);

   // Read information about the decomposition from the metatdata file.
   int distribInfoSz = 18;
   int distribInfo[18];
   pushSubDir(a_name);
   readIntegerArray("distribInfo", distribInfo, distribInfoSz);
   popSubDir();
   int lowest_rank = distribInfo[0];
   int highest_rank = distribInfo[1];
   const int* base_indices = distribInfo+2;
   const int* procDims = distribInfo+6;
   const int* dimSizesLeft = distribInfo+10;
   const int* dimNumLeft = distribInfo+14;

   // Now read the data from each rank that wrote the distributed ParallelArray.
   for (int rank = lowest_rank; rank <= highest_rank; ++rank) {
      // Get this writing processor's bulkdata file name.
      int file_num = rank*m_max_num_files/a_num_generating_processes;
      ostringstream bulkdata_file_name;
      bulkdata_file_name << m_base_name << ".g" << file_num;

      // Open the bulkdata file.
      openFileAndRoot(bulkdata_file_name.str(),
         m_bulkdata_file,
         m_bulkdata_root);
      bulkdata_file_name.str("");

      // Open the dataset for this rank's distributed ParallelArray.
      ostringstream array_name;
      array_name << m_metadata_group_names.top() << "\\" << a_name
                 << ".p" << rank;
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
      hid_t dset = H5Dopen(m_bulkdata_root, array_name.str().c_str());
#else
      hid_t dset =
         H5Dopen(m_bulkdata_root, array_name.str().c_str(), H5P_DEFAULT);
#endif
#else
      hid_t dset = H5Dopen(m_bulkdata_root, array_name.str().c_str());
#endif
      if (dset < 0) {
         LOKI_ABORT("Unable to open dataset for a ParallelArray.");
      }
      array_name.str("");

      // Get the dataspace for the ParallelArray.
      hsize_t dims[4];
      hid_t dspace = H5Dget_space(dset);
      H5Sget_simple_extent_dims(dspace, dims, NULL);

      // Describe the array in memory.
      hid_t memspace = H5Screate_simple(4, dims, NULL);

      // Allocate a buffer for this rank's data.
      double *distBuff = new double[dims[0]*dims[1]*dims[2]*dims[3]];

      // Read the ParallelArray.
      herr_t errf = H5Dread(dset,
         H5T_NATIVE_DOUBLE,
         memspace,
         dspace,
         H5P_DEFAULT,
         distBuff);
      if (errf < 0) {
         LOKI_ABORT("Unable to read to dataset for a ParallelArray.");
      }

      // We need to put the put the local distributed data into its location in
      // the serial global array.
      getArrayDataBox(rank,
         lowest_rank,
         base_indices,
         procDims,
         dimSizesLeft,
         dimNumLeft,
         data_box);

      // Now transfer the data in the local box into the serial array.
      int nd[4], istart[4], iend[4];
      for (int i = 0; i < 4; ++i) {
         nd[i] = data_box.numberCells(i);
         if (data_box.lower(i) == serial_box.lower(i)) {
            istart[i] = 0;
         }
         else {
            istart[i] = -base_indices[i];
         }
         if (data_box.upper(i) == serial_box.upper(i)) {
            iend[i] = nd[i];
         }
         else {
            iend[i] = nd[i] + base_indices[i];
         }
      }
      for (int i3 = istart[3]; i3 < iend[3]; ++i3) {
         for (int i2 = istart[2]; i2 < iend[2]; ++i2) {
            for (int i1 = istart[1]; i1 < iend[1]; ++i1) {
               for (int i0 = istart[0]; i0 < iend[0]; ++i0) {
                  int serialIdx =
                     (data_box.lower(3)-base_indices[3]+i3)*(ns0*ns1*ns2) +
                     (data_box.lower(2)-base_indices[2]+i2)*(ns0*ns1) +
                     (data_box.lower(1)-base_indices[1]+i1)*ns0 +
                     data_box.lower(0)-base_indices[0]+i0;
                  int distIdx =
                     i3*(nd[0]*nd[1]*nd[2]) + i2*(nd[0]*nd[1]) + i1*nd[0] + i0;
                  serialBuff[serialIdx] = distBuff[distIdx];
               }
            }
         }
      }

      delete [] distBuff;

      H5Sclose(memspace);
      // Close the dataset and dataspace for the double array.
      errf = H5Dclose(dset);
      if (errf < 0) {
         LOKI_ABORT("Unable to close dataset for a double array.");
      }
      errf = H5Sclose(dspace);
      if (errf < 0) {
         LOKI_ABORT("Unable to close dataspace for a double array.");
      }

      // Close the bulkdata file and its root group.
      closeFileAndRoot(m_bulkdata_file, m_bulkdata_root);
   }
}


void
RestartReader::getArrayDataBox(
   int a_p,
   int a_baseProc,
   const int* a_dimBase,
   const int* a_dimProc,
   const int* a_dimSizesLeft,
   const int* a_dimNumLeft,
   tbox::Box& box)
{
   // Processor number offset from base (starting) processor for this
   // distribution
   int p0 = a_p - a_baseProc;

   // Compute the "processor vector" (where we are located in the "grid" of
   // distributed local arrays) pv[0], pv[1], ..., pv[numDim-1]
   // Such that 
   //    a_p = a_baseProc + pv[numDim-1] +
   //        a_dimProc[numDim-1]*(pv[numDim-2] + a_dimProc[numDim-2]*(pv[numDim-3] + ... )
   int pv[4];
   for (int d = 3; d >= 0; --d) {
      // Get the number of processors allocated to this dimension and compute
      // the value of the processor vector for this dimension.
      const int thisDimProc = a_dimProc[d];
      pv[d] = p0 - (p0/thisDimProc)*thisDimProc;
      p0 = (p0-pv[d])/thisDimProc;
   }


   // Now compute the data array bounds on the on processor a_p.
   tbox::Dimension dim(4);
   tbox::IntVector lower(dim);
   tbox::IntVector upper(dim);
   for (int d = 0; d < 4; ++d) {
      // Fill in the distribution along each dimension
      // 
      //      +------+------+------+---- ...    --+------+------+
      //        left   left   left                  right  right
      const int left   = max(1, a_dimSizesLeft[d]),
                right  = max(1, a_dimSizesLeft[d]-1);

      // Start the lower end of the data box in this dimension at the base
      // index.
      lower[d] = a_dimBase[d];

      // Add the length of a left box to its lower end for each left box
      // preceeding the data box.
      lower[d] += left*min(pv[d], a_dimNumLeft[d]);

      // Add the length of a right box to its lower end for each right box
      // preceeding the data box.
      if (pv[d] > a_dimNumLeft[d]) {
         lower[d] += right*(pv[d]-a_dimNumLeft[d]);
      }

      upper[d] = lower[d] - a_dimBase[d];
      // If the data box in this dimension is in the midst of the left boxes
      // then add the length of a left box to the lower end of this data box
      // to get the upper end in this dimension.  Otherwise it is in the midst
      // of the right boxes and we add the length of a right box.  Include the
      // ghosts.
      if (pv[d] < a_dimNumLeft[d]) {
         upper[d] += left-1-a_dimBase[d];
      }
      else {
         upper[d] += right-1-a_dimBase[d];
      }
   }

   // Now set the data box bounds.
   for (int i = 0; i < 4; ++i) {
      box.lower(i) = lower[i];
      box.upper(i) = upper[i];
   }
}

} // end namespace Loki
