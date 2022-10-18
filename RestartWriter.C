/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "RestartWriter.H"
#include "Loki_Defines.H"
#include "Loki_Utilities.H"

namespace Loki {

const int RestartWriter::s_TAG_BATON = 4379;

RestartWriter::RestartWriter(
   const string& a_base_name,
   int a_max_num_files,
   bool a_post_processing)
   : RestartReaderWriterBase(a_base_name, a_max_num_files, a_post_processing)
{
   // If this is the first writer to a bulkdata file then create and immediately
   // close that file.
   // Get its bulkdata file name.
   int file_num =
      Loki_Utilities::s_my_id*m_max_num_files/Loki_Utilities::s_num_procs;
   m_bulkdata_file_name << m_base_name << ".g" << file_num;
   // Create and close the bulkdata file.
   if (m_first_reader_writer) {
      createFileAndRoot(m_bulkdata_file_name.str(),
         m_bulkdata_file,
         m_bulkdata_root);
      closeFileAndRoot(m_bulkdata_file, m_bulkdata_root);
   }

   // Create the metadata file and its root group.
   // There seems to be some kind of race condition happening here.  The
   // metadata file uses parallel IO so the HDF5 calls are all supposed to be
   // collective.  Yet without this syncronization the H5Fcreate for the
   // metadata file fails randomly.
   MPI_Barrier(MPI_COMM_WORLD);
   m_metadata_groups.push(0);
   m_metadata_group_names.push("root");
   createFileAndRoot(a_base_name,
      m_metadata_file,
      m_metadata_groups.top(),
      true);
}


RestartWriter::~RestartWriter()
{
}


void
RestartWriter::writeIntegerValue(
   const string& a_name,
   const int& a_val,
   bool a_write_data)
{
   // Create the dataspace for a single value.
   hid_t dspace = H5Screate(H5S_SCALAR);

   // Create the dataset for the integer value.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dcreate_vers) && H5Dcreate_vers == 1
   hid_t dset = H5Dcreate(m_metadata_groups.top(),
                          a_name.c_str(),
                          H5T_STD_I32BE,
                          dspace,
                          H5P_DEFAULT);
#else
   hid_t dset = H5Dcreate(m_metadata_groups.top(),
                          a_name.c_str(),
                          H5T_STD_I32BE,
                          dspace,
                          H5P_DEFAULT,
                          H5P_DEFAULT,
                          H5P_DEFAULT);
#endif
#else
   hid_t dset = H5Dcreate(m_metadata_groups.top(),
                          a_name.c_str(),
                          H5T_STD_I32BE,
                          dspace,
                          H5P_DEFAULT);
#endif
   if (dset < 0) {
      LOKI_ABORT("Unable to create dataset for a single integer.");
   }
   hid_t memspace = H5Screate(H5S_SCALAR);
   hid_t plistID = H5Pcreate(H5P_DATASET_XFER);
   if (plistID < 0) {
      LOKI_ABORT("Could not create xfer property list for an integer.");
   }
   H5Pset_dxpl_mpio(plistID, H5FD_MPIO_INDEPENDENT);

   // Write the integer value.
   herr_t errf;
   if (a_write_data) {
      errf = H5Dwrite(dset, H5T_NATIVE_INT, memspace, dspace, plistID, &a_val);
      if (errf < 0) {
         LOKI_ABORT("Unable to write to dataset for a single integer.");
      }
   }

   H5Pclose(plistID);
   H5Sclose(memspace);
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
RestartWriter::writeIntegerArray(
   const string& a_name,
   const int* a_vals,
   int a_num_vals,
   bool a_write_data)
{
   // Create the dataspace for the integer array.
   hsize_t dims[1] = {a_num_vals}, rank = 1;
   hid_t dspace = H5Screate_simple(rank, dims, NULL);

   // Create the dataset for the integer array.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dcreate_vers) && H5Dcreate_vers == 1
   hid_t dset = H5Dcreate(m_metadata_groups.top(),
                          a_name.c_str(),
                          H5T_STD_I32BE,
                          dspace,
                          H5P_DEFAULT);
#else
   hid_t dset = H5Dcreate(m_metadata_groups.top(),
                          a_name.c_str(),
                          H5T_STD_I32BE,
                          dspace,
                          H5P_DEFAULT,
                          H5P_DEFAULT,
                          H5P_DEFAULT);
#endif
#else
   hid_t dset = H5Dcreate(m_metadata_groups.top(),
                          a_name.c_str(),
                          H5T_STD_I32BE,
                          dspace,
                          H5P_DEFAULT);
#endif
   if (dset < 0) {
      LOKI_ABORT("Unable to create dataset for an integer array");
   }
   hid_t memspace = H5Screate_simple(rank, dims, NULL);
   hid_t plistID = H5Pcreate(H5P_DATASET_XFER);
   if (plistID < 0) {
      LOKI_ABORT("Could not create xfer property list for an integer array.");
   }
   H5Pset_dxpl_mpio(plistID, H5FD_MPIO_INDEPENDENT);

   // Write the integer array.
   herr_t errf;
   if (a_write_data) {
      errf = H5Dwrite(dset, H5T_NATIVE_INT, memspace, dspace, plistID, a_vals);
      if (errf < 0) {
         LOKI_ABORT("Unable to write to dataset for a integer array.");
      }
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
RestartWriter::writeDoubleValue(
   const string& a_name,
   const double& a_val,
   bool a_write_data)
{
   // Create the dataspace for a single value.
   hid_t dspace = H5Screate(H5S_SCALAR);

   // Create the dataset for the double value.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dcreate_vers) && H5Dcreate_vers == 1
   hid_t dset = H5Dcreate(m_metadata_groups.top(),
                          a_name.c_str(),
                          H5T_NATIVE_DOUBLE,
                          dspace,
                          H5P_DEFAULT);
#else
   hid_t dset = H5Dcreate(m_metadata_groups.top(),
                          a_name.c_str(),
                          H5T_NATIVE_DOUBLE,
                          dspace,
                          H5P_DEFAULT,
                          H5P_DEFAULT,
                          H5P_DEFAULT);
#endif
#else
   hid_t dset = H5Dcreate(m_metadata_groups.top(),
                          a_name.c_str(),
                          H5T_NATIVE_DOUBLE,
                          dspace,
                          H5P_DEFAULT);
#endif
   if (dset < 0) {
      LOKI_ABORT("Unable to create dataset for a single double.");
   }
   hid_t memspace = H5Screate(H5S_SCALAR);
   hid_t plistID = H5Pcreate(H5P_DATASET_XFER);
   if (plistID < 0) {
      LOKI_ABORT("Could not create xfer property list for a single double.");
   }
   H5Pset_dxpl_mpio(plistID, H5FD_MPIO_INDEPENDENT);

   // Write the double value.
   herr_t errf;
   if (a_write_data) {
      errf = H5Dwrite(dset,
         H5T_NATIVE_DOUBLE,
         memspace,
         dspace,
         plistID,
         &a_val);
      if (errf < 0) {
         LOKI_ABORT("Unable to write to dataset for a single double.");
      }
   }

   H5Pclose(plistID);
   H5Sclose(memspace);
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
RestartWriter::writeBulkDoubleValue(
   const string& a_name,
   const double& a_val)
{

   // Write the value to the bulkdata file.

   // If this processor is the first to write to its bulkdata file then it can
   // just open it (it's already been created).  If this processor is not the
   // first to write to its bulkdata file then it must wait for the previous
   // writer to complete writing to the file and then open the bulkdata file.
   if (!m_first_reader_writer) {
      MPI_Status status;
      int info;
      int prev_writer = Loki_Utilities::s_my_id-1;
      int err = MPI_Recv(&info,
         1,
         MPI_INT,
         prev_writer,
         s_TAG_BATON,
         MPI_COMM_WORLD,
         &status);
      if (err != MPI_SUCCESS) {
         LOKI_ABORT("Could not get OK from previous writer.");
      }
   }
   openFileAndRoot(m_bulkdata_file_name.str(),
      m_bulkdata_file,
      m_bulkdata_root);

   // Write the value to the bulkdata file.
   ostringstream array_name;
   array_name << m_metadata_group_names.top() << "\\" << a_name
              << ".p" << Loki_Utilities::s_my_id;
   ReaderWriterBase::writeDoubleValue(array_name.str(),
      m_bulkdata_root,
      a_val);

   // Close the bulkdata file and its root group.
   closeFileAndRoot(m_bulkdata_file, m_bulkdata_root);

   // If this processor is not the last to write to its bulkdata file then send
   // a message to the next processor so that it may proceed with opening the
   // bulkdata file and writing to it.
   if (!m_last_reader_writer) {
      int info = 0;
      int next_writer = Loki_Utilities::s_my_id+1;
      MPI_Send(&info, 1, MPI_INT, next_writer, s_TAG_BATON, MPI_COMM_WORLD);
   }

   // The next write of any distributed data to a bulkdata file may not occur
   // until all processors writing to a given bulkdata file have completed
   // writing this distributed data.  We could be really clever and create a
   // subcommunicator for the set of processors writing to each bulkdata file
   // and barrier on that communicator but there is enough other blocking
   // communication elsewhere such a scheme will not affect overall performance.
   if (m_max_num_files != Loki_Utilities::s_num_procs) {
      MPI_Barrier(MPI_COMM_WORLD);
   }
}


void
RestartWriter::writeDoubleArray(
   const string& a_name,
   const double* a_vals,
   int a_num_vals,
   bool a_write_data)
{
   // Create the dataspace for the double array.
   hsize_t dims[1] = {a_num_vals}, rank = 1;
   hid_t dspace = H5Screate_simple(rank, dims, NULL);

   // Create the dataset for the double array.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dcreate_vers) && H5Dcreate_vers == 1
   hid_t dset = H5Dcreate(m_metadata_groups.top(),
                          a_name.c_str(),
                          H5T_NATIVE_DOUBLE,
                          dspace,
                          H5P_DEFAULT);
#else
   hid_t dset = H5Dcreate(m_metadata_groups.top(),
                          a_name.c_str(),
                          H5T_NATIVE_DOUBLE,
                          dspace,
                          H5P_DEFAULT,
                          H5P_DEFAULT,
                          H5P_DEFAULT);
#endif
#else
   hid_t dset = H5Dcreate(m_metadata_groups.top(),
                          a_name.c_str(),
                          H5T_NATIVE_DOUBLE,
                          dspace,
                          H5P_DEFAULT);
#endif
   if (dset < 0) {
      LOKI_ABORT("Unable to create dataset for a double array");
   }
   hid_t memspace = H5Screate_simple(rank, dims, NULL);
   hid_t plistID = H5Pcreate(H5P_DATASET_XFER);
   if (plistID < 0) {
      LOKI_ABORT("Could not create xfer property list for a double array.");
   }
   H5Pset_dxpl_mpio(plistID, H5FD_MPIO_INDEPENDENT);

   // Write the double array.
   herr_t errf;
   if (a_write_data) {
      errf = H5Dwrite(dset,
         H5T_NATIVE_DOUBLE,
         memspace,
         dspace,
         plistID,
         a_vals);
      if (errf < 0) {
         LOKI_ABORT("Unable to write to dataset for a double array.");
      }
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
RestartWriter::writeString(
   const string& a_name,
   const string& a_val,
   bool a_write_data)
{
   // Create the dataspace for the string.
   hsize_t dims[1] = {a_val.length()+1}, rank = 1;
   hid_t dspace = H5Screate_simple(rank, dims, NULL);

   // Create the dataset for the string.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dcreate_vers) && H5Dcreate_vers == 1
   hid_t dset = H5Dcreate(m_metadata_groups.top(),
                          a_name.c_str(),
                          H5T_NATIVE_UCHAR,
                          dspace,
                          H5P_DEFAULT);
#else
   hid_t dset = H5Dcreate(m_metadata_groups.top(),
                          a_name.c_str(),
                          H5T_NATIVE_UCHAR,
                          dspace,
                          H5P_DEFAULT,
                          H5P_DEFAULT,
                          H5P_DEFAULT);
#endif
#else
   hid_t dset = H5Dcreate(m_metadata_groups.top(),
                          a_name.c_str(),
                          H5T_NATIVE_UCHAR,
                          dspace,
                          H5P_DEFAULT);
#endif
   if (dset < 0) {
      LOKI_ABORT("Unable to create dataset for a string");
   }
   hid_t memspace = H5Screate_simple(rank, dims, NULL);

   // Write the string.
   herr_t errf;
   if (a_write_data) {
      errf = H5Dwrite(dset,
         H5T_NATIVE_UCHAR,
         memspace,
         dspace,
         H5P_DEFAULT,
         a_val.c_str());
      if (errf < 0) {
         LOKI_ABORT("Unable to write to dataset for a string.");
      }
   }

   H5Sclose(memspace);
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
RestartWriter::writeParallelArray(
   const string& a_name,
   const ParallelArray& a_array,
   bool a_write_meta_data)
{
   // We should only be writing 2D, 3D (2 spatial dims + components) or 4D data.
   int dim = a_array.dim();
   if (dim != 2 && dim != 3 && dim != 4) {
      LOKI_ABORT("Trying to write distributed data that is not 2D or 4D.");
   }

   // Write information about the decomposition to the metatdata file.

   // Gather distribution info.
   int distribInfoSz = 2+4*dim;
   int distribInfo[18];
   distribInfo[0] = a_array.procLo();
   distribInfo[1] = a_array.procHi();
   const vector<int>& dim_procs = a_array.getDimPartitions();
   const vector<int>& left_partitions_size = a_array.getLeftPartitionsSize();
   const vector<int>& num_left_partitions = a_array.getNumLeftPartitions();
   for (int d = 0; d < dim; ++d) {
      distribInfo[d+2] = -a_array.numGhosts();
      distribInfo[d+2+dim] = dim_procs[d];
      distribInfo[d+2+2*dim] = left_partitions_size[d];
      distribInfo[d+2+3*dim] = num_left_partitions[d];
   }

   // Write the distribution info.
   pushSubDir(a_name);
   writeIntegerArray("distribInfo",
      distribInfo,
      distribInfoSz,
      a_write_meta_data);
   popSubDir();

   // Write the local array to the bulkdata file.

   // If this processor is the first to write to its bulkdata file then it can
   // just open it (it's already been created).  If this processor is not the
   // first to write to its bulkdata file then it must wait for the previous
   // writer to complete writing to the file and then open the bulkdata file.
   if (!m_first_reader_writer) {
      MPI_Status status;
      int info;
      int prev_writer = Loki_Utilities::s_my_id-1;
      int err = MPI_Recv(&info,
         1,
         MPI_INT,
         prev_writer,
         s_TAG_BATON,
         MPI_COMM_WORLD,
         &status);
      if (err != MPI_SUCCESS) {
         LOKI_ABORT("Could not get OK from previous writer.");
      }
   }
   openFileAndRoot(m_bulkdata_file_name.str(),
      m_bulkdata_file,
      m_bulkdata_root);

   // Create the dataspace for the double array.
   // Early versions of HDF5 do not allow writing empty datasets.  So if a_array
   // does not contain any data we write a single scalar value of 0.  The reader
   // must take this convention into accout.
   const ParallelArray::Box& data_box = a_array.dataBox();
   bool is_empty = data_box.size() == 0;
   hsize_t dims[4] = {1, 1, 1, 1};
   if (!is_empty) {
      if (dim == 2) {
         dims[0] = data_box.numberOfCells(1);
         dims[1] = data_box.numberOfCells(0);
      }
      else if (dim == 3) {
         dims[0] = data_box.numberOfCells(2);
         dims[1] = data_box.numberOfCells(1);
         dims[2] = data_box.numberOfCells(0);
      }
      else {
         dims[0] = data_box.numberOfCells(3);
         dims[1] = data_box.numberOfCells(2);
         dims[2] = data_box.numberOfCells(1);
         dims[3] = data_box.numberOfCells(0);
      }
   }
   hid_t dspace = H5Screate_simple(dim, dims, NULL);
   if (dspace < 0) {
      LOKI_ABORT("Unable to create dataspace for a double array.");
   }

   // Write the local array.
   ostringstream array_name;
   array_name << m_metadata_group_names.top() << "\\" << a_name
              << ".p" << Loki_Utilities::s_my_id;
   if (!is_empty) {
      ReaderWriterBase::writeDoubleArray(array_name.str(),
         m_bulkdata_root,
         dspace,
         a_array.getData());
   }
   else {
      double scalar = 0.0;
      ReaderWriterBase::writeDoubleArray(array_name.str(),
         m_bulkdata_root,
         dspace,
         &scalar);
   }

   // Close the dataspace for the double array.
   herr_t errf = H5Sclose(dspace);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataspace for double array.");
   }

   // Close the bulkdata file and its root group.
   closeFileAndRoot(m_bulkdata_file, m_bulkdata_root);

   // If this processor is not the last to write to its bulkdata file then send
   // a message to the next processor so that it may proceed with opening the
   // bulkdata file and writing to it.
   if (!m_last_reader_writer) {
      int info = 0;
      int next_writer = Loki_Utilities::s_my_id+1;
      MPI_Send(&info, 1, MPI_INT, next_writer, s_TAG_BATON, MPI_COMM_WORLD);
   }

   // The next write of any distributed data to a bulkdata file may not occur
   // until all processors writing to a given bulkdata file have completed
   // writing this distributed data.  We could be really clever and create a
   // subcommunicator for the set of processors writing to each bulkdata file
   // and barrier on that communicator but there is enough other blocking
   // communication elsewhere such a scheme will not affect overall performance.
   if (m_max_num_files != Loki_Utilities::s_num_procs) {
      MPI_Barrier(MPI_COMM_WORLD);
   }
}

} // end namespace Loki
