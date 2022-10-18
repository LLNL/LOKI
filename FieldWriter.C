/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "FieldWriter.H"
#include "Directions.H"
#include "Loki_Defines.H"
#include <sstream>

namespace Loki {

FieldWriter::FieldWriter(
   const string& a_base_name,
   const ProblemDomain& a_domain,
   int a_solution_order,
   int a_time_slices_per_file)
   : ReaderWriterBase(),
     m_base_name(a_base_name),
     m_current_file_index(-1),
     m_time_slices_per_file(a_time_slices_per_file),
     m_time_slices_written_to_current_file(0),
     m_total_time_slices_written(0),
     m_processing_time_slice(false),
     m_nx(a_domain.numberOfCells(0)+a_solution_order),
     m_ny(a_domain.numberOfCells(1)+a_solution_order),
     m_x_coords(a_domain.numberOfCells(0)),
     m_y_coords(a_domain.numberOfCells(1))
{
   for (int i = 0; i < a_domain.numberOfCells(0); ++i) {
      m_x_coords[i] = a_domain.x_lo(0) + (i+0.5)*a_domain.dx(0);
   }
   for (int i = 0; i < a_domain.numberOfCells(1); ++i) {
      m_y_coords[i] = a_domain.x_lo(1) + (i+0.5)*a_domain.dx(1);
   }
}


FieldWriter::FieldWriter(
   const string& a_name,
   const vector<string>& a_field_names,
   int a_nx,
   int a_ny,
   int a_time_slices_per_field)
   : ReaderWriterBase()
{
   // Create the new file and root group.
   createFileAndRoot(a_name, m_file, m_root);

   // Create but do not write to the datasets for the fields in the post
   // processed file.
   createFieldDatasets(a_field_names,
      a_nx,
      a_ny,
      a_time_slices_per_field,
      false);

   // Close the post processed file and the root group.
   closeFileAndRoot(m_file, m_root);

   // Reopen the post processsed file and the root group.
   openFileAndRoot(a_name, m_file, m_root);
}


FieldWriter::~FieldWriter()
{
}


void
FieldWriter::startTimeSlice(
   double a_time,
   double a_dt,
   const vector<string>& a_field_names,
   bool a_proc_lo)
{
   // If already processing a time slice then this is being used improperly.
   if (m_processing_time_slice) {
      LOKI_ABORT("Time slice already being processed.  Missing call to endTimeSlice.");
   }
   m_processing_time_slice = true;

   // See if we need to move to the next field file in the series.
   bool new_file;
   if (m_time_slices_written_to_current_file == 0 ||
       m_time_slices_written_to_current_file == m_time_slices_per_file) {
      ++m_current_file_index;
      m_time_slices_written_to_current_file = 0;
      new_file = true;
   }
   else {
      new_file = false;
   }

   // Form the full name of the file and open/create it.
   ostringstream field_file_name;
   field_file_name << m_base_name << ".fields_"
                   << m_current_file_index << ".hdf";
   if (new_file && a_proc_lo) {
      // Create the new file and root group.
      createFileAndRoot(field_file_name.str(), m_file, m_root);

      // The coordinates are written only once to each field file.
      writeCoords(m_x_coords, m_y_coords);
   }
   else {
      openFileAndRoot(field_file_name.str(), m_file, m_root);
   }

   // The lowest rank processor writing field data must create but not write to
   // the datasets for this time slice of each field.  It then writes the time
   // and dt.
   if (a_proc_lo) {
      createFieldDatasets(a_field_names, m_nx, m_ny, 1, true);
      ostringstream dset_name;
      dset_name << "time_slice_" << m_total_time_slices_written << "_time";
      writeDoubleValue(dset_name.str(), m_root, a_time);
      dset_name.str("");
      dset_name << "time_slice_" << m_total_time_slices_written << "_dt";
      writeDoubleValue(dset_name.str(), m_root, a_dt);
   }
}


void
FieldWriter::startTimeSlice(
   double a_time,
   double a_dt,
   int a_num_tracking_particles,
   int a_num_probes,
   const vector<vector<double> >& a_probes,
   const tbox::IntVector& a_num_cells,
   const vector<string>& a_field_names,
   bool a_proc_lo)
{
   // See if we will create a new file or append to an existing one.
   bool new_file;
   if ((m_time_slices_written_to_current_file == 0 ||
        m_time_slices_written_to_current_file == m_time_slices_per_file)) {
      new_file = true;
   }
   else {
      new_file = false;
   }

   // Call the simple version that doesn't deal with probes.
   startTimeSlice(a_time, a_dt, a_field_names, a_proc_lo);

   // If a new file was created then the lowest rank processor writing field
   // data needs to write this method's additional data to the new file.
   if (new_file && a_proc_lo) {
      // Write the number of tracking particles and number of probes.
      writeIntegerValue("numTrackingParticles",
         m_root,
         a_num_tracking_particles);
      writeIntegerValue("numProbes", m_root, a_num_probes);

      // Now write the probe locations.
      for (int ip(0); ip < a_num_probes; ++ip) {
         ostringstream name;
         // Write the probe's x location.
         name << "ix_probe" << ip;
         int xloc = int(floor(a_probes[0][ip] * a_num_cells[X1]));
         writeIntegerValue(name.str(), m_root, xloc);
         name.str("");

         // Write the probe's x location.
         name << "iy_probe" << ip;
         int yloc = int(floor(a_probes[1][ip] * a_num_cells[X2]));
         writeIntegerValue(name.str(), m_root, yloc);
      }
   }
}


void
FieldWriter::writeField(
   const string& a_field_name,
   const double* a_field,
   const ParallelArray::Box& a_in_box,
   const ParallelArray::Box& a_out_box,
   int a_n_ghosts)
{
   // If not processing a time slice then this is being used improperly.
   if (!m_processing_time_slice) {
      LOKI_ABORT("Time slice not being processed.  Missing call to startTimeSlice.");
   }

   // If the in and out boxes are the same then the field is not parallelized
   // and we can use a simpler scheme to write the field.  But first open the
   // dataset for the plot.
   bool same_boxes = a_in_box == a_out_box;
   ostringstream name;
   name << "time_slice_" << m_total_time_slices_written << "_" << a_field_name;
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
   hid_t dsetpp = H5Dopen(m_root, name.str().c_str());
#else
   hid_t dsetpp = H5Dopen(m_root, name.str().c_str(), H5P_DEFAULT);
#endif
#else
   hid_t dsetpp = H5Dopen(m_root, name.str().c_str());
#endif
   if (dsetpp < 0) {
      LOKI_ABORT("Unable to open field dataset.");
   }
   if (same_boxes) {
      // Write the field.
      herr_t errf = H5Dwrite(dsetpp,
         H5T_NATIVE_DOUBLE,
         H5S_ALL,
         H5S_ALL,
         H5P_DEFAULT,
         a_field);
   }
   else {
      // Specify the size and shape of this part of the field.
      hsize_t field_part_dim[2];
      field_part_dim[0] = a_out_box.numberOfCells(1);
      field_part_dim[1] = a_out_box.numberOfCells(0);
      hsize_t offset_pp[2];
      offset_pp[0] = a_out_box.lower(1) + a_n_ghosts;
      offset_pp[1] = a_out_box.lower(0) + a_n_ghosts;
      hsize_t stride[2];
      stride[0] = 1;
      stride[1] = 1;
      hsize_t dim_pp[2];
      dim_pp[0] = 1;
      dim_pp[1] = 1;
      hsize_t block_pp[2];
      block_pp[0] = a_out_box.numberOfCells(1);
      block_pp[1] = a_out_box.numberOfCells(0);

      // Create memory space with size of this part of the field.  Get the
      // dataspace and select this part of the field from the dataspace.
      hid_t memspace_id = H5Screate_simple(2, field_part_dim, NULL);
      if (memspace_id < 0) {
         LOKI_ABORT("Unable to create field memory space.");
      }
      hid_t dataspace_id = H5Dget_space(dsetpp);
      if (dataspace_id < 0) {
         LOKI_ABORT("Unable to get field dataspace.");
      }

      // Select the hyperslab in the dataspace.
      herr_t errf = H5Sselect_hyperslab(dataspace_id,
         H5S_SELECT_SET,
         offset_pp,
         stride,
         dim_pp,
         block_pp);
      if (errf < 0) {
         LOKI_ABORT("Unable to select field hyperslab.");
      }

      // Create a buffer containing the data specified by a_out_box from the
      // data in a_field described by a_in_box.
      double* data_out =
         new double [a_out_box.numberOfCells(0)*a_out_box.numberOfCells(1)];
      int nx_in = a_in_box.numberOfCells(0);
      int x_start = a_out_box.lower(0) - a_in_box.lower(0);
      int x_end = x_start + a_out_box.numberOfCells(0);
      int y_start = a_out_box.lower(1) - a_in_box.lower(1);
      int y_end = y_start + a_out_box.numberOfCells(1);
      int idx = 0;
      for (int j = y_start; j < y_end; ++j) {
         for (int i = x_start; i < x_end; ++i) {
            data_out[idx++] = a_field[j*nx_in+i];
         }
      }

      // Write this part of the field to the hyperslab in the dataspace.
      errf = H5Dwrite(dsetpp,
         H5T_NATIVE_DOUBLE,
         memspace_id,
         dataspace_id,
         H5P_DEFAULT,
         data_out);
      if (!same_boxes) {
         delete [] data_out;
      }
      if (errf < 0) {
         LOKI_ABORT("Unable to write to dataset for field.");
      }

      // Release resources.
      errf = H5Sclose(memspace_id);
      if (errf < 0) {
         LOKI_ABORT("Unable to close field memory space.");
      }
      errf = H5Sclose(dataspace_id);
      if (errf < 0) {
         LOKI_ABORT("Unable to close field dataspace.");
      }
      errf = H5Dclose(dsetpp);
      if (errf < 0) {
         LOKI_ABORT("Unable to close field dataset.");
      }
   }
}


void
FieldWriter::endTimeSlice(
   bool a_proc_lo)
{
   // If not processing a time slice then this is being used improperly.
   if (!m_processing_time_slice) {
      LOKI_ABORT("Time slice not being processed.  Missing call to startTimeSlice.");
   }

   ++m_time_slices_written_to_current_file;
   ++m_total_time_slices_written;
   m_processing_time_slice = false;

   // Write the total number of time slices to the 1st field file.
   if (a_proc_lo) {
      writeNumTimeSlices();
   }

   // Close the file and root group.
   closeFileAndRoot(m_file, m_root);
}


void
FieldWriter::writePlotTimes(
   const vector<double>& a_plot_times)
{
   // Create the dataspace for the plot times.
   hsize_t dims[1];
   dims[0] = static_cast<int>(a_plot_times.size());
   hid_t dspace = H5Screate_simple(1, dims, NULL);
   if (dspace < 0) {
      LOKI_ABORT("Unable to create dataspace for the plot times.");
   }

   // Write the plot times.
   writeDoubleArray("time", m_root, dspace, &a_plot_times[0]);

   // Close the dataspace for the plot times.
   herr_t errf = H5Sclose(dspace);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataspace for the plot times.");
   }
}


void
FieldWriter::writeCoords(
   vector<double>& a_x_coords,
   vector<double>& a_y_coords)
{
   // Create the dataspace for the x coordinates.
   hsize_t dims[1];
   dims[0] = static_cast<int>(a_x_coords.size());
   hid_t dspace = H5Screate_simple(1, dims, NULL);
   if (dspace < 0) {
      LOKI_ABORT("Unable to create dataspace for the x coordinates.");
   }

   // Write the x coordinates.
   writeDoubleArray("x", m_root, dspace, &a_x_coords[0]);

   // Close the dataspace for the x coordinates.
   herr_t errf = H5Sclose(dspace);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataspace for the x coordinates.");
   }

   // Write the y coordinates.
   // Create the dataspace for the y coordinates.
   dims[0] = static_cast<int>(a_y_coords.size());
   dspace = H5Screate_simple(1, dims, NULL);
   if (dspace < 0) {
      LOKI_ABORT("Unable to create dataspace for the y coordinates.");
   }

   // Write the y coordinates.
   writeDoubleArray("y", m_root, dspace, &a_y_coords[0]);

   // Close the dataspace for the y coordinates.
   errf = H5Sclose(dspace);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataspace for the y coordinates.");
   }
}


void
FieldWriter::writeFieldTimeSlice(
   const string& a_field_name,
   int a_nx,
   int a_ny,
   int a_which_time_slice,
   const double* a_field_data)
{
   // Open the dataset for the plot.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
   hid_t dsetpp = H5Dopen(m_root, a_field_name.c_str());
#else
   hid_t dsetpp = H5Dopen(m_root, a_field_name.c_str(), H5P_DEFAULT);
#endif
#else
   hid_t dsetpp = H5Dopen(m_root, a_field_name.c_str());
#endif
   if (dsetpp < 0) {
      LOKI_ABORT("Unable to open dataset for post processed plot.");
   }

   // Specify the size and shape of each time slice's field data.
   hsize_t time_slice_dim[2];
   time_slice_dim[0] = a_ny;
   time_slice_dim[1] = a_nx;
   hsize_t offset_pp[3];
   offset_pp[0] = a_which_time_slice;
   offset_pp[1] = 0;
   offset_pp[2] = 0;
   hsize_t stride[3];
   stride[0] = 1;
   stride[1] = 1;
   stride[2] = 1;
   hsize_t dim_pp[3];
   dim_pp[0] = 1;
   dim_pp[1] = 1;
   dim_pp[2] = 1;
   hsize_t block_pp[3];
   block_pp[0] = 1;
   block_pp[1] = a_ny;
   block_pp[2] = a_nx;

   // Create memory space with size of the field's time slice.  Get the
   // dataspace and select this time slice's subset from the dataspace.
   hid_t memspace_id = H5Screate_simple(2, time_slice_dim, NULL);
   if (memspace_id < 0) {
      LOKI_ABORT("Unable to create post processed field time slice memory space.");
   }
   hid_t dataspace_id = H5Dget_space(dsetpp);
   if (dataspace_id < 0) {
      LOKI_ABORT("Unable to get post processed field dataspace.");
   }

   // Select the hyperslab in the post processed dataspace and write this time
   // slice of this field to it.
   herr_t errf = H5Sselect_hyperslab(dataspace_id,
                                     H5S_SELECT_SET,
                                     offset_pp,
                                     stride,
                                     dim_pp,
                                     block_pp);
   if (errf < 0) {
      LOKI_ABORT("Unable to select post processed field hyperslab.");
   }
   errf = H5Dwrite(dsetpp,
                   H5T_NATIVE_DOUBLE,
                   memspace_id,
                   dataspace_id,
                   H5P_DEFAULT,
                   a_field_data);
   if (errf < 0) {
      LOKI_ABORT("Unable to write to dataset for post processed field.");
   }

   // Release resources.
   errf = H5Sclose(memspace_id);
   if (errf < 0) {
      LOKI_ABORT("Unable to close post processed field memory space.");
   }
   errf = H5Sclose(dataspace_id);
   if (errf < 0) {
      LOKI_ABORT("Unable to close post processsed field dataspace.");
   }
   errf = H5Dclose(dsetpp);
   if (errf < 0) {
      LOKI_ABORT("Unable to close post processed field dataset.");
   }
}


void
FieldWriter::writeNumTimeSlices()
{
   // Write the total number of time slices to the 1st field file.

   // Open the 1st field file.
   ostringstream field_file_name;
   field_file_name << m_base_name << ".fields_" << 0 << ".hdf";
   hid_t field_file = H5Fopen(field_file_name.str().c_str(),
      H5F_ACC_RDWR,
      H5P_DEFAULT);
   if (field_file < 0) {
      ostringstream os;
      os << "Unable to open field file "
         << field_file_name.str().c_str() << ".";
      LOKI_ABORT(os.str().c_str());
   }

   // Open the root group.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Gopen_vers) && H5Gopen_vers == 1
   hid_t group = H5Gopen(field_file, "root");
#else
   hid_t group = H5Gopen(field_file, "root", H5P_DEFAULT);
#endif
#else
   hid_t group = H5Gopen(field_file, "root");
#endif
   if (group < 0) {
      LOKI_ABORT("Can not open root group of field file.");
   }

   // Get the datset for the number of time slices.
   hid_t dspace;
   hid_t dset;
   if (m_total_time_slices_written == 1) {

      // Create the dataspace for the total number of time slices.
      hsize_t dims = 1;
      dspace = H5Screate_simple(1, &dims, NULL);
      if (dspace < 0) {
         LOKI_ABORT("Unable to create dataspace for total number of time slices.");
      }
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dcreate_vers) && H5Dcreate_vers == 1
      dset = H5Dcreate(group,
                       "total_num_time_slices",
                       H5T_STD_I32BE,
                       dspace,
                       H5P_DEFAULT);
#else
      dset = H5Dcreate(group,
                       "total_num_time_slices",
                       H5T_STD_I32BE,
                       dspace,
                       H5P_DEFAULT,
                       H5P_DEFAULT,
                       H5P_DEFAULT);
#endif
#else
      dset = H5Dcreate(group,
                       "total_num_time_slices",
                       H5T_STD_I32BE,
                       dspace,
                       H5P_DEFAULT);
#endif
   }
   else {
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
      dset = H5Dopen(group, "total_num_time_slices");
#else
      dset = H5Dopen(group, "total_num_time_slices", H5P_DEFAULT);
#endif
#else
      dset = H5Dopen(group, "total_num_time_slices");
#endif
   }
   if (dset < 0) {
      LOKI_ABORT("Unable to create dataset for total number of time slices.");
   }

   // Write the current total number of time slices.
   herr_t errf = H5Dwrite(dset,
      H5T_NATIVE_INT,
      H5S_ALL,
      H5S_ALL,
      H5P_DEFAULT,
      &m_total_time_slices_written);
   if (errf < 0) {
      LOKI_ABORT("Unable to write to dataset for total number of time slices.");
   }

   // Close everything.
   errf = H5Dclose(dset);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataset for total number of time slices.");
   }
   if (m_total_time_slices_written == 1) {
      errf = H5Sclose(dspace);
      if (errf < 0) {
         LOKI_ABORT("Unable to close dataspace for total number of time slices.");
      }
   }

   // Close the file and root group.
   closeFileAndRoot(field_file, group);

   // Write the number of time slices written to the current field file to that
   // file.

   // Get the datset for the number of time slices written to the current file.
   if (m_time_slices_written_to_current_file == 1) {

      // Create the dataspace for the number of time slices written to the
      // current file.
      hsize_t dims = 1;
      dspace = H5Screate_simple(1, &dims, NULL);
      if (dspace < 0) {
         LOKI_ABORT("Unable to create dataspace for number of time slices in thiis file.");
      }
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dcreate_vers) && H5Dcreate_vers == 1
      dset = H5Dcreate(m_root,
                       "num_time_slices_in_this_file",
                       H5T_STD_I32BE,
                       dspace,
                       H5P_DEFAULT);
#else
      dset = H5Dcreate(m_root,
                       "num_time_slices_in_this_file",
                       H5T_STD_I32BE,
                       dspace,
                       H5P_DEFAULT,
                       H5P_DEFAULT,
                       H5P_DEFAULT);
#endif
#else
      dset = H5Dcreate(m_root,
                       "num_time_slices_in_this_file",
                       H5T_STD_I32BE,
                       dspace,
                       H5P_DEFAULT);
#endif
   }
   else {
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dopen_vers) && H5Dopen_vers == 1
      dset = H5Dopen(m_root, "num_time_slices_in_this_file");
#else
      dset = H5Dopen(m_root, "num_time_slices_in_this_file", H5P_DEFAULT);
#endif
#else
      dset = H5Dopen(m_root, "num_time_slices_in_this_file");
#endif
   }
   if (dset < 0) {
      LOKI_ABORT("Unable to create dataset for number of time slices in this file.");
   }

   // Write the current number of time slices in this file.
   errf = H5Dwrite(dset,
      H5T_NATIVE_INT,
      H5S_ALL,
      H5S_ALL,
      H5P_DEFAULT,
      &m_time_slices_written_to_current_file);
   if (errf < 0) {
      LOKI_ABORT("Unable to write to dataset for number of time slices in this file.");
   }

   // Close everything.
   errf = H5Dclose(dset);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataset for number of time slices in this file.");
   }
   if (m_time_slices_written_to_current_file == 1) {
      errf = H5Sclose(dspace);
      if (errf < 0) {
         LOKI_ABORT("Unable to close dataspace for number of time slices in this file.");
      }
   }
}


void
FieldWriter::createFieldDatasets(
   const vector<string>& a_field_names,
   int a_nx,
   int a_ny,
   int a_time_slices_per_field,
   bool prepend_slice)
{
   int num_plots = static_cast<int>(a_field_names.size());
   for (int plot = 0; plot < num_plots; ++plot) {
      ostringstream name;
      if (prepend_slice) {
         name << "time_slice_" << m_total_time_slices_written << "_"
              << a_field_names[plot];
      }
      else {
         name << a_field_names[plot];
      }
      hid_t dspacepp;
      if (a_time_slices_per_field == 1) {
         hsize_t dimspp[2];
         dimspp[0] = a_ny;
         dimspp[1] = a_nx;
         dspacepp = H5Screate_simple(2, dimspp, NULL);
      }
      else {
         hsize_t dimspp[3];
         dimspp[0] = a_time_slices_per_field;
         dimspp[1] = a_ny;
         dimspp[2] = a_nx;
         dspacepp = H5Screate_simple(3, dimspp, NULL);
      }
      if (dspacepp < 0) {
         LOKI_ABORT("Unable to create dataspace for plot.");
      }

      // Create the dataset for the plot.
#if (((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8) && (H5_VERS_RELEASE>=0)) || \
     ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR>8)) || \
     (H5_VERS_MAJOR>1))
#if defined (H5Dcreate_vers) && H5Dcreate_vers == 1
      hid_t dsetpp = H5Dcreate(m_root,
                               name.str().c_str(),
                               H5T_NATIVE_DOUBLE,
                               dspacepp,
                               H5P_DEFAULT);
#else
      hid_t dsetpp = H5Dcreate(m_root,
                               name.str().c_str(),
                               H5T_NATIVE_DOUBLE,
                               dspacepp,
                               H5P_DEFAULT,
                               H5P_DEFAULT,
                               H5P_DEFAULT);
#endif
#else
      hid_t dsetpp = H5Dcreate(m_root,
                               name.str().c_str(),
                               H5T_NATIVE_DOUBLE,
                               dspacepp,
                               H5P_DEFAULT);
#endif
      if (dsetpp < 0) {
         LOKI_ABORT("Unable to create dataset for plot.");
      }

      // Close the dataset and dataspace for the plot.
      herr_t errf = H5Dclose(dsetpp);
      if (errf < 0) {
         LOKI_ABORT("Unable to close dataset for plot.");
      }
      errf = H5Sclose(dspacepp);
      if (errf < 0) {
         LOKI_ABORT("Unable to close dataspace for plot.");
      }
   }
}
} // end namespace Loki
