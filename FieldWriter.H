/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _FIELD_WRITER_H_
#define _FIELD_WRITER_H_

#include "ReaderWriterBase.H"
#include "ProblemDomain.H"
#include "tbox/IntVector.H"

namespace Loki {

/**
 * A simple class to write Loki 2D field data.  The EM processor will write the
 * field data.  All field data will be written to a series of files named
 * base_name.fields_n.hdf where n is the index in the series of files.  It is
 * possible to write multiple time slices of field data to a single field file.
 * Once a given field file has this many time slices written to it the next time
 * slice will be written to the next field file in the series.  This class uses
 * the serial HDF5 interface as only the EM processor or the serial
 * postprocessor write field data.
 */
class FieldWriter : public ReaderWriterBase
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_base_name The base name common to all field files.
    * @param[in] a_domain The problem domain.
    * @param[in] a_solution_order The spatial solution order.
    * @param[in] a_time_slices_per_file The maximum number of time slices to be
    *                                   written to a given field file.
    */
   FieldWriter(
      const string& a_base_name,
      const ProblemDomain& a_domain,
      int a_solution_order,
      int a_time_slices_per_file);

   /**
    * @brief Constructor for post processing.
    *
    * @param[in] a_name The name of the post processed field file.
    * @param[in] a_field_names The names of all the fields.
    * @param[in] a_nx The number of cells in x.
    * @param[in] a_ny The number of cells in y.
    * @param[in] a_time_slices_per_field The number of field time slices.
    */
   FieldWriter(
      const string& a_name,
      const vector<string>& a_field_names,
      int a_nx,
      int a_ny,
      int a_time_slices_per_field);

   /**
    * @brief Destructor.
    */
   virtual
   ~FieldWriter();

   /**
    * @brief Start a time slice in the current field file for collision
    *        diagnostic fields.
    *
    * @param[in] a_time Simulation time associated with this time slice.
    * @param[in] a_dt Simulation time step associated with this time slice.
    * @param[in] a_field_names The names of all the fields.
    * @param[in] a_proc_lo True if used by the lowest rank RM processor.
    */
   void
   startTimeSlice(
      double a_time,
      double a_dt,
      const vector<string>& a_field_names,
      bool a_proc_lo);

   /**
    * @brief Start a time slice in the current field file.
    *
    * @param[in] a_time Simulation time associated with this time slice.
    * @param[in] a_dt Simulation time step associated with this time slice.
    * @param[in] a_num_tracking_particles Number of tracking particles in the
    *                                     problem.
    * @param[in] a_num_probes Number of probes in the problem.
    * @param[in] a_probes The probes.
    * @param[in] a_num_cells The number of cells in each dimension.
    * @param[in] a_field_names The names of all the fields.
    * @param[in] a_proc_lo True if used by the lowest rank EM processor.
    */
   void
   startTimeSlice(
      double a_time,
      double a_dt,
      int a_num_tracking_particles,
      int a_num_probes,
      const vector<vector<double> >& a_probes,
      const tbox::IntVector& a_num_cells,
      const vector<string>& a_field_names,
      bool a_proc_lo);

   /**
    * @brief Write a field for the current time slice.
    *
    * @param[in] a_field_name Name of field being written.
    * @param[in] a_field Field being written.
    * @param[in] a_in_box The box describing the extend of the field.
    * @param[in] a_out_box The box describing the part of the field to be
    *                      written.
    * @param[in] a_n_ghosts The number of ghost cells.
    */
   void
   writeField(
      const string& a_field_name,
      const double* a_field,
      const ParallelArray::Box& a_in_box,
      const ParallelArray::Box& a_out_box,
      int a_n_ghosts);

   /**
    * @brief End a time slice.
    *
    * @param[in] a_proc_lo True if used by the lowest rank EM processor.
    */
   void
   endTimeSlice(
      bool a_proc_lo);

   /**
    * @brief Write the plot times to the post processed field file.
    *
    * @param[in] a_plot_times The plot times to write.
    */
   void
   writePlotTimes(
      const vector<double>& a_plot_times);

   /**
    * @brief Write the coordinates to the post processed field file.
    *
    * @param[in] a_x_coords The x coordinates.
    * @param[in] a_y_coords The y coordinates.
    */
   void
   writeCoords(
      vector<double>& a_x_coords,
      vector<double>& a_y_coords);

   /**
    * @brief Write a time slice of a field to the post processed field file.
    */
   void
   writeFieldTimeSlice(
      const string& a_field_name,
      int a_nx,
      int a_ny,
      int a_which_time_slice,
      const double* a_field_data);

private:
   // Unimplemented default constructor.
   FieldWriter();

   // Unimplemented copy constructor.
   FieldWriter(
      const FieldWriter& a_other);

   // Unimplemented assignment operator.
   FieldWriter&
   operator = (
      const FieldWriter& a_rhs);

   // Write the total number of time slices that have been written to field
   // file(s) so far.
   void
   writeNumTimeSlices();

   // Create but do not write to the datasets for the fields.
   void
   createFieldDatasets(
      const vector<string>& a_field_names,
      int a_nx,
      int a_ny,
      int a_time_slices_per_field,
      bool prepend_slice);

   // The file being written to.
   hid_t m_file;

   // The root group in the file being written to.
   hid_t m_root;

   // The base name common to all field files.
   string m_base_name;

   // The index of the current file being written to.
   int m_current_file_index;

   // The number of time slices written to a single field file.
   int m_time_slices_per_file;

   // The number of time slices already written to the current file.
   int m_time_slices_written_to_current_file;

   // The total number of time slices written to all files.
   int m_total_time_slices_written;

   // True if currently processing a time slice.
   bool m_processing_time_slice;

   // The x dimension of fields being written.
   int m_nx;

   // The y dimension of fields being written.
   int m_ny;

   // The x coordinates.
   vector<double> m_x_coords;

   // The y coordinates.
   vector<double> m_y_coords;
};

} // end namespace Loki

#endif
