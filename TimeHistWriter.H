/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _TIME_HIST_WRITER_H_
#define _TIME_HIST_WRITER_H_

#include "ReaderWriterBase.H"

namespace Loki {

/**
 * A simple class to write Loki 0D time history data.  The EM processor will
 * write the time history data.  All time history data will be written to a
 * series of files named base_name.time_hists_n.hdf where n is the index in the
 * series of files.  This class uses the serial HDF5 interface as only the EM
 * processor or the serial postprocessor write time history data.
 */
class TimeHistWriter : public ReaderWriterBase
{
public:
   /**
    * @brief Constructor of time history writer for collision diagnostics.
    *
    * @param[in] a_name The name of the time history file written by this
    *                   writer.
    * @param[in] a_num_seq_times The number of sample times.
    * @param[in] a_seq_times The sample times of each time history.
    */
   TimeHistWriter(
      const string& a_name,
      int a_num_seq_times,
      const vector<double>& a_seq_times);

   /**
    * @brief Constructor of a time history writer.
    *
    * @param[in] a_name The name of the time history file written by this
    *                   writer.
    * @param[in] a_num_probes The number of probes in the problem.
    * @param[in] a_num_tracking_particles The number of tracking particles.
    * @param[in] a_num_seq_times The number of sample times.
    * @param[in] a_seq_times The sample times of each time history.
    */
   TimeHistWriter(
      const string& a_name,
      int a_num_probes,
      int a_num_tracking_particles,
      int a_num_seq_times,
      const vector<double>& a_seq_times);

   /**
    * @brief Constructor of time history writer for post processing.
    *
    * @param[in] a_name The name of the time history file written by this
    *                   writer.
    * @param[in] a_seq_times The sample times of each time history.
    */
   TimeHistWriter(
      const string& a_name,
      const vector<double>& a_seq_times);

   /**
    * @brief Destructor.
    */
   virtual
   ~TimeHistWriter();

   /**
    * @brief Write a time history.
    *
    * @param[in] a_time_hist_name Name of time history being written.
    * @param[in] a_num_seq_times Number of time points in time history.
    * @param[in] a_time_hist Time history being written.
    */
   void
   writeTimeHistory(
      const string& a_time_hist_name,
      int a_num_seq_times,
      const vector<double>& a_time_hist);

   /**
    * @brief Write a time history.
    *
    * @param[in] a_time_hist_name Name of time history being written.
    * @param[in] a_time_hist Time history being written.
    */
   void
   writeTimeHistory(
      const string& a_time_hist_name,
      const vector<double>& a_time_hist);

private:
   // Unimplemented default constructor.
   TimeHistWriter();

   // Unimplemented copy constructor.
   TimeHistWriter(
      const TimeHistWriter& a_other);

   // Unimplemented assignment operator.
   TimeHistWriter&
   operator = (
      const TimeHistWriter& a_rhs);

   // Does all work common to both constructors.
   void
   doCommonConstruction(
      const string& a_name,
      int a_num_seq_times,
      const double* a_seq_times,
      bool a_for_post_proc);

   // The file being written to.
   hid_t m_file;

   // The root group in the file being written to.
   hid_t m_root;
};

} // end namespace Loki

#endif
