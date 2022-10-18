/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "TimeHistWriter.H"
#include "Loki_Defines.H"

namespace Loki {

TimeHistWriter::TimeHistWriter(
   const string& a_name,
   int a_num_seq_times,
   const vector<double>& a_seq_times)
   : ReaderWriterBase()
{
   doCommonConstruction(a_name, a_num_seq_times, &a_seq_times[0], false);
}


TimeHistWriter::TimeHistWriter(
   const string& a_name,
   int a_num_probes,
   int a_num_tracking_particles,
   int a_num_seq_times,
   const vector<double>& a_seq_times)
   : ReaderWriterBase()
{
   doCommonConstruction(a_name, a_num_seq_times, &a_seq_times[0], false);

   // Write the number of probes and the number of tracking particles.
   writeIntegerValue("numProbes", m_root, a_num_probes);
   writeIntegerValue("numTrackingParticles", m_root, a_num_tracking_particles);
}

TimeHistWriter::TimeHistWriter(
   const string& a_name,
   const vector<double>& a_seq_times)
   : ReaderWriterBase()
{
   doCommonConstruction(a_name,
      static_cast<int>(a_seq_times.size()),
      &a_seq_times[0],
      true);
}


TimeHistWriter::~TimeHistWriter()
{
   // Close the file and root group.
   closeFileAndRoot(m_file, m_root);
}


void
TimeHistWriter::writeTimeHistory(
   const string& a_time_hist_name,
   int a_num_seq_times,
   const vector<double>& a_time_hist)
{
   // Create the dataspace for the time history.
   hsize_t dims[1];
   dims[0] = a_num_seq_times;
   hid_t dspace = H5Screate_simple(1, dims, NULL);
   if (dspace < 0) {
      LOKI_ABORT("Unable to create dataspace for the time history.");
   }

   writeDoubleArray(a_time_hist_name, m_root, dspace, &a_time_hist[0]);

   // Close the dataspace for the time history.
   herr_t errf = H5Sclose(dspace);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataspace for time history.");
   }
}


void
TimeHistWriter::writeTimeHistory(
   const string& a_time_hist_name,
   const vector<double>& a_time_hist)
{
   // Create the dataspace for the time history.
   hsize_t dims[1];
   dims[0] = static_cast<int>(a_time_hist.size());
   hid_t dspace = H5Screate_simple(1, dims, NULL);
   if (dspace < 0) {
      LOKI_ABORT("Unable to create dataspace for the time history.");
   }

   writeDoubleArray(a_time_hist_name,
      m_root,
      dspace,
      &a_time_hist[0]);

   // Close the dataspace for the time history.
   herr_t errf = H5Sclose(dspace);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataspace for time history.");
   }
}


void
TimeHistWriter::doCommonConstruction(
   const string& a_name,
   int a_num_seq_times,
   const double* a_seq_times,
   bool a_for_post_proc)
{
   // Create the file and root.
   createFileAndRoot(a_name, m_file, m_root);

   // Create the dataspace for the sequence times.
   hsize_t dims[1];
   dims[0] = a_num_seq_times;
   hid_t dspace = H5Screate_simple(1, dims, NULL);
   if (dspace < 0) {
      LOKI_ABORT("Unable to create dataspace for the sequence times.");
   }

   // Write the sequence times to the file.
   if (a_for_post_proc) {
      writeDoubleArray("series_time", m_root, dspace, a_seq_times);
   }
   else {
      writeDoubleArray("sequence_times", m_root, dspace, a_seq_times);
   }

   // Close the dataspace for the time history.
   herr_t errf = H5Sclose(dspace);
   if (errf < 0) {
      LOKI_ABORT("Unable to close dataspace for time history.");
   }
}
} // end namespace Loki
