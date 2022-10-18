/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _TIME_HIST_READER_H_
#define _TIME_HIST_READER_H_

#include "ReaderWriterBase.H"
#include <string>
#include <vector>

using namespace std;

namespace Loki {

/**
 * A simple class to read Loki 0D time history data from a time history file
 * produced by a TimeHistWriter.  It uses the serial HDF5 interface as only
 * the serial postprocessor and test harness read time history data.
 */
class TimeHistReader : public ReaderWriterBase
{
public:
   /**
    * @brief Constructor of a time history reader.
    *
    * @param[in] a_name The name of the time history file read by this reader.
    */
   TimeHistReader(
      const string& a_name);

   /**
    * @brief Destructor.
    */
   virtual
   ~TimeHistReader();

   /**
    * @brief Read a time history from the time history file.
    *
    * @param[in] a_name Name of time history.
    * @param[out] a_vals Repository for time history.
    */
   void
   readTimeHistory(
      const string& a_name,
      vector<double>& a_vals);

   /**
    * @brief Read the number of probes.
    *
    * @param[out] a_num_probes The number of probes.
    */
   void
   readNumProbes(
      int& a_num_probes);

   /**
    * @brief Read the number of tracking particles.
    *
    * @param[out] a_num_tracking_particles The number of tracking particles.
    */
   void
   readNumTrackingParticles(
      int& a_num_tracking_particles);

private:
   // Unimplemented default constructor.
   TimeHistReader();

   // Unimplemented copy constructor.
   TimeHistReader(
      const TimeHistReader& a_other);

   // Unimplemented assignment operator.
   TimeHistReader&
   operator = (
      const TimeHistReader& a_rhs);

   // The file being read from.
   hid_t m_file;

   // The root group in the file being read from.
   hid_t m_root;
};

} // end namespace Loki

#endif
