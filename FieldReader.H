/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _FIELD_READER_H_
#define _FIELD_READER_H_

#include "ReaderWriterBase.H"

using namespace std;

namespace Loki {

/**
 * A simple class to read Loki 2D field data as written by FieldWriter.  It
 * uses the serial HDF5 interface as only the serial postprocessor and test
 * harness read field data.
 */
class FieldReader : public ReaderWriterBase
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_name The name of the field file.
    */
   FieldReader(
      const string& a_name);

   /**
    * @brief Destructor.
    */
   virtual
   ~FieldReader();

   /**
    * @brief Read the time for the named time slice.
    *
    * @param[in] a_name Name of the time for this time slice.
    * @param[out] a_time The time for the named time slice.
    */
   void
   readTime(
      const string& a_name,
      double& a_time);

   /**
    * @brief Read the coordinates.
    *
    * @param[out] a_x_coords The x coordinates.
    * @param[out] a_y_coords The y coordinates.
    */
   void
   readCoords(
      vector<double>& a_x_coords,
      vector<double>& a_y_coords);

   /**
    * @brief Read a field for the current time slice.
    *
    * @param[in] a_field_name Name of field being read.
    * @param[out] a_field Field being read.
    */
   void
   readField(
      const string& a_field_name,
      vector<double>& a_field);

   /**
    * @brief Read the total number of time slices.
    *
    * @param[out] a_total_num_time_slices The total number of time slices.
    */
   void
   readTotalNumTimeSlices(
      int& a_total_num_time_slices);

   /**
    * @brief Read the number of time slices in this field file.
    *
    * @param[out] a_num_time_slices_in_file The number of time slices in this
    *                                       field file.
    */
   void
   readNumTimeSlicesInFile(
      int& a_num_time_slices_in_file);

private:
   // Unimplemented default constructor.
   FieldReader();

   // Unimplemented copy constructor.
   FieldReader(
      const FieldReader& a_other);

   // Unimplemented assignment operator.
   FieldReader&
   operator = (
      const FieldReader& a_rhs);

   // The file being read from.
   hid_t m_file;

   // The root group in the file being read from.
   hid_t m_root;
};

} // end namespace Loki

#endif
