/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _READER_WRITER_BASE_H_
#define _READER_WRITER_BASE_H_

#include "hdf5.h"
#include <string>
#include <vector>

using namespace std;

namespace Loki {

/**
 * A class to encapsulate implementation details common to Loki file readers and
 * writers.  All derived classes are specialized readers and writers for
 * specific types of Loki data: fields, time histories, and restart data.
 */
class ReaderWriterBase
{
public:
protected:
   // Default constructor.
   ReaderWriterBase();

   // Destructor.
   virtual
   ~ReaderWriterBase();

   // Create the file to be read from/written to and its root group.
   void
   createFileAndRoot(
      const string& a_file_name,
      hid_t& a_file,
      hid_t& a_root,
      bool a_parallel_io = false);

   // Open an existing file to be read from/written to and its root group.
   void
   openFileAndRoot(
      const string& a_file_name,
      hid_t& a_file,
      hid_t& a_root);

   // Close the file to be read from/written to and its root group.
   void
   closeFileAndRoot(
      hid_t& a_file,
      hid_t& a_root);

   // Create the named group embedded in the supplied home group or file.
   void
   createGroup(
      const string& a_group_name,
      hid_t& a_group,
      hid_t& a_home);

   // Open the named group embedded in the supplied home group or file.
   void
   openGroup(
      const string& a_group_name,
      hid_t& a_group,
      hid_t& a_home);

   // Close the supplied group.
   void
   closeGroup(
      hid_t& a_group);

   // Write an array of doubles to the supplied home group or file.
   void
   writeDoubleArray(
      const string& a_name,
      hid_t a_home,
      hid_t a_dspace,
      const double* a_vals);

   // Write a single double to the supplied home group or file.
   void
   writeDoubleValue(
      const string& a_name,
      hid_t a_home,
      double a_val);

   // Write an array of integers to the supplied home group or file.
   void
   writeIntegerArray(
      const string& a_name,
      hid_t a_home,
      hid_t a_dspace,
      const int* a_vals);

   // Write a single integer to the supplied home group or file.
   void
   writeIntegerValue(
      const string& a_name,
      hid_t a_home,
      int a_val);

   // Read an array of doubles from the supplied home group or file.
   void
   readDoubleArray(
      const string& a_name,
      hid_t a_home,
      vector<double>& a_vals);

   // Read a single double from the supplied home group or file.
   void
   readDoubleValue(
      const string& a_name,
      hid_t a_home,
      double& a_val);

   // Read a single integer from the supplied home group or file.
   void
   readIntegerValue(
      const string& a_name,
      hid_t a_home,
      int& a_val);

private:
   // Unimplemented copy constructor.
   ReaderWriterBase(
      const ReaderWriterBase& a_other);

   // Unimplemented assignment operator.
   ReaderWriterBase&
   operator = (
      const ReaderWriterBase& a_rhs);
};

} // end namespace Loki

#endif
