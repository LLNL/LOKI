/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _RESTART_READER_WRITER_BASE_H_
#define _RESTART_READER_WRITER_BASE_H_

#include "ReaderWriterBase.H"

#include <stack>

using namespace std;

namespace Loki {

/**
 * A class to encapsulate implementation details common to Loki restart file
 * readers and writers.
 */
class RestartReaderWriterBase : public ReaderWriterBase
{
public:
   /**
    * @brief Close the most nested group and return the file write context to
    *        the next less nested group.
    */
   void
   popSubDir();

protected:
   // Constructor.
   RestartReaderWriterBase(
      const string& a_base_name,
      int a_max_num_files,
      bool a_post_processing);

   // Destructor.
   virtual
   ~RestartReaderWriterBase();

   // Create or open the named nested group within the metadata file.  Any data
   // write is placed in the currently most nested group.
   void
   pushSubDir(
      const string& a_name,
      bool a_create_group);

   // The metadata file.
   hid_t m_metadata_file;

   // The nested groups of the metadata file.
   stack<hid_t> m_metadata_groups;

   // The names of the nested groups of the metadata file.
   stack<string> m_metadata_group_names;

   // The bulkdata file.
   hid_t m_bulkdata_file;

   // The bulkdata root group.
   hid_t m_bulkdata_root;

   // The maximum number of bulkdata files to read.
   int m_max_num_files;

   // The base name of the bulk and metadata files.
   string m_base_name;

   // If true this is the first processor to read from or write to its bulkdata
   // file.
   bool m_first_reader_writer;

   // If true this is the last processor to read from or write to its bulkdata
   // file.
   bool m_last_reader_writer;

   // If true then this is used in the context of the serial post processor.
   bool m_post_processing;

private:
   // Unimplemented default constructor.
   RestartReaderWriterBase();

   // Unimplemented copy constructor.
   RestartReaderWriterBase(
      const RestartReaderWriterBase& a_other);

   // Unimplemented assignment operator.
   RestartReaderWriterBase&
   operator = (
      const RestartReaderWriterBase& a_rhs);
};

} // end namespace Loki

#endif
