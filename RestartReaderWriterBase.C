/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "RestartReaderWriterBase.H"
#include "Loki_Defines.H"
#include "Loki_Utilities.H"

namespace Loki {

RestartReaderWriterBase::RestartReaderWriterBase(
   const string& a_base_name,
   int a_max_num_files,
   bool a_post_processing)
   : ReaderWriterBase(),
     m_base_name(a_base_name),
     m_post_processing(a_post_processing)
{
   // Set some data members.
   if (m_post_processing) {
      if (Loki_Utilities::s_num_procs != 1) {
         LOKI_ABORT("Post processor must be run serially.");
      }
      m_max_num_files = a_max_num_files;
   }
   else {
      m_max_num_files = min(a_max_num_files, Loki_Utilities::s_num_procs);
   }

   // Figure out if this processor is the first or last to read from its
   // bulkdata file.
   if (Loki_Utilities::s_my_id == 0 ||
       (Loki_Utilities::s_my_id-1)*m_max_num_files/Loki_Utilities::s_num_procs !=
        Loki_Utilities::s_my_id*m_max_num_files/Loki_Utilities::s_num_procs) {
      m_first_reader_writer = true;
   }
   else {
      m_first_reader_writer = false;
   }
   if (Loki_Utilities::s_my_id == Loki_Utilities::s_num_procs-1 ||
       Loki_Utilities::s_my_id*m_max_num_files/Loki_Utilities::s_num_procs !=
       (Loki_Utilities::s_my_id+1)*m_max_num_files/Loki_Utilities::s_num_procs) {
      m_last_reader_writer = true;
   }
   else {
      m_last_reader_writer = false;
   }
}


RestartReaderWriterBase::~RestartReaderWriterBase()
{
   while (m_metadata_groups.size() != 1) {
      closeGroup(m_metadata_groups.top());
      m_metadata_groups.pop();
      m_metadata_group_names.pop();
   }
   closeFileAndRoot(m_metadata_file, m_metadata_groups.top());
   m_metadata_groups.pop();
   m_metadata_group_names.pop();
}


void
RestartReaderWriterBase::pushSubDir(
   const string& a_name,
   bool a_create_group)
{
   // Open the named group as a group in the metadata file's currently most
   // deeply nested group.
   hid_t& cur_group = m_metadata_groups.top();
   m_metadata_groups.push(0);
   m_metadata_group_names.push(a_name);
   if (a_create_group) {
      createGroup(a_name, m_metadata_groups.top(), cur_group);
   }
   else {
      openGroup(a_name, m_metadata_groups.top(), cur_group);
   }
}


void
RestartReaderWriterBase::popSubDir()
{
   // Free the resources associated with the most nested metadata file group.
   // The root group may not be popped until the metadata file is closed.
   if (m_metadata_groups.empty()) {
      LOKI_ABORT("No groups to pop.");
   }
   else if (m_metadata_groups.size() == 1) {
      LOKI_ABORT("Attempting to pop the root group.");
   }
   else {
      closeGroup(m_metadata_groups.top());
      m_metadata_groups.pop();
      m_metadata_group_names.pop();
   }
}
} // end namespace Loki
