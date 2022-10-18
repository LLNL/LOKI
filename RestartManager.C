/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "RestartManager.H"
#include "RestartReader.H"
#include "RestartWriter.H"
#include "Loki_Defines.H"
#include "Loki_Utilities.H"
// Allows use of 'mkdir'
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>

namespace Loki {

RestartManager* RestartManager::s_restart_manager_instance = 0;

RestartManager*
RestartManager::getManager()
{
   // See if the singleton exists.  If not then build it.
   if (s_restart_manager_instance == 0) {
      s_restart_manager_instance = new RestartManager();
   }
   return s_restart_manager_instance;
}


RestartManager::RestartManager()
   :  m_step_interval(-1),
      m_time_interval(-1.0),
      m_next_write_time(0.0),
      m_write_on_time(false),
      m_restart_write_path("."),
      m_restart_read_path("."),
      m_is_from_restart(false),
      m_restart_index(0),
      m_max_restart_files(16)
{
   // Parse to check if restart is required.  The "start_from_restart" input is
   // not part of the "restart" sub-database for some reason.
   string start_from_restart("false");
   LokiInputParser pp;
   pp.query("start_from_restart", start_from_restart);
   m_is_from_restart = start_from_restart.compare("true") == 0 ? true : false;
   pp.query("max_files_for_write", m_max_restart_files);

   // Now read the other restart related input from the "restart" sub-database.
   LokiInputParser restart_pp("restart");
   parseParameters(restart_pp);

   // If we're running a restarted simulation then we must find the index of the
   // last restart file written.  We'll pick up the simulation from there.
   if (m_is_from_restart) {
      findRestartIndex();
   }
}


void
RestartManager::resetNextWriteTime(
   double a_time)
{
   if (m_write_on_time) {
      m_next_write_time = a_time;
   }
}


void
RestartManager::write(
   double a_time)
{
   // Glue together the info we already have to form the name of the file we're
   // going to write to and open it.
   ostringstream restart_filename;
   createFileName(restart_filename, m_restart_write_path, m_restart_index++);
   Loki_Utilities::printF("Writing to %s\n", restart_filename.str().c_str());

   // See if the path to the restart file exists.  If not, make it.
   if (access(m_restart_write_path.c_str(), F_OK) != 0) {
      if (mkdir(m_restart_write_path.c_str(), S_IRWXU|S_IRGRP|S_IXGRP) != 0) {
         // Don't bother with file already exists error which is not an error
         // and can occur nproc-1 times.
         if (errno != EEXIST) {
            perror("mkdir() error");
         }
      }
      else {
         puts("Created new restart directory!");
      }
      // TODO: Are we concerned with portability??
   }

   // Open the restart file.
   RestartWriter db_writer(restart_filename.str(), m_max_restart_files);

   // Loop over vector of Serializable items and have each save their state to
   // the restart file.
   for (SerializableList::iterator iter(m_items.begin()); 
        iter != m_items.end();
        ++iter) {
     (*iter)->putToRestart(db_writer, a_time);
   }

   // Update next time to write a restart file if we're doing time based
   // checkpointing.
   if (m_write_on_time) { // time
      m_next_write_time += m_time_interval;
   }
}


void
RestartManager::print()
{
   // Print info on restart frequency.
   Loki_Utilities::printF("\n********************************************************\n");
   Loki_Utilities::printF("[ Restart Manager status ]:\n");

   if (m_is_from_restart) {
      Loki_Utilities::printF("\tRunning from restart\n");
      Loki_Utilities::printF("\tRestart read directory path: %s\n",
         m_restart_read_path.c_str());
      Loki_Utilities::printF("\tRestart file index: %i\n", m_restart_index);
   }
   else {
      Loki_Utilities::printF("\tRunning from initial conditions\n");
   }
   Loki_Utilities::printF("\tRestart write directory path: %s\n",
          m_restart_write_path.c_str());

   if (m_write_on_time) {
      Loki_Utilities::printF("\t  Restart every %e time units\n",
         m_time_interval);
      Loki_Utilities::printF("\t  Next write time: %e\n", m_next_write_time);
   }
   else {
      Loki_Utilities::printF("\t  Restart every %i time steps\n",
         m_step_interval);
   }
}


void
RestartManager::parseParameters(
   LokiInputParser& a_pp)
{
   // Logic to take either time or step interval.  You can specify a restart
   // write frequency based on time steps or simulation time but not both.  One
   // or the other must be specified.
   bool has_step_interval = a_pp.contains("step_interval");
   bool has_time_interval = a_pp.contains("time_interval");
   if (has_step_interval && has_time_interval) {
      LOKI_ABORT("Must set restart frequency for only one of steps or time, not both.");
   }
   else if (!has_step_interval && !has_time_interval) {
      LOKI_ABORT("Must set restart frequency for one of steps or time.");
   }
   else {
      a_pp.query("step_interval", m_step_interval);
      a_pp.query("time_interval", m_time_interval);
      // Get the actual intervals and check for validity.  There is an
      // additional flag to indicate which scheme is being used.
      if (has_step_interval) {
         if (m_step_interval <= 0) {
            LOKI_ABORT("Input step_interval must be > 0");
         }
      }
      else {
         if (m_time_interval <= 0.0) {
            LOKI_ABORT("Input time_interval must be > 0.0");
         }
         else {
            m_write_on_time = true;
         }
      }
   }

   string tmp1;
   a_pp.query("write_directory", tmp1);
   m_restart_write_path = tmp1;
   string tmp2;
   a_pp.query("read_directory", tmp2);
   m_restart_read_path = tmp2;
}


void
RestartManager::restore()
{
   if (m_restart_index == 0) {
      LOKI_ABORT("No distributions found from which to restart ... quitting");
   }

   // Glue together the info we already have to form the name of the file we'll
   // read and open it.
   ostringstream restart_filename;
   createFileName(restart_filename, m_restart_read_path, m_restart_index);
   Loki_Utilities::printF("Restarting from %s\n",
      restart_filename.str().c_str());

   RestartReader db_reader(restart_filename.str(), m_max_restart_files);

   // Loop over vector of Serializable items and have each restore their state
   // from the restart file.
   for (SerializableList::iterator iter(m_items.begin());
        iter != m_items.end();
        ++iter) {
      (*iter)->getFromRestart(db_reader);
   }

   // The restart index does double duty.  When restoring from a restart it
   // is the index of the restart file to restore from.  After that it is the
   // index of the next restart file to write.  So reset it here, after
   // restoring from a restart.
   m_restart_index = 0;
}


void
RestartManager::findRestartIndex()
{
   // Given the restart directory read path supplied by the user, look for the
   // last restart file that has been written and set m_restart_index to that
   // file's index.
   int i = 1;
   while (true) {
      ostringstream restart_filename;
      createFileName(restart_filename, m_restart_read_path, i);
      struct stat stat_buf;
      if (stat(restart_filename.str().c_str(), &stat_buf) == 0) {
         m_restart_index = i;
         ++i;
         restart_filename.str("");
      }
      else {
         break;
      }
   }
}


RestartManager::~RestartManager()
{
}

} // end namespace Loki
