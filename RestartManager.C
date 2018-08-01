/*************************************************************************
 *
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * Written by Jeffrey Banks banksj3@rpi.edu (Rensselaer Polytechnic Institute,
 * Amos Eaton 301, 110 8th St., Troy, NY 12180); Jeffrey Hittinger
 * hittinger1@llnl.gov, William Arrighi arrighi2@llnl.gov, Richard Berger
 * berger5@llnl.gov, Thomas Chapman chapman29@llnl.gov (LLNL, P.O Box 808,
 * Livermore, CA 94551); Stephan Brunner stephan.brunner@epfl.ch (Ecole
 * Polytechnique Federale de Lausanne, EPFL SB SPC-TH, PPB 312, Station 13,
 * CH-1015 Lausanne, Switzerland).
 * CODE-744849
 *
 * All rights reserved.
 *
 * This file is part of Loki.  For details, see.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 ************************************************************************/
#include "RestartManager.H"
#include "Loki_Utilities.H"
// Allows use of 'mkdir'
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
      m_restart_index(0)
{
  // Parse to check if restart is required.  The "start_from_restart" input is
  // not part of the "restart" sub-database for some reason.
  aString start_from_restart("false");
  ParmParse pp;
  pp.query("start_from_restart", start_from_restart);
  m_is_from_restart = start_from_restart.matches("true") ? true : false;

  // Now read the other restart related input from the "restart" sub-database.
  ParmParse restart_pp("restart");
  parseParameters(restart_pp);

  // If we're running a restarted simulation then we must find the index of the
  // last restart file written.  We'll pick up the simulation from there.
  if (m_is_from_restart) {
    findRestartIndex();
  }
}


void
RestartManager::resetNextWriteTime(
   real a_time)
{
  if (m_write_on_time) {
    m_next_write_time = a_time;
  }
}


void
RestartManager::write(real a_time)
{
   // Glue together the info we already have to form the name of the file we're
   // going to write to and open it.
   char restart_filename[100];
   createFileName(restart_filename,
      m_restart_write_path.c_str(),
      m_restart_index++);
   printF("Writing to %s\n", restart_filename);

   // See if the path to the restart file exists.  If not, make it.
   if (access(m_restart_write_path.c_str(), F_OK) != 0) {
      if (mkdir(m_restart_write_path, S_IRWXU|S_IRGRP|S_IXGRP) != 0) {
         perror("mkdir() error");
      }
      else {
         puts("Created new restart directory!");
      }
      // TODO: Are we concerned with portability??
   }

   // Open the restart file.
   HDF_DataBase db_restart;
   db_restart.mount(restart_filename, "I"); // mount, I=Initialize

   // Loop over vector of Serializable items and have each save their state to
   // the restart file.
   for (SerializableList::iterator iter(m_items.begin()); 
        iter != m_items.end();
        ++iter) {
     (*iter)->putToRestart(db_restart, a_time);
   }
   db_restart.unmount();

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
   printF("\n********************************************************\n");
   printF("[ Restart Manager status ]:\n");

   if (m_is_from_restart) {
      printF("\tRunning from restart\n");
      printF("\tRestart read directory path: %s\n",
             m_restart_read_path.c_str());
      printF("\tRestart file index: %i\n", m_restart_index);
   }
   else {
      printF("\tRunning from initial conditions\n");
   }
   printF("\tRestart write directory path: %s\n",
          m_restart_write_path.c_str());

   if (m_write_on_time) {
      printF("\t  Restart every %e time units\n", m_time_interval);
      printF("\t  Next write time: %e\n", m_next_write_time);
   }
   else {
      printF("\t  Restart every %i time steps\n", m_step_interval);
   }
}


void
RestartManager::parseParameters(
   ParmParse& a_pp)
{
   // Logic to take either time or step interval.  You can specify a restart
   // write frequency based on time steps or simulation time but not both.  One
   // or the other must be specified.
   bool has_step_interval = a_pp.contains("step_interval");
   bool has_time_interval = a_pp.contains("time_interval");
   if (has_step_interval && has_time_interval) {
      OV_ABORT("Must set restart frequency for only one of steps or time, not both.");
   }
   else if (!has_step_interval && !has_time_interval) {
      OV_ABORT("Must set restart frequency for one of steps or time.");
   }
   else {
      a_pp.query("step_interval", m_step_interval);
      a_pp.query("time_interval", m_time_interval);
      // Get the actual intervals and check for validity.  There is an
      // additional flag to indicate which scheme is being used.
      if (has_step_interval) {
         if (m_step_interval <= 0) {
            OV_ABORT("Input step_interval must be > 0");
         }
      }
      else {
         if (m_time_interval <= 0.0) {
            OV_ABORT("Input time_interval must be > 0.0");
         }
         else {
            m_write_on_time = true;
         }
      }
   }

   a_pp.query("write_directory", m_restart_write_path);
   a_pp.query("read_directory", m_restart_read_path);
}


void
RestartManager::restore()
{
   if (m_restart_index == 0) {
      OV_ABORT("No distributions found from which to restart ... quitting");
   }

   // Glue together the info we already have to form the name of the file we'll
   // read and open it.
   char restart_filename[100];
   createFileName(restart_filename,
      m_restart_read_path.c_str(),
      m_restart_index );
   printF("Restarting from %s\n", restart_filename);

   HDF_DataBase db_restart;
   db_restart.mount(restart_filename, "R"); // mount, R=read-only

   // Loop over vector of Serializable items and have each restore their state
   // from the restart file.
   for (SerializableList::iterator iter(m_items.begin());
        iter != m_items.end();
        ++iter) {
      (*iter)->getFromRestart(db_restart);
   }

   db_restart.unmount();

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
   char restart_filename[100];
   int i_max(500); // largest possible number of distribution functions
   for (int i(1); i <= i_max; ++i) {
      createFileName(restart_filename, m_restart_read_path.c_str(), i);
      HDF_DataBase db_restart;
      if (db_restart.mount(restart_filename, "R") == 0) { // mount, R=read-only
         m_restart_index = i;
         db_restart.unmount();
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
