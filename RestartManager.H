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
#ifndef _RESTART_MANAGER_H_
#define _RESTART_MANAGER_H_

#include "Overture.h"
#include "ParmParse.H"

#include "tbox/Pointer.H"
#include "Serializable.H"

#include <string>
#include <vector>

namespace Loki {

/**
 * Singleton class allowing classes derived from Serializable to register
 * themselves to write to and be restored from restart files.
 */
class RestartManager
{
public:
   /**
    * @brief Singleton accessor method.
    *
    * Static method that returns the pointer to the singleton object.  If the
    * object has not been created, this method creates the singleton object.
    *
    * @returns pointer to the singleton object
    */
   static
   RestartManager*
   getManager();

   /**
    * @brief Registers an object adhering to the Serializable interface with
    * the manager.
    *
    * @param[in] a_item Pointer to an object of Serializable type
    */
   void
   registerRestart(
      Serializable* a_item)
   {
      m_items.push_back(a_item);
   }

   /**
    * @brief Returns true if there is an action to be performed at this time or
    * time step, false otherwise.
    *
    * @param[in] a_time Current simulation time.
    * @param[in] a_step Current simulation step.
    *
    * @returns boolean result
    */
   bool
   requiresAction(
      real a_time,
      int a_step)
   {
      if (m_write_on_time) { // time
         if (a_time >= m_next_write_time) {
            return true;
         }
      }
      else { // step
         if (a_step % m_step_interval == 0) {
            return true;
         }
      }
      return false;
   }

   /**
    * @brief Restore from last restart database.
    */
   void
   restore();

   /**
    * @brief Reset next write time based on restarted simulation time
    *
    * @param[in] a_time Next time to write.
    */
   void
   resetNextWriteTime(
      real a_time);

   /**
    * @brief Dumps a restart database.
    *
    * @param[in] a_time Current simulation time.
    */
   void
   write(
      real a_time);

   /**
    * @brief Returns true if current run was started from a restart.
    */
   bool
   isFromRestart()
   {
      return m_is_from_restart;
   }

   /**
    * @brief Prints the current internal state of the RestartManager.
    */
   void
   print();

   /**
    * @brief Delete singleton instance of this class when system is shut down.
    */
   static
   void
   shutdown()
   {
      if (s_restart_manager_instance) {
         delete s_restart_manager_instance;
         s_restart_manager_instance = 0;
      }
   }

private:
   // Unimplemented copy constructor.
   RestartManager(
      const RestartManager& other);

   // Unimplemented assignment operator.
   RestartManager&
   operator = (
      const RestartManager& rhs);

   // Obtain control parameters from "maxwell" section of database.
   void
   parseParameters(
      ParmParse& a_pp);

   // Find last restart file for last saved time.
   void
   findRestartIndex();

   // Glue together directory and file index to form file name.
   void
   createFileName(
      char* a_name,
      const char* a_dir,
      int a_index)
   {
      sprintf(a_name, "%s/dist_%i.hdf", a_dir, a_index);
   }

   // The singleton instance.
   static RestartManager* s_restart_manager_instance;

   // Default constructor.
   RestartManager();

   // Destructor.
   ~RestartManager();

   typedef std::vector<Serializable*> SerializableList;

   // The items to write to restore from restart.
   SerializableList m_items;

   // Number of time steps between each write.
   int m_step_interval;

   // Time between each write.
   real m_time_interval;

   // Next time to write.
   real m_next_write_time;

   // If true write base on time.  Otherwise write on time step.
   bool m_write_on_time;

   // Path to where restart files will be written.
   aString m_restart_write_path;

   // Path to where restart files will be read.
   aString m_restart_read_path;

   // Restart file prefix.
   aString m_file_name_prefix;

   // If true problem is being run from restart.
   bool m_is_from_restart;

   // This does double duty.  If we're starting a simulation from a restart file
   // it is initially the index of the restart file to restore from.  Once the
   // restore process (if any) has completed it is the index of the next new
   // restart file to write.
   int m_restart_index;
};

} // end namespace Loki

#endif
