/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _RESTART_MANAGER_H_
#define _RESTART_MANAGER_H_

#include "LokiInputParser.H"
#include "tbox/Pointer.H"
#include "Serializable.H"

#include <sstream>
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
    *        the manager.
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
    *        time step, false otherwise.
    *
    * @param[in] a_time Current simulation time.
    * @param[in] a_step Current simulation step.
    *
    * @returns boolean result
    */
   bool
   requiresAction(
      double a_time,
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
      double a_time);

   /**
    * @brief Dumps a restart database.
    *
    * @param[in] a_time Current simulation time.
    */
   void
   write(
      double a_time);

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

   /**
    * @brief Returns the current restart file index.
    */
   int
   restartIndex() const
   {
      return m_restart_index;
   }

   /**
    * @brief Returns the path where restart files will be written.
    */
   string
   restartWritePath()
   {
      return m_restart_write_path;
   }

   /**
    * @brief Returns the path where restart files will be read.
    */
   string
   restartReadPath()
   {
      return m_restart_read_path;
   }

private:
   // Unimplemented copy constructor.
   RestartManager(
      const RestartManager& other);

   // Unimplemented assignment operator.
   RestartManager&
   operator = (
      const RestartManager& rhs);

   // Obtain control parameters from the "restart" section of database.
   void
   parseParameters(
      LokiInputParser& a_pp);

   // Find last restart file for last saved time.
   void
   findRestartIndex();

   // Glue together directory and file index to form file name.
   void
   createFileName(
      ostringstream& a_name,
      const string& a_dir,
      int a_index)
   {
      a_name << a_dir << "/dist_" << a_index << ".hdf";
   }

   // The singleton instance.
   static RestartManager* s_restart_manager_instance;

   // Default constructor.
   RestartManager();

   // Destructor.
   ~RestartManager();

   typedef vector<Serializable*> SerializableList;

   // The items to write to restore from restart.
   SerializableList m_items;

   // Number of time steps between each write.
   int m_step_interval;

   // Time between each write.
   double m_time_interval;

   // Next time to write.
   double m_next_write_time;

   // If true write base on time.  Otherwise write on time step.
   bool m_write_on_time;

   // Path to where restart files will be written.
   string m_restart_write_path;

   // Path to where restart files will be read.
   string m_restart_read_path;

   // If true problem is being run from restart.
   bool m_is_from_restart;

   // This does double duty.  If we're starting a simulation from a restart file
   // it is initially the index of the restart file to restore from.  Once the
   // restore process (if any) has completed it is the index of the next new
   // restart file to write.
   int m_restart_index;

   // Maximum number of bulkdata files for each restart dump.
   int m_max_restart_files;
};

} // end namespace Loki

#endif
