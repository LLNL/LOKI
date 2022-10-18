/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _TIMER_MANAGER_H_
#define _TIMER_MANAGER_H_

#include "Loki_Defines.H"
#include "tbox/Pointer.H"

#include <mpi.h>
#include <string>
#include <vector>

using namespace std;

namespace Loki {

/**
 * TimerManager is a singleton class for the management of timers.
 */
class TimerManager
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
   TimerManager*
   getManager();

   /**
    * @brief Returns true if the named timer exists.
    *
    * @param[in] a_name Name of time to check for.
    */
   bool
   timerExists(
      const string& a_name) const;

   /**
    * @brief Create the named timer.
    *
    * Attempting to create an already existing timer results in a warning.
    *
    * @param[in] a_name Name of time to create.
    */
   void
   createTimer(
      const string& a_name);

   /**
    * @brief Start the named timer.
    *
    * Attempting to start an already running or non-existant timer results in
    * a warning.
    *
    * @param[in] a_name Name of time to start.
    */
   void
   startTimer(
      const string& a_name);

   /**
    * @brief Stop the named timer.
    *
    * Attempting to stop an already stopped or non-existant timer results in
    * a warning.
    *
    * @param[in] a_name Name of time to stop.
    */
   void
   stopTimer(
      const string& a_name);

   /**
    * @brief Stop all running timers, print elapsed time info for each existing
    *        timer, then restart all timers that were running.
    *
    * @param[in] a_step_number Simulation step.
    * @param[in] a_number_of_cells Number of cells in problem.
    */
   void
   print(
      int a_step_number,
      long int a_number_of_cells);

   /**
    * @brief Stop all running timers, print exclusive time info for each
    *        existing timer, then restart all timers that were running.
    *
    * @param[in] a_step_number Simulation step.
    * @param[in] a_number_of_cells Number of cells in problem.
    */
   void
   print_exclusive(
      int a_step_number,
      long int a_number_of_cells);

   /**
    * @brief Delete singleton instance of this class when system is shut down.
    */
   static
   void
   shutdown()
   {
      if (s_timer_manager_instance) {
         delete s_timer_manager_instance;
         s_timer_manager_instance = 0;
      }
   }

private:
   /**
    * The actual timers.
    */
   class Timer
   {
   friend class TimerManager;

   public:
      /**
       * @brief Destructor.
       */
      ~Timer();

      /**
       * @brief Returns the name of the Timer.
       */
      const string&
      name() const
      {
         return m_name;
      }

      /**
       * @brief Starts the timer.
       */
      void
      start()
      {
         if (m_is_started) {
            LOKI_ABORT("Timer " + m_name + " already running!");
         }
         m_is_started = true;
         m_is_running = true;
         m_start = MPI_Wtime();
         m_start_exclusive = m_start;
      }

      /**
       * @brief Starts the timer.
       */
      void
      stop()
      {
         if (!m_is_started) {
            LOKI_ABORT("Timer " + m_name + " not running!");
         }
         m_is_started = false;
         m_is_running = false;
         double stop = MPI_Wtime();
         m_time += double(stop - m_start);
         m_time_exclusive += double(stop - m_start_exclusive);
      }

   protected:
      // Restart the timer to time exclusively--it it now no longer embedded in
      // another timer.
      void
      restart_exclusive()
      {
         if (m_is_running) {
            LOKI_ABORT("Timer " + m_name + " already running!");
         }
         m_is_running = true;
         m_start_exclusive = MPI_Wtime();
      }

      // Stop exclusive timing of this timer.
      void
      stop_exclusive()
      {
         if (!m_is_running) {
            LOKI_ABORT("Timer " + m_name + " not running!");
         }
         m_is_running = false;
         double stop = MPI_Wtime();
         m_time_exclusive += double(stop - m_start_exclusive);
      }

      // Constructor.  May only be constructed by TimerManager.
      explicit Timer(
         const string& a_name);

      // Returns the elapsed time accumulated by this timer--the difference
      // between when the timer stopped and when it started.
      double
      elapsedTime() const
      {
         return m_time;
      }

      // Returns the exclusive time accumulated by this timer.
      double
      exclusiveTime() const
      {
         return m_time_exclusive;
      }

      // Returns true if the timer is running (i.e.is accumulating exclusive
      // time).
      bool
      isRunning() const
      {
         return m_is_running;
      }

   private:
      // Unimplemented default constructor.
      Timer();

      // Unimplemented copy constructor.
      Timer(
         const Timer& other);

      // Unimplemented assignment operator.
      Timer&
      operator = (
         const Timer& rhs);

      // Name of timer.
      string m_name;

      // True if timer is running and therefore accumulating exclusive time.
      bool m_is_running;

      // True if timer has been started.  The timer may be started but not
      // running (i.e. not accumulating exclusive time).
      bool m_is_started;

      // Accumulated elapsed time.
      double m_time;

      // Accumulated exclusive time.
      double m_time_exclusive;

      // Time when the time was started.
      double m_start;

      // Time when exclusive timing started.
      double m_start_exclusive;
   };

   // The singleton instance.
   static TimerManager* s_timer_manager_instance;

   // Constructor.
   TimerManager();

   // Destructor.
   ~TimerManager();

   // Unimplemented copy constructor.
   TimerManager(
      const TimerManager& other);

   // Unimplemented assignment operator.
   TimerManager&
   operator = (
      const TimerManager& rhs);

   // The last timer started.
   Timer* m_current_timer;

   typedef vector<Timer*> TimerList;

   // All timers that have been created.
   TimerList m_timers;

   // All running timers.  The first member is the first timer started and the
   // last member is the last one started.
   vector<Timer*> m_running_timers;
};

} // end namespace Loki

#endif
