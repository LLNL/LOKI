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
#ifndef _TIMER_MANAGER_H_
#define _TIMER_MANAGER_H_

#include "Overture.h"
#include "tbox/Pointer.H"

#include <string>
#include <vector>

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
      const std::string& a_name) const;

   /**
    * @brief Create the named timer.
    *
    * Attempting to create an already existing timer results in a warning.
    *
    * @param[in] a_name Name of time to create.
    */
   void
   createTimer(
      const std::string& a_name);

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
      const std::string& a_name);

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
      const std::string& a_name);

   /**
    * @brief Stop all running timers, print elapsed time info for each existing
    * timer, then restart all timers that were running.
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
    * existing timer, then restart all timers that were running.
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
      const std::string&
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
            OV_ABORT("Timer " + m_name + " already running!");
         }
         m_is_started = true;
         m_is_running = true;
         m_start = getCPU();
         m_start_exclusive = m_start;
      }

      /**
       * @brief Starts the timer.
       */
      void
      stop()
      {
         if (!m_is_started) {
            OV_ABORT("Timer " + m_name + " not running!");
         }
         m_is_started = false;
         m_is_running = false;
         real stop = getCPU();
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
            OV_ABORT("Timer " + m_name + " already running!");
         }
         m_is_running = true;
         m_start_exclusive = getCPU();
      }

      // Stop exclusive timing of this timer.
      void
      stop_exclusive()
      {
         if (!m_is_running) {
            OV_ABORT("Timer " + m_name + " not running!");
         }
         m_is_running = false;
         real stop = getCPU();
         m_time_exclusive += double(stop - m_start_exclusive);
      }

      // Constructor.  May only be constructed by TimerManager.
      explicit Timer(
         const std::string& a_name);

      // Returns the elapsed time accumulated by this timer--the difference
      // between when the timer stopped and when it started.
      real
      elapsedTime() const
      {
         return m_time;
      }

      // Returns the exclusive time accumulated by this timer.
      real
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
      std::string m_name;

      // True if timer is running and therefore accumulating exclusive time.
      bool m_is_running;

      // True if timer has been started.  The timer may be started but not
      // running (i.e. not accumulating exclusive time).
      bool m_is_started;

      // Accumulated elapsed time.
      real m_time;

      // Accumulated exclusive time.
      real m_time_exclusive;

      // Time when the time was started.
      real m_start;

      // Time when exclusive timing started.
      real m_start_exclusive;
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

   typedef std::vector<Timer*> TimerList;

   // All timers that have been created.
   TimerList m_timers;

   // All running timers.  The first member is the first timer started and the
   // last member is the last one started.
   std::vector<Timer*> m_running_timers;
};

} // end namespace Loki

#endif
