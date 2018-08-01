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
#include "TimerManager.H"

#include "ParallelUtility.h"

#include <numeric>

using namespace std;

namespace Loki {

TimerManager* TimerManager::s_timer_manager_instance = 0;

TimerManager*
TimerManager::getManager()
{
   // If the singleton instance does not exist then build it.
   if (!s_timer_manager_instance) {
      s_timer_manager_instance = new TimerManager();
   }
   return s_timer_manager_instance;
}


bool
TimerManager::timerExists(
   const std::string& a_name) const
{
   // Not especially efficient but we only do this once for each timer so
   // it's not a big deal.
   bool timer_found(false);
   if (!a_name.empty()) {
      for (TimerList::const_iterator iter(m_timers.begin());
           iter != m_timers.end();
           ++iter) {
         if ((*iter)->name() == a_name) {
            timer_found = true;
            break;
         }
      }
   }
   return timer_found;
}


void
TimerManager::createTimer(
   const std::string& a_name)
{
   if (!timerExists(a_name)) {
      Timer* T_tmp = new Timer(a_name);
      m_timers.push_back(T_tmp);
   }
   else {
      cout << "Warning: Cannot create timer " << a_name << ".";
      cout << "Timer by this name alreadt exists!" << endl;
   }
}


void
TimerManager::startTimer(
   const std::string& a_name)
{
   // Every time we start the timer we must do a linear search for it which
   // isn't very efficient.  We do this alot but the list is hopefully small.
   bool timer_found(false);
   for (TimerList::iterator iter(m_timers.begin());
        iter != m_timers.end();
        ++iter) {
      if ((*iter)->name() == a_name) {
         if (!(*iter)->isRunning()) {
            if (m_current_timer != 0) {
               m_current_timer->stop_exclusive();
               m_running_timers.push_back(m_current_timer);
            }
            m_current_timer = ((*iter));
            m_current_timer->start();
         }
         else {
            cout << "Warning: startTimer called for " << a_name << ", which is already running!" << endl;
         }
         timer_found = true;
         break;
      }
   }
   if (!timer_found) {
      cout << "Warning: Requested timer " << a_name << " not found!" << endl;
   }
}


void
TimerManager::stopTimer(
   const std::string& a_name)
{
   // Every time we stop the timer we must do a linear search for it which
   // isn't very efficient.  We do this alot but the list is hopefully small.
   bool timer_found(false);
   for (TimerList::iterator iter(m_timers.begin());
        iter != m_timers.end();
        ++iter) {
      if ((*iter)->name() == a_name) {
         if ((*iter)->isRunning()) {
            (*iter)->stop();
            if (m_current_timer == (*iter)) {
               m_current_timer = 0;
               if (m_running_timers.size() != 0) {
                  m_current_timer = m_running_timers.back();
                  m_running_timers.pop_back();
                  m_current_timer->restart_exclusive();
               }
            }
            else {
               for (size_t j(0); j < m_running_timers.size(); ++j) {
                  if (m_running_timers[j]->name() == a_name) {
                     for (size_t k(j+1); k < m_running_timers.size(); ++k) {
                        m_running_timers[k-1] =  m_running_timers[k];
                     }
                     m_running_timers.pop_back();
                  }
               }
            }
         }
         else {
            cout << "Warning: stopTimer called for " << a_name << ", which is not running!" << endl;
         }
         timer_found = true;
         break;
      }
   }
   if (!timer_found) {
      cout << "Warning: Requested timer " << a_name << " not found!" << endl;
   }
}


void
TimerManager::print(
   int a_step_number,
   long int a_number_of_cells)
{
   // stop all timers and record those running
   std::vector<size_t> restart_list(0);
   int k(0);
   for (TimerList::iterator iter(m_timers.begin());
        iter != m_timers.end();
        ++iter) {
      if ((*iter)->isRunning()) {
         (*iter)->stop();
         restart_list.push_back(k);
      }
      ++k;
   }

   // Get the total time.
   vector<real> net_times(m_timers.size());
   int i(0);
   for (TimerList::iterator iter(m_timers.begin());
        iter != m_timers.end();
        ++iter) {
      net_times[i] = ParallelUtility::getSum((*iter)->elapsedTime(), 0);
      ++i;
   }
   real total_time(std::accumulate(net_times.begin(), net_times.end(), 0.0));

   // Write each timer's non-exclusive stats and the stats on the total time.
   real per_cell(1.0 / static_cast<double>(a_number_of_cells));
   real per_cell_per_step(per_cell / a_step_number);
   real to_percent_time(100.0 / total_time);
   printF("   Time for                  total     per cell    per(cell*timestep)    %Total\n");
   int j(0);
   for (TimerList::iterator iter(m_timers.begin());
        iter != m_timers.end();
        ++iter) {
      std::string name((*iter)->name());
      printF("      %-19s: %3.3e   %3.3e      %3.3e          %3.2f\n",
         name.c_str(),
         net_times[j],
         net_times[j] * per_cell,
         net_times[j] * per_cell_per_step,
         net_times[j] * to_percent_time);
      ++j;
   }
   printF("      ---------------------------------------------------------\n");
   printF("      total:               %3.3e   %3.3e      %3.3e\n",
      total_time,
      total_time * per_cell,
      total_time * per_cell_per_step);

   // restart previously running timers
   for (size_t i(0); i < restart_list.size(); ++i) {
      m_timers[restart_list[i]]->start();
   }
}


void
TimerManager::print_exclusive(
   int a_step_number,
   long int a_number_of_cells)
{
   // stop all timers and record those running
   std::vector<size_t> restart_list(0);
   int k(0);
   for (TimerList::iterator iter(m_timers.begin());
        iter != m_timers.end();
        ++iter) {
      if ((*iter)->isRunning()) {
         (*iter)->stop();
         restart_list.push_back(k);
      }
      ++k;
   }

   // Get the total time.
   vector<real> net_times(m_timers.size());
   int i(0);
   for (TimerList::iterator iter(m_timers.begin());
        iter != m_timers.end();
        ++iter) {
      net_times[i] = ParallelUtility::getSum((*iter)->exclusiveTime(), 0);
      ++i;
   }
   real total_time(std::accumulate(net_times.begin(), net_times.end(), 0.0));

   // Write each timer's exclusive stats and the stats on the total time.
   real per_cell(1.0 / static_cast<double>(a_number_of_cells));
   real per_cell_per_step(per_cell / a_step_number);
   real to_percent_time(100.0 / total_time);
   printF("   Time for                  total     per cell    per(cell*timestep)    %Total\n");
   int j(0);
   for (TimerList::iterator iter(m_timers.begin());
        iter != m_timers.end();
        ++iter) {
      std::string name((*iter)->name());
      printF("      %-19s: %3.3e   %3.3e      %3.3e          %3.2f\n",
         name.c_str(),
         net_times[j],
         net_times[j] * per_cell,
         net_times[j] * per_cell_per_step,
         net_times[j] * to_percent_time);
      ++j;
   }
   printF("      ---------------------------------------------------------\n");
   printF("      total:               %3.3e   %3.3e      %3.3e\n",
      total_time,
      total_time * per_cell,
      total_time * per_cell_per_step);

   // restart previously running timers
   for (size_t i(0); i < restart_list.size(); ++i) {
      m_timers[restart_list[i]]->start();
   }
}

TimerManager::TimerManager()
   : m_current_timer(0)
{
}

TimerManager::~TimerManager()
{
   while (m_running_timers.size() > 0) {
      m_running_timers.back()->stop();
      m_running_timers.pop_back();
   }

   for (TimerList::iterator iter(m_timers.begin());
        iter != m_timers.end();
        ) {
      delete *iter;
      iter = m_timers.erase(iter);
   }
}

TimerManager::Timer::Timer(
   const std::string& a_name)
   : m_name(a_name),
     m_is_running(false),
     m_is_started(false),
     m_time(0.0),
     m_time_exclusive(0.0),
     m_start(0.0),
     m_start_exclusive(0.0)
{
}

TimerManager::Timer::~Timer()
{
}

} // end namespace Loki
