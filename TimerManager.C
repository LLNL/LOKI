/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "TimerManager.H"

#include "Loki_Utilities.H"

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
   const string& a_name) const
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
   const string& a_name)
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
   const string& a_name)
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
            cout << "Warning: startTimer called for " << a_name
                 << ", which is already running!" << endl;
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
   const string& a_name)
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
            cout << "Warning: stopTimer called for " << a_name
                 << ", which is not running!" << endl;
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
   vector<size_t> restart_list(0);
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
   vector<double> net_times(m_timers.size());
   int i(0);
   for (TimerList::iterator iter(m_timers.begin());
        iter != m_timers.end();
        ++iter) {
      net_times[i] = Loki_Utilities::getSum((*iter)->elapsedTime(), 0);
      ++i;
   }
   double total_time(accumulate(net_times.begin(), net_times.end(), 0.0));

   // Write each timer's non-exclusive stats and the stats on the total time.
   double per_cell(1.0 / static_cast<double>(a_number_of_cells));
   double per_cell_per_step(per_cell / a_step_number);
   double to_percent_time(100.0 / total_time);
   Loki_Utilities::printF("   Time for                  total     per cell    per(cell*timestep)    %Total\n");
   int j(0);
   for (TimerList::iterator iter(m_timers.begin());
        iter != m_timers.end();
        ++iter) {
      string name((*iter)->name());
      Loki_Utilities::printF("      %-19s: %3.3e   %3.3e      %3.3e          %3.2f\n",
         name.c_str(),
         net_times[j],
         net_times[j] * per_cell,
         net_times[j] * per_cell_per_step,
         net_times[j] * to_percent_time);
      ++j;
   }
   Loki_Utilities::printF("      ---------------------------------------------------------\n");
   Loki_Utilities::printF("      total:               %3.3e   %3.3e      %3.3e\n",
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
   vector<size_t> restart_list(0);
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
   vector<double> net_times(m_timers.size());
   int i(0);
   for (TimerList::iterator iter(m_timers.begin());
        iter != m_timers.end();
        ++iter) {
      net_times[i] = Loki_Utilities::getSum((*iter)->exclusiveTime(), 0);
      ++i;
   }
   double total_time(accumulate(net_times.begin(), net_times.end(), 0.0));

   // Write each timer's exclusive stats and the stats on the total time.
   double per_cell(1.0 / static_cast<double>(a_number_of_cells));
   double per_cell_per_step(per_cell / a_step_number);
   double to_percent_time(100.0 / total_time);
   Loki_Utilities::printF("   Time for                  total     per cell    per(cell*timestep)    %Total\n");
   int j(0);
   for (TimerList::iterator iter(m_timers.begin());
        iter != m_timers.end();
        ++iter) {
      string name((*iter)->name());
      Loki_Utilities::printF("      %-19s: %3.3e   %3.3e      %3.3e          %3.2f\n",
         name.c_str(),
         net_times[j],
         net_times[j] * per_cell,
         net_times[j] * per_cell_per_step,
         net_times[j] * to_percent_time);
      ++j;
   }
   Loki_Utilities::printF("      ---------------------------------------------------------\n");
   Loki_Utilities::printF("      total:               %3.3e   %3.3e      %3.3e\n",
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
   const string& a_name)
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
