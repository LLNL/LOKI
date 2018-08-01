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
#include "Simulation.H"

#include <limits>
#include <sstream>

#include "TimerManager.H"
#include "RestartManager.H"
#include "Loki_Utilities.H"
#include "VMSystem.H"
#include "VPSystem.H"

namespace Loki {

const int Simulation::LOKI_VERSION_MAJOR = 2;
const int Simulation::LOKI_VERSION_MINOR = 3;
const int Simulation::LOKI_VERSION_PATCHLEVEL = 0;

void
Simulation::printParameters() const
{
   // Print the general simulation parameters.
   printF("\n#*#*# General Simulation Parameters #*#*#\n");
   printF("Loki version %d.%d.%d\n",
      LOKI_VERSION_MAJOR,
      LOKI_VERSION_MINOR,
      LOKI_VERSION_PATCHLEVEL);
   printF("spatial solution order  = %d\n", m_spatial_solution_order);
   printF("temporal solution order = %d\n", m_temporal_solution_order);
   ////////////////////////////////////
   // A time step manager will eventually write these out
   printF("current time            = %e\n", m_time);
   printF("time from start of run  = %e\n", m_time_from_start_of_run);
   printF("final time              = %e\n", m_max_time);
   printF("cfl                     = %e\n", m_cfl);
   ////////////////////////////////////

   ////////////////////////////////////
   // A diagnostic manager will eventually write these out
   if (m_save_data) {
      printF("saving data to %s\n", (const char*) m_show_file_name);
      printF("  save times = %e\n", m_save_times);
   }
   else {
      printF("Not saving data\n");
   }
   ///////////////////////////////////
}


void
Simulation::writePlotFile()
{
   // If the user didn't turn off writing of plot output then write plot files.
   if (m_save_data) {
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("output");

      // Delegate actual writing to the System.
      m_system->plot(m_time,
         m_dt,
         m_probes, 
         m_num_probes,
         m_saved_seq,
         m_saved_save,
         m_show); 

      // Except for the write at the first time step, write out all timer info.
      if (m_step != 0) {
         printF("\n\n************** PRINT TIMING INFO **************\n");
         int numProcs(Communication_Manager::numberOfProcessors());
         printF("t=%3.3e, dt=%3.3e: %i time steps using %i processors\n",
            m_time,
            m_dt,
            m_step,
            numProcs);
         timers->print(m_step, m_system->problemSize());

         printF("****************** Exclusive Timers *******************\n");
         printF("t=%3.3e, dt=%3.3e: %i time steps using %i processors\n",
            m_time,
            m_dt,
            m_step,
            numProcs);
         timers->print_exclusive(m_step, m_system->problemSize());
         printF("*************** END TIMING INFO ***************\n\n\n");

         fflush(0);
         Communication_Manager::Sync();
      }

      timers->stopTimer("output");
   }
}


void
Simulation::accumulateSequences()
{
   // If the user didn't turn off writing of plot output then collect the next
   // time history.
   if (m_save_data) {
      TimerManager* timers(TimerManager::getManager());
      timers->startTimer("output");

      // Delegate actual collection to the System.
      m_system->accumulateSequences(m_time,
         m_probes,
         m_num_probes,
         m_saved_seq);

      timers->stopTimer("output");
   }
}


void
Simulation::writeCheckpointFile()
{
   TimerManager* timers(TimerManager::getManager());
   timers->startTimer("output");

   // The restart manager knows when the next checkpoint (restart) file(s) need
   // to be written.
   RestartManager* restart_manager(RestartManager::getManager());
   if (restart_manager->requiresAction(m_time, m_step)) {
     time_t wall_time = time(NULL);
     printF("\n****************************************************\n");
     printF("Wall clock at restart dump: %s", asctime(localtime(&wall_time)));
     printF("****************************************************\n\n");

     // The ghosts are updated at the top of the time step loop so at this
     // point they are out of sync and need to be updated so the restart file
     // is correct.
     m_system->updateGhosts();

     // All entities that participate in writing restart data have been
     // registered with the restart manager so we delegate this task it it.
     restart_manager->write(m_time);

     // Get info on next restart time/time step.
     restart_manager->print();
   }
   timers->stopTimer("output");
}


Simulation::Simulation(
   ParmParse& a_pp)
  :   m_spatial_solution_order(4),
      m_temporal_solution_order(4),
      m_verbosity(1),
      m_step(0),
      m_max_step(0),
      m_time(0.0),
      m_max_time(1.0),
      m_dt(-1.0),
      m_cfl(0.9),
   ////////////////////////////////////
   // Belong in a diagnostic manager
      m_save_times(1.0),
      m_save_data(true),
      m_show_file_name("ex.show"),
      m_seq_times(1.0),
      m_num_probes(-1),
      m_saved_seq(0),
      m_saved_save(0),
      m_last_seq(0),
      m_last_save(0),
   ////////////////////////////////////
      m_time_from_start_of_run(0.0),
      m_system(0)
{
   TimerManager* timers(TimerManager::getManager());
   timers->createTimer("setup");
   timers->createTimer("output");

   timers->startTimer("setup");

   // Parse all general simulation related input.
   parseParameters(a_pp);

   // This could probably have been folded into parseParameters.
   // See which system we're running and construct it.
   aString test_str("poisson");
   a_pp.query("sys_type", test_str);
   if (test_str.matches("poisson")) {
      m_system = new VPSystem(a_pp,
         m_num_probes,
         m_spatial_solution_order,
         m_temporal_solution_order);
   }
   else if (test_str.matches("maxwell")) {
      m_system = new VMSystem(a_pp,
         m_num_probes,
         m_spatial_solution_order,
         m_temporal_solution_order);
   }
   else {
      OV_ABORT("Unknown sys_type input.");
   }

   // The Simulation writes restart data so register this object with the
   // RestartManager which will use the putToRestart/getFromRestart callbacks to
   // get the restart data written/read.
   RestartManager* restart_manager(RestartManager::getManager());
   restart_manager->registerRestart(this);

   // Virtaully all the objects in play have now been built.  If we're running
   // from a restart file, ask the RestartManager to restore the state data of
   // all objects that have registered themselves with the RestartManager to
   // what is in the restart file.  As m_time is set from the restart file we
   // need to tell the RestartManager that this time is the effective simulation
   // start time and hence the time of the next (first) restart write.
   if (restart_manager->isFromRestart()) {
     restart_manager->restore();
     restart_manager->resetNextWriteTime(m_time);
   }
   restart_manager->print();

   // Initialize the System.  The name is kind of odd given that we might have
   // just restored everything from a restart file which initializes a buch of
   // stuff.
   m_system->initialize(restart_manager->isFromRestart(), m_time, m_dt);


   ////////////////////////////////////
   // Belong in a diagnostic manager
   // This sets up some general parameters of the Overture show file, what the
   // plot data is written to.
   if (m_save_data) {
      bool use_stream_mode(true);
      m_show.open(m_show_file_name, ".", use_stream_mode);
      m_show.saveGeneralComment("cgpp: Vlassov Poisson 4D");
      // flush frequency must be at least the number of frame series.
      // to write sequences it must be even bigger.
      m_show.setFlushFrequency(m_system->getNumFrameSeries());
   }
   ////////////////////////////////////


   // Everything is now built and tidied up so collecte a time history at the
   // initial time, write the plot data, and write the restart data.
   accumulateSequences();
   writePlotFile();

   writeCheckpointFile();

   Communication_Manager::Sync();
   timers->stopTimer("setup");
}

Simulation::~Simulation()
{
}


void
Simulation::advance()
{
   // Figure out a hopefully safe dt.
   selectTimeStep();

   ++m_step;

   // Advance the System to the next time.
   m_time = m_system->advance(m_time, m_dt);
   m_time_from_start_of_run += m_dt;

   if (m_verbosity >= 1) {
      printF("Step %d completed, simulation time is %f, dt is %f\n",
         m_step, m_time, m_dt);
   }

   // If it is the next time to save the time histories then do it.
   if (m_time >= (m_last_seq+1)*m_seq_times) {
      accumulateSequences();
      m_last_seq++;
   }
   // If it is the next time to write a plot then do it.
   if (m_time >= (m_last_save+1)*m_save_times) {
      writePlotFile();
      m_last_save++;
   }

   // The check for whether it time to write the next checkpoint file is
   // in this function as opposed to the other 2 immediately above.
   writeCheckpointFile();
}


void
Simulation::finalize()
{
   // We always want a final time step checkpoint file.
   writeCheckpointFile();
   Communication_Manager::Sync();
}


void
Simulation::parseParameters(
   ParmParse& a_pp)
{
   // Get the spatial and temporal solution orders.
   a_pp.query("spatial_solution_order", m_spatial_solution_order);
   if (m_spatial_solution_order != 4 && m_spatial_solution_order != 6) {
      OV_ABORT("spatial_solution_order must be 4 or 6");
   }
   a_pp.query("temporal_solution_order", m_temporal_solution_order);
   if (m_temporal_solution_order != 4 && m_temporal_solution_order != 6) {
      OV_ABORT("temporal_solution_order must be 4 or 6");
   }

   // This determines the amount of diagnositic output generated
   a_pp.query("verbosity", m_verbosity);

   // set cfl number for the case of dynamic timestep selection
   a_pp.query("cfl", m_cfl);

   // Stop when the simulation time gets here
   a_pp.query("final_time", m_max_time);

   // Stop after this number of steps
   a_pp.query("max_step", m_max_step);

   ////////////////////////////////////
   // A diagnostic manager will eventually read these in
   // This default, single probe is pretty arbitrary.  Read in if the user has
   // a better idea.
   m_num_probes = 1;
   m_probes.resize(2, m_num_probes);
   m_probes(0, 0) = 0.5;
   m_probes(1, 0) = 0.5;
   if (a_pp.query("number_of_probes", m_num_probes)) {
      m_probes.resize(2, m_num_probes);
      for (int i(0); i < m_num_probes; ++i) {
         std::ostringstream oss;
         oss << "probe." << i+1;
         ParmParse prpp(oss.str().c_str());

         Array<double> loc(2);
         prpp.getarr("location", loc, 0, 2);
         m_probes(X1, i) = loc[X1];
         m_probes(X2, i) = loc[X2];
      }
   }

   // Get the frequency in simulation time for the collection of time history
   // points.
   a_pp.query("sequence_write_times", m_seq_times);

   // Get the frequency in simulation time for writting plot data.
   a_pp.query("save_times", m_save_times);

   // User can turn off collection/writing of plot data.
   aString sdata("true");
   a_pp.query("save_data", sdata);
   m_save_data = sdata.matches("true") ? true : false;

   // Get the base part of the names of all plot files.  Overture will append
   // to this according to its own conventions.
   a_pp.query("show_file_name", m_show_file_name);
   ////////////////////////////////////

   printParameters();
}


void
Simulation::selectTimeStep()
{
   // We fiddle with dt so that it honors the requested plot frequency.  This is
   // essential to allow Fourier transforms of the plot data.
  real next_time = (m_last_save + 1) * m_save_times;
  //(floor(m_time/m_save_times+100*std::numeric_limits<real>::epsilon()) + 1) *m_save_times;
  
   const real final_time(std::min(next_time, m_max_time));

   const real dt_stable(m_system->stableDt() * m_cfl);

   const real time_remaining(final_time - m_time);

   if (dt_stable <= time_remaining) {
      const int n(static_cast<int>(ceil(time_remaining / dt_stable)));
      m_dt = time_remaining / n;
   }
   else {
      m_dt = time_remaining * (1.0 + 10*std::numeric_limits<real>::epsilon());
   }
}


void
Simulation::getFromRestart(
   const HDF_DataBase& a_db)
{
   // Check for compatibility of both Loki and Overture with the versions of
   // each used to create the restart file.
   int major_version;
   if (a_db.get(major_version, "major version") == 0) {
      int minor_version, patch_level;
      a_db.get(minor_version, "minor version");
      a_db.get(patch_level, "patch level");
   }
   else {
      aString overture_version;
      a_db.get(overture_version, "OvertureVersion");
      if (!overture_version.matches("Overture.snapshot")) {
         OV_ABORT("This version of Loki is not compatible with the version that created the restart files.");
      }
   }

   // Our initial time is the final time of the run we restart from.
   a_db.get(m_time, "time");
   m_last_save = m_time / m_save_times;  // Cast to integer
   m_last_seq = m_time / m_seq_times;  // Cast to integer
   if (m_time == m_max_time) {
      OV_ABORT("Time read from restart file matches final time. Reached end of simulation.");
   }
}


void
Simulation::putToRestart(
   HDF_DataBase& a_db,
   real a_time)
{
   NULL_USE(a_time);

   // Write the necessary restart data.  Nobody seems to actually need
   // time step, CFL, or tf.
   a_db.put(LOKI_VERSION_MAJOR, "major version");
   a_db.put(LOKI_VERSION_MINOR, "minor version");
   a_db.put(LOKI_VERSION_PATCHLEVEL, "patch level");
   a_db.put(m_time, "time");
   a_db.put(m_dt, "time step");
   a_db.put(m_cfl, "CFL");
   a_db.put(m_max_time, "tf");
}

} // end namespace Loki
