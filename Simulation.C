/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "Simulation.H"

#include <limits>
#include <sstream>

#include "TimerManager.H"
#include "RestartManager.H"
#include "VMSystem.H"
#include "VPSystem.H"

namespace Loki {

const int Simulation::LOKI_VERSION_MAJOR = 3;
const int Simulation::LOKI_VERSION_MINOR = 0;
const int Simulation::LOKI_VERSION_PATCHLEVEL = 1;
double Simulation::s_LIGHT_SPEED;
bool Simulation::s_DO_RELATIVITY;

void
Simulation::printParameters() const
{
   // Print the general simulation parameters.
   Loki_Utilities::printF("\n#*#*# General Simulation Parameters #*#*#\n");
   Loki_Utilities::printF("Loki version %d.%d.%d\n",
      LOKI_VERSION_MAJOR,
      LOKI_VERSION_MINOR,
      LOKI_VERSION_PATCHLEVEL);
   Loki_Utilities::printF("spatial solution order  = %d\n",
      m_spatial_solution_order);
   Loki_Utilities::printF("temporal solution order = %d\n",
      m_temporal_solution_order);
   ////////////////////////////////////
   // A time step manager will eventually write these out
   Loki_Utilities::printF("current time            = %e\n", m_time);
   Loki_Utilities::printF("time from start of run  = %e\n",
      m_time_from_start_of_run);
   Loki_Utilities::printF("final time              = %e\n", m_max_time);
   Loki_Utilities::printF("cfl                     = %e\n", m_cfl);
   if (m_do_maxwell || s_DO_RELATIVITY) {
      Loki_Utilities::printF("light speed             = %e\n", s_LIGHT_SPEED);
   }
   ////////////////////////////////////

   ////////////////////////////////////
   // A diagnostic manager will eventually write these out
   if (m_save_data) {
      Loki_Utilities::printF("Save times              = %e\n", m_save_times);
   }
   else {
      Loki_Utilities::printF("Not saving data\n");
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
         m_time_hist_file_name);

      // Except for the write at the first time step, write out all timer info.
      if (m_step != 0) {
         Loki_Utilities::printF("\n\n************** PRINT TIMING INFO **************\n");
         Loki_Utilities::printF("t=%3.3e, dt=%3.3e: %i time steps using %i processors\n",
            m_time,
            m_dt,
            m_step,
            Loki_Utilities::s_num_procs);
         timers->print(m_step, m_system->problemSize());

         Loki_Utilities::printF("****************** Exclusive Timers *******************\n");
         Loki_Utilities::printF("t=%3.3e, dt=%3.3e: %i time steps using %i processors\n",
            m_time,
            m_dt,
            m_step,
            Loki_Utilities::s_num_procs);
         timers->print_exclusive(m_step, m_system->problemSize());
         Loki_Utilities::printF("*************** END TIMING INFO ***************\n\n\n");
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

      // The particle boundary conditions are set at the top of the time step
      // loop so at this point they are not set and must be for the proper
      // particle time histories to be output.
      m_system->updateGhosts(true);

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
      Loki_Utilities::printF("\n****************************************************\n");
      Loki_Utilities::printF("Wall clock at restart dump: %s",
         asctime(localtime(&wall_time)));
      Loki_Utilities::printF("****************************************************\n\n");

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
   LokiInputParser& a_pp)
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
      m_save_coll_data(false),
      m_probes(2),
      m_seq_times(1.0),
      m_num_probes(-1),
      m_saved_seq(0),
      m_saved_seq_coll(0),
      m_saved_save(0),
      m_saved_save_coll(0),
      m_last_seq(0),
      m_last_save(0),
   ////////////////////////////////////
      m_time_from_start_of_run(0.0),
      m_system(0)
{
   // Create all timers that always apply.
   TimerManager* timers(TimerManager::getManager());
   timers->createTimer("setup");
   timers->createTimer("output");
   timers->createTimer("Vlasov");
   timers->createTimer("phys to phase");
   timers->createTimer("blowout");
   timers->createTimer("field driver");
   timers->createTimer("contraction");
   timers->createTimer("reduction");
   timers->createTimer("reduction to phase");
   timers->createTimer("summation");
   timers->createTimer("BC (Vlasov)");
   timers->createTimer("parallel ghost");
   timers->createTimer("collisions");
   timers->createTimer("collision kernel");
   timers->createTimer("collision comm");
   timers->createTimer("krook");
   timers->createTimer("other");
   timers->createTimer("AddSolnData");
   timers->createTimer("CopySolnData");
   timers->createTimer("ZeroSolnData");
   timers->createTimer("Tracking Particles");
   timers->createTimer("Noisy Particles");

   timers->startTimer("setup");

   // Parse all general simulation related input.
   parseParameters(a_pp);

   // See which system we're running and construct it.
   if (m_do_maxwell) {
      m_system = new VMSystem(a_pp,
         m_num_probes,
         m_spatial_solution_order,
         m_temporal_solution_order,
         m_coll_diag.on());
   }
   else {
      m_system = new VPSystem(a_pp,
         m_num_probes,
         m_spatial_solution_order,
         m_temporal_solution_order,
         m_coll_diag.on());
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

   // Initialize the System.  The name is kind of odd given that we might have
   // just restored everything from a restart file which initializes a buch of
   // stuff.
   m_system->initialize(restart_manager->isFromRestart(),
      m_time,
      m_dt,
      m_coll_diag.on());

   // Now that everything is set up, print the parameters of the problem being
   // run.
   m_system->printParameters();
   printParameters();
   restart_manager->print();


   ////////////////////////////////////
   // Belongs in a diagnostic manager
   if (m_save_data) {
      m_time_hist_file_name =
         restart_manager->restartWritePath() + ".time_hists";
   }

   if (m_save_coll_data) {
      m_coll_time_hist_file_name =
         restart_manager->restartWritePath() + "_coll.time_hists";
   }
   ////////////////////////////////////


   // Everything is now built and tidied up so collect a time history at the
   // initial time, write the plot data, and write the restart data.
   accumulateSequences();
   writePlotFile();

   writeCheckpointFile();

   fflush(0);
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
      Loki_Utilities::printF("Step %d completed, simulation time is %f, dt is %f\n",
         m_step, m_time, m_dt);
      fflush(0);
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

   if (m_save_coll_data) {
      if (m_coll_diag.doCollSeqAccum(m_time, m_step)) {
         m_system->accumulateCollisionSequences(m_time,
            m_dt,
            0,
            m_saved_seq_coll);
      }
      if (m_coll_diag.doCollPlot(m_time, m_step)) {
         m_system->plotColl(m_time,
            m_dt,
            m_saved_seq_coll,
            m_saved_save_coll,
            m_coll_time_hist_file_name);
      }
   }
   
   // The check for whether it is time to write the next checkpoint file is
   // in this function as opposed to the others immediately above.
   writeCheckpointFile();

   fflush(0);
}


void
Simulation::finalize()
{
   // We always want a final time step checkpoint file.
   writeCheckpointFile();
   fflush(0);
}


void
Simulation::parseParameters(
   LokiInputParser& a_pp)
{
   // Get the spatial and temporal solution orders.
   a_pp.query("spatial_solution_order", m_spatial_solution_order);
   if (m_spatial_solution_order != 4 && m_spatial_solution_order != 6) {
      LOKI_ABORT("spatial_solution_order must be 4 or 6");
   }
   a_pp.query("temporal_solution_order", m_temporal_solution_order);
   if (m_temporal_solution_order != 4 && m_temporal_solution_order != 6) {
      LOKI_ABORT("temporal_solution_order must be 4 or 6");
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
   if (a_pp.contains("number_of_probes")) {
      a_pp.query("number_of_probes", m_num_probes);
      m_probes[X1].resize(m_num_probes);
      m_probes[X2].resize(m_num_probes);
      for (int i(0); i < m_num_probes; ++i) {
         ostringstream oss;
         oss << "probe." << i+1;
         LokiInputParser prpp(oss.str().c_str());

         vector<double> loc(2);
         prpp.getarr("location", loc, 0, 2);
         m_probes[X1][i] = loc[X1];
         m_probes[X2][i] = loc[X2];
      }
   }
   else {
      m_num_probes = 1;
      m_probes[X1].resize(m_num_probes);
      m_probes[X2].resize(m_num_probes);
      m_probes[X1][0] = 0.5;
      m_probes[X1][0] = 0.5;
   }

   // Get the frequency in simulation time for the collection of time history
   // points.
   a_pp.query("sequence_write_times", m_seq_times);

   // Get the frequency in simulation time for writing plot data.
   a_pp.query("save_times", m_save_times);

   // User can turn off collection/writing of plot data.
   string sdata("true");
   a_pp.query("save_data", sdata);
   m_save_data = sdata.compare("true") == 0 ? true : false;

   // User can turn on collision diagnostics.
   string s_coll_data = "false";
   a_pp.query("save_coll_data", s_coll_data);
   m_save_coll_data = s_coll_data.compare("true") == 0 ? true : false;

   m_coll_diag.parseParameters(a_pp);

   // Find out if we're running a VP or VM problem.
   string sys_str("poisson");
   a_pp.query("sys_type", sys_str);
   if (sys_str.compare("poisson") == 0) {
      m_do_maxwell = false;
   }
   else if (sys_str.compare("maxwell") == 0) {
      m_do_maxwell = true;
   }
   else {
      LOKI_ABORT("Unknown sys_type input.");
   }

   // Find out if we're running the problem relativistically.
   string rel_str("false");
   a_pp.query("do_relativity", rel_str);
   s_DO_RELATIVITY = rel_str.compare("false") == 0 ? false : true;

   // Get the speed of light.
   if (a_pp.contains("light_speed")) {
      a_pp.query("light_speed", s_LIGHT_SPEED);
   }
   else if (m_do_maxwell || s_DO_RELATIVITY) {
      LOKI_ABORT("Must supply light speed.");
   }
   ////////////////////////////////////
}


void
Simulation::selectTimeStep()
{
   // We fiddle with dt so that it honors the requested plot frequency.  This is
   // essential to allow Fourier transforms of the plot data.
  double next_time = (m_last_save + 1) * m_save_times;
  //(floor(m_time/m_save_times+100*numeric_limits<double>::epsilon()) + 1) *m_save_times;
  
   const double final_time(min(next_time, m_max_time));

   const double dt_stable(m_system->stableDt() * m_cfl);

   const double time_remaining(final_time - m_time);

   if (dt_stable <= time_remaining) {
      const int n(static_cast<int>(ceil(time_remaining / dt_stable)));
      m_dt = time_remaining / n;
   }
   else {
      m_dt = time_remaining * (1.0 + 10*numeric_limits<double>::epsilon());
   }
}


void
Simulation::getFromRestart(
   RestartReader& a_reader)
{
   // Check that the number of processes that generated the data is the same as
   // the number of processors trying to read it.
   int generating_processors;
   a_reader.readIntegerValue("generating processes", generating_processors);
   if (Loki_Utilities::s_num_procs != generating_processors) {
      LOKI_ABORT("Attemping to restart a problem on a different number of processors.");
   }

   // Check for compatibility of Loki with the version used to create the
   // restart file.
   int major_version, minor_version, patch_level;
   a_reader.readIntegerValue("major version", major_version);
   a_reader.readIntegerValue("minor version", minor_version);
   a_reader.readIntegerValue("patch level", patch_level);
   if (major_version < 3 || minor_version < 0 || patch_level < 0) {
      LOKI_ABORT("Restart data generated with incompatible version of Loki.");
   }

   // Our initial time is the final time of the run we restart from.
   a_reader.readDoubleValue("time", m_time);
   m_last_save = m_time / m_save_times;  // Cast to integer
   m_last_seq = m_time / m_seq_times;  // Cast to integer
   if (m_time == m_max_time) {
      LOKI_ABORT("Time read from restart file matches final time. Reached end of simulation.");
   }
}


void
Simulation::putToRestart(
   RestartWriter& a_writer,
   double a_time)
{
   NULL_USE(a_time);

   // Write the necessary restart data.  Nobody seems to actually need
   // time step, CFL, or tf.
   bool write_data = Loki_Utilities::s_my_id == 0;
   a_writer.writeIntegerValue("generating processes",
      Loki_Utilities::s_num_procs,
      write_data);
   a_writer.writeIntegerValue("major version", LOKI_VERSION_MAJOR, write_data);
   a_writer.writeIntegerValue("minor version", LOKI_VERSION_MINOR, write_data);
   a_writer.writeIntegerValue("patch level",
      LOKI_VERSION_PATCHLEVEL,
      write_data);
   a_writer.writeDoubleValue("time", m_time, write_data);
   a_writer.writeDoubleValue("time step", m_dt, write_data);
   a_writer.writeDoubleValue("CFL", m_cfl, write_data);
   a_writer.writeDoubleValue("tf", m_max_time, write_data);
}

} // end namespace Loki
