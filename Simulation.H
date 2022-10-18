/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <string>
#include <vector>

#include "LokiInputParser.H"
#include "Loki_Defines.H"
#include "tbox/Pointer.H"
#include "System.H"

#include "Serializable.H"

namespace Loki {

/**
 * Simulation driver class.
 *
 * Responsible for generic initialization, advancement, and finalization of
 * simulation.  This class specifically controls the monitoring of the number
 * of steps and of the simulation time in order to execute plotfile and
 * checkpoint file creation as well as simulation termination.
 *
 * Required input keys and data types:
 *
 *    -\b max_step
 *      positiveinteger value maximum number of simulation steps
 *
 *    -\b max_time
 *      positive double value final time at which simulation terminates
 *
 *    -\b fixed_dt
 *      positive double value, if set, for the time step to be used
 *
 *    -\b max_dt_grow
 *      positive double value of the multiple of the current time step to use
 *      for the next time step.  Must be grater than 1.  Multually exclusive
 *      with fixed_dt
 *
 *    -\b initial_dt_fraction
 *      positive double value fraction of the initial stable time step to use.
 *      Must be less than 1.  Multually exclusive with fixed_dt
 *
 *    -\b checkpoint_interval
 *      integer value specifying the number of steps between checkpoint dumps
 *
 *    -\b checkpoint_prefix
 *      string used as prefix for checkpoint file names ["chk"]
 *
 *    -\b plot_interval
 *      integer value specifying the number of steps between plot dumps
 *
 *    -\b plot_prefix
 *      string used as prefix for plot file names ["plt"]
 *
 *    -\b verbosity
 *      integer flag specifying the verbosity of logging output; zero turns
 *      of output, and increasing values produce more detailed logs
 *
 * The following represents a sample input entry:
 *
 * @verbatim
 *  verbosity = 1
 *  cfl = .8
 *  final_time = 100
 *  max_step = 1000
 * @endverbatim
*/
class Simulation : public Serializable
{
public:
   /**
    * @brief The speed of light for the simulation.
    */
   static double s_LIGHT_SPEED;

   /**
    * @brief If true use relativistic calculations.
    */
   static bool s_DO_RELATIVITY;

   /**
    * @brief Loki major version number.
    */
   static const int LOKI_VERSION_MAJOR;

   /**
    * @brief Loki minor version number.
    */
   static const int LOKI_VERSION_MINOR;

   /**
    * @brief Loki patch level number.
    */
   static const int LOKI_VERSION_PATCHLEVEL;

   /**
    * @brief Constructor.
    *
    * @param[in] a_pp Input database from which to intialize.
    */
   explicit Simulation(
      LokiInputParser& a_pp);

   /**
    * @brief Destructor.
    */
   virtual
   ~Simulation();

   /**
    * @brief Returns true if simulation is incomplete, false otherwise.
    *
    * @return True if simulation is incomplete.
    */
   bool
   notDone()
   {
      return ((m_step < m_max_step) && (m_time < m_max_time));
   }

   /**
    * @brief Advances the simulation a single timestep and writes
    *        restart/diagnostics as necessary.
    */
   void
   advance();

   /**
    * @brief Cleans up and closes open resources at completion of simulation.
    */
   void
   finalize();

   //@{

   //! @name Implementation of Serializable interface.

   /**
    * @brief Read the simulation parameters from the restart database.
    *
    * @param[in] a_reader The restart file reader.
    */
   virtual
   void
   getFromRestart(
      RestartReader& a_reader);

   /**
    * @brief Write the simulation parameters to the restart database.
    *
    * @param[in] a_writer The restart file writer.
    * @param[in] a_time Current simulation time.
    */
   virtual
   void
   putToRestart(
      RestartWriter& a_writer,
      double a_time);

  //@}

private:
   // Embedded class to hold collision diagnostic parameters.
   class CollisionDiagnosticParameters
   {
   public:
      // Default constructor.
      CollisionDiagnosticParameters() {
         m_saveSteps = m_writeSteps = m_lastSave = m_lastSaveStep =
            m_lastWrite = m_lastWriteStep = -1;
         m_saveTimes = m_writeTimes = -1.0;
         m_on = false;
      }

      // Destructor
      ~CollisionDiagnosticParameters() {
      }

      // Returns true if collision sequences should be accumulated at this
      // timestep.
      bool
      doCollSeqAccum(double a_time,
         int a_step) {
         if (m_on) {
            if (m_writeTimes > 0.0) {
               if (a_time < (m_lastWrite+1)*m_writeTimes) {
                  return false;
               }
               else {
                  ++m_lastWrite;
                  return true;
               }
            }
            else {
               if (m_lastWriteStep >= 0 &&
                   a_step - m_lastWriteStep < m_writeSteps) {
                  return false;
               }
               else {
                  m_lastWriteStep = a_step;
                  return true;
               }
            }
         }
      }

      // Returns true if collision plot data should be saved at this timestep.
      bool
      doCollPlot(double a_time,
         int a_step) {
         if (m_on) {
            if (m_saveTimes > 0.0) {
               if (a_time < (m_lastSave+1)*m_saveTimes) {
                  return false;
               }
               else {
                  ++m_lastSave;
                  return true;
               }
            }
            else {
               if (m_lastSaveStep >= 0 &&
                   a_step - m_lastSaveStep < m_saveSteps) {
                  return false;
               }
               else {
                  m_lastSaveStep = a_step;
                  return true;
               }
            }
         }
      }

      // Read input parameters from the top level of the input database.
      void
      parseParameters(
         LokiInputParser& a_pp) {
         if (a_pp.contains("coll_save_steps")) {
            a_pp.get("coll_save_steps", m_saveSteps);
            m_on = true;
         }

         if (a_pp.contains("coll_save_times")) {
            a_pp.get("coll_save_times", m_saveTimes);
            m_on = true;
         }

         if (a_pp.contains("coll_seq_write_steps")) {
            a_pp.get("coll_seq_write_steps", m_writeSteps);
            m_on = true;
         }

         if (a_pp.contains("coll_seq_write_times")) {
            a_pp.get("coll_seq_write_times", m_writeTimes);
            m_on = true;
         }

         if (m_on) {
            if ((m_saveSteps > 0 && m_saveTimes > 0.0) ||
                (m_saveSteps < 0 && m_saveTimes < 0.0)) {
               LOKI_ABORT("Invalid parameters for collision 2D field saves");
            }

            if ((m_writeSteps > 0 && m_writeTimes > 0.0) ||
                (m_writeSteps < 0 && m_writeTimes < 0.0)) {
               LOKI_ABORT("Invalid parameters for collision sequence writes");
            }
         }
      }

      bool
      on() {
         return m_on;
      }

   private:
      // "write" is for history saves; "save" is for 2D spatial field saves.
      int m_saveSteps, m_writeSteps;
      double m_saveTimes, m_writeTimes;
      bool m_on;

      int m_lastSave, m_lastSaveStep, m_lastWrite, m_lastWriteStep;
   };

   // Unimplemented default constructor.
   Simulation();

   // Unimplemented copy constructor.
   Simulation(
      const Simulation& other);

   // Unimplemented assignment operator.
   Simulation&
   operator = (
      const Simulation& rhs);

   // Obtain control parameters from the top level of the input database.
   void
   parseParameters(
      LokiInputParser& a_ppsim);

   // Print class' parameters.
   void
   printParameters() const;

   // Write simulation data to plot file.
   void
   writePlotFile();

   // Accumulate 1D sequence data.
   void
   accumulateSequences();

   // Write data to checkpoint file.
   void
   writeCheckpointFile();

   // Select time step based on stable dt and final simulation time.
   void
   selectTimeStep();

   // Spatial order of solution
   int m_spatial_solution_order;

   // Temporal order of solution.
   int m_temporal_solution_order;

   // If non-zero, prints out more information.
   int m_verbosity;

   // The current simulation step.
   int m_step;

   // The maximum number of simulation steps to take.
   int m_max_step;

   // The current simulation time.
   double m_time;

   // The maximum simulation time.
   double m_max_time;

   // The current stable time step.
   double m_dt;

   // Time step safety factor.
   double m_cfl;

   ////////////////////////////////////
   // Belong in a diagnostic manager
   //
   // The frequency in simulation time for writing plot data.
   double    m_save_times;
   //
   // If false, no time history or plot data will be collected or written no
   // matter what frequencies have been set.
   bool    m_save_data;
   bool m_save_coll_data;
   //
   // The base name of the time history files.
   string m_time_hist_file_name;
   //
   // The base name of the collision time history files.
   string m_coll_time_hist_file_name;
   //
   // Probe locations as a fraction of the domain in (x,y)
   vector<vector<double> > m_probes;
   //
   // The frequency in simulation time for the collection of time history
   // points.
   double m_seq_times;
   //
   // The number of probes.
   int m_num_probes;
   //
   // The number of time history points accumulated so far.
   int m_saved_seq, m_saved_seq_coll;
   // The last time history point index, including restarts.
   int m_last_seq;
   //
   // The number of plot cycles that have occurred so far.
   int m_saved_save, m_saved_save_coll;
   // The last plot index, including restarts.
   int m_last_save;

   // Sequences are saved every m_seq_times starting from 0.  This is true even
   // on restart when m_time does not start from 0.  So we need a 0 based time
   // to properly trigger sequence saving.  If m_time is used, the sequences
   // are saved far too often on restart and consequently the code runs MUCH
   // slower.
   double m_time_from_start_of_run;

   // If true then this is a VM problem.
   bool m_do_maxwell;

   // The system we are solving.
   tbox::Pointer<System> m_system;

   // Information about requested collision operator diagnostics.
   CollisionDiagnosticParameters m_coll_diag;
};

} // end namespace Loki

#endif
