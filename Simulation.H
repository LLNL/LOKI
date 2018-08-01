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
#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <string>
#include <vector>

#include "Overture.h"
#include "ParmParse.H"
#include "Ogshow.h"
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
 *      positive real value final time at which simulation terminates
 *
 *    -\b fixed_dt
 *      positive real value, if set, for the time step to be used
 *
 *    -\b max_dt_grow
 *      positive real value of the multiple of the current time step to use for
 *      the next time step.  Must be grater than 1.  Multually exclusive with
 *      fixed_dt
 *
 *    -\b initial_dt_fraction
 *      positive real value fraction of the initial stable time step to use.
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
    * @brief Constructor.
    *
    * @param[in] a_pp Input database from which to intialize.
    */
   explicit Simulation(
      ParmParse& a_pp);

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
    * restart/diagnostics as necessary.
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
    * @param[in] a_db The restart database.
    */
   virtual
   void
   getFromRestart(
      const HDF_DataBase& a_db);

   /**
    * @brief Write the simulation parameters to the restart database.
    *
    * @param[in] a_db The restart database.
    * @param[in] a_time Current simulation time.
    */
   virtual
   void
   putToRestart(
      HDF_DataBase& a_db,
      real a_time);

  //@}

private:
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
      ParmParse& a_ppsim);

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
   real m_time;

   // The maximum simulation time.
   real m_max_time;

   // The current stable time step.
   real m_dt;

   // Time step safety factor.
   real m_cfl;

   ////////////////////////////////////
   // Belong in a diagnostic manager
   //
   // The frequency in simulation time for writting plot data.
   real    m_save_times;
   //
   // If false, no time history or plot data will be collected or written no
   // matter what frequencies have been set.
   bool    m_save_data;
   //
   // The base name for plot files.  Overture appends its own info to this.
   aString m_show_file_name;
   //
   // Probe locations as a fraction of the domain in (x,y)
   RealArray m_probes;
   //
   // The frequency in simulation time for the collection of time history
   // points.
   real m_seq_times;
   //
   // The number of probes.
   int m_num_probes;
   //
   // The number of time history points accumulated so far.
   int m_saved_seq;
   // The last time history point index, including restarts.
   int m_last_seq;
   //
   // The number of plot cycles that have occurred so far.
   int m_saved_save;
   // The last plot index, including restarts.
   int m_last_save;
   //
   // The Overture show file object for plot data.
   Ogshow m_show;
   ////////////////////////////////////

   // Sequences are saved every m_seq_times starting from 0.  This is true even
   // on restart when m_time does not start from 0.  So we need a 0 based time
   // to properly trigger sequence saving.  If m_time is used, the sequences
   // are saved far too often on restart and consequently the code runs MUCH
   // slower.
   real m_time_from_start_of_run;

   // The system we are solving.
   tbox::Pointer<System> m_system;

   // Loki versioning info.
   static const int LOKI_VERSION_MAJOR;
   static const int LOKI_VERSION_MINOR;
   static const int LOKI_VERSION_PATCHLEVEL;
};

} // end namespace Loki

#endif
