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
#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include "Overture.h"
#include "Ogshow.h"

namespace Loki {

/**
 * Abstract base class that defines the interface of the system we are solving
 * either Vlasov-Poisson or Vlasov-Maxwell.
 */
class System
{
public:
   /**
    * @brief Default Constructor
    */
   System();

   /**
    * @brief Destructor
    */
   virtual
   ~System();

   /**
    * @brief Computes and returns a stable time step estimate.
    */
   virtual
   Real
   stableDt() = 0;

   /**
    * @brief Advance the solution by a single time step.  Returns the
    * simulation time advanced to.
    *
    * @param[in] a_time Current simulation time.
    * @param[in] a_dt Current simulation time step.
    */
   virtual
   Real
   advance(
      Real a_time,
      Real a_dt) = 0;

   /**
    * @brief Initialize the system given the current simulation time.
    *
    * @param[in] a_is_from_restart If true simulation is being run from
    *                              restart.
    * @param[in] a_time Current simulation time.
    * @param[in] a_dt Current simulation time step.
    */
   virtual
   void
   initialize(
      bool a_is_from_restart,
      real a_time,
      real a_dt) = 0;

   /**
    * @brief Returns number of fields to be plotted.
    */
   virtual
   int
   getNumFrameSeries() = 0;

   /**
    * @brief Write the field and accumulated sequences to the Overture show
    * file.
    *
    * @param[in] a_time Current simulation time.
    * @param[in] a_dt Current simulation time step.
    * @param[in] a_probes Coordinates of probes.
    * @param[in] a_num_probes Number of probes.
    * @param[in] a_saved_seq Number of sequence time steps saved.
    * @param[in] a_saved_save Number of plots saved.
    * @param[out] a_show The Overture show file.
    */
   virtual
   void
   plot(
      real a_time,
      real a_dt,
      const RealArray& a_probes, 
      int a_num_probes,
      int a_saved_seq,
      int& a_saved_save,
      Ogshow& a_show) = 0;

   /**
    * @brief Accumulate another time point for the sequences.
    *
    * @param[in] a_time Current simulation time.
    * @param[in] a_probes Coordinates of probes.
    * @param[in] a_num_probes Number of probes.
    * @param[out] a_saved_seq Number of sequence time steps saved.
    */
   virtual
   void
   accumulateSequences(
      real a_time,
      const RealArray& a_probes, 
      int a_num_probes,
      int& a_saved_seq) = 0;

   /**
    * @brief Returns the number of cells in the simulation.
    */
   virtual
   long int
   problemSize() const = 0;

   /**
    * @brief Updates ghost cells in KineticSpecies distributions.
    */
   virtual
   void
   updateGhosts() = 0;
};

}  // end namespace Loki

#endif
