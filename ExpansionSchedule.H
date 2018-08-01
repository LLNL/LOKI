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
#ifndef _EXPANSION_SCHEDULE_H_
#define _EXPANSION_SCHEDULE_H_

#include "Overture.h"
#include "ParallelUtility.h"

#include "tbox/Box.H"

#include <vector>

namespace Loki {

/**
 * Performs the communication necessary to export a configuration space quantity
 * to phase space.  The configuration space processor(s) communicate to a
 * single phase space processor covering a given piece of configuration space.
 * These phase space processors then communicate the info to the other phase
 * space processors covering the same piece of configuration space.
 */
class ExpansionSchedule
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_box The box representing the piece of configuration space
    *                  associated with this object's processor.
    * @param[in] a_number_of_species The total number of KineticSpecies.
    * @param[in] a_config_space_id The id for the piece of configuration space
    *                              associated with this object's processor.
    * @param[in] a_comm Global communicator.
    */
   ExpansionSchedule(
      const tbox::Box& a_box,
      int a_number_of_species,
      int a_config_space_id,
      const MPI_Comm& a_comm);

   /**
    * @brief Destructor.
    */
   ~ExpansionSchedule();

   /**
    * @brief Perform the expansion.
    *
    * @param[in] a_global_config_array The global array defined on the
    *                                  configuration space processors.
    * @param[out] a_local_phase_array Part of a_global_config_array owned by
    *                                 a phase space processor.
    */
   void
   execute(
      const realArray& a_global_config_array,
      RealArray& a_local_phase_array);

private:
   // Unimplemented default constructor.
   ExpansionSchedule();

   // Unimplemented copy constructor.
   ExpansionSchedule(
      const ExpansionSchedule& other);

   // Unimplemented assignment operator.
   ExpansionSchedule&
   operator = (
      const ExpansionSchedule& a_rhs);

   // Returns true if this processor only deals with quantities defined in
   // phase space.
   bool isPhaseSpaceProcessor()
   {
      return m_is_phase_space_processor;
   }

   // True if this processor is a phase space processor.
   bool m_is_phase_space_processor;

   // Phase space extents of each processor.
   Array<IndexBox> m_phase_copy_to;

   // The communicator for a phase space processor assoicated with a specific
   // piece of configuration space.
   MPI_Comm m_expansion_comm;
};

} // end namespace Loki

#endif
