/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _EXPANSION_SCHEDULE_H_
#define _EXPANSION_SCHEDULE_H_

#include "ParallelArray.H"
#include "ProblemDomain.H"
#include <mpi.h>

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
    * @param[in] a_phase_space_dist A ParallelArray distributed over phase space
    * @param[in] a_domain The problem domain for the species performing the
    *                     expansion.
    * @param[in] a_comm The communicator for the species performing the
    *                   expansion.
    */
   ExpansionSchedule(
      const ParallelArray& a_phase_space_dist,
      const ProblemDomain& a_domain,
      MPI_Comm a_comm);

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
      const ParallelArray& a_global_config_array,
      ParallelArray& a_local_phase_array);

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

   // Precomputes the phase space head processors a given configuration space
   // processor sends data to.
   void
   computeHeadProcTargets(
      const ParallelArray& a_global_config_array);

   // Precomputes the configuration space processors a given phase space
   // head processor receives data from.
   void
   computeConfigRecvTargets(
      const ParallelArray& a_global_config_array,
      const ParallelArray& a_local_phase_array);

   // Send data from configuration space processors to the phase space head
   // processors.
   void
   sendConfigSpaceDataToHeadProcTargets(
      const ParallelArray& a_global_config_array);

   // A phase space processor receives data from configuration space processors.
   void
   headProcsRecvConfigSpaceData(
      const ParallelArray& a_global_config_array,
      ParallelArray& a_local_phase_array);

   // Construct the communicators for phase space processors covering the same
   // piece of configuration space.
   void
   constructExpansionCommunicator(
      const ProblemDomain& a_domain,
      MPI_Comm a_comm);

   // A ParallelArray distributed over phase space. 
   const ParallelArray& m_phase_space_dist;

   // If true, the communication pattern must be precompute.
   bool m_need_communication_pattern;

   // The phase space head processors that a given configuration space processor
   // sends data to.
   vector<int> m_send_targets;

   // The box describing the part of the data sent from a given configuration
   // space processor to a phase space head processor.
   vector<ParallelArray::Box> m_send_boxes;

   // The configuration space processors that the a given phase space processor
   // receives data from.
   vector<int> m_recv_targets;

   // The box describing the part of the data sent from each configuration space
   // processor.
   vector<ParallelArray::Box> m_recv_boxes;

   // Communicator for phase space processors covering the same piece of
   // configuration space.
   MPI_Comm m_expansion_comm;

   // The processor id on m_expansion_comm.
   int m_expansion_comm_id;

   // Tag for send of configuration space data to phase space processors.
   static const int s_TAG;
};

} // end namespace Loki

#endif
