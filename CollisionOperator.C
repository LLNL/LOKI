/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "CollisionOperator.H"
#include "KineticSpecies.H"

namespace Loki {

const int CollisionOperator::s_DIAGNOSTIC_WORK_SIZE = 4;

CollisionOperator::CollisionOperator(
   const ProblemDomain& a_domain,
   int a_num_dparams,
   int a_num_iparams)
   : m_domain(a_domain),
     m_dparameters(a_num_dparams),
     m_iparameters(a_num_iparams)
{
}


CollisionOperator::~CollisionOperator()
{
}


void
CollisionOperator::initialize(
   KineticSpecies* a_species)
{
   m_proc_lo = a_species->procLo();
   m_proc_hi = a_species->procHi();
   m_comm = a_species->communicator();
   m_data_box = a_species->dataBox();
   m_interior_box = a_species->interiorBox();
   m_velocities = a_species->velocities();
}

} // end namespace Loki
