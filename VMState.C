/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "VMState.H"

namespace Loki {

VMState::VMState()
{
}

VMState::VMState(
   const VMState& other)
{
   // Clone other into this.
   m_kspv.resize(other.m_kspv.size());
   for (int s(0); s < static_cast<int>(other.m_kspv.size()); ++s) {
      m_kspv[s] = other.m_kspv[s]->clone();
   }
   m_maxwell = other.m_maxwell->clone();
}

VMState::~VMState()
{
}

} // end namespace Loki
