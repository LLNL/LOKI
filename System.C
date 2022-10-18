/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "System.H"
#include "KineticSpeciesPtrVect.H"

namespace Loki {

System::System(
   LokiInputParser& a_pp,
   int a_spatial_solution_order,
   int a_temporal_solution_order)
   : m_cdim(CDIM),
     m_pdim(PDIM),
     m_cfg_domain(new ProblemDomain(m_cdim, a_spatial_solution_order, a_pp)),
     m_length_seq(8),
     m_length_seq_coll(8),
     m_spatial_solution_order(a_spatial_solution_order),
     m_temporal_solution_order(a_temporal_solution_order)
{
}


System::~System()
{
}


void
System::PopulateInterspeciesData(
   KineticSpeciesPtrVect& ks)
{
   const int numKS = static_cast<int>(ks.size());

   if (numKS < 1) {
      LOKI_ABORT("No species.");
   }

   int speciesHead[numKS];
   double speciesMass[numKS];

   // It is assumed here that the species are in the same order on all
   // processes.  To check this, root could broadcast its species name list for
   // all processes to compare.

   for (int s = 0; s < numKS; ++s) {
      speciesHead[s] = ks[s]->headRank();
      speciesMass[s] = ks[s]->mass();

      //cout << rank << ": species " << s << " head " << speciesHead[s] << endl;
   }

   for (int s = 0; s < numKS; ++s) {
      ks[s]->SetSpeciesHeadsAndMass(numKS, speciesHead, speciesMass);
   }
}

} // end namespace Loki
