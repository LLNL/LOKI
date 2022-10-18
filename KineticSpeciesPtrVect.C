/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "KineticSpeciesPtrVect.H"

namespace Loki {

KineticSpeciesPtrVect::KineticSpeciesPtrVect()
{
}

KineticSpeciesPtrVect::KineticSpeciesPtrVect(
   const KineticSpeciesPtrVect& other)
{
   // Add clones of the contents of other to this.
   resize(other.size());
   for (int s(0); s < static_cast<int>(other.size()); ++s) {
      (*this)[s] = other[s]->clone();
   }
}

KineticSpeciesPtrVect::~KineticSpeciesPtrVect()
{
}

KineticSpeciesPtrVect::Iterator::Iterator(
   KineticSpeciesPtrVect& a_in,
   bool a_from_end)
{
   // Store in m_local_indices the indices of the all the species in a_in that
   // are local to this processor.  These are the species over which this
   // iterator will operate.
   for (int s(0); s < static_cast<int>(a_in.size()); ++s) {
      if (a_in[s]->isInRange(Loki_Utilities::s_my_id)) {
         m_local_indicies.push_back(s);
      }
   }

   // Point to the proper end of m_local_indicies.
   if (a_from_end) {
      m_current = m_local_indicies.end();
   }
   else {
      m_current = m_local_indicies.begin();
   }
   m_from_end = a_from_end;

   // Set up pointer to the iterated container.
   m_base = &a_in;
}


KineticSpeciesPtrVect::Iterator::Iterator(
   const KineticSpeciesPtrVect::Iterator& a_other)
{
   // Just copy the data members.
   m_local_indicies = a_other.m_local_indicies;
   m_current = a_other.m_current;
   m_from_end = a_other.m_from_end;
   m_base = a_other.m_base;
}

KineticSpeciesPtrVect::ConstIterator::ConstIterator(
   const KineticSpeciesPtrVect& a_in,
   bool a_from_end)
{
   // Store in m_local_indices the indices of the all the species in a_in that
   // are local to this processor.  These are the species over which this
   // iterator will operate.
   for (int s(0); s < static_cast<int>(a_in.size()); ++s) {
      if (a_in[s]->isInRange(Loki_Utilities::s_my_id)) {
         m_local_indicies.push_back(s);
      }
   }

   // Point to the proper end of m_local_indicies.
   if (a_from_end) {
      m_current = m_local_indicies.end();
   }
   else {
      m_current = m_local_indicies.begin();
   }
   m_from_end = a_from_end;

   // Set up pointer to the iterated container.
   m_base = &a_in;
}


KineticSpeciesPtrVect::ConstIterator::ConstIterator(
   const KineticSpeciesPtrVect::ConstIterator& a_other)
{
   // Just copy the data members.
   m_local_indicies = a_other.m_local_indicies;
   m_current = a_other.m_current;
   m_from_end = a_other.m_from_end;
   m_base = a_other.m_base;
}


void
KineticSpeciesPtrVect::FindSpeciesVelocityBoundingBox(
   double *vlo,
   double *vhi) const
{
   const int num_kinetic_species = static_cast<int>(size());
   if (num_kinetic_species < 1) {
      LOKI_ABORT("No species in FindSpeciesVelocityBoundingBox.");
   }

   for (int s=0; s<num_kinetic_species; ++s) {
      const vector<double>& lower = (*this)[s]->GetGlobalDomainLower();
      const vector<double>& upper = (*this)[s]->GetGlobalDomainUpper();

      if (s == 0) {
         vlo[0] = lower[V1];
         vlo[1] = lower[V2];
         vhi[0] = upper[V1];
         vhi[1] = upper[V2];
      }
      else {
         vlo[0] = min(vlo[0], lower[V1]);
         vlo[1] = min(vlo[1], lower[V2]);
         vhi[0] = max(vhi[0], upper[V1]);
         vhi[1] = max(vhi[1], upper[V2]);
      }
   }
}

} // end namespace Loki
