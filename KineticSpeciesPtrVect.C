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
   int my_id(Communication_Manager::My_Process_Number);
   my_id = std::max(0, my_id);
   for (int s(0); s < static_cast<int>(a_in.size()); ++s) {
      if (a_in[s]->isInRange(my_id)) {
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
   int my_id(Communication_Manager::My_Process_Number);
   my_id = std::max(0, my_id);
   for (int s(0); s < static_cast<int>(a_in.size()); ++s) {
      if (a_in[s]->isInRange(my_id)) {
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

} // end namespace Loki
