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
#ifndef _VM_STATE_H_
#define _VM_STATE_H_

#include "KineticSpeciesPtrVect.H"
#include "Maxwell.H"
#include "tbox/Pointer.H"

namespace Loki {

/**
 * Glorified struct holding the state variables (KineticSpecies and Maxwell)
 * of a VMSystem.
 */
class VMState
{
public:
   /**
    * @brief Default constructor.
    */
   VMState();

   /**
    * @brief Copy constructor.
    *
    * @param[in] other The thing being copied.
    */
   VMState(
      const VMState& other);

   /**
    * @brief Destructor.
    */
   ~VMState();

   /**
    * @brief Accessor to Kinetic Species.
    */
   KineticSpeciesPtrVect&
   kineticSpecies()
   {
      return m_kspv;
   }

   /**
    * @brief Accessor to Kinetic Species.
    */
   const KineticSpeciesPtrVect&
   kineticSpecies() const
   {
      return m_kspv;
   }

   /**
    * @brief Accessor to Maxwell.
    */
   tbox::Pointer<Maxwell>&
   maxwell()
   {
      return m_maxwell;
   }

   /**
    * @brief Accessor to Maxwell.
    */
   const tbox::Pointer<Maxwell>&
   maxwell() const
   {
      return m_maxwell;
   }

private:
   // Unimplemented assignment operator.
   VMState&
   operator = (
      const VMState& rhs);

   // Holder for KineticSpecies data.
   KineticSpeciesPtrVect m_kspv;

   // Holder for Maxwell data.
   tbox::Pointer<Maxwell> m_maxwell;
};

} // end namespace Loki

#endif
