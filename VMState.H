/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
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
