/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _VP_STATE_H_
#define _VP_STATE_H_

#include "KineticSpeciesPtrVect.H"
#include "Poisson.H"
#include "tbox/Pointer.H"

namespace Loki {

/**
 * Glorified struct holding the state variables (KineticSpecies and Poisson)
 * of a VPSystem.
 */
class VPState
{
public:
   /**
    * @brief Default constructor.
    */
   VPState();

   /**
    * @brief Copy constructor.
    *
    * @param[in] other The thing being copied.
    */
   VPState(
      const VPState& other);

   /**
    * @brief Destructor.
    */
   ~VPState();

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
    * @brief Accessor to Poisson.
    */
   tbox::Pointer<Poisson>&
   poisson()
   {
      return m_poisson;
   }

   /**
    * @brief Accessor to Poisson.
    */
   const tbox::Pointer<Poisson>&
   poisson() const
   {
      return m_poisson;
   }

private:
   // Unimplemented assignment operator.
   VPState&
   operator = (
      const VPState& rhs);

   // Holder for KineticSpecies data.
   KineticSpeciesPtrVect m_kspv;

   // Holder for Poisson data.
   tbox::Pointer<Poisson> m_poisson;
};

} // end namespace Loki

#endif
