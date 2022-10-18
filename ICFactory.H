/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _IC_FACTORY_H_
#define _IC_FACTORY_H_

#include "LokiInputParser.H"

namespace Loki {

class KineticSpecies;
class ICInterface;

/**
 * Factory for creating initial condition objects.  Only contains static create
 * method.
 */
class ICFactory
{
public:
   /**
    * @brief Creates the requested ICInterface.
    *
    * @param[in] a_pp Where the initial condition's input parameters will be
    *                 found.
    * @param[in] a_species The species the initial condition works on.
    *
    * @return The requested ICInterface.
    */
   static
   ICInterface*
   create(
      LokiInputParser& a_pp,
      KineticSpecies* a_species);

private:
   // Unimplemented default constructor.
   ICFactory();

   // Unimplemented destructor.
   ~ICFactory();

   // Unimplemented copy constructor.
   ICFactory(const ICFactory& other);

   // Unimplemented assignment operator.
   ICFactory&
   operator = (
      const ICFactory& rhs);
};

} // end namespace Loki

#endif
