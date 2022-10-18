/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "ICFactory.H"
#include "KineticSpecies.H"

// Add new derived ICInterface headers here
#include "PerturbedMaxwellianIC.H"
#include "External2DIC.H"
#include "InterpenetratingStreamIC.H"

namespace Loki {

ICInterface*
ICFactory::create(
   LokiInputParser& a_pp,
   KineticSpecies* a_species)
{
   // Get the name of initial condition.  This is required.
   string ic_name("Perturbed Maxwellian");
   a_pp.query("name", ic_name);
   ICInterface* ic_interface = 0;

   // If the user forgot the name, error.  Otherwise try to build what they
   // asked for.
   if (PerturbedMaxwellianIC::isType(ic_name)) {
      ic_interface = new PerturbedMaxwellianIC(a_pp,
         a_species->domain(),
         a_species->mass());
   }
   else if (External2DIC::isType(ic_name)) {
      ic_interface = new External2DIC(a_pp,
         a_species->domain(),
         a_species->mass(),
         a_species->numGhosts(),
         a_species->globalBox());
   }
   else if (InterpenetratingStreamIC::isType(ic_name)) {
      ic_interface = new InterpenetratingStreamIC(a_pp,
         a_species->domain(),
         a_species->mass(),
         a_species->numGhosts());
   }
   // Add new cases here in this form:
   // else if (NewIC::isType(ic_name)) {
   //    ic_interface = new NewIC(a_pp);
   // }
   else {
      ostringstream error_msg;
      error_msg << "Unknown initial condition \"" << ic_name
                << "\" ... quitting";
      LOKI_ABORT(error_msg);
   }

   return ic_interface;
}

} // end namespace Loki
