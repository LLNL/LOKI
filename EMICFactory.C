/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "EMICFactory.H"
#include "SimpleEMIC.H"
#include "Loki_Defines.H"

// Add new derived EMICInterface headers here

namespace Loki {

EMICInterface*
EMICFactory::create(
   LokiInputParser& a_pp)
{
  // Get the name of the EM initial condition.  This is required.
   string name("none");
   a_pp.query("name", name);
   EMICInterface* em_ic_interface = 0;

   // If the user forgot the name, error.  Otherwise try to build what they
   // asked for.
   if (name.compare("none") == 0) {
      LOKI_ABORT("Missing electromagnetic initial condition ... quitting");
   }
   else {
      if (SimpleEMIC::isType(name)) {
         em_ic_interface = new SimpleEMIC(a_pp);
      }
      // Add new cases here in this form:
      // else if (NewEMIC::isType(name)) {
      //    em_ic_interface = new NewEMIC(a_pp);
      // }
      else {
         ostringstream error_msg;
         error_msg << "Unknown electromagnetic initial condition \"" << name
                   << "\" ... quitting";
      }
   }

   return em_ic_interface;
}

} // end namespace Loki
