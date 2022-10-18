/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "VELICFactory.H"
#include "Loki_Defines.H"
#include "SimpleVELIC.H"

// Add new derived VELICInterface headers here

namespace Loki {

VELICInterface*
VELICFactory::create(
   LokiInputParser& a_pp)
{
   // Get the name of the transverse drift velocity initializer.  This is
   // required.
   string name("none");
   a_pp.query("name", name);
   VELICInterface* vel_ic_interface = 0;

   // If the user forgot the name, error.  Otherwise try to build what they
   // asked for.
   if (name.compare("none") == 0) {
      LOKI_ABORT("Missing transverse drift velocity initial condition ... quitting");
   }
   else {
      if (SimpleVELIC::isType(name)) {
         vel_ic_interface = new SimpleVELIC(a_pp);
      }
      // Add new cases here in this form:
      // else if (NewVELIC::isType(name)) {
      //    vel_ic_interface = new NewVELIC(a_pp);
      // }
      else {
         ostringstream error_msg;
         error_msg << "Unknown transverse drift velocity initial condition \""
                   << name << "\" ... quitting";
         LOKI_ABORT(error_msg);
      }
   }

   return vel_ic_interface;
}

} // end namespace Loki
