/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "ElectricPotentialDriverFactory.H"
#include "LokiInputParser.H"
#include "Loki_Defines.H"

// Add new derived ElectricPotentialDriver headers here
#include "ShapedRampedCosinePotentialDriver.H"

namespace Loki {

ElectricPotentialDriver*
ElectricPotentialDriverFactory::create(
   LokiInputParser& a_pp)
{
   // Get the name of the electric potential driver.  This is required.
   string name("none");
   a_pp.query("name", name);
   ElectricPotentialDriver* driver = 0;

   // If the user forgot the name, error.  Otherwise try to build what they
   // asked for.
   if (name.compare("none") == 0) {
      LOKI_ABORT("Missing external electric potential driver ... quitting");
   }
   else {
      if (ShapedRampedCosinePotentialDriver::isType(name)) {
         driver = new ShapedRampedCosinePotentialDriver(a_pp);
      }
      // Add new cases here in this form:
      // else if (NewDriver::isType(name)) {
      //    driver = new NewDriver(a_pp);
      // }
      else {
         ostringstream error_msg;
         error_msg << "Unknown external electric potential \"" << name
                   << "\" ... quitting";
         LOKI_ABORT(error_msg);
      }
   }
   return driver;
}

} // end namespace Loki
