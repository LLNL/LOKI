/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "ElectricFieldDriverFactory.H"
#include "LokiInputParser.H"
#include "Loki_Defines.H"

// Add new derived ElectricFieldDriver headers here
#include "ShapedRampedCosineDriver.H"

namespace Loki {

ElectricFieldDriver*
ElectricFieldDriverFactory::create(
   LokiInputParser& a_pp,
   int a_driver_num,
   int a_species_num)
{
   // Get the name of the electric field driver.  This is required.
   string name("none");
   a_pp.query("name", name);
   ElectricFieldDriver* driver = 0;

   // If the user forgot the name, error.  Otherwise try to build what they
   // asked for.
   if (name.compare("none") == 0) {
      ostringstream msg;
      msg << "Missing external electric field driver number " << a_driver_num
          << " ... quitting";
      LOKI_ABORT(msg.str());
   }
   else {
      if (ShapedRampedCosineDriver::isType(name)) {
         driver =
            new ShapedRampedCosineDriver(a_pp, a_driver_num, a_species_num);
      }
      // Add new cases here in this form:
      // else if (NewDriver::isType(name)) {
      //    driver = new NewDriver(a_pp);
      // }
      else {
         ostringstream error_msg;
         error_msg << "Unknown external electric field \"" << name
                   << "\" ... quitting";
         LOKI_ABORT(error_msg);
      }
   }
   return driver;
}

} // end namespace Loki
