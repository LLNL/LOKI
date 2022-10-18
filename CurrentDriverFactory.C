/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "CurrentDriverFactory.H"
#include "LokiInputParser.H"
#include "Loki_Defines.H"

// Add new derived CurrentDriver headers here
#include "GaussianCurrentDriver.H"
#include "ShapedRampedCosineCurrentDriver.H"
#include "UnidirectionalCurrentDriver.H"

namespace Loki {

CurrentDriver*
CurrentDriverFactory::create(
   LokiInputParser& a_pp,
   const ProblemDomain& a_domain,
   int a_driver_num)
{
   string name("none");
   a_pp.query("name", name);
   CurrentDriver* driver = 0;

   if (name.compare("none") == 0) {
      ostringstream msg;
      msg << "Missing current driver number " << a_driver_num
          << " ... quitting";
      LOKI_ABORT(msg.str());
   }
   else {
      if (ShapedRampedCosineCurrentDriver::isType(name)) {
         driver = new ShapedRampedCosineCurrentDriver(a_pp);
      }
      else if (UnidirectionalCurrentDriver::isType(name)) {
         driver = new UnidirectionalCurrentDriver(a_pp);
      }
      else if (GaussianCurrentDriver::isType(name)) {
         driver = new GaussianCurrentDriver(a_pp, a_domain);
      }
      // Add new cases here in this form:
      // else if (NewDriver::isType(name)) {
      //    driver = new NewDriver(a_pp);
      // }
      else {
         ostringstream error_msg;
         error_msg << "Unknown external current \"" << name
                   << "\" ... quitting";
         LOKI_ABORT(error_msg);
      }
   }
   return driver;
}

} // end namespace Loki
