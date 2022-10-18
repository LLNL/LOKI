/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "TZSourceFactory.H"
#include "Loki_Defines.H"
#include "TrigTZSource.H"
#include "ElectronTrigTZSource.H"
#include "TwoSpecies_ElectronTrigTZSource.H"
#include "TwoSpecies_IonTrigTZSource.H"

// Add new derived TZSourceInterface headers here

namespace Loki {

TZSourceInterface*
TZSourceFactory::create(
   LokiInputParser& a_pp)
{
   string name("none");
   a_pp.query("name", name);
   TZSourceInterface* tz_source_interface = 0;

   if (name.compare("none") != 0) {
      if (TrigTZSource::isType(name)) {
         tz_source_interface = new TrigTZSource(a_pp);
      }
      // Add new cases here in this form:
      // else if (NewTZSource::isType(name)) {
      //    tz_source_interface = new NewTZSource(a_pp);
      // }

       else if (ElectronTrigTZSource::isType(name)) {
          tz_source_interface = new ElectronTrigTZSource(a_pp);
       }
       else if (TwoSpecies_ElectronTrigTZSource::isType(name)) {
          tz_source_interface = new TwoSpecies_ElectronTrigTZSource(a_pp);
       }
       else if (TwoSpecies_IonTrigTZSource::isType(name)) {
          tz_source_interface = new TwoSpecies_IonTrigTZSource(a_pp);
       }
      else {
         ostringstream error_msg;
         error_msg << "Unknown twilight zone source \"" << name
                   << "\" ... quitting";
         LOKI_ABORT(error_msg);
      }
   }

   return tz_source_interface;
}

} // end namespace Loki
