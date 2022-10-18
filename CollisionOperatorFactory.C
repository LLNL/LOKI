/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "CollisionOperatorFactory.H"
#include "LokiInputParser.H"
#include "KineticSpecies.H"

// Add new derived CollisionOperator headers here
#include "PitchAngleCollisionOperator.H"
#include "RosenbluthCollisionOperator.H"
#include "OriginalRosenbluthCollisionOperator.H"

namespace Loki {

CollisionOperator*
CollisionOperatorFactory::create(
   LokiInputParser& a_pp,
   const KineticSpecies* a_species)
{
   // Get the name of the collision operator.  This is required.
   string name("none");
   a_pp.query("name", name);
   CollisionOperator* coll_op = 0;

   // If the user forgot the name, error.  Otherwise try to build what they
   // asked for.
   if (name.compare("none") == 0) {
      LOKI_ABORT("Missing collision operator ... quitting");
   }
   else {
      if (PitchAngleCollisionOperator::isType(name)) {
         coll_op = new PitchAngleCollisionOperator(a_pp, a_species);
      }
      else if (RosenbluthCollisionOperator::isType(name)) {
         coll_op = new RosenbluthCollisionOperator(a_pp, a_species);
      }
      else if (OriginalRosenbluthCollisionOperator::isType(name)) {
         coll_op = new OriginalRosenbluthCollisionOperator(a_pp, a_species);
      }
      // Add new cases here in this form:
      // else if (NewOperator::isType(name)) {
      //    coll_op = new NewOperator(a_pp);
      // }
      else {
         ostringstream error_msg;
         error_msg << "Unknown collision operator \"" << name
                   << "\" ... quitting";
         LOKI_ABORT(error_msg);
      }
   }
   return coll_op;
}

} // end namespace Loki
