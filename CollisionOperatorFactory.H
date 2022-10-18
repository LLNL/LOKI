/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _COLLISION_OPERATOR_FACTORY_H_
#define _COLLISION_OPERATOR_FACTORY_H_

namespace Loki {

class LokiInputParser;
class CollisionOperator;
class KineticSpecies;

/**
 * Factory for creating collision operators.  Only contains static create
 * method.
 */
class CollisionOperatorFactory
{
public:
   /**
    * @brief Creates the requested CollisionOperator.
    *
    * @param[in] a_pp Where the operator's input parameters will be found.
    * @param[in] a_species The species this operator is associated with.
    *
    * @return The requested CollisionOperator.
    */
   static
   CollisionOperator*
   create(
      LokiInputParser& a_pp,
      const KineticSpecies* a_species);

private:
   // Unimplemented default constructor.
   CollisionOperatorFactory();

   // Unimplemented destructor.
   ~CollisionOperatorFactory();

   // Unimplemented copy constructor.
   CollisionOperatorFactory(
      const CollisionOperatorFactory& other);

   // Unimplemented assignment operator.
   CollisionOperatorFactory&
   operator = (
      const CollisionOperatorFactory& rhs);
};

} // end namespace Loki

#endif
