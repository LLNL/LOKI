/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _VEL_IC_FACTORY_H_
#define _VEL_IC_FACTORY_H_

namespace Loki {

class LokiInputParser;
class VELICInterface;

/**
 * Factory for creating transverse drift velocity initial condition objects.
 * Only contains static create method.
 */
class VELICFactory
{
public:
   /**
    * @brief Creates the requested VELICInterface.
    *
    * @param[in] a_pp Where the transverse drift velocity initial condition's
    *                 input parameters will be found.
    *
    * @return The requested VELICInterface.
    */
   static
   VELICInterface*
   create(
      LokiInputParser& a_pp);

private:
   // Unimplemented default constructor.
   VELICFactory();

   // Unimplemented destructor.
   ~VELICFactory();

   // Unimplemented copy constructor.
   VELICFactory(const VELICFactory& other);

   // Unimplemented assignment operator.
   VELICFactory&
   operator = (
      const VELICFactory& rhs);
};

} // end namespace Loki

#endif
