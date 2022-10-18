/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _EM_IC_FACTORY_H_
#define _EM_IC_FACTORY_H_

namespace Loki {

class LokiInputParser;
class EMICInterface;

/**
 * Factory for creating electromagnetic initial condition objects.  Only
 * contains static create method.
 */
class EMICFactory
{
public:
   /**
    * @brief Creates the requested EMICInterface.
    *
    * @param[in] a_pp Where the electromagnetic initial condition's input
    *                 parameters will be found.
    *
    * @return The requested EMICInterface.
    */
   static
   EMICInterface*
   create(
      LokiInputParser& a_pp);

private:
   // Unimplemented default constructor.
   EMICFactory();

   // Unimplemented destructor.
   ~EMICFactory();

   // Unimplemented copy constructor.
   EMICFactory(const EMICFactory& other);

   // Unimplemented assignment operator.
   EMICFactory&
   operator = (
      const EMICFactory& rhs);
};

} // end namespace Loki

#endif
