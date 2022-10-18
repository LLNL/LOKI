/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _CURRENT_DRIVER_FACTORY_H_
#define _CURRENT_DRIVER_FACTORY_H_

namespace Loki {

class LokiInputParser;
class ProblemDomain;
class CurrentDriver;

/**
 * Factory for creating current drivers.  Only contains static create method.
 */
class CurrentDriverFactory
{
public:
   /**
    * @brief Creates the requested CurrentDriver.
    *
    * @param[in] a_pp Where the driver's input parameters will be found.
    * @param[in] a_domain Information about the problem domain.
    * @param[in] a_driver_num One based index of this driver.
    *
    * @return The requested CurrentDriver.
    */
   static
   CurrentDriver*
   create(
      LokiInputParser& a_pp,
      const ProblemDomain& a_domain,
      int a_driver_num);

private:
   // Unimplemented default constructor.
   CurrentDriverFactory();

   // Unimplemented destructor.
   ~CurrentDriverFactory();

   // Unimplemented copy constructor.
   CurrentDriverFactory(
      const CurrentDriverFactory& other);

   // Unimplemented assignment operator.
   CurrentDriverFactory&
   operator = (
      const CurrentDriverFactory& rhs);
};

} // end namespace Loki

#endif
