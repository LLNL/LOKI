/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _ELECTRIC_POTENTIAL_DRIVER_FACTORY_H_
#define _ELECTRIC_POTENTIAL_DRIVER_FACTORY_H_

namespace Loki {

class LokiInputParser;
class ElectricPotentialDriver;
/**
 * Factory for creating electric potential drivers.  Only contains static
 * create method.
 */
class ElectricPotentialDriverFactory
{
public:
   /**
    * @brief Creates the requested ElectricPotentialDriver.
    *
    * @param[in] a_pp Where the driver's input parameters will be found.
    *
    * @return The requested ElectricPotentialDriver.
    */
   static
   ElectricPotentialDriver*
   create(
      LokiInputParser& a_pp);

private:
   // Unimplemented default constructor.
   ElectricPotentialDriverFactory();

   // Unimplemented destructor.
   ~ElectricPotentialDriverFactory();

   // Unimplemented copy constructor.
   ElectricPotentialDriverFactory(
      const ElectricPotentialDriverFactory& other);

   // Unimplemented assignment operator.
   ElectricPotentialDriverFactory&
   operator = (
      const ElectricPotentialDriverFactory& rhs);
};

} // end namespace Loki

#endif
