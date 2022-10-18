/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _ELECTRIC_FIELD_DRIVER_FACTORY_H_
#define _ELECTRIC_FIELD_DRIVER_FACTORY_H_

namespace Loki {

class LokiInputParser;
class ElectricFieldDriver;

/**
 * Factory for creating electric field drivers.  Only contains static create
 * method.
 */
class ElectricFieldDriverFactory
{
public:
   /**
    * @brief Creates the requested ElectricFieldDriver.
    *
    * @param[in] a_pp Where the driver's input parameters will be found.
    * @param[in] a_driver_num One based index of this driver.
    * @param[in] a_species_num One based index of the KineticSpecies this
    *                          driver is associated with.
    *
    * @return The requested ElectricFieldDriver.
    */
   static
   ElectricFieldDriver*
   create(
      LokiInputParser& a_pp,
      int a_driver_num,
      int a_species_num);

private:
   // Unimplemented default constructor.
   ElectricFieldDriverFactory();

   // Unimplemented destructor.
   ~ElectricFieldDriverFactory();

   // Unimplemented copy constructor.
   ElectricFieldDriverFactory(
      const ElectricFieldDriverFactory& other);

   // Unimplemented assignment operator.
   ElectricFieldDriverFactory&
   operator = (
      const ElectricFieldDriverFactory& rhs);
};

} // end namespace Loki

#endif
