/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _ELECTRIC_POTENTIAL_DRIVER_H_
#define _ELECTRIC_POTENTIAL_DRIVER_H_

#include "ParallelArray.H"
#include "ProblemDomain.H"  

namespace Loki {

/**
 * ElectricPotentialDriver is an abstract base class for all
 * ElectricPotentialDriver implementations.
 */
class ElectricPotentialDriver
{
public:
   /**
    * @brief Default Constructor.
    */
   ElectricPotentialDriver();

   /**
    * @brief Destructor.
    */
   virtual
   ~ElectricPotentialDriver();

   /**
    * @brief Compute the value of the electric potential at the supplied time.
    *
    * @param[out] a_phi The electric potential.
    * @param[in] a_domain Information about the problem domain.
    * @param[in] a_time The time at which the potential is to be computed.
    */
   virtual
   void
   evaluate(
      ParallelArray& a_phi,
      const ProblemDomain& a_domain,
      double a_time) const = 0;

   /**
    * @brief Print the driver's parameters.
    */
   virtual
   void
   printParameters() const = 0;

private:
   // Unimplemented copy constructor.
   ElectricPotentialDriver(
      const ElectricPotentialDriver& other);

   // Unimplemented assignment operator.
   ElectricPotentialDriver&
   operator = (
      const ElectricPotentialDriver& rhs);
};

} // end namespace Loki

#endif
