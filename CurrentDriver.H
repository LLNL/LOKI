/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _CURRENT_DRIVER_H_
#define _CURRENT_DRIVER_H_

#include "tbox/Box.H"
#include "ParallelArray.H"
#include "ProblemDomain.H"

namespace Loki {

/**
 * CurrentDriver is an abstract base class for all CurrentDriver
 * implementations.
 */
class CurrentDriver
{
public:
   /**
    * @brief Default Constructor.
    */
   CurrentDriver();

   /**
    * @brief Destructor.
    */
   virtual
   ~CurrentDriver();

   /**
    * @brief Compute and accumulate the antenna source at the supplied time.
    *
    * @param[out] a_antenna_source The antenna source accumulator.
    * @param[in] a_omega_eff2 Square of the effective plasma frequency.
    * @param[in] a_domain Information about the problem domain.
    * @param[in] a_time The time at which the current is to be computed.
    */
   virtual
   void
   evaluate(
      ParallelArray& a_antenna_source,
      const ParallelArray& a_omega_eff2,
      const ProblemDomain& a_domain,
      double a_time) const = 0;

   /**
    * @brief Returns true if the dispersion relation is needed.
    */
   virtual
   bool
   needsDispersionRelation() = 0;

   /**
    * @brief Print the driver's parameters.
    */
   virtual
   void
   printParameters() const = 0;

private:
   // Unimplemented copy constructor.
   CurrentDriver(
      const CurrentDriver& other);

   // Unimplemented assignment operator.
   CurrentDriver&
   operator = (
      const CurrentDriver& rhs);
};

} // end namespace Loki

#endif
