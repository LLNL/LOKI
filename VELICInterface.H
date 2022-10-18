/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _VEL_IC_INTERFACE_H_
#define _VEL_IC_INTERFACE_H_

#include "ParallelArray.H"
#include "ProblemDomain.H"

namespace Loki {

/**
 * VELICInterface is an abstract base class for transverse drift velocity
 * initial condition implementations needed for a Vlasov-Maxwell system.
 */
class VELICInterface
{
public:
   /**
    * @brief Default Constructor.
    */
   VELICInterface();

   /**
    * @brief Destructor.
    */
   virtual
   ~VELICInterface();

   /**
    * @brief Compute the initial conditions of the transverse drift velocity
    *        that this VELICInterface is concerned with.
    *
    * @param[out] a_vz Transverse drift velocity.
    * @param[in] a_domain Information about the problem domain.
    */
   virtual
   void
   set(
      ParallelArray& a_vz,
      const ProblemDomain& a_domain) const = 0;

   /**
    * @brief Print the parameters defining the initial condition.
    */
   virtual
   void
   printParameters() const = 0;

   /**
    * @brief Returns the number of EM waves initializing the transverse drift
    *        velocity.
    */
   virtual
   int
   numWaves() = 0;

private:
   // Unimplemented copy constructor.
   VELICInterface(
      const VELICInterface& other);

   // Unimplemented assignment operator.
   VELICInterface&
   operator = (
      const VELICInterface& rhs);
};

} // end namespace Loki

#endif
