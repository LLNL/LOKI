/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _EM_IC_INTERFACE_H_
#define _EM_IC_INTERFACE_H_

#include "ParallelArray.H"
#include "ProblemDomain.H"

namespace Loki {

/**
 * EMICInterface is an abstract base class for all electromagnetic field
 * initial condition implementations needed for a Vlasov-Maxwell system.
 */
class EMICInterface
{
public:
   /**
    * @brief Default Constructor.
    */
   EMICInterface();

   /**
    * @brief Destructor.
    */
   virtual
   ~EMICInterface();

   /**
    * @brief Compute the initial conditions of the electromagnetic field, a_u,
    *        that this EMICInterface is concerned with.
    *
    * @param[out] a_u All components of the electromagnetic field.
    * @param[in] a_domain Information about the problem domain.
    */
   virtual
   void
   set(
      ParallelArray& a_u,
      const ProblemDomain& a_domain) const = 0;

   /**
    * @brief Print the parameters defining the initial condition.
    */
   virtual
   void
   printParameters() const = 0;

   /**
    * @brief Returns true if this EM initializer handles the E field.
    */
   virtual
   bool
   initializesE() = 0;

   /**
    * @brief Returns true if this EM initializer handles the B field.
    */
   virtual
   bool
   initializesB() = 0;

   /**
    * @brief Returns the number of EM waves being initialized.
    */
   virtual
   int
   numWaves() = 0;

private:
   // Unimplemented copy constructor.
   EMICInterface(
      const EMICInterface& other);

   // Unimplemented assignment operator.
   EMICInterface&
   operator = (
      const EMICInterface& rhs);
};

} // end namespace Loki

#endif
