/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _TZ_SOURCE_INTERFACE_H_
#define _TZ_SOURCE_INTERFACE_H_

#include "ParallelArray.H"
#include "ProblemDomain.H"

namespace Loki {

/**
 * TZSourceInterface is an abstract base class for twilight zone Vlasov RHS
 * source implementations.
 */
class TZSourceInterface
{
public:
   /**
    * @brief Destructor.
    */
   virtual
   ~TZSourceInterface();

   /**
    * @brief Compute the twilight zone source and add it to a_dist_func.
    *
    * @param[out] a_dist_func Distribution function.
    * @param[in] a_domain Information about the problem domain.
    * @param[in] a_time Current simulation time.
    * @param[in] a_velocities Velocities, relativistic or otherwise.
    */
   virtual
   void
   set(
      ParallelArray& a_dist_func,
      const ProblemDomain& a_domain,
      double a_time,
      const ParallelArray& a_velocities) const = 0;

   /**
    * @brief Compute the error in the twilight zone source solution.
    *
    * @param[out] a_tz_error Error in the twilight zone source solution.
    * @param[in] a_dist_func Distribution function.
    * @param[in] a_domain Information about the problem domain.
    * @param[in] a_time Current simulation time.
    * @param[in] a_velocities Velocities, relativistic or otherwise.
    */
   virtual
   void
   computeError(
      ParallelArray& a_tz_error,
      const ParallelArray& a_dist_func,
      const ProblemDomain& a_domain,
      double a_time,
      const ParallelArray& a_velocities) const = 0;

   /**
    * @brief Print the parameters defining the twilight zone source.
    */
   virtual
   void
   printParameters() const = 0;

protected:
   /**
    * @brief Default Constructor.
    */
   TZSourceInterface();

private:
   // Unimplemented copy constructor.
   TZSourceInterface(
      const TZSourceInterface& other);

   // Unimplemented assignment operator.
   TZSourceInterface&
   operator = (
      const TZSourceInterface& rhs);
};

} // end namespace Loki

#endif
