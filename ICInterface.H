/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _IC_INTERFACE_H_
#define _IC_INTERFACE_H_

#include "ProblemDomain.H"
#include "ThermalBehavior.H"
#include "ParallelArray.H"

namespace Loki {

/**
 * ICInterface is an abstract base class for all initial condition
 * implementations.
 */
class ICInterface
{
public:
   /**
    * @brief Default Constructor.
    *
    * @param[in] a_domain Information about the problem domain.
    */
   ICInterface(
      const ProblemDomain& a_domain);

   /**
    * @brief Destructor.
    */
   virtual
   ~ICInterface();

   /**
    * @brief Retrieve cached initial condition and set a_u to it.
    *
    * @param[out] a_u The KineticSpecies distribution function to set.
    */
   virtual
   void
   set(
      ParallelArray& a_u) const = 0;

   /**
    * @brief Compute and cache the initial condition for the points in the
    *        supplied array.
    *
    * @param[in] a_u Domain over which initial condition is cached.
    */
   virtual
   void
   cache(
      const ParallelArray& a_u) = 0;

   /**
    * @brief Retrieve the cached initial condition at the requested point.
    *
    * @param[in] a_i1 x index
    * @param[in] a_i2 y index
    * @param[in] a_i3 vx index
    * @param[in] a_i4 vy index
    */
   virtual
   double
   getIC_At_Pt(
      int a_i1,
      int a_i2,
      int a_i3,
      int a_i4) const = 0;

   /**
    * @brief Print the driver's parameters.
    */
   virtual
   void
   printParameters() const = 0;

   /**
    * @brief Returns the initial x flow velocity.
    */
   double
   vflowinitx()
   {
      return m_vflowinitx;
   }

   /**
    * @brief Returns the initial y flow velocity.
    */
   double
   vflowinity()
   {
      return m_vflowinity;
   }

protected:
   // ProblemDomain for the species associated with this initial condition.
   const ProblemDomain& m_domain;

   // Initial x flow velocity.
   double m_vflowinitx;

   // Initial y flow velocity.
   double m_vflowinity;

   // Relative abundance of the species.
   double m_frac;

   // If true the cached initial condition may be factored into 2 2D arrays,
   // m_fx(x,y) and m_fv(vx,vy).  Otherwise it is not factorable and the entire
   // 4D array, m_f, must be cached.
   bool m_factorable;

   // When factorable, the configuration space term.
   ParallelArray m_fx;

   // When factorable, the velocity space term.
   ParallelArray m_fv;

   // When not factorable, the entire initial condition.
   ParallelArray m_f;

private:
   // Unimplemented copy constructor.
   ICInterface(
      const ICInterface& other);

   // Unimplemented assignment operator.
   ICInterface&
   operator = (
      const ICInterface& rhs);
};

} // end namespace Loki

#endif
