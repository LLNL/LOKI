/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _THERMAL_BEHAVIOR_H_
#define _THERMAL_BEHAVIOR_H_

#include "tbox/DescribedClass.H"

namespace Loki {

/**
 * Base class for different plasma initial condition thermal behaviors.
 */
class ThermalBehavior : public tbox::DescribedClass {
public:
   /** @brief Destructor.
    */
   virtual
   ~ThermalBehavior();

   /** @brief Computes the contribution to the distribution function due to the
    *         thermal behavior of the plasma.
    *
    * @param[in] a_spatial_factor Factor for spatial perturbation of
    *                             distribution function.
    * @param[in] a_x3 vx/px coordinate
    * @param[in] a_x4 vy/py coordinate
    */
   virtual
   double
   thermalFactor(
      double a_spatial_factor,
      double a_x3,
      double a_x4) const = 0;

   double
   fnorm()
   {
      return m_fnorm;
   }

protected:
   // Constructor.
   ThermalBehavior(
      double a_vflowinitx,
      double a_vflowinity,
      double a_vx0,
      double a_vy0,
      double a_mass,
      double a_tx,
      double a_ty);

   // For a non-relativistic problem this is the initial x flow velocity.
   // For a relativistic problem this is the initial x flow momentum.
   double m_flowinitx;

   // For a non-relativistic problem this is the initial y flow velocity.
   // For a relativistic problem this is the initial y flow momentum.
   double m_flowinity;

   // For a non-relativistic problem this is the spatial perturbation of the
   // distribution function in vx.
   // For a relativistic problem this is the spatial perturbation of the
   // distribution function in px.
   double m_x0;

   // For a non-relativistic problem this is the spatial perturbation of the
   // distribution function in vy.
   // For a relativistic problem this is the spatial perturbation of the
   // distribution function in py.
   double m_y0;

   // For a non-relativistic problem this is the square of the x thermal
   // velocity.
   // For a relativistric problem this is the x thermal momentum.
   double m_thermalx;

   // For a non-relativistic problem this is the square of the y thermal
   // velocity.
   // For a relativistric problem this is the y thermal momentum.
   double m_thermaly;

   // Normalization factor for thermal behavior.
   double m_fnorm;

private:
   // Unimplemented default constructor.
   ThermalBehavior();

   // Unimplemented copy constructor.
   ThermalBehavior(
      const ThermalBehavior& a_other);

   // Unimplemented assignment operator.
   ThermalBehavior&
   operator = (
      const ThermalBehavior& a_rhs);
};

} // end namespace Loki

#endif
