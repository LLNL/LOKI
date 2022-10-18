/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _MAXWELLIAN_THERMAL_H_
#define _MAXWELLIAN_THERMAL_H_

#include "ThermalBehavior.H"

namespace Loki {

/**
 * Class to implement Maxwellian thermal behavior.  For a relativistic problem
 * the Maxwellian is WRT momentum space.  For a non-relativistic problem the
 * Maxwellian is WRT velocity space.
 */
class MaxwellianThermal : public ThermalBehavior
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_vflowinitx The x flow velocity.
    * @param[in] a_vflowinity The y flow velocity.
    * @param[in] a_vx0 The spatial perturbation of the distribution function in
    *                  vx.
    * @param[in] a_vy0 The spatial perturbation of the distribution function in
    *                  vy.
    * @param[in] a_mass The species mass.
    * @param[in] a_tx The temperature in x.
    * @param[in] a_ty The temperature in y.
    */
   MaxwellianThermal(
      double a_vflowinitx,
      double a_vflowinity,
      double a_vx0,
      double a_vy0,
      double a_mass,
      double a_tx,
      double a_ty);

   /**
    * @brief Destructor.
    */
   ~MaxwellianThermal();

   /**
    * @brief Computes the contribution to the distribution function due to the
    *        thermal behavior of the plasma.
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
      double a_x4) const;

private:
   // Unimplemented default constructor.
   MaxwellianThermal();

   // Unimplemented copy constructor.
   MaxwellianThermal(
      const MaxwellianThermal& a_other);

   // Unimplemented assignment operator;
   MaxwellianThermal&
   operator = (
      const MaxwellianThermal& a_rhs);
};

} // end namespace Loki

#endif
