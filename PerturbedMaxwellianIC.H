/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _PERTURBED_MAXWELLIAN_IC_H_
#define _PERTURBED_MAXWELLIAN_IC_H_

#include "ICInterface.H"
#include "tbox/Pointer.H"
#include "LokiInputParser.H"
#include <vector>

namespace Loki {

/**
 * Implements interface of ICInterface for a perturbed maxwellian initial
 * condition.  This class actually covers 3 such initial conditions:
 * 1) "built in"
 * 2) Landau damping
 * 3) Maxwellian with noise
 */
class PerturbedMaxwellianIC : public ICInterface
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_pp Where the input parameters are to be found.
    * @param[in] a_domain Information about the problem domain.
    * @param[in] a_mass Mass of the species being initialized.
    */
   PerturbedMaxwellianIC(
      LokiInputParser& a_pp,
      const ProblemDomain& a_domain,
      double a_mass);

   /**
    * @brief Destructor.
    */
   virtual
   ~PerturbedMaxwellianIC();

   //@{

   //! @name Implementation of ICInterface interface.

   /**
    * @brief Retrieve cached initial condition and set a_u to it.
    *
    * @param[out] a_u The KineticSpecies distribution function to set.
    */
   virtual
   void
   set(
      ParallelArray& a_u) const;

   /**
    * @brief Compute and cache the initial condition for the points in the
    *        supplied array.
    *
    * @param[in] a_u Domain over which initial condition is cached.
    */
   virtual
   void
   cache(
      const ParallelArray& a_u);

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
      int a_i4) const;

   /**
    * @brief Print the driver's parameters.
    */
   virtual
   void
   printParameters() const;

   //@}

   /**
    * @brief Returns true if name describes this class' type.
    *
    * @param[in] a_name Name of a type of initial condition.
    *
    * @return true if name describes this class' type.
    */
   static
   bool
   isType(
      const string& a_name);

private:
   // Unimplemented default constructor.
   PerturbedMaxwellianIC();

   // Unimplemented copy constructor.
   PerturbedMaxwellianIC(
      const PerturbedMaxwellianIC& other);

   // Unimplemented assignment operator.
   PerturbedMaxwellianIC&
   operator = (
      const PerturbedMaxwellianIC& rhs);

   // Parse the input parameters.
   void
   parseParameters(
      LokiInputParser& a_pp);

   // Type of initial condition:
   // 1 = built in
   // 2 = Landau damping
   // 3 = Maxwellian with noise
   int m_ic_option;

   // If true the thermal behavior is maxwellian.  Otherwise it is Juttner.
   bool m_maxwellian_thermal;

   // The implementation of the requested thermal behavior.
   tbox::Pointer<ThermalBehavior> m_thermal;

   // The parameters necessary for the evaluation of this initial condition.
   double m_tx;
   double m_ty;
   double m_vx0;
   double m_vy0;
   double m_x_wave_number;
   double m_y_wave_number;
   double m_flow_vel_phi;
   double m_a;
   double m_kx1;
   double m_ky1;
   double m_b;
   double m_kx2;
   double m_c;
   double m_ky2;
   double m_spatial_phi;

   // Domain length in x.
   double m_Lx;

   // Noise amplitudes for Maxwellian initial condition
   vector<double> m_noise_amp;

   // Noise phases for Maxwellian initial condition
   vector<double> m_noise_phase;
};

} // end namespace Loki

#endif
