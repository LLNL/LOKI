/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _GAUSSIAN_CURRENT_DRIVER_H_
#define _GAUSSIAN_CURRENT_DRIVER_H_

#include "CurrentDriver.H"
#include "LokiInputParser.H"

namespace Loki {

/**
 * Implementation of the CurrentDriver interface for a current driver to
 * generate a Gaussian beam.
 */
class GaussianCurrentDriver : public CurrentDriver
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_pp Where the input parameters are to be found.
    * @param[in] a_domain Information about the problem domain.
    */
   GaussianCurrentDriver(
      LokiInputParser& a_pp,
      const ProblemDomain& a_domain);

   /**
    * @brief Destructor.
    */
   virtual
   ~GaussianCurrentDriver();

   //@{

   //! @name Implementation of CurrentDriver interface.

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
      double a_time) const;

   /**
    * @brief Print the driver's parameters.
    */
   virtual
   void
   printParameters() const;

   /**
    * @brief Returns true if the dispersion relation is needed.
    */
   virtual
   bool
   needsDispersionRelation();

   //@}

   /**
    * @brief Returns true if name describes this class' type.
    *
    * @param[in] a_name Name of a type of current driver.
    *
    * @return true if name describes this class' type.
    */
   static
   bool
   isType(
      const string& a_name);

private:
   // Unimplemented default constructor.
   GaussianCurrentDriver();

   // Unimplemented copy constructor.
   GaussianCurrentDriver(
      const GaussianCurrentDriver& other);

   // Unimplemented assignment operator.
   GaussianCurrentDriver&
   operator = (
      const GaussianCurrentDriver& rhs);

   // Parse the input parameters.
   void
   parseParameters(
      LokiInputParser& a_pp);

   // The name of this type of current driver.
   static const string s_CLASS_NAME;

   // Indices of specific parameters in m_parameters.
   enum {
      PLANE,
      X0,
      FOCAL_PT,
      FWHM,
      LAMBDA,
      OMEGA,
      NENC,
      J0,
      T0,
      TRAMPUP,
      THOLD,
      TRAMPDOWN,
      NUM_PARAMS
   };

   // The parameters necessary for the evaluation of this driver.
   vector<double> m_parameters;

   // The in-plane current shaping envelope perpendicular to the current
   // direction.
   vector<double> m_shaping;

   // The in-plane phase of the driver.
   vector<double> m_phase;

   // Low index of m_shaping.  Essentially the domain y low for a driver in
   // the y-z plane or the domain x low for a driver in the x-z plane.
   int m_shape_lo;

   // High index of m_shaping.  Essentially the domain y high for a driver in
   // the y-z plane or the domain x high for a driver in the x-z plane.
   int m_shape_hi;

   // The x or y index at which the current source is applied.
   int m_x_apply;
};

} // end namespace Loki

#endif
