/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _UNIDIRECTIONAL_CURRENT_DRIVER_H_
#define _UNIDIRECTIONAL_CURRENT_DRIVER_H_

#include "CurrentDriver.H"
#include "LokiInputParser.H"

namespace Loki {

/**
 * Implementation of the CurrentDriver interface for a current driver
 * propagating a wave in a single direction.
 */
class UnidirectionalCurrentDriver : public CurrentDriver
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_pp Where the input parameters are to be found.
    */
   UnidirectionalCurrentDriver(
      LokiInputParser& a_pp);

   /**
    * @brief Destructor.
    */
   virtual
   ~UnidirectionalCurrentDriver();

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
   UnidirectionalCurrentDriver();

   // Unimplemented copy constructor.
   UnidirectionalCurrentDriver(
      const UnidirectionalCurrentDriver& other);

   // Unimplemented assignment operator.
   UnidirectionalCurrentDriver&
   operator = (
      const UnidirectionalCurrentDriver& rhs);

   // Parse the input parameters.
   void
   parseParameters(
      LokiInputParser& a_pp);

   // The name of this type of current driver.
   static const string s_CLASS_NAME;

   // Indices of specific parameters in m_parameters.
   enum {
      OMEGA,
      J0,
      T0,
      TRAMPUP,
      THOLD,
      TRAMPDOWN,
      POLAR_ANGLE,
      X0X,
      X0Y,
      X1X,
      X1Y,
      BETA,
      LIGHT_SPEED,
      SHAPE,
      WIDTH,
      J0_INPLANE,
      J0_OUTPLANE,
      L,
      PHI,
      NUM_PARAMS
   };

   // Matrix of constant geometric factors,
   vector<double> m_R;

   // The parameters necessary for the evaluation of this driver.
   vector<double> m_parameters;
};

} // end namespace Loki

#endif
