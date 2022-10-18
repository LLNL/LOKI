/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _TWOELECTRON_TRIG_TZ_SOURCE_H_
#define _TWOELECTRON_TRIG_TZ_SOURCE_H_

#include "TZSourceInterface.H"
#include "LokiInputParser.H"
#include <vector>

namespace Loki {

/**
 * Implements the interface of TZSourceInterface for a trigonometric twilight
 * zone source having the following functional form:
 *
 * @verbatim
 *  f = f + amp*cos(x)*cos(y)*cos(t)
 * @endverbatim
 */
class TwoSpecies_ElectronTrigTZSource : public TZSourceInterface
{
public:
   /**
    * @brief Default Constructor.
    *
    * @param[in] a_pp Where the input parameters are to be found.
    */
   TwoSpecies_ElectronTrigTZSource(
      LokiInputParser& a_pp);

   /**
    * @brief Destructor.
    */
   virtual
   ~TwoSpecies_ElectronTrigTZSource();

   //@{

   //! @name Implementation of TZSourceInterface interface.

   /**
    * @brief Compute the twilight zone source for this class and add it to the
    *        distribution function.
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
      const ParallelArray& a_velocities) const;

   /**
    * @brief Compute the error in the twilight zone source solution.
    *
    * @param[out] a_tz_error Errir in the twilight zone source solution.
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
      const ParallelArray& a_velocities) const;

   /**
    * @brief Print the parameters defining the initial condition.
    */
   virtual
   void
   printParameters() const;

   //@}

   /**
    * @brief Returns true if name describes this class' type.
    *
    * @param[in] a_name Name of a type of twilight zone source.
    *
    * @return true if name describes this class' type.
    */
   static
   bool
   isType(
      const string& a_name);

private:
   // Unimplemented copy constructor.
   TwoSpecies_ElectronTrigTZSource(
      const TwoSpecies_ElectronTrigTZSource& other);

   // Unimplemented assignment operator.
   TwoSpecies_ElectronTrigTZSource&
   operator = (
      const TwoSpecies_ElectronTrigTZSource& rhs);

   // The name of this type of twilight zone source.
   static const string s_CLASS_NAME;

   // Indices of specific parameters in m_dparameters.
   enum {
      AMP,
      EMASS,
      IMASS,
      NUM_DPARAMS};

   // The double valued parameters.
   vector<double> m_dparameters;
};

} // end namespace Loki
#endif
