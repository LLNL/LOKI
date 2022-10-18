/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _SIMPLE_VEL_IC_H_
#define _SIMPLE_VEL_IC_H_

#include "VELICInterface.H"
#include "LokiInputParser.H"

namespace Loki {

/**
 * VELICInterface implements the VELICInterface interface for a simple wave
 * propagating in a single direction and having the following functional form:
 *
 * @verbatim
 *  vz = vz + zamp*cos(k*propagation_dir + phi)
 * @endverbatim
 */
class SimpleVELIC : public VELICInterface
{
public:
   /**
    * @brief Default Constructor.
    *
    * @param[in] a_pp Where the input parameters are to be found.
    */
   SimpleVELIC(
      LokiInputParser& a_pp);

   /**
    * @brief Destructor.
    */
   virtual
   ~SimpleVELIC();

   //@{

   //! @name Implementation of VELICInterface interface.

   /**
    * @brief Compute the initial conditions of the transverse drift velocity
    *        that this SimpleVELIC is concerned with.
    *
    * @param[out] a_vz Transverse drift velocity.
    * @param[in] a_domain Information about the problem domain.
    */
   virtual
   void
   set(
      ParallelArray& a_vz,
      const ProblemDomain& a_domain) const;

   /**
    * @brief Print the parameters defining the initial condition.
    */
   virtual
   void
   printParameters() const;

   /**
    * @brief Returns the number of EM waves initializing the transverse drift
    *        velocity.
    */
   virtual
   int
   numWaves();

   //@}

   /**
    * @brief Returns true if name describes this class' type.
    *
    * @param[in] a_name Name of a type of electromagnetic initial condition.
    *
    * @return true if name describes this class' type.
    */
   static
   bool
   isType(
      const string& a_name);

private:
   // Unimplemented copy constructor.
   SimpleVELIC(
      const SimpleVELIC& other);

   // Unimplemented assignment operator.
   SimpleVELIC&
   operator = (
      const SimpleVELIC& rhs);

   // The name of this type of electromagnetic initial condition.
   static const string s_CLASS_NAME;

   // The number of wave initializers.
   int m_num_waves;

   // The amplitude of each wave's contribution to the transverse drift
   // velocity.
   vector<double> m_amp;

   // The x wave number of each wave.
   vector<double> m_x_wave_number;

   // The y wave number of each wave.
   vector<double> m_y_wave_number;

   // The phase of each wave's contribution to the transverse drift velocity.
   vector<double> m_phi;
};

} // end namespace Loki
#endif
