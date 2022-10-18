/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _SIMPLE_EM_IC_H_
#define _SIMPLE_EM_IC_H_

#include "EMICInterface.H"
#include "LokiInputParser.H"
#include <vector>

namespace Loki {

/**
 * SimpleEMICInterface implements the EMICInterface interface for a simple wave
 * propagating in a single direction and having the following functional form:
 *
 * @verbatim
 *  E = (Ex, Ey, Ez)cos(k*propagation_dir + phi)
 * @endverbatim
 */
class SimpleEMIC : public EMICInterface
{
public:
   /**
    * @brief Default Constructor.
    *
    * @param[in] a_pp Where the input parameters are to be found.
    */
   SimpleEMIC(
      LokiInputParser& a_pp);

   /**
    * @brief Destructor.
    */
   virtual
   ~SimpleEMIC();

   //@{

   //! @name Implementation of EMICInterface interface.

   /**
    * @brief Compute the initial conditions of the component of the
    *        electromagnetic field, a_u, that this SimpleEMIC is concerned
    *        with.
    *
    * @param[out] a_u All components of the electromagnetic filed.
    * @param[in] a_domain Information about the problem domain.
    */
   virtual
   void
   set(
      ParallelArray& a_u,
      const ProblemDomain& a_domain) const;

   /**
    * @brief Print the parameters defining the initial condition.
    */
   virtual
   void
   printParameters() const;

   /**
    * @brief Returns true if this EM initializer handles the E field.
    */
   virtual
   bool
   initializesE();

   /**
    * @brief Returns true if this EM initializer handles the B field.
    */
   virtual
   bool
   initializesB();

   /**
    * @brief Returns the number of EM waves being initialized.
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
   SimpleEMIC(
      const SimpleEMIC& other);

   // Unimplemented assignment operator.
   SimpleEMIC&
   operator = (
      const SimpleEMIC& rhs);

   // The name of this type of electromagnetic initial condition.
   static const string s_CLASS_NAME;

   // The 2 fields.
   enum EM_Var{
      E = 1,
      B = 4};

   // The 2 propagation directions.
   enum Prop_Dir{
      X = 1,
      Y = 2};

   // Indices of specific parameters in m_iparameters.
   enum {
      FIELD,
      NUM_IPARAMS};

   // The number of wave initializers.
   int m_num_waves;

   // The x amplitude of each wave.
   vector<double> m_x_amp;

   // The y amplitude of each wave.
   vector<double> m_y_amp;

   // The z amplitude of each wave.
   vector<double> m_z_amp;

   // The x wave number of each wave.
   vector<double> m_x_wave_number;

   // The y wave number of each wave.
   vector<double> m_y_wave_number;

   // The phase of each wave.
   vector<double> m_phi;

   // The integer valued parameters.
   vector<int> m_iparameters;
};

} // end namespace Loki
#endif
