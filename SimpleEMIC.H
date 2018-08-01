/*************************************************************************
 *
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * Written by Jeffrey Banks banksj3@rpi.edu (Rensselaer Polytechnic Institute,
 * Amos Eaton 301, 110 8th St., Troy, NY 12180); Jeffrey Hittinger
 * hittinger1@llnl.gov, William Arrighi arrighi2@llnl.gov, Richard Berger
 * berger5@llnl.gov, Thomas Chapman chapman29@llnl.gov (LLNL, P.O Box 808,
 * Livermore, CA 94551); Stephan Brunner stephan.brunner@epfl.ch (Ecole
 * Polytechnique Federale de Lausanne, EPFL SB SPC-TH, PPB 312, Station 13,
 * CH-1015 Lausanne, Switzerland).
 * CODE-744849
 *
 * All rights reserved.
 *
 * This file is part of Loki.  For details, see.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 ************************************************************************/
#ifndef _SIMPLE_EM_IC_H_
#define _SIMPLE_EM_IC_H_

#include "EMICInterface.H"
#include "ParmParse.H"

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
      ParmParse& a_pp);

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
      RealArray& a_u,
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
      const aString& a_name);

private:
   // Unimplemented copy constructor.
   SimpleEMIC(
      const SimpleEMIC& other);

   // Unimplemented assignment operator.
   SimpleEMIC&
   operator = (
      const SimpleEMIC& rhs);

   // The name of this type of electromagnetic initial condition.
   static const aString s_CLASS_NAME;

   // The 2 fields.
   enum EM_Var{
      E = 1,
      B = 4};

   // The 2 propagation directions.
   enum Prop_Dir{
      X = 1,
      Y = 2};

   // Indices of specific parameters in m_dparameters.
   enum {
      XAMP,
      YAMP,
      ZAMP,
      X_WAVE_NUMBER,
      Y_WAVE_NUMBER,
      PHI,
      NUM_DPARAMS};

   // Indices of specific parameters in m_iparameters.
   enum {
      FIELD,
      NUM_IPARAMS};

   // The double valued parameters.
   RealArray m_dparameters;

   // The integer valued parameters.
   IntegerArray m_iparameters;
};

} // end namespace Loki
#endif
