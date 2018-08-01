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
#ifndef _EXTERNAL_2D_IC_H_
#define _EXTERNAL_2D_IC_H_

#include "ICInterface.H"
#include "Overture.h"
#include "ParmParse.H"

namespace Loki {

/**
 * Implements interface of ICInterface by adding the velocity modulation to
 * an externally supplied 2D distribution.
 */
class External2DIC : public ICInterface
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_pp Where the input parameters are to be found.
    * @param[in] a_vflowinitx The initial flow velocity in x.
    * @param[in] a_vflowinity The initial flow velocity in y.
    */
   External2DIC(
      ParmParse& a_pp,
      real a_vflowinitx,
      real a_vflowinity);

   /**
    * @brief Destructor.
    */
   virtual
   ~External2DIC();

   //@{

   //! @name Implementation of ICInterface interface.

   /**
    * @brief Compute the initial conditions of u at the supplied time.
    *
    * @param[out] a_u The thing to compute.
    * @param[in] a_domain Information about the problem domain.
    * @param[in] a_grown_global_box The global box for the entire problem with
    *                               an additional set of ghosts because that's
    *                               how Overture does things.
    * @param[in] a_time The time at which u is to be computed.
    */
   virtual
   void
   set(
      RealArray& a_u,
      const ProblemDomain& a_domain,
      const tbox::Box& a_grown_global_box,
      real a_time) const;

   /**
    * @brief Print the driver's parameters.
    */
   virtual
   void
   printParameters() const;

   // FIXME: When setting boundary conditions the KineticSpecies passes the
   // initial condition parameters to its fortran routines which "know"
   // something about how initial conditions are set.  This seems very wrong.
   /**
    * @brief Returns object's parameters.
    */
   virtual
   const RealArray&
   parameters() const
   {
      return m_parameters;
   }

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
      const aString& a_name);

private:
   // Unimplemented default constructor.
   External2DIC();

   // Unimplemented copy constructor.
   External2DIC(
      const External2DIC& other);

   // Unimplemented assignment operator.
   External2DIC&
   operator = (
      const External2DIC& rhs);

   // Parse the input parameters.
   void
   parseParameters(
      ParmParse& a_pp);

   // Indices of specific parameters in m_parameters.
   enum {
      ALPHA,
      BETA,
      VX0,
      VY0,
      VFLOWINITX,
      VFLOWINITY,
      X_WAVE_NUMBER,
      Y_WAVE_NUMBER,
      FLOW_VEL_PHI,
      NUM_PARAMS};

   // The parameters necessary for the evaluation of this initial condition.
   RealArray m_parameters;

   // The name of the 2D external distribution.
   aString m_ext_2d_file;
};

} // end namespace Loki

#endif
