/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _INTERPENETRATING_STREAM_IC_H_
#define _INTERPENETRATING_STREAM_IC_H_

#include "ICInterface.H"
#include "tbox/Pointer.H"
#include "ParallelArray.H"
#include "LokiInputParser.H"

namespace Loki {

/**
 * Implements interface of ICInterface for the initial condition of an
 * interpenetrating stream.
 */
class InterpenetratingStreamIC : public ICInterface
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_pp Where the input parameters are to be found.
    * @param[in] a_domain Information about the problem domain.
    * @param[in] a_mass The species mass.
    * @param[in] a_n_ghosts The number of ghosts.
    */
   InterpenetratingStreamIC(
      LokiInputParser& a_pp,
      const ProblemDomain& a_domain,
      double a_mass,
      int a_n_ghosts);

   /**
    * @brief Destructor.
    */
   virtual
   ~InterpenetratingStreamIC();

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
    * @brief Print the initial condition's parameters.
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
   InterpenetratingStreamIC();

   // Unimplemented copy constructor.
   InterpenetratingStreamIC(
      const InterpenetratingStreamIC& other);

   // Unimplemented assignment operator.
   InterpenetratingStreamIC&
   operator = (
      const InterpenetratingStreamIC& rhs);

   // Print the initial condition's box syntax parameters.
   void
   printParametersBox() const;

   // Print the initial condition's half plane syntax parameters.
   void
   printParametersHalfPlane() const;

   // Parse the input parameters.
   void
   parseParameters(
      LokiInputParser& a_pp);

   // Parse the input parameters for the box syntax.
   void
   parseParametersBox(
      LokiInputParser& a_pp);

   // Parse the input parameters for the half plane syntax.
   void
   parseParametersHalfPlane(
      LokiInputParser& a_pp);

   // If true the thermal behavior is maxwellian.  Otherwise it is Juttner.
   bool m_maxwellian_thermal;

   // The implementation of the requested thermal behavior.
   tbox::Pointer<ThermalBehavior> m_thermal;

   // The second configuration space term for half-plane syntax, and two sided
   // or centered.
   ParallelArray m_fx2;

   // The second velocity space term for half-plane syntax and two sided.
   ParallelArray m_fv2;

   // The second thermal behavior for half-plane syntax and two sided.
   tbox::Pointer<ThermalBehavior> m_thermal2;

   // Flag to indicate if the original rotated half-plane oriented syntax or the
   // box oriented syntax is used.
   bool m_box_syntax;

   // The parameters necessary for the evaluation of this initial condition.
   // Parameters common to both box and half-plane sytaxes.
   double m_beta;

   // Parameters that only pertain to the half-plane syntax.
   double m_tl;
   double m_tt;
   double m_vl0;
   double m_vt0;
   double m_theta;
   double m_d;
   double m_frac2;
   double m_floor;
   bool m_two_sided;
   bool m_centered;

   // Parameters that only pertain to the box syntax.
   double m_tx;
   double m_ty;
   int m_num_boxes;
   vector<double> m_box_x_hi;
   vector<double> m_box_x_lo;
   vector<double> m_box_y_hi;
   vector<double> m_box_y_lo;
   vector<double> m_box_frac;
};

} // end namespace Loki

#endif
