/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _EXTERNAL_2D_IC_H_
#define _EXTERNAL_2D_IC_H_

#include "ICInterface.H"
#include "tbox/Pointer.H"
#include "LokiInputParser.H"

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
    * @param[in] a_domain Information about the problem domain.
    * @param[in] a_mass Mass of the species being initialized.
    * @param[in] a_num_ghosts The number of ghosts in each direction.
    * @param[in] a_global_box The problem domain + boundary.
    */
   External2DIC(
      LokiInputParser& a_pp,
      const ProblemDomain& a_domain,
      double a_mass,
      int a_num_ghosts,
      const tbox::Box& a_global_box);

   /**
    * @brief Destructor.
    */
   virtual
   ~External2DIC();

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
      LokiInputParser& a_pp);

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

   // The name of the 2D external distribution file.
   string m_ext_2d_file;

   // The data read in from the external distribution file.
   vector<double> m_data;

   // Number of data points in x.
   int m_nx;

   // Number of ghosts in each direction.
   int m_num_ghosts;
};

} // end namespace Loki

#endif
