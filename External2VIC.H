/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _EXTERNAL_2V_IC_H_
#define _EXTERNAL_2V_IC_H_

#include "ICInterface.H"
#include <map>

namespace Loki {

/**
 * Implements interface of ICInterface by applying the externally supplied
 * velocity space behavior everywhere in configuration space.
 */
class External2VIC : public ICInterface
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_file File to read initial condition from.
    * @param[in] a_domain Information about the problem domain.
    * @param[in] a_num_ghosts The number of ghosts in each direction.
    */
   External2VIC(
      const string& a_file,
      const ProblemDomain& a_domain,
      int a_num_ghosts);

   /**
    * @brief Destructor.
    */
   virtual
   ~External2VIC();

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
   External2VIC();

   // Unimplemented copy constructor.
   External2VIC(
      const External2VIC& other);

   // Unimplemented assignment operator.
   External2VIC&
   operator = (
      const External2VIC& rhs);

   // The velocity space data.
   vector<double> m_vspace_data;

   // Number of data points in vx.
   int m_nvx;

   // Number of ghosts in each direction.
   int m_num_ghosts;
};

} // end namespace Loki

#endif
