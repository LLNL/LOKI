/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _INTERPOLATOR_H_
#define _INTERPOLATOR_H_

#include "ParallelArray.H"
#include "ProblemDomain.H"

namespace Loki {

/**
 * Abstract base class defining the interface to evaluate a 2D field at a point.
 * The evaluation is via a set of basis functions of a given order and the
 * corresponding basis function coefficients.  An Nth order interpolator has
 * N*N basis functions defined on an NxN interpolatin grid.
 */
class Interpolator
{
public:
   /**
    * @brief Destructor.
    */
   virtual
   ~Interpolator()
   {
   }

   /**
    * @brief Returns interpolation of a_interp_coefs at the evaluation point
    *        set in the last call to setEvaluationPoint.
    *
    * @param[in] a_interp_coefs Coefficients for the basis functions.
    *
    * @return The interpolation of a_interp_coefs at the last evaluation point
    *         set.
    */
   virtual
   double
   interpolate(
      double** a_interp_coefs) = 0;

   /**
    * @brief Computes the basis function coefficients that interpolate to the
    *        supplied value at the evaluation point set in the last call to
    *        setEvaluationPoint.
    *
    * @param[in] a_evaluated_value The value at the evaluation point.
    * @param[out] a_interp_coefs Coefficients for the basis functions.
    */
   virtual
   void
   deterpolate(
      double a_evaluated_value,
      double** a_interp_coefs) = 0;

   /**
    * @brief Set up the interpolator for a specific evaluation point.
    *
    * Once the interpolator has been set for a specific evaluation point by a
    * call to this function, interpolate or deterpolate may be called repeatedly
    * to interpolate or deterpolate different fields at that evaluation point.
    *
    * @param[in] a_eval_x x coordinate of the evaluated result.
    * @param[in] a_eval_y y coordinate of the evaluated result.
    */
   void
   setEvaluationPoint(
      double a_eval_x,
      double a_eval_y)
   {
      // Find the "fractional" grid index of the evaluation point.
      double x_fractional_idx = (a_eval_x-m_domain_x_lo) / m_domain_dx;
      double y_fractional_idx = (a_eval_y-m_domain_y_lo) / m_domain_dy;

      // Now compute the interpolating grid indices.
      computeInterpolatingGridIndices(x_fractional_idx, y_fractional_idx);

      // Evaluate the basis functions at the evaluation point.
      evalBasisFuncs(x_fractional_idx, y_fractional_idx);
   }

   /**
    * @brief Return the requested x interpolating grid index.
    *
    * @param[in] a_i The requested grid index.
    */
   virtual
   int
   xInterpolatingGridIndex(
      int a_i) = 0;

   /**
    * @brief Return the requested y interpolating grid index.
    *
    * @param[in] a_i The requested grid index.
    */
   virtual
   int
   yInterpolatingGridIndex(
      int a_i) = 0;

protected:
   // Constructor.
   Interpolator(
      const ProblemDomain& a_domain)
   {
      m_domain_x_lo = a_domain.lower(0);
      m_domain_y_lo = a_domain.lower(1);
      m_domain_dx = a_domain.dx(0);
      m_domain_dy = a_domain.dx(1);
   }

   // Computes the indices of the interpolating grid for an evaluation point at
   // the supplied "fractional" grid indices.
   virtual
   void
   computeInterpolatingGridIndices(
      double a_x_fractional_idx,
      double a_y_fractional_idx) = 0;

   // Evaluate the basis functions at the supplied "fractional" grid indices.
   virtual
   void
   evalBasisFuncs(
      double a_x_fractional_idx,
      double a_y_fractional_idx) = 0;

protected:
   // The x position of the domain's left edge.  Does not include ghosts.
   double m_domain_x_lo;

   // The y position of the domain's bottom edge.  Does not include ghosts.
   double m_domain_y_lo;

   // The domain's grid spacing in x.
   double m_domain_dx;

   // The domain's grid spacing in y.
   double m_domain_dy;

private:
   // Unimplemented default constructor.
   Interpolator();
};

} // end namespace Loki

#endif
