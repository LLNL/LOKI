/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _INTERPOLATOR6_H_
#define _INTERPOLATOR6_H_

#include "Interpolator.H"

namespace Loki {

/**
 * Concrete class implementing the interface of Interpolator to evaluate a 2D
 * field at a point using a 6th order interpolation.
 */
class Interpolator6 : public Interpolator
{
public:
   /**
    * @brief Constructor
    *
    * @param[in] a_domain The ProblemDomain the field is associated with.
    */
   Interpolator6(
      const ProblemDomain& a_domain)
      : Interpolator(a_domain)
   {
   }

   /**
    * @brief Destructor
    */
   virtual
   ~Interpolator6()
   {
   }

   /**
    * @brief Returns interpolatation of a_interp_coefs at the evaluation point
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
      double** a_interp_coefs)
   {
      double result = 0.0;
      for (int xidx = 0; xidx < 6; ++xidx) {
         for (int yidx = 0; yidx < 6; ++yidx) {
            result += a_interp_coefs[xidx][yidx]*m_f[xidx][yidx];
         }
      }
      return result;
   }

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
      double** a_interp_coefs)
   {
      for (int xidx = 0; xidx < 6; ++xidx) {
         for (int yidx = 0; yidx < 6; ++yidx) {
            a_interp_coefs[xidx][yidx] =
               a_evaluated_value*m_f[xidx][yidx]/m_f_scale;
         }
      }
   }

   /**
    * @brief Return the requested x interpolating grid index.
    *
    * @param[in] a_i The requested grid index.
    */
   virtual
   int
   xInterpolatingGridIndex(
      int a_i)
   {
      return m_x_grid_indices[a_i];
   }

   /**
    * @brief Return the requested y interpolating grid index.
    *
    * @param[in] a_i The requested grid index.
    */
   virtual
   int
   yInterpolatingGridIndex(
      int a_i)
   {
      return m_y_grid_indices[a_i];
   }

private:
   // Unimplemented default constructor.
   Interpolator6();

   // Unimplemented copy constructor.
   Interpolator6(
      const Interpolator6& a_other);

   // Unimplemented assignment operator.
   Interpolator6&
   operator = (
      const Interpolator6& a_rhs);

   // Computes the indices of the interpolating grid for an evaluation point at
   // the supplied "fractional" grid indices.
   virtual
   void
   computeInterpolatingGridIndices(
      double a_x_fractional_idx,
      double a_y_fractional_idx)
   {
      if (a_x_fractional_idx >= 0.0) {
         m_x_grid_indices[3] = (int)a_x_fractional_idx;
      }
      else {
         m_x_grid_indices[3] = (int)a_x_fractional_idx -1;
      }
      m_x_grid_indices[0] = m_x_grid_indices[3] - 3;
      m_x_grid_indices[1] = m_x_grid_indices[3] - 2;
      m_x_grid_indices[2] = m_x_grid_indices[3] - 1;
      m_x_grid_indices[4] = m_x_grid_indices[3] + 1;
      m_x_grid_indices[5] = m_x_grid_indices[3] + 2;
      if (a_y_fractional_idx >= 0.0) {
         m_y_grid_indices[3] = (int)a_y_fractional_idx;
      }
      else {
         m_y_grid_indices[3] = (int)a_y_fractional_idx -1;
      }
      m_y_grid_indices[0] = m_y_grid_indices[3] - 3;
      m_y_grid_indices[1] = m_y_grid_indices[3] - 2;
      m_y_grid_indices[2] = m_y_grid_indices[3] - 1;
      m_y_grid_indices[4] = m_y_grid_indices[3] + 1;
      m_y_grid_indices[5] = m_y_grid_indices[3] + 2;
   }

   // Evaluate the basis functions at the supplied "fractional" grid indices.
   virtual
   void
   evalBasisFuncs(
      double a_x_fractional_idx,
      double a_y_fractional_idx)
   {
      // Compute a parametric x and y position relative to the interpolating
      // points.  Then evaluate the basis functions at that parametric
      // coordinate.
      double x = a_x_fractional_idx - m_x_grid_indices[2];
      double x2 = x*x;
      double x3 = x*x2;
      double x4 = x*x3;
      double x5 = x*x4;
      double y = a_y_fractional_idx - m_y_grid_indices[2];
      double y2 = y*y;
      double y3 = y*y2;
      double y4 = y*y3;
      double y5 = y*y4;
      double fx[6], fy[6];
      fx[0] = -x5/120.0 +           x3/24.0                  - x/30.0;
      fx[1] =  x5/24.0  + x4/24.0 - x3*7.0/24.0 - x2/24.0    + x/4.0;
      fx[2] = -x5/12.0  - x4/6.0  + x3*7.0/12.0 + x2*2.0/3.0 - x;
      fx[3] =  x5/12.0  + x4/4.0  - x3*5.0/12.0 - x2*5.0/4.0 + x/3.0  + 1.0;
      fx[4] = -x5/24.0  - x4/6.0  + x3/24.0     + x2*2.0/3.0 + x/2.0;
      fx[5] =  x5/120.0 + x4/24.0 + x3/24.0     - x2/24.0    - x/20.0;
      fy[0] = -y5/120.0 +           y3/24.0                  - y/30.0;
      fy[1] =  y5/24.0  + y4/24.0 - y3*7.0/24.0 - y2/24.0    + y/4.0;
      fy[2] = -y5/12.0  - y4/6.0  + y3*7.0/12.0 + y2*2.0/3.0 - y;
      fy[3] =  y5/12.0  + y4/4.0  - y3*5.0/12.0 - y2*5.0/4.0 + y/3.0  + 1.0;
      fy[4] = -y5/24.0  - y4/6.0  + y3/24.0     + y2*2.0/3.0 + y/2.0;
      fy[5] =  y5/120.0 + y4/24.0 + y3/24.0     - y2/24.0    - y/20.0;
      m_f_scale = 0.0;
      for (int xidx = 0; xidx < 6; ++xidx) {
         for (int yidx = 0; yidx < 6; ++yidx) {
            double this_f = fx[xidx]*fy[yidx];
            m_f_scale += this_f*this_f;
            m_f[xidx][yidx] = this_f;
         }
      }
   }

   // The x indices of the interpolating grid.
   int m_x_grid_indices[6];

   // The y indices of the interpolating grid.
   int m_y_grid_indices[6];

   // Basis functions evaluated at the evaluation point.
   double m_f[6][6];

   // Squared norm of basis functions evaluated at the evaluation point.
   double m_f_scale;
};

} // end namespace Loki

#endif
