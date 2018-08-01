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
#ifndef _INTERPOLATOR_H_
#define _INTERPOLATOR_H_

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
      double x_fractional_idx =
         (a_eval_x-m_domain_x_lo) / m_domain_dx + 2*m_num_x_ghosts;
      double y_fractional_idx =
         (a_eval_y-m_domain_y_lo) / m_domain_dy + 2*m_num_y_ghosts;

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

   /**
    * @brief Accessor for number of x grid points in the local box.
    *
    * Useful info for client to index into field being interpolated.
    */
   int
   nx()
   {
      return m_nx;
   }

   /**
    * @brief Accessor for number of y grid points in the local box.
    *
    * Useful info for client to index into field being interpolated.
    */
   int
   ny()
   {
      return m_ny;
   }

   /**
    * @brief Accessor for total number of grid points in the local box.
    *
    * Useful info for client to index into field being interpolated.
    */
   int
   nxny()
   {
      return m_nxny;
   }

protected:
   // Constructor.
   Interpolator(
      const tbox::Box& a_local_box,
      const ProblemDomain& a_domain,
      int a_num_x_ghosts,
      int a_num_y_ghosts)
   {
      m_domain_x_lo = a_domain.lower(0);
      m_domain_y_lo = a_domain.lower(1);
      m_domain_dx = a_domain.dx(0);
      m_domain_dy = a_domain.dx(1);
      m_nx = a_local_box.upper(0) - a_local_box.lower(0) + 1;
      m_ny = a_local_box.upper(1) - a_local_box.lower(1) + 1;
      m_nxny = m_nx*m_ny;
      m_num_x_ghosts = a_num_x_ghosts;
      m_num_y_ghosts = a_num_y_ghosts;
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

   // The number of ghosts in x.
   int m_num_x_ghosts;

   // The number of ghosts in y.
   int m_num_y_ghosts;

   // The number of x grid points in the local box.
   int m_nx;

   // The number of y grid points in the local box.
   int m_ny;

   // The total number of grid points in the local box.
   int m_nxny;

private:
   // Unimplemented default constructor.
   Interpolator();
};

} // end namespace Loki

#endif
