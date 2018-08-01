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
#ifndef _KROOK_LAYER_H_
#define _KROOK_LAYER_H_

#include "Overture.h"
#include "OvertureTypes.h"
#include "ParmParse.H"
#include "HDF5_DataBase.h"

#include "ProblemDomain.H"
#include "Directions.H"
#include "tbox/Dimension.H"

namespace Loki {

/**
 * Class describing a Krook layer used to attenuate KineticSpecies distribution
 * functions in defined layers at the domain boundaries.
 */
class KrookLayer
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_dim The configuration space dimensions of the problem.
    * @param[in] a_pp Where the input parameter will be found
    * @param[in] a_domain The global problem domain.
    */
   KrookLayer(
      const tbox::Dimension& a_dim,
      ParmParse& a_pp,
      const ProblemDomain& a_domain);

   /**
    * @brief Destructor.
    */
   ~KrookLayer();

   /**
    * @brief Returns true if a Krook layer in any direction has been defined.
    */
   bool
   hasKrookLayer() const
   {
      return m_has_layer;
   }

   /**
    * @brief Returns true if a Krook layer on the specified side of the
    * specified direction has been defined.
    *
    * @param[in] a_dir Direction of interest.
    * @param[in] a_side Side of interest.
    */
   bool
   hasKrookLayer(
      int a_dir,
      int a_side) const
   {
      if (a_side==LO) {
         return m_has_layer_lo[a_dir];
      }
      return m_has_layer_hi[a_dir];
   }

   /**
    * @brief Returns 1 if a Krook layer in any direction has been defined, 0
    * otherwise.
    */
   int
   hasKrookLayerAsInteger() const
   {
      return ((m_has_layer) ? 1 : 0);
   }

   /**
    * @brief Print this class' parameters.
    */
   void
   printParameters() const;

   /**
    * @brief Write information about Krook layers to restart database.
    *
    * @param[out] a_db Restart database.
    */
   void
   putToDatabase(
      HDF_DataBase& a_db) const;

   /**
    * @brief Read information about Krook layers to restart database.
    *
    * @param[in] a_db Restart database.
    */
   void
   getFromDatabase(
      const HDF_DataBase& a_db);

   /**
    * @brief Returns the coordinates of the lower end of each Krook layer.
    */
   const Array<double>&
   lower() const
   {
      return m_x_lo;
   }

   /**
    * @brief Returns the coordinates of the upper end of each Krook layer.
    */
   const Array<double>&
   upper() const
   {
      return m_x_hi;
   }

   /**
    * @brief Returns the coordinates of the lower end of the Krook layer in the
    * requested direction.
    *
    * @param[in] a_dir Direction of interest.
    */
   double
   lower(
      int a_dir) const
   {
      return m_x_lo[a_dir];
   }

   /**
    * @brief Returns the coordinates of the upper end of the Krook layer in the
    * requested direction.
    *
    * @param[in] a_dir Direction of interest.
    */
   double
   upper(
      int a_dir) const
   {
      return m_x_hi[a_dir];
   }

   /**
    * @brief Returns the exponent in the evaluation of each Krook layer.
    */
   double
   power() const
   {
      return m_power;
   }

   /**
    * @brief Returns the overall scaling in the evaluation of each Krook layer.
    */
   double
   coefficient() const
   {
      return m_coefficient;
   }

private:
   // Unimplemented default constructor.
   KrookLayer();

   // Unimplemented copy constructor.
   KrookLayer(
      const KrookLayer& other);

   // Unimplemented assignment operator.
   KrookLayer&
   operator = (
      const KrookLayer& rhs);

   // Parse this class' parameters.
   void
   parseParameters(
      ParmParse& a_pp);

   // Check the validity of this class' parameters.
   void
   checkParameters(
      const ProblemDomain& a_domain);

   // Configuration space dimensions of the problem.
   tbox::Dimension m_dim;

   // Coordinates of the lower end of each Krook layer.
   Array<double> m_x_lo;

   // Coordinates of the upper end of each Krook layer.
   Array<double> m_x_hi;

   // For each dimension, is true if there is a Krook layer on the lower side
   // of that dimension.
   Array<bool> m_has_layer_lo;

   // For each dimension, is true if there is a Krook layer on the upper side
   // of that dimension.
   Array<bool> m_has_layer_hi;

   // Is true if a Krook layer in any direction has been defined.
   bool m_has_layer;

   // The exponent in the evaluation of each Krook layer.
   real m_power;

   // The overall scaling in the evaluation of each Krook layer.
   real m_coefficient;
};

} // end namespace Loki

#endif
