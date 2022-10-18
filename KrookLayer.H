/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _KROOK_LAYER_H_
#define _KROOK_LAYER_H_

#include "LokiInputParser.H"
#include "RestartReader.H"
#include "RestartWriter.H"
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
      LokiInputParser& a_pp,
      const ProblemDomain& a_domain);

   /**
    * @brief Destructor.
    */
   ~KrookLayer();

   /**
    * @brief Compute the damping coefficient nu for each point of the
    *        distribution function which overlaps the layer.
    *
    * @param[in] a_dist_func The partitioned distribution function.
    * @param[in] a_spatial_order The problem's spatial order.
    */
   void
   initialize(
      const ParallelArray& a_dist_func,
      int a_spatial_order);

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
    *        specified direction has been defined.
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
    *        otherwise.
    */
   int
   hasKrookLayerAsInteger() const
   {
      return ((m_has_layer) ? 1 : 0);
   }

   /**
    * @brief Returns true if the Krook layer overlaps this processor's extent.
    */
   bool
   overlaps() const
   {
      return m_overlaps;
   }

   /**
    * @brief Print this class' parameters.
    */
   void
   printParameters() const;

   /**
    * @brief Write information about Krook layers to restart database.
    *
    * @param[out] a_writer Restart file writer.
    * @param[in] a_write_data If true this processor actually writes its data.
    */
   void
   putToDatabase(
      RestartWriter& a_writer,
      bool a_write_data) const;

   /**
    * @brief Read information about Krook layers to restart database.
    *
    * @param[in] a_reader Restart file reader.
    */
   void
   getFromDatabase(
      RestartReader& a_reader);

   /**
    * @brief Returns the coordinates of the lower end of each Krook layer.
    */
   const vector<double>&
   lower() const
   {
      return m_x_lo;
   }

   /**
    * @brief Returns the coordinates of the upper end of each Krook layer.
    */
   const vector<double>&
   upper() const
   {
      return m_x_hi;
   }

   /**
    * @brief Returns the coordinates of the lower end of the Krook layer in the
    *        requested direction.
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
    *        requested direction.
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
    * @brief Returns the damping coefficient.
    */
   const ParallelArray&
   nu() const
   {
      return m_nu;
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
      LokiInputParser& a_pp);

   // Check the validity of this class' parameters.
   void
   checkParameters();

   // Configuration space dimensions of the problem.
   tbox::Dimension m_dim;

   // Coordinates of the lower end of each Krook layer.
   vector<double> m_x_lo;

   // Coordinates of the upper end of each Krook layer.
   vector<double> m_x_hi;

   // For each dimension, is true if there is a Krook layer on the lower side
   // of that dimension.
   deque<bool> m_has_layer_lo;

   // For each dimension, is true if there is a Krook layer on the upper side
   // of that dimension.
   deque<bool> m_has_layer_hi;

   // Is true if a Krook layer in any direction has been defined.
   bool m_has_layer;

   // If true the Krook layer overlaps this processor's extent.
   bool m_overlaps;

   // The exponent in the evaluation of each Krook layer.
   double m_power;

   // The overall scaling in the evaluation of each Krook layer.
   double m_coefficient;

   // The problem domain.
   ProblemDomain m_domain;

   // The damping coefficient.
   ParallelArray m_nu;
};

} // end namespace Loki

#endif
