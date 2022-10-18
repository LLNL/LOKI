/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _EXTERNAL_DIST_KROOK_LAYER_H_
#define _EXTERNAL_DIST_KROOK_LAYER_H_

#include "LokiInputParser.H"
#include "RestartReader.H"
#include "RestartWriter.H"
#include "ProblemDomain.H"
#include "Directions.H"
#include "External2VIC.H"
#include "tbox/Dimension.H"

namespace Loki {

/**
 * Class describing something similar to a Krook layer used to force the
 * distribution function in a spatial layer toward a user supplied distribution
 * dependent only on velocity.
 */
class ExternalDistKrookLayer
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_dim The configuration space dimensions of the problem.
    * @param[in] a_pp Where the input parameter will be found
    * @param[in] a_domain The global problem domain.
    */
   ExternalDistKrookLayer(
      const tbox::Dimension& a_dim,
      LokiInputParser& a_pp,
      const ProblemDomain& a_domain);

   /**
    * @brief Destructor.
    */
   ~ExternalDistKrookLayer();

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
    * @brief Returns true if an ExternalDistKrookLayer in any direction has
    *        been defined.
    */
   bool
   hasKrookLayer() const
   {
      return m_has_layer;
   }

   /**
    * @brief Returns true if an ExternalDistKrookLayer on the specified side
    *        of the specified direction has been defined.
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
    * @brief Returns 1 if an ExternalDistKrookLayer in any direction has been
    *        defined, 0 otherwise.
    */
   int
   hasKrookLayerAsInteger() const
   {
      return ((m_has_layer) ? 1 : 0);
   }

   /**
    * brief Returns true if the ExternalDistKrookLayer overlaps this processor's
    *       extent.
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
    * @brief Write information about ExternalDistKrookLayer to restart
    *        database.
    *
    * @param[out] a_writer Restart file writer.
    * @param[in] a_write_data If true this processor actually writes its data.
    */
   void
   putToDatabase(
      RestartWriter& a_writer,
      bool a_write_data) const;

   /**
    * @brief Read information about ExternalDistKrookLayer from restart
    *        database.
    *
    * @param[in] a_reader Restart file reader.
    */
   void
   getFromDatabase(
      RestartReader& a_reader);

   /**
    * @brief Returns the coordinates of the lower end of the lower layers in
    *        each directrion.
    */
   const vector<double>&
   lolo() const
   {
      return m_x_lo_lo;
   }

   /**
    * @brief Returns the coordinates of the upper end of the lower layers in
    *        each direction.
    */
   const vector<double>&
   hilo() const
   {
      return m_x_hi_lo;
   }

   /**
    * @brief Returns the coordinates of the lower end of the upper layers in
    *        each direction.
    */
   const vector<double>&
   lohi() const
   {
      return m_x_lo_hi;
   }

   /**
    * @brief Returns the coordinates of the upper end of the upper layers in
    *        each direction.
    */
   const vector<double>&
   hihi() const
   {
      return m_x_hi_hi;
   }

   /**
    * @brief Returns the coordinates of the lower end of the lower layer in the
    *        requested direction.
    *
    * @param[in] a_dir Direction of interest.
    */
   double
   lolo(
      int a_dir) const
   {
      return m_x_lo_lo[a_dir];
   }

   /**
    * @brief Returns the coordinates of the upper end of the lower layer in the
    *        requested direction.
    *
    * @param[in] a_dir Direction of interest.
    */
   double
   hilo(
      int a_dir) const
   {
      return m_x_hi_lo[a_dir];
   }

   /**
    * @brief Returns the coordinates of the lower end of the upper layer in the
    *        requested direction.
    *
    * @param[in] a_dir Direction of interest.
    */
   double
   lohi(
      int a_dir) const
   {
      return m_x_lo_hi[a_dir];
   }

   /**
    * @brief Returns the coordinates of the upper end of the upper layer in the
    *        requested direction.
    *
    * @param[in] a_dir Direction of interest.
    */
   double
   hihi(
      int a_dir) const
   {
      return m_x_hi_hi[a_dir];
   }

   /**
    * @brief Returns the damping coefficient.
    */
   const ParallelArray&
   nu() const
   {
      return m_nu;
   }

   /**
    * @brief Returns the external distribution initial condition.
    */
   External2VIC*
   externalDistIC() const
   {
      return m_external_distribution_ic;
   }

private:
   // Unimplemented default constructor.
   ExternalDistKrookLayer();

   // Unimplemented copy constructor.
   ExternalDistKrookLayer(
      const ExternalDistKrookLayer& other);

   // Unimplemented assignment operator.
   ExternalDistKrookLayer&
   operator = (
      const ExternalDistKrookLayer& rhs);

   // Parse this class' parameters.
   void
   parseParameters(
      LokiInputParser& a_pp);

   // Check the validity of this class' parameters.
   void
   checkParameters();

   // Configuration space dimensions of the problem.
   tbox::Dimension m_dim;

   // Coordinates of the lower end of the lower ExternalDistKrookLayer in each
   // dimension.
   vector<double> m_x_lo_lo;

   // Coordinates of the upper end of the lower ExternalDistKrookLayer in each
   // dimension.
   vector<double> m_x_hi_lo;

   // Coordinates of the lower end of the upper ExternalDistKrookLayer in each
   // dimension.
   vector<double> m_x_lo_hi;

   // Coordinates of the upper end of the upper ExternalDistKrookLayer in each
   // dimension.
   vector<double> m_x_hi_hi;

   // For each dimension, is true if there is an ExternalDistKrookLayer on the
   // lower side of that dimension.
   deque<bool> m_has_layer_lo;

   // For each dimension, is true if there is an ExternalDistKrookLayer on the
   // upper side of that dimension.
   deque<bool> m_has_layer_hi;

   // Is true if an ExternalDistKrookLayer in any direction has been defined.
   bool m_has_layer;

   // If true the ExternalDistKrookLayer overlaps this processor's extent.
   bool m_overlaps;

   // The overall scaling in the evaluation of each ExternalDistKrookLayer.
   double m_coefficient;

   // The problem domain.
   ProblemDomain m_domain;

   // The damping coefficient.
   ParallelArray m_nu;

   // Name of file holding the data for the External2VIC.
   string m_ext_2v_file;

   // The external distribution initial condition.
   External2VIC* m_external_distribution_ic;
};

} // end namespace Loki

#endif
