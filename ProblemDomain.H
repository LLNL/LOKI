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
#ifndef _PROBLEM_DOMAIN_H_
#define _PROBLEM_DOMAIN_H_

#include "Overture.h"
#include "OvertureTypes.h"
#include "ParmParse.H"
#include "HDF5_DataBase.h"

#include "tbox/Dimension.H"
#include "tbox/IntVector.H"
#include "tbox/Box.H"
#include "tbox/DescribedClass.H"

namespace Loki {

/**
 * Description of the global Cartesian problem domain.  The problem domain will
 * be 2D for configuration space entities (Poisson & Maxwell), or 4D for phase
 * space entites (KineticSpecies).
 */
class ProblemDomain : public tbox::DescribedClass
{
public:
   /**
    * @brief Constructor to create a ProblemDomain from a restart file.
    *        This is used by the post-processor.
    *
    * @param[in] a_dim Dimension of domain.
    * @param[in] a_db Restart database
    */
   ProblemDomain(
      const tbox::Dimension& a_dim,
      HDF_DataBase& a_db);

   /**
    * @brief Constructor to create a ProblemDomain from an input file.
    *        This is used by the System to construct the configuration space
    *        ProblemDomain.
    *
    * @param[in] a_dim Dimension of domain.
    * @param[in] a_spatial_solution_order Spatial order of solution.
    * @param[in] a_pp Where input parameters will be found.
    */
   ProblemDomain(
      const tbox::Dimension& a_dim,
      int a_spatial_solution_order,
      ParmParse& a_pp);

   /**
    * @brief Constructor from caller's data.
    *        This is used by the KineticSpecies which takes the data from a
    *        configuration space ProblemDomain, reads the info for the other
    *        2 dimensions, glues it all together and calls this constructor.
    *
    * @param[in] a_dim Dimension of domain.
    * @param[in] a_n_cells Number of cells in each direction.
    * @param[in] a_x_lo Coordinate of lower end of each direction.
    * @param[in] a_x_hi Coordinate of upper end of each direction.
    * @param[in] a_is_periodic Is true in any direction that is periodic.
    */
   ProblemDomain(
      const tbox::Dimension& a_dim,
      const tbox::IntVector& a_n_cells,
      const Array<double>&   a_x_lo,
      const Array<double>&   a_x_hi,
      const Array<bool>&     a_is_periodic);

   /**
    * @brief Copy constructor.
    *
    * @param[in] a_other Thing to be copied.
    */
   ProblemDomain(
      const ProblemDomain& a_other);

   /**
    * @brief Copy destructor.
    */
   ~ProblemDomain();

   /**
    * @brief Returns the lower coordinate of requested direction.
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
    * @brief Returns the upper coordinate of requested direction.
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
    * @brief Returns the lower coordinate of all directions.
    */
   const Array<double>&
   lower() const
   {
      return m_x_lo;
   }

   /**
    * @brief Returns the upper coordinate of all directions.
    */
   const Array<double>&
   upper() const
   {
      return m_x_hi;
   }

   /**
    * @brief Returns the cell size in requested direction.
    *
    * @param[in] a_dir Direction of interest.
    */
   double
   dx(
      int a_dir) const
   {
      return m_dx[a_dir];
   }

   /**
    * @brief Returns the lower coordinate of requested direction.
    *
    * @param[in] a_dir Direction of interest.
    */
   double
   x_lo(
      int a_dir) const
   {
      return m_x_lo[a_dir];
   }

   /**
    * @brief Returns the upper coordinate of requested direction.
    *
    * @param[in] a_dir Direction of interest.
    */
   double
   x_hi(
      int a_dir) const
   {
      return m_x_hi[a_dir];
   }

   /**
    * @brief Returns the cell size in all directions.
    */
   const Array<double>&
   dx() const
   {
      return m_dx;
   }

   /**
    * @brief Returns true if the domain is periodic in the requested direction.
    *
    * @param[in] a_dir Direction of interest.
    */
   bool
   isPeriodic(
      int a_dir) const
   {
      return m_is_periodic[a_dir];
   }

   /**
    * @brief Returns 1 if the domain is periodic in the requested direction.
    *
    * @param[in] a_dir Direction of interest.
    */
   int
   isPeriodicAsInt(
      int a_dir) const
   {
      return (m_is_periodic[a_dir] ? 1 : 0);
   }

   /**
    * @brief Returns boolean describing periodicity in all directions.
    */
   const Array<bool>&
   isPeriodic() const
   {
      return m_is_periodic;
   }

   /**
    * @brief Returns the number of cells in the requested direction.
    *
    * @param[in] a_dir Direction of interest.
    */
   int
   numberOfCells(
      int a_dir) const
   {
      return m_n_cells[a_dir];
   }

   /**
    * @brief Returns the number of cells in all directions.
    */
   const tbox::IntVector&
   numberOfCells() const
   {
      return m_n_cells;
   }

   /**
    * @brief Returns box describing global problem domain.
    */
   const tbox::Box&
   box() const
   {
      return m_box;
   }

   /**
    * @brief Returns dimension of global problem domain.
    */
   const tbox::Dimension&
   dim() const
   {
      return m_dim;
   }

   /**
    * @brief Print this class' parameters.
    */
   void
   printParameters() const;

   /**
    * @brief Write this class to a restart database.
    *
    * @param[out] a_db Restart database.
    */
   void
   putToDatabase(
      HDF_DataBase& a_db) const;

   /**
    * @brief Read global problem domain parameters from a restart database.
    *
    * @param[in] a_db Restart database.
    */
   void
   getFromDatabase(
      const HDF_DataBase& a_db);

private:
   // Unimplemented default constructor.
   ProblemDomain();

   // Unimplemented assignment operator.
   ProblemDomain&
   operator = (
      const ProblemDomain& rhs);

   // Read input parameters from the top level of the input database.
   void
   parseParameters(
      int a_spatial_solution_order,
      ParmParse& a_pp);

   // Set up m_dx and m_box from input parameters.
   void
   define();

   // Dimension of global problem domain.
   tbox::Dimension m_dim;

   // Lower coordinates in each direction.
   Array<double> m_x_lo;

   // Upper coordinates in each direction.
   Array<double> m_x_hi;

   // Cell size in each direction.
   Array<double> m_dx;

   // Number of cells in each direction.
   tbox::IntVector m_n_cells;

   // Box describing global problem domain.
   tbox::Box m_box;

   // True for each direction that is periodic.
   Array<bool> m_is_periodic;
};

/**
 * @brief Macro to extract information about extent of ProblemDomain for
 * fortran calls.
 */
#define PROBLEMDOMAIN_TO_FORT(__pd) \
   *(__pd.lower()).dataPtr(), \
   *(__pd.upper()).dataPtr(), \
   *(__pd.dx()).dataPtr()

/**
 * @brief Macro to define information about extent of ProblemDomain for fortran
 * calls.
 */
#define PROBLEMDOMAIN_TO_FORT_PROTO \
   const real& xlo, const real& xhi, const real& dx

} // end namespace Loki

#endif
