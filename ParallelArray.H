/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _PARALLEL_ARRAY_H_
#define _PARALLEL_ARRAY_H_

#include <deque>
#include <vector>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string.h>

using namespace std;

// ARRAY_ASSERT is the error checking macro.  As the performance of this class
// is central to the performance of applications using it error checking may
// be disabled by including -DDEBUG_CHECK_ASSERTIONS in the compilation.
#define ARRAY_ABORT(message) \
{ \
   ostringstream os; \
   os << "Error occured in file " << __FILE__ << ", function " << __func__ << ", line " << __LINE__ << "." << endl; \
   cout << os.str();                \
   cout << message << endl;         \
   cout << "I am now going to purposely abort so that you can get a traceback from a debugger" << endl;             \
   abort();                         \
}

#ifdef DEBUG_CHECK_ASSERTIONS

#define ARRAY_ASSERT(EXP) \
  do {                                                \
      if (!(EXP)) {                                   \
         ostringstream os;                            \
         os << "Failed assertion: " << # EXP << ends; \
         ARRAY_ABORT(os.str());                       \
      }                                               \
   } while (0)
#else

/*
 * No assertion checking
 */
#define ARRAY_ASSERT(EXP)

#endif

/**
 * A class for N dimensional arrays of doubles distributed over a set of
 * processors.  More precisely, ParallelArray describes a potentially
 * multi-component field defined on an N dimensional base space.  Currently the
 * base space may have a maximum dimension of 6.
 *
 * The base space is represented by a 0 based index space.  The lower end of
 * each dimension in the base space has index 0 while the upper end has index
 * number of cells -1.  The base space contains a layer of cells for boundary
 * conditions.  This boundary layer is number of ghost zones wide where the
 * number of ghost zones depends on the spatial order of accuracy.  Thus the
 * index space of each dimension of the global ParallelArray spans
 * [-nghost, ncells - 1 + nghost].  Each domain's index space is a subset of
 * the global index space.  An interior box which represents the index space of
 * zones stricly on the interior of a domain is maintained as is a local box
 * which is the interior box plus boundary zones and a data box which is the
 * interior box plus boundary and ghost zones.

 * At the time that the array is partitioned among its processors the details of
 * all the communication necessary to update ghost data and periodic boundary
 * conditions is determined.  This is done in such a way that the resulting
 * communication is deterministic and needs no synchronization points
 * (MPI_Barrier).  This scheme will result in more ghost/periodic boundary data
 * packets being sent than in other approaches.  However, the total amount of
 * ghost/periodic boundary data transferred is the minimum required and no
 * "meta-data" needs to be communicated.  The periodic boundary condition
 * communication information is set up in a specific manner.  Only periodic
 * boundary "faces" are communicated.  Periodic boundaries are communicated
 * sequentially from the lowest to the highest dimension.  The communicated
 * data has been determined in such a way that communication of a higher
 * dimension periodic boundary face also handles the communication of all
 * periodic boundary "edges" or "corners" formed with lower dimension periodic
 * boundaries.  This greatly simplifies what needs to be determined in order to
 * handle periodic data.
 *
 * The public interface is intentionally thin.  There are dimension specific
 * "()" operators to optimize performance as these operators are likely to
 * appear in the innermost of a set of loops over dimensions.
 */
class ParallelArray
{
public:
   // A simple class to describe the index space of different segments of the
   // array such as the interior zones or interior+ghost+boundary zones.
   class Box
   {
      public:
         // Default constructor.  Creates and undimensioned Box.
         Box();

         // Creates an empty Box of a given dimension.
         Box(
            int a_dim);

         // Destructor.
         ~Box();

         // Copy constructor.
         Box(
            const Box& a_other);

         // Assignment operator.
         Box&
         operator = (
            const Box& a_rhs);

         // Redimensions a Box.
         void
         redimension(
            int a_dim)
         {
            m_dim = a_dim;
            m_upper.assign(a_dim, 0);
            m_lower.assign(a_dim, -1);
         }

         // Returns the dimension of the Box.
         int
         dim() const
         {
            ARRAY_ASSERT(m_dim != -1);
            return m_dim;
         }

         // Returns a const upper index space bound for the requested dimension.
         const int&
         upper(
            int a_dim) const
         {
            ARRAY_ASSERT(m_dim != -1);
            return m_upper[a_dim];
         }

         // Returns a mutable upper index space bound for the requested
         // dimension.
         int&
         upper(
            int a_dim)
         {
            ARRAY_ASSERT(m_dim != -1);
            return m_upper[a_dim];
         }

         // Returns a const lower index space bound for the requested dimension.
         const int&
         lower(
            int a_dim) const
         {
            ARRAY_ASSERT(m_dim != -1);
            return m_lower[a_dim];
         }

         // Returns a mutable lower index space bound for the requested
         // dimension.
         int&
         lower(
            int a_dim)
         {
            ARRAY_ASSERT(m_dim != -1);
            return m_lower[a_dim];
         }

         // Returns the number of cells in the requested dimension.
         int
         numberOfCells(
            int a_dim) const
         {
            ARRAY_ASSERT(m_dim != -1);
            return upper(a_dim) - lower(a_dim) + 1;
         }

         // Returns the total number of cells in the box.
         int
         size() const
         {
            ARRAY_ASSERT(m_dim != -1);
            int result = 1;
            for (int i = 0; i < m_dim; ++i) {
               result *= upper(i) - lower(i) + 1;
            }
            return result;
         }

         // Returns true if the the size is 0.
         bool
         empty() const
         {
            return size() == 0;
         }

         // Construct the intersection of this box with a_other.
         Box
         intersect(
            const Box& a_other) const
         {
            ARRAY_ASSERT(m_dim == a_other.m_dim);
            Box intersection(m_dim);
            for (int i = 0; i < m_dim; ++i) {
               int max_lo = max(lower(i), a_other.lower(i));
               int min_hi = min(upper(i), a_other.upper(i));
               if (max_lo > min_hi) {
                  for (int j = 0; j < m_dim; ++j) {
                     intersection.lower(i) = 0;
                     intersection.upper(i) = -1;
                  }
                  return intersection;
               }
               intersection.lower(i) = max_lo;
               intersection.upper(i) = min_hi;
            }
            return intersection;
         }

         bool
         operator == (
            const ParallelArray::Box& a_other) const
         {
            bool are_equal = m_dim == a_other.m_dim &&
               m_upper == a_other.m_upper && m_lower == a_other.m_lower;
         }

      private:
         // The dimension of the box.
         int m_dim;

         // The upper index space bound of each dimension.
         vector<int> m_upper;

         // The lower index space bound of each dimension.
         vector<int> m_lower;
   };

   // Each of the 3^(m_num_dist_dim) communications is one of these types.
   enum Comm_Type {
      NO_COMM = -1,      // Indicator for no communication type.
      X1_PER_COMM,       // Communication of boundary periodic only in x1.
      X2_PER_COMM,       // Communication of boundary periodic only in x2.
      X3_PER_COMM,       // Communication of boundary periodic only in x3.
      X4_PER_COMM,       // Communication of boundary periodic only in x4.
      X5_PER_COMM,       // Communication of boundary periodic only in x5.
      X6_PER_COMM,       // Communication of boundary periodic only in x6.
      X1_PER_COMM_LOCAL, // Local data copy of boundary periodic only in x1.
      X2_PER_COMM_LOCAL, // Local data copy of boundary periodic only in x2.
      X3_PER_COMM_LOCAL, // Local data copy of boundary periodic only in x3.
      X4_PER_COMM_LOCAL, // Local data copy of boundary periodic only in x4.
      X5_PER_COMM_LOCAL, // Local data copy of boundary periodic only in x5.
      X6_PER_COMM_LOCAL, // Local data copy of boundary periodic only in x6.
      CORNER_PER_COMM,   // Communication of boundary periodic in multiple
                         // dimensions which is a no-op.
      GHOST_COMM,        // Communication of ghost data.
      SELF_COMM,         // Communication with self which is a no-op.
      NUM_COMM_TYPES
   };

   /**
    * @brief Default constructor.  Creates an undimensioned and unpartitioned
    *        ParallelArray.  A call to either of the partition methods which
    *        set the dimension and number of distributed dimensions must be
    *        made before attemping to invoke any other method.
    */
   ParallelArray();

   /**
    * @brief Constructor to create a dimensioned but unpartitioned
    *        ParallelArray.  A call to either of the partition methods which do
    *        not set the dimension and number of distributed dimensions must be
    *        made before attemping to invoke any other method.
    *
    * @param[in] a_dim The dimension of the array.
    * @param[in] a_num_dist_dim Dimensions 0 through a_num_dist_dim-1 are
    *                           distributed among the processors.  The index
    *                           space of higher dimensions is common to all
    *                           processors.
    */
   ParallelArray(
      int a_dim,
      int a_num_dist_dim);

   /**
    * @brief Constructor to create a dimensioned and partitioned ParallelArray.
    *        Following construction with this constructor the array is fully
    *        functional.
    *
    * @param[in] a_dim The dimension of the array.
    * @param[in] a_num_dist_dim Dimensions 0 through a_num_dist_dim-1 are
    *                           distributed among the processors.  The index
    *                           space of higher dimensions is common to all
    *                           processors.
    * @param[in] a_rank_lo The lowest ranking process in the range of processes
    *                      the array is distribued among.
    * @param[in] a_rank_hi The highest ranking process in the range of processes
    *                      the array is distribued among.
    * @param[in] a_nghosts The number of ghost cells required.
    * @param[in] a_is_periodic True for each dimension which is periodic.
    * @param[in] a_num_cells Number of cells in the problem in each dimension.
    */
   ParallelArray(
      int a_dim,
      int a_num_dist_dim,
      int a_rank_lo,
      int a_rank_hi,
      int a_nghosts,
      const deque<bool>& a_is_periodic,
      const vector<int>& a_num_cells);

   /**
    * @brief Constructor to create a dimensioned ParallelArray on the specified
    *        base space which is partitioned on a single processor.  This is
    *        really a local rather than a parallel array.  It represents the
    *        local part of a distributed array.  Following construction with
    *        this constructor the array is fully functional.  As this array is
    *        local, the communication methods are no-ops.
    *
    * @param[in] a_base_space The index space over which the array is defined.
    * @param[in] a_num_dist_dim Dimensions 0 through a_num_dist_dim-1 are
    *                           distributed among the processors.  The index
    *                           space of higher dimensions is common to all
    *                           processors.
    * @param[in] a_nghosts The number of ghost cells required.
    * @param[in] a_num_global_cells The global extent of which this
    *                               ParallelArray is a local part.
    */
   ParallelArray(
      const Box& a_base_space,
      int a_num_dist_dim,
      int a_nghosts,
      const vector<int>& a_num_global_cells);

   /**
    * @brief Copy constructor.
    *
    * @param[in] a_other What to copy.
    */
   ParallelArray(
      const ParallelArray& a_other);

   /**
    * @brief Assignment operator.
    *
    * @param[in] a_rhs Right hand side of assignment.
    */
   ParallelArray&
   operator = (
      const ParallelArray& a_rhs);

   /**
    * @brief Destructor.
    */
   ~ParallelArray();

   /**
    * @brief To be called after construction with the default constructor.
    *        Determines the partitioning of the array given a range of
    *        processors and required number of ghosts.  This method also
    *        determines the processors and data sent and received in ghost data
    *        and periodic boundary condition communications.
    * @param[in] a_dim The dimension of the array.
    * @param[in] a_num_dist_dim Dimensions 0 through a_num_dist_dim-1 are
    *                           distributed among the processors.  The index
    *                           space of higher dimensions is common to all
    *                           processors.
    * @param[in] a_rank_lo The lowest ranking process in the range of processes
    *                      the array is distribued among.
    * @param[in] a_rank_hi The highest ranking process in the range of processes
    *                      the array is distribued among.
    * @param[in] a_nghosts The number of ghost cells required.
    * @param[in] a_is_periodic True for each dimension which is periodic.
    * @param[in] a_num_cells Number of cells in the problem in each dimension.
    */
   void
   partition(
      int a_dim,
      int a_num_dist_dim,
      int a_rank_lo,
      int a_rank_hi,
      int a_nghosts,
      const deque<bool>& a_is_periodic,
      const vector<int>& a_num_cells);

   /**
    * @brief Determines the partitioning of the array given a range of
    *        processors and required number of ghosts.  This method also
    *        determines the processors and data sent and received in ghost data
    *        and periodic boundary condition communications.
    * @param[in] a_rank_lo The lowest ranking process in the range of processes
    *                      the array is distribued among.
    * @param[in] a_rank_hi The highest ranking process in the range of processes
    *                      the array is distribued among.
    * @param[in] a_nghosts The number of ghost cells required.
    * @param[in] a_is_periodic True for each dimension which is periodic.
    * @param[in] a_num_cells Number of cells in the problem in each dimension.
    */
   void
   partition(
      int a_rank_lo,
      int a_rank_hi,
      int a_nghosts,
      const deque<bool>& a_is_periodic,
      const vector<int>& a_num_cells);

   /**
    * @brief To be called after construction with the default constructor.
    *        Places the array on the supplied base space and partitions it on a
    *        single processor.  An array partitioned with this method is
    *        effectively an undistributed local array.
    *
    * @param[in] a_base_space The index space over which the array is defined.
    * @param[in] a_num_dist_dim Dimensions 0 through a_num_dist_dim-1 are
    *                           distributed among the processors.  The index
    *                           space of higher dimensions is common to all
    *                           processors.
    * @param[in] a_nghosts The number of ghost cells required.
    * @param[in] a_num_global_cells The global extent of which this
    *                               ParallelArray is a local part.
    */
   void
   partition(
      const Box& a_base_space,
      int a_num_dist_dim,
      int a_nghosts,
      const vector<int>& a_num_global_cells);

   /**
    * @brief Places the array on the supplied base space and partitions it on a
    *        single processor.  An array partitioned with this method is
    *        effectively an undistributed local array.
    *
    * @param[in] a_base_space The index space over which the array is defined.
    * @param[in] a_nghosts The number of ghost cells required.
    * @param[in] a_num_global_cells The global extent of which this
    *                               ParallelArray is a local part.
    */
   void
   partition(
      const Box& a_base_space,
      int a_nghosts,
      const vector<int>& a_num_global_cells);

   /**
    * @brief Returns the number of partitions in each partitioned dimension.
    */
   const vector<int>&
   getDimPartitions() const
   {
      return m_dim_partitions;
   }

   /**
    * @brief Returns the number of cells in partitions to the left (lower) end
    *        of each dimension.
    */
   const vector<int>&
   getLeftPartitionsSize() const
   {
      ARRAY_ASSERT(m_partitioned);
      return m_left_partitions_size;
   }

   /**
    * @brief Returns the number of partitions to the left (lower) end of each
    *        dimension.
    */
   const vector<int>&
   getNumLeftPartitions() const
   {
      ARRAY_ASSERT(m_partitioned);
      return m_num_left_partitions;
   }

   /**
    * @brief Returns the dimension of the array.
    */
   int
   dim() const
   {
      return m_dim;
   }

   /**
    * @brief Returns the number of distributed dimensions of the array.
    */
   int
   distDim() const
   {
      return m_num_dist_dim;
   }

   /**
    * @brief Returns the number of ghost zones.
    */
   int
   numGhosts() const
   {
      ARRAY_ASSERT(m_partitioned);
      return m_nghosts;
   }

   /**
    * @brief Returns the global rank given a 0 based index rank in the
    *        partitions of each dimension.
    *
    * @param[in] a_idx_rank The 0 based index rank, one index for each
    *                       partitioned dimension, of the partition whose
    *                       global rank is to be determined.
    */
   int
   getGlobalRank(
      const vector<int> a_idx_rank) const;

   /**
    * @brief Returns the index rank of the array on the processor with the
    *        given the global rank.
    *
    * @param[in] a_global_rank The global rank.
    * @param[out] a_index_rank The 0 based index rank, one index for each
    *                          partitioned dimension, of the array on the
    *                          processor with the supplied global rank.
    */
   void
   getIndexRank(
      int a_global_rank,
      vector<int>& a_index_rank) const
   {
      ARRAY_ASSERT(m_partitioned);
      ARRAY_ASSERT(a_global_rank >= m_proc_lo);
      ARRAY_ASSERT(a_global_rank <= m_proc_hi);
      ARRAY_ASSERT(a_index_rank.size() == m_dim);
      getIndexRank(a_global_rank, m_proc_lo, m_dim_partitions, a_index_rank);
   }

   /**
    * @brief Returns the lowest rank processor this array is distributed over.
    */
   int
   procLo() const
   {
      return m_proc_lo;
   }

   /**
    * @brief Returns the highest rank processor this array is distributed over.
    */
   int
   procHi() const
   {
      return m_proc_hi;
   }

   /**
    * @brief Communicates ghost data among the processes that this array is
    *        distributed over.
    */
   void
   communicateGhostData()
   {
      ARRAY_ASSERT(m_partitioned);
      communicateData(GHOST_COMM, NO_COMM);
   }

   /**
    * @brief Communicate periodic boundaries.  As mentioned above, the periodic
    *        boundaries are communicated sequentially from lowest to highest
    *        dimension.  CORNER_PER_COMM cases are automatically handled by
    *        virtue of how the periodic communicated data is set up.
    */
   void
   communicatePeriodicBoundaries()
   {
      ARRAY_ASSERT(m_partitioned);
      for (int i = 0; i < m_num_dist_dim; ++i) {
         if (m_periodic[i]) {
            if (i == 0) {
               communicateData(X1_PER_COMM, X1_PER_COMM_LOCAL);
            }
            else if (i == 1) {
               communicateData(X2_PER_COMM, X2_PER_COMM_LOCAL);
            }
            else if (i == 2) {
               communicateData(X3_PER_COMM, X3_PER_COMM_LOCAL);
            }
            else if (i == 3) {
               communicateData(X4_PER_COMM, X4_PER_COMM_LOCAL);
            }
            else if (i == 4) {
               communicateData(X5_PER_COMM, X5_PER_COMM_LOCAL);
            }
            else if (i == 5) {
               communicateData(X6_PER_COMM, X6_PER_COMM_LOCAL);
            }
         }
      }
   }

   /**
    * @brief Returns the box describing the extent of the interior of this
    *        parallel array.  This box does not include ghost/boundary cells.
    */
   const Box&
   interiorBox() const
   {
      ARRAY_ASSERT(m_partitioned);
      return m_interior_box;
   }

   /**
    * @brief Returns the box describing the extent of the interior of this
    *        parallel array on the processor with the given rank.
    *
    * @param[in] a_rank The rank whose interior box is requested.
    */
   Box
   interiorBox(
      int a_rank) const;

   /**
    * @brief Returns the box describing the extent of the interior of this
    *        parallel array along with its boundary cells.
    */
   const Box&
   localBox() const
   {
      ARRAY_ASSERT(m_partitioned);
      return m_local_box;
   }

   /**
    * @brief Returns the box describing the extent of the interior of this
    *        parallel array along with its boundary cells on the processor
    *        with the given rank.
    *
    * @param[in] a_rank The rank whose local box is requested.
    */
   Box
   localBox(
      int a_rank) const;

   /**
    * @brief Returns the box describing the extent of the interior of this
    *        parallel array along with its ghost and boundary cells.
    */
   const Box&
   dataBox() const
   {
      ARRAY_ASSERT(m_partitioned);
      return m_data_box;
   }

   /**
    * @brief Returns the box describing the extent of the interior of this
    *        parallel array along with its ghost/boundary cells on the processor
    *        with the given rank.
    *
    * @param[in] a_rank The rank whose data box is requested.
    */
   Box
   dataBox(
      int a_rank) const;

   /**
    * @brief Returns a const pointer to the array's underlying data.
    */
   const double*
   getData() const
   {
      ARRAY_ASSERT(m_partitioned);
      return m_data;
   }

   /**
    * @brief Returns a mutable pointer to the array's underlying data.
    */
   double*
   getData()
   {
      ARRAY_ASSERT(m_partitioned);
      return m_data;
   }

   /**
    * @brief Sets all values of array to specified value.
    *
    * @param[in] a_new_val The new array value.
    */
   void
   operator = (
      double a_new_val)
   {
      ARRAY_ASSERT(m_partitioned);
      for (int i = 0; i < m_data_box.size(); ++i) {
         m_data[i] = a_new_val;
      }
   }

   /**
    * @brief Multiplies all values of the array by specified factor.
    *
    * @param[in] a_factor The factor.
    */
   void
   operator *= (
      double a_factor)
   {
      ARRAY_ASSERT(m_partitioned);
      for (int i = 0; i < m_data_box.size(); ++i) {
         m_data[i] *= a_factor;
      }
   }

   /**
    * @brief Adds a specified value to all values of the array.
    *
    * @param[in] a_addend The value to add.
    */
   void
   operator += (
      double a_addend)
   {
      ARRAY_ASSERT(m_partitioned);
      for (int i = 0; i < m_data_box.size(); ++i) {
         m_data[i] += a_addend;
      }
   }

   /**
    * @brief Adds another ParallelArray to this.
    *
    * @param[in] a_addend The array to add.
    */
   void
   operator += (
      const ParallelArray& a_addend)
   {
      ARRAY_ASSERT(m_partitioned);
      ARRAY_ASSERT(a_addend.m_partitioned);
      ARRAY_ASSERT(m_data_box == a_addend.m_data_box);
      for (int i = 0; i < m_data_box.size(); ++i) {
         m_data[i] += a_addend.m_data[i];
      }
   }

   /**
    * @brief Mutable data accessor for a 1D array.
    */
   double&
   operator () (
      int a_i0)
   {
      ARRAY_ASSERT(m_partitioned);
      ARRAY_ASSERT(m_dim == 1);
      ARRAY_ASSERT(a_i0 >= m_data_box.lower(0) && a_i0 <= m_data_box.upper(0));
      return m_data[dataArrayIndex(a_i0)];
   }

   /**
    * @brief Const data accessor for a const 1D array.
    */
   const double&
   operator () (
      int a_i0) const
   {
      ARRAY_ASSERT(m_partitioned);
      ARRAY_ASSERT(m_dim == 1);
      ARRAY_ASSERT(a_i0 >= m_data_box.lower(0) && a_i0 <= m_data_box.upper(0));
      return m_data[dataArrayIndex(a_i0)];
   }

   /**
    * @brief Mutable data accessor for a 2D array.
    */
   double&
   operator () (
      int a_i0,
      int a_i1)
   {
      ARRAY_ASSERT(m_partitioned);
      ARRAY_ASSERT(m_dim == 2);
      ARRAY_ASSERT(a_i0 >= m_data_box.lower(0) && a_i0 <= m_data_box.upper(0));
      ARRAY_ASSERT(a_i1 >= m_data_box.lower(1) && a_i1 <= m_data_box.upper(1));
      return m_data[dataArrayIndex(a_i0, a_i1)];
   }

   /**
    * @brief Const data accessor for a const 2D array.
    */
   const double&
   operator () (
      int a_i0,
      int a_i1) const
   {
      ARRAY_ASSERT(m_partitioned);
      ARRAY_ASSERT(m_dim == 2);
      ARRAY_ASSERT(a_i0 >= m_data_box.lower(0) && a_i0 <= m_data_box.upper(0));
      ARRAY_ASSERT(a_i1 >= m_data_box.lower(1) && a_i1 <= m_data_box.upper(1));
      return m_data[dataArrayIndex(a_i0, a_i1)];
   }

   /**
    * @brief Mutable data accessor for a 3D array.
    */
   double&
   operator () (
      int a_i0,
      int a_i1,
      int a_i2)
   {
      ARRAY_ASSERT(m_partitioned);
      ARRAY_ASSERT(m_dim == 3);
      ARRAY_ASSERT(a_i0 >= m_data_box.lower(0) && a_i0 <= m_data_box.upper(0));
      ARRAY_ASSERT(a_i1 >= m_data_box.lower(1) && a_i1 <= m_data_box.upper(1));
      ARRAY_ASSERT(a_i2 >= m_data_box.lower(2) && a_i2 <= m_data_box.upper(2));
      return m_data[dataArrayIndex(a_i0, a_i1, a_i2)];
   }

   /**
    * @brief Const data accessor for a const 3D array.
    */
   const double&
   operator () (
      int a_i0,
      int a_i1,
      int a_i2) const
   {
      ARRAY_ASSERT(m_partitioned);
      ARRAY_ASSERT(m_dim == 3);
      ARRAY_ASSERT(a_i0 >= m_data_box.lower(0) && a_i0 <= m_data_box.upper(0));
      ARRAY_ASSERT(a_i1 >= m_data_box.lower(1) && a_i1 <= m_data_box.upper(1));
      ARRAY_ASSERT(a_i2 >= m_data_box.lower(2) && a_i2 <= m_data_box.upper(2));
      return m_data[dataArrayIndex(a_i0, a_i1, a_i2)];
   }

   /**
    * @brief Mutable data accessor for a 4D array.
    */
   double&
   operator () (
      int a_i0,
      int a_i1,
      int a_i2,
      int a_i3)
   {
      ARRAY_ASSERT(m_partitioned);
      ARRAY_ASSERT(m_dim == 4);
      ARRAY_ASSERT(a_i0 >= m_data_box.lower(0) && a_i0 <= m_data_box.upper(0));
      ARRAY_ASSERT(a_i1 >= m_data_box.lower(1) && a_i1 <= m_data_box.upper(1));
      ARRAY_ASSERT(a_i2 >= m_data_box.lower(2) && a_i2 <= m_data_box.upper(2));
      ARRAY_ASSERT(a_i3 >= m_data_box.lower(3) && a_i3 <= m_data_box.upper(3));
      return m_data[dataArrayIndex(a_i0, a_i1, a_i2, a_i3)];
   }

   /**
    * @brief Const data accessor for a const 4D array.
    */
   const double&
   operator () (
      int a_i0,
      int a_i1,
      int a_i2,
      int a_i3) const
   {
      ARRAY_ASSERT(m_partitioned);
      ARRAY_ASSERT(m_dim == 4);
      ARRAY_ASSERT(a_i0 >= m_data_box.lower(0) && a_i0 <= m_data_box.upper(0));
      ARRAY_ASSERT(a_i1 >= m_data_box.lower(1) && a_i1 <= m_data_box.upper(1));
      ARRAY_ASSERT(a_i2 >= m_data_box.lower(2) && a_i2 <= m_data_box.upper(2));
      ARRAY_ASSERT(a_i3 >= m_data_box.lower(3) && a_i3 <= m_data_box.upper(3));
      return m_data[dataArrayIndex(a_i0, a_i1, a_i2, a_i3)];
   }

   /**
    * @brief Mutable data accessor for a 5D array.
    */
   double&
   operator () (
      int a_i0,
      int a_i1,
      int a_i2,
      int a_i3,
      int a_i4)
   {
      ARRAY_ASSERT(m_partitioned);
      ARRAY_ASSERT(m_dim == 5);
      ARRAY_ASSERT(a_i0 >= m_data_box.lower(0) && a_i0 <= m_data_box.upper(0));
      ARRAY_ASSERT(a_i1 >= m_data_box.lower(1) && a_i1 <= m_data_box.upper(1));
      ARRAY_ASSERT(a_i2 >= m_data_box.lower(2) && a_i2 <= m_data_box.upper(2));
      ARRAY_ASSERT(a_i3 >= m_data_box.lower(3) && a_i3 <= m_data_box.upper(3));
      ARRAY_ASSERT(a_i4 >= m_data_box.lower(4) && a_i4 <= m_data_box.upper(4));
      return m_data[dataArrayIndex(a_i0, a_i1, a_i2, a_i3, a_i4)];
   }

   /**
    * @brief Const data accessor for a const 5D array.
    */
   const double&
   operator () (
      int a_i0,
      int a_i1,
      int a_i2,
      int a_i3,
      int a_i4) const
   {
      ARRAY_ASSERT(m_partitioned);
      ARRAY_ASSERT(m_dim == 5);
      ARRAY_ASSERT(a_i0 >= m_data_box.lower(0) && a_i0 <= m_data_box.upper(0));
      ARRAY_ASSERT(a_i1 >= m_data_box.lower(1) && a_i1 <= m_data_box.upper(1));
      ARRAY_ASSERT(a_i2 >= m_data_box.lower(2) && a_i2 <= m_data_box.upper(2));
      ARRAY_ASSERT(a_i3 >= m_data_box.lower(3) && a_i3 <= m_data_box.upper(3));
      ARRAY_ASSERT(a_i4 >= m_data_box.lower(4) && a_i4 <= m_data_box.upper(4));
      return m_data[dataArrayIndex(a_i0, a_i1, a_i2, a_i3, a_i4)];
   }

   /**
    * @brief Mutable data accessor for a 6D array.
    */
   double&
   operator () (
      int a_i0,
      int a_i1,
      int a_i2,
      int a_i3,
      int a_i4,
      int a_i5)
   {
      ARRAY_ASSERT(m_partitioned);
      ARRAY_ASSERT(m_dim == 6);
      ARRAY_ASSERT(a_i0 >= m_data_box.lower(0) && a_i0 <= m_data_box.upper(0));
      ARRAY_ASSERT(a_i1 >= m_data_box.lower(1) && a_i1 <= m_data_box.upper(1));
      ARRAY_ASSERT(a_i2 >= m_data_box.lower(2) && a_i2 <= m_data_box.upper(2));
      ARRAY_ASSERT(a_i3 >= m_data_box.lower(3) && a_i3 <= m_data_box.upper(3));
      ARRAY_ASSERT(a_i4 >= m_data_box.lower(4) && a_i4 <= m_data_box.upper(4));
      ARRAY_ASSERT(a_i5 >= m_data_box.lower(5) && a_i5 <= m_data_box.upper(5));
      return m_data[dataArrayIndex(a_i0, a_i1, a_i2, a_i3, a_i4, a_i5)];
   }

   /**
    * @brief Const data accessor for a const 6D array.
    */
   const double&
   operator () (
      int a_i0,
      int a_i1,
      int a_i2,
      int a_i3,
      int a_i4,
      int a_i5) const
   {
      ARRAY_ASSERT(m_partitioned);
      ARRAY_ASSERT(m_dim == 6);
      ARRAY_ASSERT(a_i0 >= m_data_box.lower(0) && a_i0 <= m_data_box.upper(0));
      ARRAY_ASSERT(a_i1 >= m_data_box.lower(1) && a_i1 <= m_data_box.upper(1));
      ARRAY_ASSERT(a_i2 >= m_data_box.lower(2) && a_i2 <= m_data_box.upper(2));
      ARRAY_ASSERT(a_i3 >= m_data_box.lower(3) && a_i3 <= m_data_box.upper(3));
      ARRAY_ASSERT(a_i4 >= m_data_box.lower(4) && a_i4 <= m_data_box.upper(4));
      ARRAY_ASSERT(a_i5 >= m_data_box.lower(5) && a_i5 <= m_data_box.upper(5));
      return m_data[dataArrayIndex(a_i0, a_i1, a_i2, a_i3, a_i4, a_i5)];
   }

private:
   // A class to describe the information to be passed between this part of the
   // parallel array and the part on another processor.
   class Communication
   {
      public:
         // Constructor for the communication information of a array of a given
         // dimension.
         Communication(
            int a_dim);

         // Destructor.
         ~Communication();

         // Copy constructor.
         Communication(
            const Communication& a_other);

         // Assignment operator.
         Communication&
         operator = (
            const Communication& a_rhs);

         // The rank of the other domain in this communication.
         int m_partner_rank;

         // The classification of this commnication.
         Comm_Type m_comm_type;

         // Box describing the data to be sent.
         Box m_send_box;

         // Tag to add to s_TAG_BASE to form MPI tag for this send box.
         int m_send_tag;

         // Box describing where received data goes.
         Box m_recv_box;

         // Tag to add to s_TAG_BASE to form MPI tag for this receive box.
         int m_recv_tag;

      private:
         // Unimplemented default constructor.
         Communication();
   };

   // Fills m_dim_partitions with the number of partitions in each dimension
   // given the number of processes the array is distributed over.  The scheme
   // used tries to form the "squarest" possible domains.
   void
   assignProcs();

   // Returns the largest prime factor of a_n.
   int
   greatestPrimeFactor(
      int a_n) const;

   // Returns the smallest prime factor of a_n.
   int
   leastPrimeFactor(
      int a_n) const;

   // Given the global rank, the lowest rank the array is distributed over, and
   // a vector of the number of processors in each dimension fills in a vector
   // of index ranks in each dimension.
   void
   getIndexRank(
      int a_global_rank,
      int a_rank_lo,
      const vector<int> a_dim_partitions,
      vector<int>& a_index_rank) const;

   // Determines the communication layout, the extent of the interior, and the
   // extent of the interior with ghosts/boundary.
   void
   setupLocalDomain();

   // Returns the MPI tag for the receive box associated with the neighbor with
   // the supplied index shift.
   int
   getRecvTag(
      const vector<int>& a_index_shift)
   {
      int tag = 0;
      int shift = 1;
      for (int dim = 0; dim < m_num_dist_dim; ++dim) {
         tag += pow(3, dim)*a_index_shift[dim];
         shift *= 3;
      }
      tag += (shift-1)/2+1;
      return tag;
   }

   // Returns the MPI tag for the send box associated with the neighbor with
   // the supplied index shift.
   int
   getSendTag(
      const vector<int>& a_index_shift)
   {
      int tag = 0;
      int shift = 1;
      for (int dim = 0; dim < m_num_dist_dim; ++dim) {
         tag -= pow(3, dim)*a_index_shift[dim];
         shift *= 3;
      }
      tag += (shift-1)/2+1;
      return tag;
   }

   // Communicates the requested types of data given the communication layout
   // for this array's distribution.
   void
   communicateData(
      Comm_Type a_parallel_comm_type,
      Comm_Type a_local_comm_type);

   // Computes the offset into m_data for a 1D array given the requested index.
   int
   dataArrayIndex(
      int a_i0) const
   {
      return a_i0-m_data_box.lower(0);
   }

   // Computes the offset into m_data for a 2D array given the requested
   // indexes.
   int
   dataArrayIndex(
      int a_i0,
      int a_i1) const
   {
      return (a_i1-m_data_box.lower(1))*m_data_box.numberOfCells(0)+
         a_i0-m_data_box.lower(0);
   }

   // Computes the offset into m_data for a 3D array given the requested
   // indexes.
   int
   dataArrayIndex(
      int a_i0,
      int a_i1,
      int a_i2) const
   {
      return ((a_i2-m_data_box.lower(2))*m_data_box.numberOfCells(1)+
         a_i1-m_data_box.lower(1))*m_data_box.numberOfCells(0)+
         a_i0-m_data_box.lower(0);
   }

   // Computes the offset into m_data for a 4D array given the requested
   // indexes.
   int
   dataArrayIndex(
      int a_i0,
      int a_i1,
      int a_i2,
      int a_i3) const
   {
      return (((a_i3-m_data_box.lower(3))*m_data_box.numberOfCells(2)+
         a_i2-m_data_box.lower(2))*m_data_box.numberOfCells(1)+
         a_i1-m_data_box.lower(1))*m_data_box.numberOfCells(0)+
         a_i0-m_data_box.lower(0);
   }

   // Computes the offset into m_data for a 5D array given the requested
   // indexes.
   int
   dataArrayIndex(
      int a_i0,
      int a_i1,
      int a_i2,
      int a_i3,
      int a_i4) const
   {
      return ((((a_i4-m_data_box.lower(4))*m_data_box.numberOfCells(3)+
         a_i3-m_data_box.lower(3))*m_data_box.numberOfCells(2)+
         a_i2-m_data_box.lower(2))*m_data_box.numberOfCells(1)+
         a_i1-m_data_box.lower(1))*m_data_box.numberOfCells(0)+
         a_i0-m_data_box.lower(0);
   }

   // Computes the offset into m_data for a 6D array given the requested
   // indexes.
   int
   dataArrayIndex(
      int a_i0,
      int a_i1,
      int a_i2,
      int a_i3,
      int a_i4,
      int a_i5) const
   {
      return (((((a_i5-m_data_box.lower(5))*m_data_box.numberOfCells(4)+
         a_i4-m_data_box.lower(4))*m_data_box.numberOfCells(3)+
         a_i3-m_data_box.lower(3))*m_data_box.numberOfCells(2)+
         a_i2-m_data_box.lower(2))*m_data_box.numberOfCells(1)+
         a_i1-m_data_box.lower(1))*m_data_box.numberOfCells(0)+
         a_i0-m_data_box.lower(0);
   }

   // Transfers data from m_data to a send buffer given a send box for a 1D
   // array.
   void
   transferToSendBuffer1D(
      const Communication& a_comm,
      double* a_send_buffer) const;

   // Transfers data from m_data to a send buffer given a send box for a 2D
   // array.
   void
   transferToSendBuffer2D(
      const Communication& a_comm,
      double* a_send_buffer) const;

   // Transfers data from m_data to a send buffer given a send box for a 3D
   // array.
   void
   transferToSendBuffer3D(
      const Communication& a_comm,
      double* a_send_buffer) const;

   // Transfers data from m_data to a send buffer given a send box for a 4D
   // array.
   void
   transferToSendBuffer4D(
      const Communication& a_comm,
      double* a_send_buffer) const;

   // Transfers data from m_data to a send buffer given a send box for a 5D
   // array.
   void
   transferToSendBuffer5D(
      const Communication& a_comm,
      double* a_send_buffer) const;

   // Transfers data from m_data to a send buffer given a send box for a 6D
   // array.
   void
   transferToSendBuffer6D(
      const Communication& a_comm,
      double* a_send_buffer) const;

   // Transfers data from a receive buffer to m_data given a receive box for a
   // 1D array.
   void
   transferFromRecvBuffer1D(
      const Communication& a_comm,
      const double* a_recv_buffer);

   // Transfers data from a receive buffer to m_data given a receive box for a
   // 2D array.
   void
   transferFromRecvBuffer2D(
      const Communication& a_comm,
      const double* a_recv_buffer);

   // Transfers data from a receive buffer to m_data given a receive box for a
   // 3D array.
   void
   transferFromRecvBuffer3D(
      const Communication& a_comm,
      const double* a_recv_buffer);

   // Transfers data from a receive buffer to m_data given a receive box for a
   // 4D array.
   void
   transferFromRecvBuffer4D(
      const Communication& a_comm,
      const double* a_recv_buffer);

   // Transfers data from a receive buffer to m_data given a receive box for a
   // 5D array.
   void
   transferFromRecvBuffer5D(
      const Communication& a_comm,
      const double* a_recv_buffer);

   // Transfers data from a receive buffer to m_data given a receive box for a
   // 6D array.
   void
   transferFromRecvBuffer6D(
      const Communication& a_comm,
      const double* a_recv_buffer);

   // Transfers data from a_send_box to a_recv_box for a 1D array.
   void
   transferLocalData1D(
      const Box& a_send_box,
      const Box& a_recv_box);

   // Transfers data from a_send_box to a_recv_box for a 2D array.
   void
   transferLocalData2D(
      const Box& a_send_box,
      const Box& a_recv_box);

   // Transfers data from a_send_box to a_recv_box for a 3D array.
   void
   transferLocalData3D(
      const Box& a_send_box,
      const Box& a_recv_box);

   // Transfers data from a_send_box to a_recv_box for a 4D array.
   void
   transferLocalData4D(
      const Box& a_send_box,
      const Box& a_recv_box);

   // Transfers data from a_send_box to a_recv_box for a 5D array.
   void
   transferLocalData5D(
      const Box& a_send_box,
      const Box& a_recv_box);

   // Transfers data from a_send_box to a_recv_box for a 6D array.
   void
   transferLocalData6D(
      const Box& a_send_box,
      const Box& a_recv_box);

   // The dimension of this array.
   int m_dim;

   // Dimensions 0 through m_num_dist_dim-1 are distributed.
   int m_num_dist_dim;

   // The global processor ID.
   int m_my_id;

   // Lowest processor rank array is distributed over.
   int m_proc_lo;

   // Highest processor rank array is distributed over.
   int m_proc_hi;

   // Box describing the extent of the interior of this part of the global
   // parallel array.
   Box m_interior_box;

   // Box describing the extent of the interior and any boundary cells of this
   // part of the global parallel array.
   Box m_local_box;

   // Box describing the extent of the interior with any ghost zones and
   // boundary zones of this part of the global parallel array.
   Box m_data_box;

   // The Communication objects for all ghost and periodic boundary condition
   // communication.
   vector<Communication> m_comms;

   // The number of each type of communication.
   vector<int> m_comm_type_counts;

   // The array of data including interior, ghost, and boundary zones.
   double* m_data;

   // If true problem is periodic in that dimension.
   deque<bool> m_periodic;

   // If true the array has been partitioned.
   bool m_partitioned;

   // The number of partitions in each dimension.
   vector<int> m_dim_partitions;

   // The number of cells in the partitions to the left (lower) end of each
   // dimension.  These partitions contain an extra cell than do the other,
   // upper, partitions.
   vector<int> m_left_partitions_size;

   // The number of partitions to the left (lower) end of each dimension.
   vector<int> m_num_left_partitions;

   // The number of cells in each dimension over with the array is partitioned.
   vector<int> m_num_cells;

   // The number of ghosts in each dimension.
   int m_nghosts;

   // Base MPI tag for communication.
   // This is set to 3645.  The message tag for each point to point
   // communication is this base value + a neighbor ID.  There are
   // pow(3, m_num_dist_dim) neighbors.  Hence the message tags generated by
   // this class are in the range [s_TAG_BASE+1, s_TAG_BASE+pow(3, m_num_dist_dim)].
   // Applications should be mindful of the use of message tags in this range.
   static const int s_TAG_BASE;
};

#endif
