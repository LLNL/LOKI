/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "ParallelArray.H"
#include <mpi.h>
#include <map>

const int ParallelArray::s_TAG_BASE = 3645;

ParallelArray::ParallelArray()
   : m_dim(-1),
     m_num_dist_dim(-1),
     m_proc_lo(-1),
     m_proc_hi(-1),
     m_interior_box(),
     m_local_box(),
     m_data_box(),
     m_comm_type_counts(NUM_COMM_TYPES, 0),
     m_data(0),
     m_partitioned(false)
{
   MPI_Comm_rank(MPI_COMM_WORLD, &m_my_id);
}

ParallelArray::ParallelArray(
   int a_dim,
   int a_num_dist_dim)
   : m_dim(a_dim),
     m_num_dist_dim(a_num_dist_dim),
     m_proc_lo(-1),
     m_proc_hi(-1),
     m_interior_box(m_dim),
     m_local_box(m_dim),
     m_data_box(m_dim),
     m_comm_type_counts(NUM_COMM_TYPES, 0),
     m_data(0),
     m_periodic(m_num_dist_dim, false),
     m_partitioned(false),
     m_dim_partitions(a_num_dist_dim, 1)
{
   ARRAY_ASSERT(m_dim >= 1 && m_dim <= 6);
   ARRAY_ASSERT(m_num_dist_dim > 0 && m_num_dist_dim <= m_dim);
   MPI_Comm_rank(MPI_COMM_WORLD, &m_my_id);
}


ParallelArray::ParallelArray(
   int a_dim,
   int a_num_dist_dim,
   int a_rank_lo,
   int a_rank_hi,
   int a_nghosts,
   const deque<bool>& a_is_periodic,
   const vector<int>& a_num_cells)
   : m_dim(a_dim),
     m_num_dist_dim(a_num_dist_dim),
     m_proc_lo(-1),
     m_proc_hi(-1),
     m_interior_box(m_dim),
     m_local_box(m_dim),
     m_data_box(m_dim),
     m_comm_type_counts(NUM_COMM_TYPES, 0),
     m_data(0),
     m_periodic(m_num_dist_dim, false),
     m_partitioned(false),
     m_dim_partitions(a_num_dist_dim, 1)
{
   ARRAY_ASSERT(m_dim >= 1 && m_dim <= 6);
   ARRAY_ASSERT(m_num_dist_dim > 0 && m_num_dist_dim <= m_dim);
   MPI_Comm_rank(MPI_COMM_WORLD, &m_my_id);
   partition(a_rank_lo,
      a_rank_hi,
      a_nghosts,
      a_is_periodic,
      a_num_cells);
}


ParallelArray::ParallelArray(
   const Box& a_base_space,
   int a_num_dist_dim,
   int a_nghosts,
   const vector<int>& a_num_global_cells)
   : m_dim(a_base_space.dim()),
     m_num_dist_dim(a_num_dist_dim),
     m_proc_lo(-1),
     m_proc_hi(-1),
     m_interior_box(m_dim),
     m_local_box(m_dim),
     m_data_box(m_dim),
     m_comm_type_counts(NUM_COMM_TYPES, 0),
     m_data(0),
     m_periodic(m_num_dist_dim, false),
     m_partitioned(false),
     m_dim_partitions(a_num_dist_dim, 1)
{
   ARRAY_ASSERT(m_dim >= 1 && m_dim <= 6);
   ARRAY_ASSERT(m_num_dist_dim > 0 && m_num_dist_dim <= m_dim);
   MPI_Comm_rank(MPI_COMM_WORLD, &m_my_id);
   partition(a_base_space, a_nghosts, a_num_global_cells);
}


ParallelArray::ParallelArray(
   const ParallelArray& a_other)
   : m_dim(a_other.m_dim),
     m_num_dist_dim(a_other.m_num_dist_dim),
     m_my_id(a_other.m_my_id),
     m_proc_lo(a_other.m_proc_lo),
     m_proc_hi(a_other.m_proc_hi),
     m_interior_box(a_other.m_interior_box),
     m_local_box(a_other.m_local_box),
     m_data_box(a_other.m_data_box),
     m_comms(a_other.m_comms),
     m_comm_type_counts(a_other.m_comm_type_counts),
     m_periodic(a_other.m_periodic),
     m_partitioned(a_other.m_partitioned),
     m_dim_partitions(a_other.m_dim_partitions),
     m_left_partitions_size(a_other.m_left_partitions_size),
     m_num_left_partitions(a_other.m_num_left_partitions),
     m_num_cells(a_other.m_num_cells),
     m_nghosts(a_other.m_nghosts)
{
   // If a_other is partitioned then the data array needs to be constructed
   // for this and it must be populated with a_other's data.
   if (m_partitioned) {
      int data_size = m_data_box.size();
      m_data = new double [data_size];
      memcpy(m_data, a_other.m_data, data_size*sizeof(double));
   }
   else {
      m_data = 0;
   }
}


ParallelArray&
ParallelArray::operator = (
   const ParallelArray& a_rhs)
{
   // Delete/allocate/reallocate m_data up front while we know the states of
   // both this and the rhs.
   if (m_partitioned) {
      if (!a_rhs.m_partitioned) {
         delete [] m_data;
         m_data = 0;
      }
      else {
         int rhs_data_size = a_rhs.m_data_box.size();
         if (m_data_box.size() != rhs_data_size) {
            delete [] m_data;
            m_data = new double [rhs_data_size];
         }
         memcpy(m_data, a_rhs.m_data, rhs_data_size*sizeof(double));
      }
   }
   else if (a_rhs.m_partitioned) {
      int rhs_data_size = a_rhs.m_data_box.size();
      m_data = new double [rhs_data_size];
      memcpy(m_data, a_rhs.m_data, rhs_data_size*sizeof(double));
   }
   else {
      m_data = 0;
   }
   m_dim = a_rhs.m_dim;
   m_num_dist_dim = a_rhs.m_num_dist_dim;
   m_proc_lo = a_rhs.m_proc_lo;
   m_proc_hi = a_rhs.m_proc_hi;
   m_interior_box = a_rhs.m_interior_box;
   m_local_box = a_rhs.m_local_box;
   m_data_box = a_rhs.m_data_box;
   m_comms = a_rhs.m_comms;
   m_comm_type_counts = a_rhs.m_comm_type_counts;
   m_periodic = a_rhs.m_periodic;
   m_partitioned = a_rhs.m_partitioned;
   m_dim_partitions = a_rhs.m_dim_partitions;
   m_left_partitions_size = a_rhs.m_left_partitions_size;
   m_num_left_partitions = a_rhs.m_num_left_partitions;
   m_num_cells = a_rhs.m_num_cells;
   m_nghosts = a_rhs.m_nghosts;
}


ParallelArray::~ParallelArray()
{
   if (m_data != 0) {
      delete [] m_data;
   }
}


void
ParallelArray::partition(
   int a_dim,
   int a_num_dist_dim,
   int a_rank_lo,
   int a_rank_hi,
   int a_nghosts,
   const deque<bool>& a_is_periodic,
   const vector<int>& a_num_cells)
{
   ARRAY_ASSERT(m_dim == -1);
   ARRAY_ASSERT(!m_partitioned);
   ARRAY_ASSERT(a_dim >= 1 && a_dim <= 6);
   ARRAY_ASSERT(a_num_dist_dim > 0 && a_num_dist_dim <= a_dim);

   // Set the dimension and number of distributed dimensions.
   m_dim = a_dim;
   m_num_dist_dim = a_num_dist_dim;

   // Dimension the interior and data boxes.
   m_interior_box.redimension(m_dim);
   m_local_box.redimension(m_dim);
   m_data_box.redimension(m_dim);

   // Initialize the periodicity.
   m_periodic.resize(m_num_dist_dim, false);
   m_dim_partitions.resize(m_num_dist_dim, 1);

   // Now the array looks like it came from a non-default constructor so call
   // the partition method appropriate for that constructor.
   partition(a_rank_lo,
      a_rank_hi,
      a_nghosts,
      a_is_periodic,
      a_num_cells);
}


void
ParallelArray::partition(
   int a_rank_lo,
   int a_rank_hi,
   int a_nghosts,
   const deque<bool>& a_is_periodic,
   const vector<int>& a_num_cells)
{
   ARRAY_ASSERT(m_dim != -1);
   ARRAY_ASSERT(!m_partitioned);
   ARRAY_ASSERT(a_rank_hi >= a_rank_lo);
   ARRAY_ASSERT(a_is_periodic.size() == m_num_dist_dim);
   ARRAY_ASSERT(a_num_cells.size() == m_dim);

   // Set the low and high processors.
   m_proc_lo = a_rank_lo;
   m_proc_hi = a_rank_hi;

   // Set the number of ghost and number of cells in each dimension.
   m_nghosts = a_nghosts;
   m_num_cells = a_num_cells;

   // Store periodicity information.
   for (int dim = 0; dim < m_num_dist_dim; ++dim) {
      m_periodic[dim] = a_is_periodic[dim];
   }

   // Distribute the processors among the different dimensions.
   assignProcs();

   // If the array is not partitioned onto this processor make the interior and
   // data boxes empty to indicate this and partitioning is complete.
   if (m_my_id < m_proc_lo || m_my_id > m_proc_hi) {
      m_left_partitions_size.resize(m_dim);
      m_num_left_partitions.resize(m_dim);
      for (int i = 0; i < m_dim; ++i) {
         m_interior_box.lower(i) = 0;
         m_interior_box.upper(i) = -1;
         m_data_box.lower(i) = 0;
         m_data_box.upper(i) = -1;
         m_local_box.lower(i) = 0;
         m_local_box.upper(i) = -1;
         m_left_partitions_size[i] = 0;
         m_num_left_partitions[i] = 0;
      }
      m_partitioned = true;
      return;
   }

   // The array is paritioned onto this processor so allocate the necessary
   // number of communication neighbors, determine how many processors partition
   // each dimension, determine which part of the array is owned by each
   // processor, and the data send/received by each communication neighbor.
   m_comms.resize(static_cast<int>(pow(3, m_num_dist_dim)),
      Communication(m_dim));
   setupLocalDomain();
   m_data = new double [m_data_box.size()];
   m_partitioned = true;
   operator=(0.0);
}


void
ParallelArray::partition(
   const Box& a_base_space,
   int a_num_dist_dim,
   int a_nghosts,
   const vector<int>& a_num_global_cells)
{
   ARRAY_ASSERT(m_dim == -1);
   ARRAY_ASSERT(!m_partitioned);
   ARRAY_ASSERT(a_base_space.dim() >= 1 && a_base_space.dim() <= 6);
   ARRAY_ASSERT(a_num_dist_dim > 0 && a_num_dist_dim <= a_base_space.dim());

   // Set the dimension and number of distributed dimensions.
   m_dim = a_base_space.dim();
   m_num_dist_dim = a_num_dist_dim;

   // Dimension the interior and data boxes.
   m_interior_box.redimension(m_dim);
   m_local_box.redimension(m_dim);
   m_data_box.redimension(m_dim);

   // Initialize the periodicity.
   m_periodic.resize(m_num_dist_dim, false);
   m_dim_partitions.resize(m_num_dist_dim, 1);

   // Now the array looks like it came from a non-default constructor so call
   // the partition method appropriate for that constructor.
   partition(a_base_space, a_nghosts, a_num_global_cells);
}


void
ParallelArray::partition(
   const Box& a_base_space,
   int a_nghosts,
   const vector<int>& a_num_global_cells)
{
   ARRAY_ASSERT(m_dim != -1);
   ARRAY_ASSERT(!m_partitioned);
   ARRAY_ASSERT(a_base_space.dim() == m_dim);
   ARRAY_ASSERT(a_num_global_cells.size() == m_num_dist_dim);

   // Set the low and high processors.
   m_proc_lo = m_my_id;
   m_proc_hi = m_my_id;

   // Set the number of ghosts.
   m_nghosts = a_nghosts;

   // Store interior box and data box information.  Also store the number of
   // cells in each dimension.  If the base space box is empty then this array
   // is also empty.
   m_left_partitions_size.resize(m_dim);
   m_num_left_partitions.resize(m_dim);
   if (a_base_space.empty()) {
      for (int dim = 0; dim < m_dim; ++dim) {
         m_interior_box.lower(dim) = 0;
         m_interior_box.upper(dim) = -1;
         m_data_box.lower(dim) = 0;
         m_data_box.upper(dim) = -1;
         m_local_box.lower(dim) = 0;
         m_local_box.upper(dim) = -1;
         m_num_cells.push_back(0);
         m_left_partitions_size[dim] = 0;
         m_num_left_partitions[dim] = 0;
      }
   }
   else {
      for (int dim = 0; dim < m_num_dist_dim; ++dim) {
         m_interior_box.lower(dim) = a_base_space.lower(dim);
         m_interior_box.upper(dim) = a_base_space.upper(dim);
         m_data_box.lower(dim) = m_interior_box.lower(dim) - m_nghosts;
         m_data_box.upper(dim) = m_interior_box.upper(dim) + m_nghosts;
         if (m_data_box.lower(dim) < 0) {
            m_local_box.lower(dim) = m_data_box.lower(dim);
         }
         else {
            m_local_box.lower(dim) = m_data_box.lower(dim) + m_nghosts;
         }
         if (m_data_box.upper(dim) >= a_num_global_cells[dim]) {
            m_local_box.upper(dim) = m_data_box.upper(dim);
         }
         else {
            m_local_box.upper(dim) = m_data_box.upper(dim) - m_nghosts;
         }
         m_num_cells.push_back(a_base_space.numberOfCells(dim));
         m_left_partitions_size[dim] = m_interior_box.numberOfCells(dim);
         m_num_left_partitions[dim] = 1;
      }
      for (int dim = m_num_dist_dim; dim < m_dim; ++dim) {
         m_interior_box.lower(dim) = a_base_space.lower(dim);
         m_interior_box.upper(dim) = a_base_space.upper(dim);
         m_data_box.lower(dim) = m_interior_box.lower(dim);
         m_data_box.upper(dim) = m_interior_box.upper(dim);
         m_local_box.lower(dim) = m_data_box.lower(dim);
         m_local_box.upper(dim) = m_data_box.upper(dim);
         m_num_cells.push_back(a_base_space.numberOfCells(dim));
         m_left_partitions_size[dim] = m_interior_box.numberOfCells(dim);
         m_num_left_partitions[dim] = 1;
      }
   }

   // Allocate the data.
   m_data = new double [m_data_box.size()];
   m_partitioned = true;
   operator=(0.0);
}


int
ParallelArray::getGlobalRank(
   const vector<int> a_idx_rank) const
{
   // This function may only be called during partitioning after the array has
   // been properly initialized or after partitioning.
   ARRAY_ASSERT(m_num_dist_dim != -1);
   int rank = a_idx_rank[0];
   for (int i = 1; i < m_num_dist_dim; ++i) {
      rank = rank*m_dim_partitions[i] + a_idx_rank[i];
   }
   rank += m_proc_lo;
   return rank;
}


ParallelArray::Box
ParallelArray::interiorBox(
   int a_rank) const
{
   ARRAY_ASSERT(m_partitioned);
   if (a_rank < m_proc_lo || a_rank > m_proc_hi) {
      Box the_box(m_dim);
      for (int dim = 0; dim < m_dim; ++dim) {
         the_box.lower(dim) = 0;
         the_box.upper(dim) = -1;
      }
      return the_box;
   }
   else if (a_rank == m_my_id) {
     return m_interior_box;
   }
   else {
      Box the_box(m_dim);
      // Get the index rank of the processor of interest.
      vector<int> rank_idx(m_num_dist_dim);
      getIndexRank(a_rank, m_proc_lo, m_dim_partitions, rank_idx);

      // Assign interior grid size.  To maintain load balance, we want each
      // partition of a dimension to have as close to equal number of zones as
      // possible.  We distribute any extra zones in a dimension among the the
      // lower partitions of that dimension until all extras have been
      // distributed.
      for (int dim = 0; dim < m_num_dist_dim; ++dim) {
         int Nloc = m_num_cells[dim]/m_dim_partitions[dim];
         int Nloc_extra = m_num_cells[dim]%m_dim_partitions[dim];
         int lower, upper;
         if (rank_idx[dim] < Nloc_extra) {
            ++Nloc;
            lower = rank_idx[dim]*Nloc;
         }
         else {
            lower = rank_idx[dim]*Nloc+Nloc_extra;
         }
         upper = Nloc+lower-1;

         the_box.lower(dim) = lower;
         the_box.upper(dim) = upper;
      }
      for (int dim = m_num_dist_dim; dim < m_dim; ++dim) {
         the_box.lower(dim) = 0;
         the_box.upper(dim) = m_num_cells[dim]-1;
      }
      return the_box;
   }
}


ParallelArray::Box
ParallelArray::localBox(
   int a_rank) const
{
   Box the_box = dataBox(a_rank);
   for (int i = 0; i < m_num_dist_dim; ++i) {
      if (the_box.lower(i) >= 0) {
         the_box.lower(i) = the_box.lower(i) + m_nghosts;
      }
      if (the_box.upper(i) < m_num_cells[i]) {
         the_box.upper(i) = the_box.upper(i) - m_nghosts;
      }
   }
   return the_box;
}


ParallelArray::Box
ParallelArray::dataBox(
   int a_rank) const
{
   ARRAY_ASSERT(m_partitioned);
   if (a_rank < m_proc_lo || a_rank > m_proc_hi) {
      Box the_box(m_dim);
      for (int dim = 0; dim < m_dim; ++dim) {
         the_box.lower(dim) = 0;
         the_box.upper(dim) = -1;
      }
      return the_box;
   }
   else if (a_rank == m_my_id) {
     return m_data_box;
   }
   else {
      Box the_box(m_dim);
      // Get the index rank of the processor of interest.
      vector<int> rank_idx(m_num_dist_dim);
      getIndexRank(a_rank, m_proc_lo, m_dim_partitions, rank_idx);

      // Assign interior grid size.  To maintain load balance, we want each
      // partition of a dimension to have as close to equal number of zones as
      // possible.  We distribute any extra zones in a dimension among the the
      // lower partitions of that dimension until all extras have been
      // distributed.
      for (int dim = 0; dim < m_num_dist_dim; ++dim) {
         int Nloc = m_num_cells[dim]/m_dim_partitions[dim];
         int Nloc_extra = m_num_cells[dim]%m_dim_partitions[dim];
         int lower, upper;
         if (rank_idx[dim] < Nloc_extra) {
            ++Nloc;
            lower = rank_idx[dim]*Nloc;
         }
         else {
            lower = rank_idx[dim]*Nloc+Nloc_extra;
         }
         upper = Nloc+lower-1;

         the_box.lower(dim) = lower - m_nghosts;
         the_box.upper(dim) = upper + m_nghosts;
      }
      for (int dim = m_num_dist_dim; dim < m_dim; ++dim) {
         the_box.lower(dim) = 0;
         the_box.upper(dim) = m_num_cells[dim]-1;
      }
      return the_box;
   }
}


// Private methods.
void
ParallelArray::assignProcs()
{
   // Distribute the available processors among the distributed dimensions
   // trying to make the "squarest" domains possible.
   int num_ranks = m_proc_hi - m_proc_lo + 1;
   for (int used = 1; used < num_ranks;) {
      int bigRatioDim = 0;
      float ratio =
         ((float)m_num_cells[0])/((float)m_dim_partitions[0]);
      for (int i = 1; i < m_num_dist_dim; ++i)   {
         float nRatio =
            ((float)m_num_cells[i])/((float)m_dim_partitions[i]);
         if (nRatio > ratio ||
             (nRatio == ratio && nRatio > 1 &&
              m_dim_partitions[i] < m_dim_partitions[bigRatioDim])) {
            ratio = nRatio;
            bigRatioDim = i;
         }
      }
      int factor = greatestPrimeFactor(num_ranks/used);
      m_dim_partitions[bigRatioDim] *= factor;
      used *= factor;
   }
}


int
ParallelArray::leastPrimeFactor(
   int a_n) const
{
   if (a_n < 2) {
      return a_n;
   }
   if (a_n % 2 == 0) {
      return 2;
   }

   int sqrtn = int(sqrt(double(a_n+0.5)));
   for (int i = 3; i <= sqrtn; i += 2) {
      if (a_n % i == 0) {
         return i;
      }
   }
   return a_n;
}


int
ParallelArray::greatestPrimeFactor(
   int a_n) const
{
   // After dividing out all the smaller prime factors, what's left is the
   // largest one.
   for (int lpf = leastPrimeFactor(a_n);
        lpf != a_n;
        lpf = leastPrimeFactor(a_n)) {
      a_n = a_n/lpf;
   }

   return a_n;
}


void
ParallelArray::getIndexRank(
   int a_global_rank,
   int a_rank_lo,
   const vector<int> a_dim_partitions,
   vector<int>& a_index_rank) const
{
   int rank = a_global_rank - a_rank_lo;
   for (int i = 0; i < m_num_dist_dim; ++i) {
      int factor = 1;
      for (int j = i+1; j < m_num_dist_dim; ++j) {
         factor *= a_dim_partitions[j];
      }
      a_index_rank[i] = rank/factor;
      rank -= a_index_rank[i]*factor;
   }
}


void
ParallelArray::setupLocalDomain()
{
   // Get the index rank of this processor.
   vector<int> rank_idx(m_num_dist_dim);
   getIndexRank(m_my_id, m_proc_lo, m_dim_partitions, rank_idx);

   // Assign interior grid size.  To maintain load balance, we want each
   // partition of a dimension to have as close to equal number of zones as
   // possible.  We distribute any extra zones in a dimension among the the
   // lower partitions of that dimension until all extras have been distributed.
   m_left_partitions_size.resize(m_dim);
   m_num_left_partitions.resize(m_dim);
   for (int dim = 0; dim < m_num_dist_dim; ++dim) {
      int Nloc = m_num_cells[dim]/m_dim_partitions[dim];
      int Nloc_extra = m_num_cells[dim]%m_dim_partitions[dim];
      if (Nloc_extra == 0) {
         m_left_partitions_size[dim] = Nloc;
         m_num_left_partitions[dim] = m_dim_partitions[dim];
      }
      else {
         m_left_partitions_size[dim] = Nloc + 1;
         m_num_left_partitions[dim] = Nloc_extra;
      }
      int lower, upper;
      if (rank_idx[dim] < Nloc_extra) {
         ++Nloc;
         lower = rank_idx[dim]*Nloc;
      }
      else {
         lower = rank_idx[dim]*Nloc+Nloc_extra;
      }
      upper = Nloc+lower-1;

      m_interior_box.lower(dim) = lower;
      m_interior_box.upper(dim) = upper;
      m_data_box.lower(dim) = lower - m_nghosts;
      m_data_box.upper(dim) = upper + m_nghosts;
      if (m_data_box.lower(dim) < 0) {
         m_local_box.lower(dim) = m_data_box.lower(dim);
      }
      else {
         m_local_box.lower(dim) = m_data_box.lower(dim) + m_nghosts;
      }
      if (m_data_box.upper(dim) >= m_num_cells[dim]) {
         m_local_box.upper(dim) = m_data_box.upper(dim);
      }
      else {
         m_local_box.upper(dim) = m_data_box.upper(dim) - m_nghosts;
      }
   }
   for (int dim = m_num_dist_dim; dim < m_dim; ++dim) {
      int lower = 0, upper = m_num_cells[dim]-1;
      m_interior_box.lower(dim) = lower;
      m_interior_box.upper(dim) = upper;
      m_data_box.lower(dim) = lower;
      m_data_box.upper(dim) = upper;
      m_local_box.lower(dim) = lower;
      m_local_box.upper(dim) = upper;
   }

   // Now set up the communication objects to describe the information sent from
   // and received by this process' parallel array and the other processes that
   // the parallel array is distributed among.

   // Define a mesh with 3 domains in each direction.  These are all the
   // communication neighbors of the domain whose communication information we
   // are trying define.  The domain whose communication information we're
   // defining sits at the center of the mesh.
   const vector<int> comm_procs_dim(m_num_dist_dim, 3);
   vector<int> commNbr_rank_idx(m_num_dist_dim);
   for (int commNbr = 0;
        commNbr < static_cast<int>(pow(3, m_num_dist_dim));
        ++commNbr) {
      // Get the index rank of this communication neighbor relative to the mesh
      // of communication neighbors.  Since the domain of interest is at the
      // index rank containing all 1's the shift of this communication neighbor
      // relative to the domain of interest, idx_shift, is its index rank - 1.
      // We then use idx_shift to convert index rank into an index rank relative
      // to the mesh of all the domains this array is distributed over.
      getIndexRank(commNbr, 0, comm_procs_dim, commNbr_rank_idx);
      vector<int> idx_shift(m_num_dist_dim);
      for (int i = 0; i < m_num_dist_dim; ++i) {
         idx_shift[i] = commNbr_rank_idx[i] - 1;
         commNbr_rank_idx[i] = rank_idx[i] + idx_shift[i];
         // The communication neighbor may be a periodic domain which at this
         // point is represented by an index which is -1 or m_dim_partitions[i].
         // So we need to convert these indices into the proper value for that
         // periodic domain.  The modulus operation
         // commNbr_rank_idx[i]%m_dim_partitions[i] will do this when
         // commNbr_rank_idx[i] != -1.  The modulus operator is implementation
         // depdendent for negative dividends.  So we need to handle the -1
         // index case separately.
         if (commNbr_rank_idx[i] == -1) {
            commNbr_rank_idx[i] += m_dim_partitions[i];
         }
         else {
            commNbr_rank_idx[i] = commNbr_rank_idx[i]%m_dim_partitions[i];
         }
      }

      // Now convert the communication neighbor's valid index rank into a
      // global rank.
      int commNbr_rank = getGlobalRank(commNbr_rank_idx);

      // Check to see if this communication neighbor is in any distributed
      // dimension's periodic boundary.  If it's in more than 1 periodic
      // boundary then it's in some kind of "corner" periodic boundary.
      deque<bool> in_periodic_bdy(m_num_dist_dim, false);
      bool found_periodic_bdy = false;
      bool corner_periodic_bdy = false;
      for (int i = 0; i < m_num_dist_dim; ++i) {
         if ((rank_idx[i] == 0 && idx_shift[i] < 0) ||
             (rank_idx[i] == m_dim_partitions[i] - 1 && idx_shift[i] > 0)) {
            in_periodic_bdy[i] = true;
            if (found_periodic_bdy) {
               corner_periodic_bdy = true;
            }
            else {
               found_periodic_bdy = true;
            }
         }
      }

      // A process "communicates" with itself only to enforce periodicity.
      // Also, communication of "corner" periodic boundaries are not explicitly
      // executed.  In these cases set the Communication object with an invalid
      // partner rank and SELF_COMM/CORNER_PER_COMM communication type.
      // Otherwise fully define the Communication object.
      if (commNbr_rank == m_my_id && !found_periodic_bdy) {
         m_comms[commNbr].m_partner_rank = -1;
         m_comms[commNbr].m_comm_type = SELF_COMM;
         ++m_comm_type_counts[SELF_COMM];
      }
      else if (corner_periodic_bdy) {
         m_comms[commNbr].m_partner_rank = -1;
         m_comms[commNbr].m_comm_type = CORNER_PER_COMM;
         ++m_comm_type_counts[CORNER_PER_COMM];
      }
      else {
         m_comms[commNbr].m_partner_rank = commNbr_rank;

         // Assume that this is a communication of ghost data.  If it turns out
         // to actually be communication of a periodic dimension this will be
         // overridden.
         Comm_Type comm_type = GHOST_COMM;

         // Set the MPI tags for the send and receive boxes.
         m_comms[commNbr].m_recv_tag = getRecvTag(idx_shift);
         m_comms[commNbr].m_send_tag = getSendTag(idx_shift);

         // For each distributed dimension of this communication neighbor check
         // if it is actually a periodic boundary communication neighbor and
         // set up the send/recv boxes for the communication.
         for (int i = 0; i < m_num_dist_dim; ++i) {
            // If this communication neighbor is in a periodic boundary and
            // periodicity is to be applied then the send/recv boxes of the next
            // lower dimension in which periodicity is to be applied may need to
            // be adjusted.  This is necessary so that "corner" or "edge"
            // periodic boundaries are communicated by sweeps of periodic faces.
            if (m_periodic[i] && in_periodic_bdy[i]) {
               // Look at the lower dimensions.
               for (int j = i-1; j >= 0; --j) {
                  // If this dimension is periodic and this communication
                  // neighbor is at the upper or lower end of this dimension
                  // then we may need to adjust the send/recv boxes.
                  if (m_periodic[j] && idx_shift[j] == 0 &&
                      (rank_idx[j] == 0 ||
                       rank_idx[j] == m_dim_partitions[j] - 1)) {
                     // If all dimensions of this communication neighbor other
                     // than the pair of periodic dimensions we're dealing with
                     // (i and j) are not in their periodic boundaries then the
                     // send/recv boxes need to be adjusted 
                     bool extend_dim = true;
                     for (int k = 0; k < m_num_dist_dim; ++k) {
                        if (k == i || k == j) {
                           continue;
                        }
                        if (in_periodic_bdy[k]) {
                           extend_dim = false;
                           break;
                        }
                     }
                     // Adjust the send/recv boxes of this next lower periodic
                     // dimension for this communication neighbor.
                     if (extend_dim) {
                        if (rank_idx[j] == 0) {
                           m_comms[commNbr].m_recv_box.lower(j) -= m_nghosts;
                           m_comms[commNbr].m_send_box.lower(j) -= m_nghosts;
                        }
                        if (rank_idx[j] == m_dim_partitions[j] - 1) {
                           m_comms[commNbr].m_recv_box.upper(j) += m_nghosts;
                           m_comms[commNbr].m_send_box.upper(j) += m_nghosts;
                        }
                     }
                  }
               }
            }
            // Define the send and receive boxes for this dimension.
            if (idx_shift[i] < 0) {
               m_comms[commNbr].m_recv_box.lower(i) =
                  m_interior_box.lower(i) - m_nghosts;
               m_comms[commNbr].m_recv_box.upper(i) =
                  m_interior_box.lower(i) - 1;
               m_comms[commNbr].m_send_box.lower(i) =
                  m_interior_box.lower(i);
               m_comms[commNbr].m_send_box.upper(i) =
                  m_interior_box.lower(i) + m_nghosts - 1;
            }
            else if (idx_shift[i] > 0) {
               m_comms[commNbr].m_recv_box.lower(i) =
                  m_interior_box.upper(i) + 1;
               m_comms[commNbr].m_recv_box.upper(i) =
                  m_interior_box.upper(i) + m_nghosts;
               m_comms[commNbr].m_send_box.lower(i) =
                  m_interior_box.upper(i) - m_nghosts + 1;
               m_comms[commNbr].m_send_box.upper(i) =
                  m_interior_box.upper(i);
            }
            else {
               m_comms[commNbr].m_recv_box.lower(i) = m_interior_box.lower(i);
               m_comms[commNbr].m_recv_box.upper(i) = m_interior_box.upper(i);
               m_comms[commNbr].m_send_box.lower(i) = m_interior_box.lower(i);
               m_comms[commNbr].m_send_box.upper(i) = m_interior_box.upper(i);
            }
            // See if this is some form of periodic communication.  We will set
            // the communication type below.  We've already handled any "corner"
            // periodic communication cases earlier so only one of these
            // periodic communications is now possible.
            if (in_periodic_bdy[i]) {
               if (i == 0) {
                  if (commNbr_rank == m_my_id) {
                     comm_type = X1_PER_COMM_LOCAL;
                  }
                  else {
                     comm_type = X1_PER_COMM;
                  }
               }
               else if (i == 1) {
                  if (commNbr_rank == m_my_id) {
                     comm_type = X2_PER_COMM_LOCAL;
                  }
                  else {
                     comm_type = X2_PER_COMM;
                  }
               }
               else if (i == 2) {
                  if (commNbr_rank == m_my_id) {
                     comm_type = X3_PER_COMM_LOCAL;
                  }
                  else {
                     comm_type = X3_PER_COMM;
                  }
               }
               else if (i == 3) {
                  if (commNbr_rank == m_my_id) {
                     comm_type = X4_PER_COMM_LOCAL;
                  }
                  else {
                     comm_type = X4_PER_COMM;
                  }
               }
               else if (i == 4) {
                  if (commNbr_rank == m_my_id) {
                     comm_type = X5_PER_COMM_LOCAL;
                  }
                  else {
                     comm_type = X5_PER_COMM;
                  }
               }
               else if (i == 5) {
                  if (commNbr_rank == m_my_id) {
                     comm_type = X6_PER_COMM_LOCAL;
                  }
                  else {
                     comm_type = X6_PER_COMM;
                  }
               }
            }
         }
         // Now deal with the undistributed dimensions.  Basically, all we
         // need to do is to set the send/recv box for these dimensions.
         for (int i = m_num_dist_dim; i < m_dim; ++i) {
            m_comms[commNbr].m_recv_box.lower(i) = m_interior_box.lower(i);
            m_comms[commNbr].m_recv_box.upper(i) = m_interior_box.upper(i);
            m_comms[commNbr].m_send_box.lower(i) = m_interior_box.lower(i);
            m_comms[commNbr].m_send_box.upper(i) = m_interior_box.upper(i);
         }
         // Set up the communication type based on which if any periodic
         // boundaries this communication neighbor falls on.
         m_comms[commNbr].m_comm_type = comm_type;
         ++m_comm_type_counts[comm_type];
      }
   }
}


void
ParallelArray::communicateData(
   Comm_Type a_parallel_comm_type,
   Comm_Type a_local_comm_type)
{
   int num_parallel_comms = m_comm_type_counts[a_parallel_comm_type];
   int num_local_comms = m_comm_type_counts[a_local_comm_type];

   // If there's no communication of these types we're done.
   if (num_parallel_comms == 0 && num_local_comms == 0) {
      return;
   }

   int num_comms = static_cast<int>(m_comms.size());

   // First post the receives and sends for the parallel data communication.
   // Post all the receives.
   vector<double*> recvBuffers(num_parallel_comms);
   vector<MPI_Request> recvReqs(num_parallel_comms);
   map<int, int> recv_idx_to_comm;
   int which_comm = 0;
   if (num_parallel_comms != 0) {
      for (int i = 0; i < num_comms; ++i) {
         if (m_comms[i].m_comm_type == a_parallel_comm_type) {
            // Allocate the receive buffer for this message.
            int thisRecvBufferSize = m_comms[i].m_recv_box.size();
            double* thisRecvBuffer = new double [thisRecvBufferSize];
            recvBuffers[which_comm] = thisRecvBuffer;

            // Post the receive.
            MPI_Irecv(thisRecvBuffer,
               thisRecvBufferSize,
               MPI_DOUBLE,
               m_comms[i].m_partner_rank,
               s_TAG_BASE+m_comms[i].m_recv_tag,
               MPI_COMM_WORLD,
               &recvReqs[which_comm]);

            // Track which communication neighbor this receive message is
            // associated with.
            recv_idx_to_comm[which_comm] = i;
            ++which_comm;
         }
      }
   }

   // Post all the sends.
   vector<double*> sendBuffers(num_parallel_comms);
   vector<MPI_Request> sendReqs(num_parallel_comms);
   which_comm = 0;
   if (num_parallel_comms != 0) {
      for (int i = 0; i < num_comms; ++i) {
         if (m_comms[i].m_comm_type == a_parallel_comm_type) {
            // Allocate and fill the send buffer.
            int thisSendBufferSize = m_comms[i].m_send_box.size();
            double* thisSendBuffer = new double [thisSendBufferSize];
            sendBuffers[which_comm] = thisSendBuffer;
            if (m_dim == 1) {
               transferToSendBuffer1D(m_comms[i], thisSendBuffer);
            }
            else if (m_dim == 2) {
               transferToSendBuffer2D(m_comms[i], thisSendBuffer);
            }
            else if (m_dim == 3) {
               transferToSendBuffer3D(m_comms[i], thisSendBuffer);
            }
            else if (m_dim == 4) {
               transferToSendBuffer4D(m_comms[i], thisSendBuffer);
            }
            else if (m_dim == 5) {
               transferToSendBuffer5D(m_comms[i], thisSendBuffer);
            }
            else {
               transferToSendBuffer6D(m_comms[i], thisSendBuffer);
            }

            // Post the send.
            MPI_Isend(thisSendBuffer,
               thisSendBufferSize,
               MPI_DOUBLE,
               m_comms[i].m_partner_rank,
               s_TAG_BASE+m_comms[i].m_send_tag,
               MPI_COMM_WORLD,
               &sendReqs[which_comm]);
            ++which_comm;
         }
      }
   }

   // While waiting for the dust to settle for the parallel communication,
   // do the local data copies for any local data communication.
   if (num_local_comms != 0) {
      // If there is any local communications then there must be a pair.
      ARRAY_ASSERT(num_local_comms == 2);
      int local_comm_idx[2];
      int which_local_comm = 0;
      for (int i = 0; i < num_comms; ++i) {
         if (m_comms[i].m_comm_type == a_local_comm_type) {
            local_comm_idx[which_local_comm++] = i;
         }
      }
      ARRAY_ASSERT(m_comms[local_comm_idx[0]].m_send_tag ==
                   m_comms[local_comm_idx[1]].m_recv_tag &&
                   m_comms[local_comm_idx[0]].m_recv_tag ==
                   m_comms[local_comm_idx[1]].m_send_tag);
      if (m_dim == 1) {
         transferLocalData1D(m_comms[local_comm_idx[0]].m_send_box,
            m_comms[local_comm_idx[1]].m_recv_box);
         transferLocalData1D(m_comms[local_comm_idx[1]].m_send_box,
            m_comms[local_comm_idx[0]].m_recv_box);
      }
      else if (m_dim == 2) {
         transferLocalData2D(m_comms[local_comm_idx[0]].m_send_box,
            m_comms[local_comm_idx[1]].m_recv_box);
         transferLocalData2D(m_comms[local_comm_idx[1]].m_send_box,
            m_comms[local_comm_idx[0]].m_recv_box);
      }
      else if (m_dim == 3) {
         transferLocalData3D(m_comms[local_comm_idx[0]].m_send_box,
            m_comms[local_comm_idx[1]].m_recv_box);
         transferLocalData3D(m_comms[local_comm_idx[1]].m_send_box,
            m_comms[local_comm_idx[0]].m_recv_box);
      }
      else if (m_dim == 4) {
         transferLocalData4D(m_comms[local_comm_idx[0]].m_send_box,
            m_comms[local_comm_idx[1]].m_recv_box);
         transferLocalData4D(m_comms[local_comm_idx[1]].m_send_box,
            m_comms[local_comm_idx[0]].m_recv_box);
      }
      else if (m_dim == 5) {
         transferLocalData5D(m_comms[local_comm_idx[0]].m_send_box,
            m_comms[local_comm_idx[1]].m_recv_box);
         transferLocalData5D(m_comms[local_comm_idx[1]].m_send_box,
            m_comms[local_comm_idx[0]].m_recv_box);
      }
      else if (m_dim == 6) {
         transferLocalData6D(m_comms[local_comm_idx[0]].m_send_box,
            m_comms[local_comm_idx[1]].m_recv_box);
         transferLocalData6D(m_comms[local_comm_idx[1]].m_send_box,
            m_comms[local_comm_idx[0]].m_recv_box);
      }
   }

   // Now process the receives and ensure the sends have completed for the
   // parallel data communication.
   // Process all the receives.
   for (int i = 0; i < num_parallel_comms; ++i) {
      // Find a completed receive.
      int recvIdx;
      MPI_Status stat;
      MPI_Waitany(num_parallel_comms, &recvReqs[0], &recvIdx, &stat);

      // Get the receive buffer and the communication neighbor with which this
      // receive is associated.
      int comm = recv_idx_to_comm[recvIdx];
      double* thisRecvBuffer = recvBuffers[recvIdx];

      // Fill the data buffer with the receive buffer for this message and free
      // the buffer.
      if (m_dim == 1) {
         transferFromRecvBuffer1D(m_comms[comm], thisRecvBuffer);
      }
      else if (m_dim == 2) {
         transferFromRecvBuffer2D(m_comms[comm], thisRecvBuffer);
      }
      else if (m_dim == 3) {
         transferFromRecvBuffer3D(m_comms[comm], thisRecvBuffer);
      }
      else if (m_dim == 4) {
         transferFromRecvBuffer4D(m_comms[comm], thisRecvBuffer);
      }
      else if (m_dim == 5) {
         transferFromRecvBuffer5D(m_comms[comm], thisRecvBuffer);
      }
      else {
         transferFromRecvBuffer6D(m_comms[comm], thisRecvBuffer);
      }
      delete [] thisRecvBuffer;
   }

   // Ensure all the sends have completed and free the send buffers.
   if (num_parallel_comms != 0) {
      vector<MPI_Status> sendStatus(num_parallel_comms);
      MPI_Waitall(num_parallel_comms, &sendReqs[0], &sendStatus[0]);
      for (int i = 0; i < num_parallel_comms; ++i) {
         delete [] sendBuffers[i];
      }
   }
}


void
ParallelArray::transferToSendBuffer1D(
   const Communication& a_comm,
   double* a_send_buffer) const
{
   vector<int> istart(1, 0), iend(1, 0);
   for (int dim = 0; dim < 1; ++dim) {
      istart[dim] = a_comm.m_send_box.lower(dim);
      iend[dim] = a_comm.m_send_box.upper(dim);
   }
   int buff_idx = 0;
   for (int i0 = istart[0]; i0 <= iend[0]; ++i0) {
      int data_idx = dataArrayIndex(i0);
      a_send_buffer[buff_idx] = m_data[data_idx];
      ++buff_idx;
   }
}


void
ParallelArray::transferToSendBuffer2D(
   const Communication& a_comm,
   double* a_send_buffer) const
{
   vector<int> istart(2, 0), iend(2, 0);
   for (int dim = 0; dim < 2; ++dim) {
      istart[dim] = a_comm.m_send_box.lower(dim);
      iend[dim] = a_comm.m_send_box.upper(dim);
   }
   int buff_idx = 0;
   for (int i1 = istart[1]; i1 <= iend[1]; ++i1) {
      for (int i0 = istart[0]; i0 <= iend[0]; ++i0) {
         int data_idx = dataArrayIndex(i0, i1);
         a_send_buffer[buff_idx] = m_data[data_idx];
         ++buff_idx;
      }
   }
}


void
ParallelArray::transferToSendBuffer3D(
   const Communication& a_comm,
   double* a_send_buffer) const
{
   vector<int> istart(3, 0), iend(3, 0);
   for (int dim = 0; dim < 3; ++dim) {
      istart[dim] = a_comm.m_send_box.lower(dim);
      iend[dim] = a_comm.m_send_box.upper(dim);
   }
   int buff_idx = 0;
   for (int i2 = istart[2]; i2 <= iend[2]; ++i2) {
      for (int i1 = istart[1]; i1 <= iend[1]; ++i1) {
         for (int i0 = istart[0]; i0 <= iend[0]; ++i0) {
            int data_idx = dataArrayIndex(i0, i1, i2);
            a_send_buffer[buff_idx] = m_data[data_idx];
            ++buff_idx;
         }
      }
   }
}


void
ParallelArray::transferToSendBuffer4D(
   const Communication& a_comm,
   double* a_send_buffer) const
{
   vector<int> istart(4, 0), iend(4, 0);
   for (int dim = 0; dim < 4; ++dim) {
      istart[dim] = a_comm.m_send_box.lower(dim);
      iend[dim] = a_comm.m_send_box.upper(dim);
   }
   int buff_idx = 0;
   for (int i3 = istart[3]; i3 <= iend[3]; ++i3) {
      for (int i2 = istart[2]; i2 <= iend[2]; ++i2) {
         for (int i1 = istart[1]; i1 <= iend[1]; ++i1) {
            for (int i0 = istart[0]; i0 <= iend[0]; ++i0) {
               int data_idx = dataArrayIndex(i0, i1, i2, i3);
               a_send_buffer[buff_idx] = m_data[data_idx];
               ++buff_idx;
            }
         }
      }
   }
}


void
ParallelArray::transferToSendBuffer5D(
   const Communication& a_comm,
   double* a_send_buffer) const
{
   vector<int> istart(5, 0), iend(5, 0);
   for (int dim = 0; dim < 5; ++dim) {
      istart[dim] = a_comm.m_send_box.lower(dim);
      iend[dim] = a_comm.m_send_box.upper(dim);
   }
   int buff_idx = 0;
   for (int i4 = istart[4]; i4 <= iend[4]; ++i4) {
      for (int i3 = istart[3]; i3 <= iend[3]; ++i3) {
         for (int i2 = istart[2]; i2 <= iend[2]; ++i2) {
            for (int i1 = istart[1]; i1 <= iend[1]; ++i1) {
               for (int i0 = istart[0]; i0 <= iend[0]; ++i0) {
                  int data_idx = dataArrayIndex(i0, i1, i2, i3, i4);
                  a_send_buffer[buff_idx] = m_data[data_idx];
                  ++buff_idx;
               }
            }
         }
      }
   }
}


void
ParallelArray::transferToSendBuffer6D(
   const Communication& a_comm,
   double* a_send_buffer) const
{
   vector<int> istart(6, 0), iend(6, 0);
   for (int dim = 0; dim < 6; ++dim) {
      istart[dim] = a_comm.m_send_box.lower(dim);
      iend[dim] = a_comm.m_send_box.upper(dim);
   }
   int buff_idx = 0;
   for (int i5 = istart[5]; i5 <= iend[5]; ++i5) {
      for (int i4 = istart[4]; i4 <= iend[4]; ++i4) {
         for (int i3 = istart[3]; i3 <= iend[3]; ++i3) {
            for (int i2 = istart[2]; i2 <= iend[2]; ++i2) {
               for (int i1 = istart[1]; i1 <= iend[1]; ++i1) {
                  for (int i0 = istart[0]; i0 <= iend[0]; ++i0) {
                     int data_idx = dataArrayIndex(i0, i1, i2, i3, i4, i5);
                     a_send_buffer[buff_idx] = m_data[data_idx];
                     ++buff_idx;
                  }
               }
            }
         }
      }
   }
}


void
ParallelArray::transferFromRecvBuffer1D(
   const Communication& a_comm,
   const double* a_recv_buffer)
{
   vector<int> istart(1, 0), iend(1, 0);
   for (int dim = 0; dim < 1; ++dim) {
      istart[dim] = a_comm.m_recv_box.lower(dim);
      iend[dim] = a_comm.m_recv_box.upper(dim);
   }
   int buff_idx = 0;
   for (int i0 = istart[0]; i0 <= iend[0]; ++i0) {
      int data_idx = dataArrayIndex(i0);
      m_data[data_idx] = a_recv_buffer[buff_idx];
      ++buff_idx;
   }
}


void
ParallelArray::transferFromRecvBuffer2D(
   const Communication& a_comm,
   const double* a_recv_buffer)
{
   vector<int> istart(2, 0), iend(2, 0);
   for (int dim = 0; dim < 2; ++dim) {
      istart[dim] = a_comm.m_recv_box.lower(dim);
      iend[dim] = a_comm.m_recv_box.upper(dim);
   }
   int buff_idx = 0;
   for (int i1 = istart[1]; i1 <= iend[1]; ++i1) {
      for (int i0 = istart[0]; i0 <= iend[0]; ++i0) {
         int data_idx = dataArrayIndex(i0, i1);
         m_data[data_idx] = a_recv_buffer[buff_idx];
         ++buff_idx;
      }
   }
}


void
ParallelArray::transferFromRecvBuffer3D(
   const Communication& a_comm,
   const double* a_recv_buffer)
{
   vector<int> istart(3, 0), iend(3, 0);
   for (int dim = 0; dim < 3; ++dim) {
      istart[dim] = a_comm.m_recv_box.lower(dim);
      iend[dim] = a_comm.m_recv_box.upper(dim);
   }
   int buff_idx = 0;
   for (int i2 = istart[2]; i2 <= iend[2]; ++i2) {
      for (int i1 = istart[1]; i1 <= iend[1]; ++i1) {
         for (int i0 = istart[0]; i0 <= iend[0]; ++i0) {
            int data_idx = dataArrayIndex(i0, i1, i2);
            m_data[data_idx] = a_recv_buffer[buff_idx];
            ++buff_idx;
         }
      }
   }
}


void
ParallelArray::transferFromRecvBuffer4D(
   const Communication& a_comm,
   const double* a_recv_buffer)
{
   vector<int> istart(4, 0), iend(4, 0);
   for (int dim = 0; dim < 4; ++dim) {
      istart[dim] = a_comm.m_recv_box.lower(dim);
      iend[dim] = a_comm.m_recv_box.upper(dim);
   }
   int buff_idx = 0;
   for (int i3 = istart[3]; i3 <= iend[3]; ++i3) {
      for (int i2 = istart[2]; i2 <= iend[2]; ++i2) {
         for (int i1 = istart[1]; i1 <= iend[1]; ++i1) {
            for (int i0 = istart[0]; i0 <= iend[0]; ++i0) {
               int data_idx = dataArrayIndex(i0, i1, i2, i3);
               m_data[data_idx] = a_recv_buffer[buff_idx];
               ++buff_idx;
            }
         }
      }
   }
}


void
ParallelArray::transferFromRecvBuffer5D(
   const Communication& a_comm,
   const double* a_recv_buffer)
{
   vector<int> istart(5, 0), iend(5, 0);
   for (int dim = 0; dim < 5; ++dim) {
      istart[dim] = a_comm.m_recv_box.lower(dim);
      iend[dim] = a_comm.m_recv_box.upper(dim);
   }
   int buff_idx = 0;
   for (int i4 = istart[4]; i4 <= iend[4]; ++i4) {
      for (int i3 = istart[3]; i3 <= iend[3]; ++i3) {
         for (int i2 = istart[2]; i2 <= iend[2]; ++i2) {
            for (int i1 = istart[1]; i1 <= iend[1]; ++i1) {
               for (int i0 = istart[0]; i0 <= iend[0]; ++i0) {
                  int data_idx = dataArrayIndex(i0, i1, i2, i3, i4);
                  m_data[data_idx] = a_recv_buffer[buff_idx];
                  ++buff_idx;
               }
            }
         }
      }
   }
}


void
ParallelArray::transferFromRecvBuffer6D(
   const Communication& a_comm,
   const double* a_recv_buffer)
{
   vector<int> istart(6, 0), iend(6, 0);
   for (int dim = 0; dim < 6; ++dim) {
      istart[dim] = a_comm.m_recv_box.lower(dim);
      iend[dim] = a_comm.m_recv_box.upper(dim);
   }
   int buff_idx = 0;
   for (int i5 = istart[5]; i5 <= iend[5]; ++i5) {
      for (int i4 = istart[4]; i4 <= iend[4]; ++i4) {
         for (int i3 = istart[3]; i3 <= iend[3]; ++i3) {
            for (int i2 = istart[2]; i2 <= iend[2]; ++i2) {
               for (int i1 = istart[1]; i1 <= iend[1]; ++i1) {
                  for (int i0 = istart[0]; i0 <= iend[0]; ++i0) {
                     int data_idx = dataArrayIndex(i0, i1, i2, i3, i4, i5);
                     m_data[data_idx] = a_recv_buffer[buff_idx];
                     ++buff_idx;
                  }
               }
            }
         }
      }
   }
}


void
ParallelArray::transferLocalData1D(
   const Box& a_send_box,
   const Box& a_recv_box)
{
   ARRAY_ASSERT(a_send_box.dim() == 1 && a_recv_box.dim() == 1);
   ARRAY_ASSERT(a_send_box.numberOfCells(0) == a_recv_box.numberOfCells(0));
   int i0end = a_send_box.numberOfCells(0);
   for (int i0 = 0; i0 < i0end; ++i0) {
      int recv_data_idx = dataArrayIndex(a_recv_box.lower(0) + i0);
      int send_data_idx = dataArrayIndex(a_send_box.lower(0) + i0);
      m_data[recv_data_idx] = m_data[send_data_idx];
   }
}


void
ParallelArray::transferLocalData2D(
   const Box& a_send_box,
   const Box& a_recv_box)
{
   ARRAY_ASSERT(a_send_box.dim() == 2 && a_recv_box.dim() == 2);
   ARRAY_ASSERT(a_send_box.numberOfCells(0) == a_recv_box.numberOfCells(0));
   ARRAY_ASSERT(a_send_box.numberOfCells(1) == a_recv_box.numberOfCells(1));
   int i0end = a_send_box.numberOfCells(0);
   int i1end = a_send_box.numberOfCells(1);
   for (int i1 = 0; i1 < i1end; ++i1) {
      for (int i0 = 0; i0 < i0end; ++i0) {
         int recv_data_idx = dataArrayIndex(a_recv_box.lower(0) + i0,
            a_recv_box.lower(1) + i1);
         int send_data_idx = dataArrayIndex(a_send_box.lower(0) + i0,
            a_send_box.lower(1) + i1);
         m_data[recv_data_idx] = m_data[send_data_idx];
      }
   }
}


void
ParallelArray::transferLocalData3D(
   const Box& a_send_box,
   const Box& a_recv_box)
{
   ARRAY_ASSERT(a_send_box.dim() == 3 && a_recv_box.dim() == 3);
   ARRAY_ASSERT(a_send_box.numberOfCells(0) == a_recv_box.numberOfCells(0));
   ARRAY_ASSERT(a_send_box.numberOfCells(1) == a_recv_box.numberOfCells(1));
   ARRAY_ASSERT(a_send_box.numberOfCells(2) == a_recv_box.numberOfCells(2));
   int i0end = a_send_box.numberOfCells(0);
   int i1end = a_send_box.numberOfCells(1);
   int i2end = a_send_box.numberOfCells(2);
   for (int i2 = 0; i2 < i2end; ++i2) {
      for (int i1 = 0; i1 < i1end; ++i1) {
         for (int i0 = 0; i0 < i0end; ++i0) {
            int recv_data_idx = dataArrayIndex(a_recv_box.lower(0) + i0,
               a_recv_box.lower(1) + i1,
               a_recv_box.lower(2) + i2);
            int send_data_idx = dataArrayIndex(a_send_box.lower(0) + i0,
               a_send_box.lower(1) + i1,
               a_send_box.lower(2) + i2);
            m_data[recv_data_idx] = m_data[send_data_idx];
         }
      }
   }
}


void
ParallelArray::transferLocalData4D(
   const Box& a_send_box,
   const Box& a_recv_box)
{
   ARRAY_ASSERT(a_send_box.dim() == 4 && a_recv_box.dim() == 4);
   ARRAY_ASSERT(a_send_box.numberOfCells(0) == a_recv_box.numberOfCells(0));
   ARRAY_ASSERT(a_send_box.numberOfCells(1) == a_recv_box.numberOfCells(1));
   ARRAY_ASSERT(a_send_box.numberOfCells(2) == a_recv_box.numberOfCells(2));
   ARRAY_ASSERT(a_send_box.numberOfCells(3) == a_recv_box.numberOfCells(3));
   int i0end = a_send_box.numberOfCells(0);
   int i1end = a_send_box.numberOfCells(1);
   int i2end = a_send_box.numberOfCells(2);
   int i3end = a_send_box.numberOfCells(3);
   for (int i3 = 0; i3 < i3end; ++i3) {
      for (int i2 = 0; i2 < i2end; ++i2) {
         for (int i1 = 0; i1 < i1end; ++i1) {
            for (int i0 = 0; i0 < i0end; ++i0) {
               int recv_data_idx = dataArrayIndex(a_recv_box.lower(0) + i0,
                  a_recv_box.lower(1) + i1,
                  a_recv_box.lower(2) + i2,
                  a_recv_box.lower(3) + i3);
               int send_data_idx = dataArrayIndex(a_send_box.lower(0) + i0,
                  a_send_box.lower(1) + i1,
                  a_send_box.lower(2) + i2,
                  a_send_box.lower(3) + i3);
               m_data[recv_data_idx] = m_data[send_data_idx];
            }
         }
      }
   }
}


void
ParallelArray::transferLocalData5D(
   const Box& a_send_box,
   const Box& a_recv_box)
{
   ARRAY_ASSERT(a_send_box.dim() == 5 && a_recv_box.dim() == 5);
   ARRAY_ASSERT(a_send_box.numberOfCells(0) == a_recv_box.numberOfCells(0));
   ARRAY_ASSERT(a_send_box.numberOfCells(1) == a_recv_box.numberOfCells(1));
   ARRAY_ASSERT(a_send_box.numberOfCells(2) == a_recv_box.numberOfCells(2));
   ARRAY_ASSERT(a_send_box.numberOfCells(3) == a_recv_box.numberOfCells(3));
   ARRAY_ASSERT(a_send_box.numberOfCells(4) == a_recv_box.numberOfCells(4));
   int i0end = a_send_box.numberOfCells(0);
   int i1end = a_send_box.numberOfCells(1);
   int i2end = a_send_box.numberOfCells(2);
   int i3end = a_send_box.numberOfCells(3);
   int i4end = a_send_box.numberOfCells(4);
   for (int i4 = 0; i4 < i4end; ++i4) {
      for (int i3 = 0; i3 < i3end; ++i3) {
         for (int i2 = 0; i2 < i2end; ++i2) {
            for (int i1 = 0; i1 < i1end; ++i1) {
               for (int i0 = 0; i0 < i0end; ++i0) {
                  int recv_data_idx = dataArrayIndex(a_recv_box.lower(0) + i0,
                     a_recv_box.lower(1) + i1,
                     a_recv_box.lower(2) + i2,
                     a_recv_box.lower(3) + i3,
                     a_recv_box.lower(4) + i4);
                  int send_data_idx = dataArrayIndex(a_send_box.lower(0) + i0,
                     a_send_box.lower(1) + i1,
                     a_send_box.lower(2) + i2,
                     a_send_box.lower(3) + i3,
                     a_send_box.lower(4) + i4);
                  m_data[recv_data_idx] = m_data[send_data_idx];
               }
            }
         }
      }
   }
}


void
ParallelArray::transferLocalData6D(
   const Box& a_send_box,
   const Box& a_recv_box)
{
   ARRAY_ASSERT(a_send_box.dim() == 6 && a_recv_box.dim() == 6);
   ARRAY_ASSERT(a_send_box.numberOfCells(0) == a_recv_box.numberOfCells(0));
   ARRAY_ASSERT(a_send_box.numberOfCells(1) == a_recv_box.numberOfCells(1));
   ARRAY_ASSERT(a_send_box.numberOfCells(2) == a_recv_box.numberOfCells(2));
   ARRAY_ASSERT(a_send_box.numberOfCells(3) == a_recv_box.numberOfCells(3));
   ARRAY_ASSERT(a_send_box.numberOfCells(4) == a_recv_box.numberOfCells(4));
   ARRAY_ASSERT(a_send_box.numberOfCells(5) == a_recv_box.numberOfCells(5));
   int i0end = a_send_box.numberOfCells(0);
   int i1end = a_send_box.numberOfCells(1);
   int i2end = a_send_box.numberOfCells(2);
   int i3end = a_send_box.numberOfCells(3);
   int i4end = a_send_box.numberOfCells(4);
   int i5end = a_send_box.numberOfCells(5);
   for (int i5 = 0; i5 < i5end; ++i5) {
      for (int i4 = 0; i4 < i4end; ++i4) {
         for (int i3 = 0; i3 < i3end; ++i3) {
            for (int i2 = 0; i2 < i2end; ++i2) {
               for (int i1 = 0; i1 < i1end; ++i1) {
                  for (int i0 = 0; i0 < i0end; ++i0) {
                     int recv_data_idx = dataArrayIndex(a_recv_box.lower(0) + i0,
                        a_recv_box.lower(1) + i1,
                        a_recv_box.lower(2) + i2,
                        a_recv_box.lower(3) + i3,
                        a_recv_box.lower(4) + i4,
                        a_recv_box.lower(5) + i5);
                     int send_data_idx = dataArrayIndex(a_send_box.lower(0) + i0,
                        a_send_box.lower(1) + i1,
                        a_send_box.lower(2) + i2,
                        a_send_box.lower(3) + i3,
                        a_send_box.lower(4) + i4,
                        a_send_box.lower(5) + i5);
                     m_data[recv_data_idx] = m_data[send_data_idx];
                  }
               }
            }
         }
      }
   }
}


ParallelArray::Communication::Communication(
   int a_dim)
   : m_send_box(a_dim),
     m_recv_box(a_dim)
{
}


ParallelArray::Communication::~Communication()
{
}


ParallelArray::Communication::Communication(
   const ParallelArray::Communication& a_other)
   : m_partner_rank(a_other.m_partner_rank),
     m_comm_type(a_other.m_comm_type),
     m_send_box(a_other.m_send_box),
     m_send_tag(a_other.m_send_tag),
     m_recv_box(a_other.m_recv_box),
     m_recv_tag(a_other.m_recv_tag)
{
}


ParallelArray::Communication&
ParallelArray::Communication::operator = (
   const ParallelArray::Communication& a_rhs)
{
   m_partner_rank = a_rhs.m_partner_rank;
   m_comm_type = a_rhs.m_comm_type;
   m_send_box = a_rhs.m_send_box;
   m_send_tag = a_rhs.m_send_tag;
   m_recv_box = a_rhs.m_recv_box;
   m_recv_tag = a_rhs.m_recv_tag;
}


ParallelArray::Box::Box()
   : m_dim(-1)
{
}


ParallelArray::Box::Box(
   int a_dim)
   : m_dim(a_dim),
     m_upper(a_dim, 0),
     m_lower(a_dim, -1)
{
}


ParallelArray::Box::~Box()
{
}


ParallelArray::Box::Box(
   const ParallelArray::Box& a_other)
   : m_dim(a_other.m_dim),
     m_upper(a_other.m_upper),
     m_lower(a_other.m_lower)
{
}


ParallelArray::Box&
ParallelArray::Box::operator = (
   const ParallelArray::Box& a_rhs)
{
   m_dim = a_rhs.m_dim;
   m_upper = a_rhs.m_upper;
   m_lower = a_rhs.m_lower;
}
