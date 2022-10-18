/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "Loki_Utilities.H"

namespace Loki {

int Loki_Utilities::s_my_id = -1;
int Loki_Utilities::s_num_procs = -1;

bool
Loki_Utilities::reduceBoolean(
   bool a_bool)
{
   // Booleans are some implementation dependent form of an int.  Assign the
   // boolean's value to a literal int and reduce it as there is no MPI type
   // for a bool.
   int local(a_bool ? 1 : 0);
   int global(1);
   MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
   return (global == 1);
}


double
Loki_Utilities::getSum(
   double a_val,
   int a_processor,
   MPI_Comm a_comm)
{
   double sum = 0.0;
   if (a_processor == -1) {
      MPI_Allreduce(&a_val, &sum, 1, MPI_DOUBLE, MPI_SUM, a_comm);
   }
   else {
      MPI_Reduce(&a_val, &sum, 1, MPI_DOUBLE, MPI_SUM, a_processor, a_comm);
   }
   return sum;
}


void
Loki_Utilities::getSums(
   const double* a_vals,
   double* a_sums,
   int a_size,
   int a_processor,
   MPI_Comm a_comm)
{
   if (a_processor == -1) {
      MPI_Allreduce(a_vals, a_sums, a_size, MPI_DOUBLE, MPI_SUM, a_comm);
   }
   else {
      MPI_Reduce(a_vals,
         a_sums,
         a_size,
         MPI_DOUBLE,
         MPI_SUM,
         a_processor,
         a_comm);
   }
}


double
Loki_Utilities::getMaxValue(
   double a_val,
   int a_processor,
   MPI_Comm a_comm)
{
   double maxVal = a_val;
   if (a_processor == -1) {
      MPI_Allreduce(&a_val, &maxVal, 1, MPI_DOUBLE, MPI_MAX, a_comm);
   }
   else {
      MPI_Reduce(&a_val,
         &maxVal,
         1,
         MPI_DOUBLE,
         MPI_MAX,
         a_processor,
         a_comm);
   }
   return maxVal;
}


void
Loki_Utilities::getMaxValues(
   const double* a_vals,
   double* a_maxs,
   int a_size,
   int a_processor,
   MPI_Comm a_comm)
{
   if (a_processor == -1) {
      MPI_Allreduce(a_vals, a_maxs, a_size, MPI_DOUBLE, MPI_MAX, a_comm);
   }
   else {
      MPI_Reduce(a_vals,
         a_maxs,
         a_size,
         MPI_DOUBLE,
         MPI_MAX,
         a_processor,
         a_comm);
   }
}


double
Loki_Utilities::getMinValue(
   double a_val,
   int a_processor,
   MPI_Comm a_comm)
{
   double minVal = a_val;
   if (a_processor == -1) {
      MPI_Allreduce(&a_val, &minVal, 1, MPI_DOUBLE, MPI_MIN, a_comm);
   }
   else {
      MPI_Reduce(&a_val, &minVal, 1, MPI_DOUBLE, MPI_MIN, a_processor, a_comm);
   }
   return minVal;
}


void
Loki_Utilities::getMinValues(
   const double* a_vals,
   double* a_mins,
   int a_size,
   int a_processor,
   MPI_Comm a_comm)
{
   if (a_processor == -1) {
      MPI_Allreduce(a_vals, a_mins, a_size, MPI_DOUBLE, MPI_MIN, a_comm);
   }
   else {
      MPI_Reduce(a_vals,
         a_mins,
         a_size,
         MPI_DOUBLE,
         MPI_MIN,
         a_processor,
         a_comm);
   }
}

} // end namespace Loki
