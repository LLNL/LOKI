/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _LOKI_UTILITIES_H
#define _LOKI_UTILITIES_H

#include "ProblemDomain.H"
#include <mpi.h>
#include <stdio.h>
#include <stdarg.h>

namespace Loki {

/**
 * Application wide utility functions.
 */
class Loki_Utilities
{
public:
   /**
    * @brief Turn a boolean into an int, perform a logical and reduction on it,
    *        and return the result as a boolean.
    *
    * @param[in] a_bool The boolean to be reduced.
    */
   static
   bool
   reduceBoolean(
      bool a_bool);

   /**
    * @brief Sums a value over all processors in a communicator and returns
    *        the result to the specified processor.
    *
    * @param[in] a_val Value from a processor to be summed.
    * @param[in] a_processor Processor to receive result.  A value of -1
    *                        means all processors get the result.
    * @param[in] a_comm The communicator the summation occurs over.
    */
   static
   double
   getSum(
      double a_val,
      int a_processor = -1,
      MPI_Comm a_comm = MPI_COMM_WORLD);

   /**
    * @brief Sums each entry from an array of values over all processors in a
    *        communicator and places the results in an array on the specified
    *        processor.
    *
    * @param[in] a_vals Values from a processor to be summed.
    * @param[out] a_sums Sums for processor a_processor.  Not used on other
    *                    processors.
    * @param[in] a_size The size of both a_vals and a_sums.
    * @param[in] a_processor Processor to receive result.  A value of -1
    *                        means all processors get the result.
    * @param[in] a_comm The communicator the summation occurs over.
    */
   static
   void
   getSums(
      const double* a_vals,
      double* a_sums,
      int a_size,
      int a_processor = -1,
      MPI_Comm a_comm = MPI_COMM_WORLD);

   /**
    * @brief Gets the maximum of a value over all processors in a communicator
    *        and returns the result to the specified processor.
    *
    * @param[in] a_val Value from a processor to be maxed.
    * @param[in] a_processor Processor to receive result.  A value of -1
    *                        means all processors get the result.
    * @param[in] a_comm The communicator the operation occurs over.
    */
   static
   double
   getMaxValue(
      double a_val,
      int a_processor = -1,
      MPI_Comm a_comm = MPI_COMM_WORLD);

   /**
    * @brief Gets the maximum of each entry from an array of values over all
    *        processors in a communicator and places the results in an array on
    *        the specified processor.
    *
    * @param[in] a_vals Values from a processor to be maxed.
    * @param[out] a_maxs Maxima for processor a_processor.  Not used on other
    *                    processors.
    * @param[in] a_size The size of both a_vals and a_maxs.
    * @param[in] a_processor Processor to receive result.  A value of -1
    *                        means all processors get the result.
    * @param[in] a_comm The communicator the operation occurs over.
    */
   static
   void
   getMaxValues(
      const double* a_vals,
      double* a_maxs,
      int a_size,
      int a_processor = -1,
      MPI_Comm a_comm = MPI_COMM_WORLD);

   /**
    * @brief Gets the minimum of a value over all processors in a communicator
    *        and returns the result to the specified processor.
    *
    * @param[in] a_val Value from a processor to be mined.
    * @param[in] a_processor Processor to receive result.  A value of -1
    *                        means all processors get the result.
    * @param[in] a_comm The communicator the operation occurs over.
    */
   static
   double
   getMinValue(
      double a_val,
      int a_processor = -1,
      MPI_Comm a_comm = MPI_COMM_WORLD);

   /**
    * @brief Gets the minimum of each entry from an array of values over all
    *        processors in a communicator and places the results in an array on
    *        the specified processor.
    *
    * @param[in] a_vals Values from a processor to be mined.
    * @param[out] a_mins Minima for processor a_processor.  Not used on other
    *                    processors.
    * @param[in] a_size The size of both a_vals and a_mins.
    * @param[in] a_processor Processor to receive result.  A value of -1
    *                        means all processors get the result.
    * @param[in] a_comm The communicator the operation occurs over.
    */
   static
   void
   getMinValues(
      const double* a_vals,
      double* a_mins,
      int a_size,
      int a_processor = -1,
      MPI_Comm a_comm = MPI_COMM_WORLD);

   // Global processor ID.
   static int s_my_id;

   // Global number of processors.
   static int s_num_procs;

   /**
    * @brief Sets the values of s_my_id and s_num_procs.
    */
   static
   void
   initialize()
   {
      MPI_Comm_rank(MPI_COMM_WORLD, &s_my_id);
      MPI_Comm_size(MPI_COMM_WORLD, &s_num_procs);
   }

   /**
    * @brief Prints a message only on process 0.  Essentially printf on
    *        process 0.
    *
    * @param[in] a_format A printf output format.
    */
   static
   void
   printF(
     const char* a_format,
      ...)
   {
      if (s_my_id >0) {
         return;
      }
      va_list args;
      va_start(args, a_format);
      vprintf(a_format, args);
      va_end(args);
   }

private:
   // Unimplemented default constructor.
   Loki_Utilities();

   // Unimplemented destructor.
   ~Loki_Utilities();

   // Unimplemented copy constructor.
   Loki_Utilities(
      const Loki_Utilities& other);

   // Unimplemented assignment operator.
   Loki_Utilities&
   operator = (
      const Loki_Utilities& rhs);

};

} // end namespace Loki

#endif
