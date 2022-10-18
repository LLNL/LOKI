/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _RESTART_READER_H_
#define _RESTART_READER_H_

#include "RestartReaderWriterBase.H"
#include "ParallelArray.H"
#include "tbox/Box.H"

namespace Loki {

/**
 * A class to read Loki restart files.
 * Each restart dump consists of a collection of individual files.  One is
 * designated for metadata.  The others are a user supplied maximum number of
 * bulkdata files which contain each processor's local 4D distribution data
 * and the local 2D field data from the electromagnetic processors.  The
 * metadata file is read using HDF5's parallel IO capabilities.  Specifically,
 * H5FD_MPIO_INDEPENDENT is used as the data transfer protocol.  As opposed to
 * the writing of the metadata each processor must read the metatdata.  Bulkdata
 * is read in the same way that it was written, using a baton passing scheme.
 * A group of processors is assigned a bulkdata file to read from.  Each
 * processor reads its local data and then signals the next processor in the
 * group that it may proceed with its reads.  Hence, reading bulkdata is
 * parallelized but increased parallelism comes at the cost of more bulkdata
 * files.
 */
class RestartReader : public RestartReaderWriterBase
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_base_name The base name common to all restart files.
    * @param[in] a_max_num_files The maximum number of bulk data files to read.
    * @param[in] a_post_processing If true this is being called from the serial
    *                              post processor.
    */
   RestartReader(
      const string& a_base_name,
      int a_max_num_files,
      bool a_post_processing = false);

   /**
    * @brief Destructor.
    */
   virtual
   ~RestartReader();

   /**
    * @brief Create the named nested group within the metadata file.  Any data
    *        write is placed in the currently most nested group.
    *
    * @param[in] a_name The name of the group.
    */
   void
   pushSubDir(
      const string& a_name)
   {
      RestartReaderWriterBase::pushSubDir(a_name, false);
   }

   /**
    * @brief Read a single integer from the group with read context in the
    *        metadata file.  This should be called collectively.  If not called
    *        collectively then the metadata file can not be closed as it will be
    *        structured differently on different processors.
    *
    * @param[in] a_name Name of the value.
    * @param[in] a_val The value.
    */
   void
   readIntegerValue(
      const string& a_name,
      int& a_val);

   /**
    * @brief Read an array of integers from the group with read context in the
    *        metadata file.  This should be called collectively.  If not called
    *        collectively then the metadata file can not be closed as it will be
    *        structured differently on different processors.
    *
    * @param[in] a_name Name of the values.
    * @param[in] a_vals The values.
    * @param[in] a_num_vals The number of integers in the array.
    */
   void
   readIntegerArray(
      const string& a_name,
      int* a_vals,
      int a_num_vals);

   /**
    * @brief Read a single double from the group with read context in the
    *        metadata file.  This should be called collectively.  If not called
    *        collectively then the metadata file can not be closed as it will
    *        be structured differently on different processors.
    *
    * @param[in] a_name Name of the value.
    * @param[in] a_val The value.
    */
   void
   readDoubleValue(
      const string& a_name,
      double& a_val);

   /**
    * @brief This differs from writeDoubleValue above in that the value being
    *        read is not metadata.  It is bulk data which is unique to each
    *        processor.  Thus is is similar to writeParallelArray except for
    *        the nature of what is being read.
    *
    * @param[in] a_name Name of the value.
    * @param[in] a_val The value.
    */
   void
   readBulkDoubleValue(
      const string& a_name,
      double& a_val);

   /**
    * @brief Read an array of doubles from the group with read context in the
    *        metadata file.  This should be called collectively.  If not called
    *        collectively then the metadata file can not be closed as it will
    *        be structured differently on different processors.
    *
    * @param[in] a_name Name of the values.
    * @param[in] a_vals The values.
    * @param[in] a_num_vals The number of doubles in the array.
    */
   void
   readDoubleArray(
      const string& a_name,
      double* a_vals,
      int a_num_vals);

   /**
    * @brief Read a string frp, the group with read context in the metadata
    *        file.  This should be called collectively.  If not called
    *        collectively then the metadata file can not be closed as it will
    *        be structured differently on different processors.
    *
    * @param[in] a_name Name of the string.
    * @param[in] a_val The string.
    */
   void
   readString(
      const string& a_name,
      string& a_val);

   /**
    * @brief Read a ParallelArray.  The bulk data is read from the root group of
    *        a bulkdata file.  The metadata is read from the root group of the
    *        metadata file.  This method is only for the use of the main code.
    *
    * @param[in] a_name Name of the ParallelArray.
    * @param[in] a_array The ParallelArray.
    */
   void
   readParallelArray(
      const string& a_name,
      ParallelArray& a_array);

   /**
    * @brief Read a ParallelArray.  The bulk data is read from the root group of
    *        a bulkdata file.  The metadata is read from the root group of the
    *        metadata file.  This method is only for the use of the post
    *        processor.
    *
    * @param[in] a_name Name of the ParallelArray.
    * @param[in] a_array The ParallelArray.
    * @param[in] a_num_generating_processes The number of processors that
    *                                       generated the restart dump.
    */
   void
   readParallelArray(
      const string& a_name,
      ParallelArray& a_array,
      int a_num_generating_processes);

private:
   // Unimplemented default constructor.
   RestartReader();

   // Unimplemented copy constructor.
   RestartReader(
      const RestartReader& a_other);

   // Unimplemented assignment operator.
   RestartReader&
   operator = (
      const RestartReader& a_rhs);

   // Get the data box for a processor owning a distributed ParallelArray.
   void
   getArrayDataBox(
      int a_p,
      int a_baseProc,
      const int* a_dimBase,
      const int* a_dimProc,
      const int* a_dimSizesLeft,
      const int* a_dimNumLeft,
      tbox::Box& a_box);

   // Tag for baton passing communication.
   static const int s_TAG_BATON;
};

} // end namespace Loki

#endif

