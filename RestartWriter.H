/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _RESTART_WRITER_H_
#define _RESTART_WRITER_H_

#include "RestartReaderWriterBase.H"
#include "ParallelArray.H"
#include <sstream>

namespace Loki {

/**
 * A class to write Loki restart files.
 * Each restart dump consists of a collection of individual files.  One is
 * designated for metadata.  The others are a user supplied maximum number of
 * bulkdata files which contain each processor's local 4D distribution data
 * and the local 2D field data from the electromagnetic processors.  The
 * metadata file is written using HDF5's parallel IO capabilities.
 * Specifically, H5FD_MPIO_INDEPENDENT is used as the data transfer protocol.
 * As this file contains metadata which is common to groups of processors, only
 * a subset of all the processors need to actually perform the data writes.
 * Metadata relating to the electromagnetics is written only by the
 * electromagnetic processor.  Metadata relating to each KineticSpecies is
 * written only by the lowest rank processor of the processors that species is
 * distributed over.  Bulkdata is written using a baton passing scheme.  A group
 * of processors is assigned a bulkdata file to write to.  Each processor
 * writes its local data and then signals the next processor in the group that
 * it may proceed with its writes.  Hence, writing bulkdata is parallelized
 * but increased parallelism comes at the cost of more bulkdata files.
 */
class RestartWriter : public RestartReaderWriterBase
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_base_name The base name common to all restart files.
    * @param[in] a_max_num_files The maximum number of bulk data files to write.
    * @param[in] a_post_processing If true this is being called from the serial
    *                              post processor.
    */
   RestartWriter(
      const string& a_base_name,
      int a_max_num_files,
      bool a_post_processing = false);

   /**
    * @brief Destructor.
    */
   virtual
   ~RestartWriter();

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
      RestartReaderWriterBase::pushSubDir(a_name, true);
   }

   /**
    * @brief Write a single integer to the currently most nested group in the
    *        metadata file.  This should be called collectively.  Only one
    *        processor has and writes the appropriate metadata.  The other
    *        processors only create all the necessary database entities.  If
    *        not called collectively then the metadata file can not be closed
    *        as it will be structured differently on different processors.
    *
    * @param[in] a_name Name of the value.
    * @param[in] a_val The value.
    * @param[in] a_write_data If true this processor actually writes its data.
    */
   void
   writeIntegerValue(
      const string& a_name,
      const int& a_val,
      bool a_write_data);

   /**
    * @brief Write an array of integers to the currently most nested group in
    *        the metadata file.  This should be called collectively.  Only one
    *        processor has and writes the appropriate metadata.  The other
    *        processors only create all the necessary database entities.  If
    *        not called collectively then the metadata file can not be closed
    *        as it will be structured differently on different processors.
    *
    * @param[in] a_name Name of the values.
    * @param[in] a_vals The values.
    * @param[in] a_num_vals The number of integers in the array.
    * @param[in] a_write_data If true this processor actually writes its data.
    */
   void
   writeIntegerArray(
      const string& a_name,
      const int* a_vals,
      int a_num_vals,
      bool a_write_data);

   /**
    * @brief Write a single double to the currently most nested group in the
    *        metadata file.  This should be called collectively.  Only one
    *        processor has and writes the appropriate metadata.  The other
    *        processors only create all the necessary database entities.  If
    *        not called collectively then the metadata file can not be closed
    *        as it will be structured differently on different processors.
    *
    * @param[in] a_name Name of the value.
    * @param[in] a_val The value.
    * @param[in] a_write_data If true this processor actually writes its data.
    */
   void
   writeDoubleValue(
      const string& a_name,
      const double& a_val,
      bool a_write_data);

   /**
    * @brief This differs from writeDoubleValue above in that the value being
    *        written is not metadata.  It is bulk data which is unique to each
    *        processor.  This is similar to writeParallelArray except for the
    *        nature of what is being written.
    *
    * @param[in] a_name Name of the value.
    * @param[in] a_val The value.
    */
   void
   writeBulkDoubleValue(
      const string& a_name,
      const double& a_val);

   /**
    * @brief Write an array of doubles to the currently most nested group in the
    *        metadata file.  This should be called collectively.  Only one
    *        processor has and writes the appropriate metadata.  The other
    *        processors only create all the necessary database entities.  If
    *        not called collectively then the metadata file can not be closed
    *        as it will be structured differently on different processors.
    *
    * @param[in] a_name Name of the values.
    * @param[in] a_vals The values.
    * @param[in] a_num_vals The number of doubles in the array.
    * @param[in] a_write_data If true this processor actually writes its data.
    */
   void
   writeDoubleArray(
      const string& a_name,
      const double* a_vals,
      int a_num_vals,
      bool a_write_data);

   /**
    * @brief Write a string to the currently most nested group in the metadata
    *        file.  This should be called collectively.  Only one processor
    *        has and writes the appropriate metadata.  The other processors
    *        only create all the necessary database entities.  If not called
    *        collectively then the metadata file can not be closed as it will
    *        be structured differently on different processors.
    *
    * @param[in] a_name Name of the string.
    * @param[in] a_val The string.
    * @param[in] a_write_data If true this processor actually writes its data.
    */
   void
   writeString(
      const string& a_name,
      const string& a_val,
      bool a_write_data);

   /**
    * @brief Write a ParallelArray.  The bulk data is written to the root group
    *        of a bulkdata file.  The metadata is written to the metadata file.
    *
    * @param[in] a_name Name of the ParallelArray.
    * @param[in] a_array The ParallelArray.
    * @param[in] a_write_meta_data If true write the associated metadata.
    */
   void
   writeParallelArray(
      const string& a_name,
      const ParallelArray& a_array,
      bool a_write_meta_data);

private:
   // Unimplemented default constructor.
   RestartWriter();

   // Unimplemented copy constructor.
   RestartWriter(
      const RestartWriter& a_other);

   // Unimplemented assignment operator.
   RestartWriter&
   operator = (
      const RestartWriter& a_rhs);

   // The name of this processor's bulkddata file.
   ostringstream m_bulkdata_file_name;

   // Tag for baton passing communication.
   static const int s_TAG_BATON;
};

} // end namespace Loki

#endif
