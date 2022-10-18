/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _SERIALIZABLE_H_
#define _SERIALIZABLE_H_

#include "RestartReader.H"
#include "RestartWriter.H"

namespace Loki {

/**
 * Serializable is an abstract base class defining the interface to any class
 * that can write to and be restored from restart files.
 */
class Serializable
{
public:
   /**
    * @brief Constructor.
    */
   Serializable();

   /**
    * @brief Destructor.
    */
   virtual
   ~Serializable();

   /**
    * @brief Read saved information from the restart database.
    *
    * @param[in] a_reader The restart file reader.
    */
   virtual
   void
   getFromRestart(
      RestartReader& a_reader) = 0;

   /**
    * @brief Write information to save to the restart database.
    *
    * @param[out] a_writer The restart file writer.
    * @param[in] a_time Current simulation time.
    */
   virtual
   void
   putToRestart(
      RestartWriter& a_writer,
      double a_time) = 0;

private:
   // Unimplemented copy constructor.
   Serializable(
      const Serializable& other);

   // Unimplemented assignment operator.
   Serializable&
   operator = (
      const Serializable& rhs);

};

} // end namespace Loki

#endif
