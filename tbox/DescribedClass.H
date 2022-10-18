/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef included_tbox_DescribedClass
#define included_tbox_DescribedClass

namespace Loki {
namespace tbox {

/**
 * @brief
 * Required base class for all classes to be used with smart pointers.
 *
 * All of the classes pointed to by a Loki smart pointer need to be
 * derived from this class.  Loki smart pointers require that the
 * C++ run-time type identification (RTTI) dynamic casting operator (
 * dynamic_cast() ) work on the objects being pointed to.  The dynamic
 * cast operator requires the classes are in the same class hierarchy.
 * Deriving from DesribedClasss guarantees this hence the reason for
 * this unattractive requirement.
 *
 * Notes:
 *
 * The Loki developers recognize this requirement is a not ideal and
 * that the dynamic_cast that is being done implicitly in
 * the smart pointer assignment is considered to be a poor choice by
 * many other reference counting pointer implementations.  With
 * hindsight a different approach would have been used, likely more
 * consistent with the Boost reference counting pointers.
 *
 * @see Pointer
 */
class DescribedClass
{
public:
   /**
    * The default constructor for DescribedClass does nothing interesting.
    */
   DescribedClass();

   /**
    * The virtual destructor for DescribedClass does nothing interesting.
    */
   virtual ~DescribedClass();
};

}
}

#endif
