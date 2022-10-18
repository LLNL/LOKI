/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef included_tbox_ConstPointerBase
#define included_tbox_ConstPointerBase

#include "ReferenceCounter.H"
#include "DescribedClass.H"

namespace Loki {
namespace tbox {

/**
 * Class ConstPointerBase is an abstract base class used by template
 * class ConstPointer<TYPE> for type-safe conversion between various
 * pointer types.  It forms the base of the RTTI conversion hierarchy for
 * pointers.  Both the non-const pointer base class and the const pointer
 * class are subclasses of the const pointer base class.  This structure
 * ensures that RTTI conversion of pointers from const to const, non-const
 * to non-const, and non-const to const work as expected but that conversion
 * from const to non-const fail at compile-time.
 *
 * @see ConstPointer
 * @see PointerBase
 * @see Pointer
 */
class ConstPointerBase
{
public:
   ConstPointerBase();

   virtual ~ConstPointerBase();

   virtual ReferenceCounter *
   getSubclassReferenceCounter() const = 0;

   virtual const DescribedClass *
   getSubclassPointer() const = 0;
};

}
}

#endif
