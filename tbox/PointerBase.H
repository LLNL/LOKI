/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef included_tbox_PointerBase
#define included_tbox_PointerBase

#include "ConstPointerBase.H"
#include "ReferenceCounter.H"

namespace Loki {
namespace tbox {

/**
 * Class PointerBase is a base class used by template class
 * Pointer<TYPE> for type-safe conversion between non-const
 * pointer types.  It is a subclass of ConstPointerBase.  Since
 * the non-const pointer class only takes this as a base class (and not
 * the const pointer base class), const pointers cannot be converted
 * into non-const pointers.
 *
 * @see ConstPointerBase
 * @see Pointer
 */
class PointerBase:public ConstPointerBase
{
public:
   PointerBase();

   virtual ~PointerBase();

   virtual ReferenceCounter *
   getSubclassReferenceCounter() const = 0;

   virtual const DescribedClass *
   getSubclassPointer() const = 0;
};

}
}

#endif
