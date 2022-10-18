/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef included_tbox_ReferenceCounter
#define included_tbox_ReferenceCounter

#include <cstddef>

namespace Loki {
namespace tbox {

/**
 * Class ReferenceCounter manages the shared reference counter and arena
 * resources used by Pointer.  It uses a local free pool of objects to speed
 * memory allocation and deallocation.  The locally cached free pool can be
 * freed by calling freeCachedCopies().
 *
 * {\b Do not subclass!}  Changing the size of a ReferenceCounter
 * object will cause my simple memory allocation mechanism to break in
 * horrible and disgusting ways.
 *
 * @see Pointer
 */
class ReferenceCounter
{
public:
   /**
    * Create a ReferenceCounter.
    * The number of references is set to one.
    */
   ReferenceCounter();

   /**
    * Destructor for ReferenceCounter.  The destructor releases
    * the managed memory arena if its count has gone to zero.
    */
   ~ReferenceCounter();

   /**
    * Decrement the number of references.  True is returned if the
    * reference count has gone to zero; false otherwise.
    */
   bool
   deleteReference()
   {
      return --d_references == 0;
   }

   /**
    * Increment the number of references.
    */
   void
   addReference()
   {
      ++d_references;
   }

   /**
    * Class-specific operator new.  Data is allocated off of an
    * internal free list to speed memory allocation.
    */
   void *
   operator new (
      size_t bytes);

   /**
    * Class-specific operator delete.  Freed data is returned to
    * an internal free list for re-use by operator new.
    */
   void
   operator delete (
      void* what);

   /**
    * Returns true if ReferenceCounter class has been be stopped.
    * This method is used only for a debugging check assert in
    * the pointer class and should not normally be used.
    */
   static bool isFinalized() {
      return s_is_finalized;
   };
  

   /**
    * Release the memory for all currently cached ReferenceCounter
    * copies.
    */
   static void
   finalizeCallback();

private:
   ReferenceCounter(
      const ReferenceCounter&);                 // not implemented
   void
   operator = (
      const ReferenceCounter&);                 // not implemented

   int d_references;

   // Next element on free list
   ReferenceCounter* d_next;

   // Free list of ReferenceCounter objects
   static ReferenceCounter* s_free_list;

   static bool s_is_finalized;
};

}
}

#endif
