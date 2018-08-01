/*************************************************************************
 *
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * Written by Jeffrey Banks banksj3@rpi.edu (Rensselaer Polytechnic Institute,
 * Amos Eaton 301, 110 8th St., Troy, NY 12180); Jeffrey Hittinger
 * hittinger1@llnl.gov, William Arrighi arrighi2@llnl.gov, Richard Berger
 * berger5@llnl.gov, Thomas Chapman chapman29@llnl.gov (LLNL, P.O Box 808,
 * Livermore, CA 94551); Stephan Brunner stephan.brunner@epfl.ch (Ecole
 * Polytechnique Federale de Lausanne, EPFL SB SPC-TH, PPB 312, Station 13,
 * CH-1015 Lausanne, Switzerland).
 * CODE-744849
 *
 * All rights reserved.
 *
 * This file is part of Loki.  For details, see.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
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
