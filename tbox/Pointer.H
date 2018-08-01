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
#ifndef included_tbox_Pointer
#define included_tbox_Pointer

#include "ReferenceCounter.H"
#include "PointerBase.H"

#include <typeinfo>

namespace Loki {
namespace tbox {

/**
 * Class Pointer<TYPE> defines a smart pointer to TYPE.  It frees the
 * user from explicitly deleting and tracking aliases for object pointers.
 * It manages all reference counting and deallocation of the pointer (even
 * if the data was originally allocated from a memory arena).  When the
 * reference count on a Pointer<TYPE> object goes to zero, the object
 * is automatically deallocated.  A block with a references count and arena
 * pointer is allocated for all non-NULL pointers.  These reference counted
 * blocks are freed at the end of the lifetime of the pointer.
 *
 * Non-const pointers can be created only from other non-const pointers.
 * The non-const and const pointer classes have been designed so that an
 * attempted conversion from a const pointer into a non-const pointer causes
 * a compile-time error.  The const convention for this class is that any
 * member function that does not change the pointer object itself may be
 * const.  Note that this means that a const pointer does not mean that the
 * pointed-to object is const; these semantics require the const pointer class.
 *
 * Class Pointer<TYPE> performs type-checking when assigning pointers
 * of different TYPEs.  If a bad type conversion is performed, then the
 * destination pointer is set to NULL.
 *
 * @see DescribedClass
 * @see ConstPointer
 * @see PointerBase
 * @see ReferenceCounter
 */
template<class TYPE>
class Pointer:public PointerBase
{
public:
   /**
    * The default constructor creates a null pointer.
    */
   Pointer()
      : d_object(0),
        d_counter(0)
   {
   }

   /**
    * Create a smart pointer with value ptr.  If managed is true, then
    * deallocation of the object pointed to by ptr will be taken care of
    * by the smart pointer.  This form assumes the pointer was allocated
    * using the standard new operator.
    */
   explicit Pointer(
      TYPE* ptr,
      const bool managed = true)
      : d_object(ptr)
   {
      if (d_object && managed) {
         d_counter = new ReferenceCounter;
      } else {
         d_counter = 0;
      }
   }

   /**
    * The pointer const constructor creates a smart pointer reference
    * aliased to the argument.
    */
   Pointer(
      const Pointer<TYPE>& ptr)
      : PointerBase(),
        d_object(ptr.d_object),
        d_counter(ptr.d_counter)
   {
      if (d_counter) d_counter->addReference();
   }

   /**
    * Create a pointer by attempting to type-cast the argument to TYPE.
    * If the type-cast fails, then the destination pointer will be set
    * to NULL.
    */
   Pointer(
      const PointerBase& ptr)
   {
      const DescribedClass* sub_ptr = ptr.getSubclassPointer();
      if (sub_ptr) {
         d_object = (TYPE *)dynamic_cast<const TYPE *>(sub_ptr);
      } else {
         d_object = 0;
      }

      if (d_object) {
         d_counter = ptr.getSubclassReferenceCounter();
         if (d_counter) d_counter->addReference();
      } else {
         d_counter = 0;
      }
   }

   /**
    * The pointer destructor frees the pointer data if the reference
    * count drops to zero.  The object is deallocated from the memory
    * pool (if it was specified in the constructor call).
    */
   ~Pointer()
   {
      if (d_counter && d_counter->deleteReference()) deleteObject();
   }

   /**
    * Smart pointer assignment.  The left hand side points to the
    * right hand side and the reference count is incremented by one.
    */
   Pointer<TYPE>&
   operator = (
      const Pointer<TYPE>& ptr)
   {
      if (this != &ptr) {
         if (d_counter && d_counter->deleteReference()) deleteObject();
         d_object = ptr.d_object;
         d_counter = ptr.d_counter;
         if (d_counter) d_counter->addReference();
      }
      return *this;
   }

   /**
    * Create a managed smart pointer with value ptr.  The object pointed
    * to by ptr will be deallocated via delete when the reference count
    * goes to zero.
    */
   Pointer<TYPE>&
   operator = (
      TYPE* ptr)
   {
      if (d_counter && d_counter->deleteReference()) deleteObject();
      d_object = ptr;
      if (d_object) {
         d_counter = new ReferenceCounter;
      } else {
         d_counter = 0;
      }
      return *this;
   }

   /**
    * Attempt to convert the argument pointer to a Pointer<TYPE>.
    * If the type conversion fails, then the destination pointer will be
    * set to NULL.
    */
   Pointer<TYPE>&
   operator = (
      const PointerBase& ptr)
   {
      if (this != &ptr) {
         if (d_counter && d_counter->deleteReference()) deleteObject();

         const DescribedClass* sub_ptr = ptr.getSubclassPointer();
         if (sub_ptr) {
            d_object = (TYPE *)dynamic_cast<const TYPE *>(sub_ptr);
         } else {
            d_object = 0;
         }

         if (d_object) {
            d_counter = ptr.getSubclassReferenceCounter();
            if (d_counter) d_counter->addReference();
         } else {
            d_counter = 0;
         }
      }
      return *this;
   }

   /**
    * Check whether two smart pointers point to the same object.
    */
   bool
   operator == (
      const Pointer<TYPE>& rhs) const
   {
      return d_object == rhs.d_object;
   }

   /**
    * Check whether two smart pointers point to different objects.
    */
   bool
   operator != (
      const Pointer<TYPE>& rhs) const
   {
      return d_object != rhs.d_object;
   }

   /**
    * Delegate member operations to the pointed-to object.  C++ defines
    * the ``->'' operator in a funny way to support delegation.  The
    * statement ptr->foo() acts as if ptr where actually a pointer
    * to an object with member function foo() instead of a class that
    * holds that pointer.  This member function is const since it cannot
    * change the pointer, although the pointed-to object may change.
    */
   TYPE *
   operator -> () const
   {
      return d_object;
   }

   /**
    * Dereference the smart pointer.  This member function is const since
    * it cannot change the pointer, although the pointed-to object may
    * change.
    */
   TYPE&
   operator * () const
   {
      return *d_object;
   }

   /**
    * Implicit conversion of the smart pointer to the pointed-to object.
    * This conversion operator is const since it cannot change the pointer,
    * although the pointed-to object may change.
    */
   operator TYPE
   * () const
   {
      return d_object;
   }

   /**
    * Explicitly convert the smart pointer to the pointed-to object.
    * This member function is const since it cannot change the pointer,
    * although the pointed-to object may change.
    */
   TYPE *
   getPointer() const
   {
      return d_object;
   }

   /**
    * Set the pointed-to object to ptr.
    *
    * NOTE: reset name comes from boost library.
    */
   void
   reset(
      TYPE* ptr)
   {
      if (d_counter && d_counter->deleteReference()) deleteObject();
      d_object = ptr;
      if (d_object) {
         d_counter = new ReferenceCounter;
      } else {
         d_counter = 0;
      }
   }

   /**
    * Check whether the smart pointer points to NULL.
    */
   bool
   isNull() const
   {
      return !d_object;
   }

   /**
    * Return true if the pointer is non-NULL.
    */
   operator bool () const
   {
     return d_object ? true : false;
   }

   /**
    * Return true if the pointer is NULL and false otherwise.  This operator
    * mimics the semantics of !p applied to a (regular) pointer p.
    */
   bool
   operator ! () const
   {
      return !d_object;
   }

   /**
    * Set the smart pointer to NULL.
    */
   void
   setNull()
   {
      if (d_counter && d_counter->deleteReference()) deleteObject();
      d_object = 0;
      d_counter = 0;
   }

   /**
    * Return a pointer to the internal reference counter.  This routine
    * should not be called by the casual user.
    */
   ReferenceCounter *
   getReferenceCounter() const
   {
      return d_counter;
   }

private:
   /*!
    * @brief Unused conversion, meant to prevent unintended conversion
    * of Pointer to int, which can cause serious hidden bugs.
    */
   operator int () const;

   void
   deleteObject()
   {
      delete d_object;
      delete d_counter;

#ifdef DEBUG_INITIALIZE_UNDEFINED
      // For a little performance gain only do this when debugging.
      d_object = 0;
      d_counter = 0;
#endif
   }

   ReferenceCounter *
   getSubclassReferenceCounter() const
   {
      return (ReferenceCounter *)d_counter;
   }

   const DescribedClass *
   getSubclassPointer() const
   {
      // the explicit cast is needed by the brain-damaged SGI C++ compiler
      return (DescribedClass *)d_object;
   }

   TYPE* d_object;
   ReferenceCounter* d_counter;
};

}
}

#endif
