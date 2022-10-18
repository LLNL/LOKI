/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef included_tbox_ConstPointer
#define included_tbox_ConstPointer

#include "ConstPointerBase.H"
#include "ReferenceCounter.H"

namespace Loki {
namespace tbox {

template<class TYPE>
class Pointer;

/**
 * Class ConstPointer<TYPE> defines a smart const pointer to TYPE.
 * It frees the user from explicitly deleting and tracking aliases for
 * object pointers.  It manages all reference counting and deallocation
 * of the pointer (even if the data was originally allocated from a memory
 * arena).  When the reference count on a ConstPointer<TYPE> object
 * goes to zero, the object is automatically deallocated.  A block with a
 * references count and arena pointer is allocated for all non-NULL pointers.
 * These reference counted blocks are freed at the end of the lifetime of
 * the pointer.
 *
 * Const pointers can be created from either const or non-const pointers.
 * The non-const and const pointer classes have been designed so that an
 * attempted conversion from a const pointer into a non-const pointer causes
 * a compile-time error.
 *
 * Class ConstPointer<TYPE> performs type-checking when assigning
 * pointers of different TYPEs.  If a bad type conversion is performed,
 * then the destination pointer is set to NULL.
 *
 * @see ConstPointerBase
 * @see ReferenceCounter
 * @see PointerBase
 * @see Pointer
 */
template<class TYPE>
class ConstPointer:public ConstPointerBase
{
public:
   /**
    * The default constructor creates a null pointer.
    */
   ConstPointer()
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
   ConstPointer(
      const TYPE* ptr,
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
    * The pointer copy constructor creates a smart pointer reference
    * aliased to the argument.
    */
   ConstPointer(
      const ConstPointer<TYPE>& ptr)
      : ConstPointerBase(),
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
   ConstPointer(
      const ConstPointerBase& ptr)
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
   ~ConstPointer()
   {
      if (d_counter && d_counter->deleteReference()) deleteObject();
   }

   /**
    * Smart pointer assignment.  The left hand side points to the
    * right hand side and the reference count is incremented by one.
    */
   ConstPointer<TYPE>&
   operator = (
      const ConstPointer<TYPE>& ptr)
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
   ConstPointer<TYPE>&
   operator = (
      const TYPE* ptr)
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
    * Attempt to convert the argument pointer to a ConstPointer<TYPE>.
    * If the type conversion fails, then the destination pointer will be set
    * to NULL.
    */
   ConstPointer<TYPE>&
   operator = (
      const ConstPointerBase& ptr)
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
      const ConstPointer<TYPE>& rhs) const
   {
      return d_object == rhs.d_object;
   }

   /**
    * Check whether two smart pointers point to different objects.
    */
   bool
   operator != (
      const ConstPointer<TYPE>& rhs) const
   {
      return d_object != rhs.d_object;
   }

   /**
    * Delegate member operations to the pointed-to object.  C++ defines
    * the ``->'' operator in a funny way to support delegation.  The
    * statement ptr->foo() acts as if ptr where actually a pointer
    * to an object with member function foo() instead of a class that
    * holds that pointer.
    */
   const TYPE *
   operator -> () const
   {
      return d_object;
   }

   /**
    * Dereference the smart pointer.  The pointer returned is a const
    * pointer to the object.
    */
   const TYPE&
   operator * () const
   {
      return *d_object;
   }

   /**
    * Implicit conversion of the smart pointer to the pointed-to object.
    * The pointer returned is a const pointer to the object.
    */
   operator const TYPE
   * () const
   {
      return d_object;
   }

   /**
    * Explicitly convert the smart pointer to the pointed-to object.
    * The pointer returned is a const pointer to the object.
    */
   const TYPE *
   getPointer() const
   {
      return d_object;
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
   void
   deleteObject()
   {
      delete d_object;
      delete d_counter;

      d_object = 0;
      d_counter = 0;
   }

   ReferenceCounter *
   getSubclassReferenceCounter() const
   {
      return d_counter;
   }

   const DescribedClass *
   getSubclassPointer() const
   {
      return (const DescribedClass*)d_object;
   }

   const TYPE* d_object;
   ReferenceCounter* d_counter;
};

}
}

#endif
