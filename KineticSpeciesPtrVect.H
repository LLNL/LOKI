/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _KINETIC_SPECIES_PTR_VECT_H_
#define _KINETIC_SPECIES_PTR_VECT_H_

#include <list>
#include "tbox/Pointer.H"
#include "KineticSpecies.H"

namespace Loki {

typedef tbox::Pointer<KineticSpecies> KineticSpeciesPtr;

/**
 * A container for Pointer<KineticSpecies> and the means to iterate over the
 * species in the container that are distributed to the processor.
 */
class KineticSpeciesPtrVect : public vector<KineticSpeciesPtr>
{
public:
   /**
    * @brief Default constructor.
    */
   KineticSpeciesPtrVect();

   /**
    * @brief Copy constructor.
    *
    * @param[in] other The thing being copied.
    */
   KineticSpeciesPtrVect(
      const KineticSpeciesPtrVect& other);

   /**
    * @brief Destructor.
    */
   ~KineticSpeciesPtrVect();

   /**
    * @brief Iterator over the species in a KineticSpeciesPtrVect that are
    *        distributed to the processor.
    */
   class Iterator
   {
   public:
      /**
       * @brief Copy constructor.
       *
       * @param[in] a_other The thing to be copied.
       */
      Iterator(
         const Iterator& a_other);

      /**
       * @brief Assignment operator.
       *
       * @param[in] a_rhs Right hand side of assignment.
       */
      Iterator&
      operator = (
         const Iterator& a_rhs)
      {
         // If this and a_rhs are different, copy their internals.
         if (this != &a_rhs) {
            m_local_indicies = a_rhs.m_local_indicies;
            m_current = a_rhs.m_current;
            m_from_end = a_rhs.m_from_end;
            m_base = a_rhs.m_base;
         }
         return *this;
      }

      /**
       * @brief Reset iteration from its start.
       */
      void
      reset()
      {
         // Reset the current iteration index to the proper end of
         // m_local_indices.
         if (m_from_end) {
            m_current = m_local_indicies.end();
         }
         else {
            m_current = m_local_indicies.begin();
         }
      }

      /**
       * @brief Pre-increment operator.
       */
      Iterator&
      operator ++ ()
      {
         ++m_current;
         return *this;
      }

      /**
       * @brief Post-increment operator.
       */
      Iterator
      operator ++ (int)
      {
         Iterator tmp(*this);
         ++m_current;
         return tmp;
      }

      /**
       * @brief Dereference operator.
       */
      KineticSpeciesPtr&
      operator * () const
      {
         return (*m_base)[*m_current];
      }

      /**
       * @brief Equality operator.
       *
       * @param[in] a_other Iterator equality is being checked against.
       */
      bool
      operator == (
         const Iterator& a_other) const
      {
         bool same_bases = m_base == a_other.m_base;
         bool this_at_end = m_current == m_local_indicies.end();
         bool other_at_end =
            a_other.m_current == a_other.m_local_indicies.end();
         bool both_at_end = this_at_end && other_at_end;
         bool neither_at_end = !this_at_end && !other_at_end;
         return (same_bases &&
            (both_at_end || (neither_at_end && (*(*this) == *a_other))));
      }

      /**
       * @brief Inequality operator.
       *
       * @param[in] a_other Iterator inequality is being checked against.
       */
      bool
      operator != (
         const Iterator& a_other) const
      {
         bool same_bases = m_base == a_other.m_base;
         bool this_at_end = m_current == m_local_indicies.end();
         bool other_at_end =
            a_other.m_current == a_other.m_local_indicies.end();
         bool both_at_end = this_at_end && other_at_end;
         bool neither_at_end = !this_at_end && !other_at_end;
         return (!same_bases ||
            (!both_at_end && (!neither_at_end || (*(*this) != *a_other))));
      }

   private:
      friend class KineticSpeciesPtrVect;

      // Unimplemented default constructor.
      Iterator();

      // Constructor only accessed by KineticSpeciesPtrVect.
      Iterator(
         KineticSpeciesPtrVect& a_in,
         bool a_from_end);

      // The indices of the KineticSpecies local to this processor.
      list<int> m_local_indicies;

      // The fundamental iterator.
      list<int>::iterator m_current;

      // If true, iterate from the end.
      bool m_from_end;

      // Pointer to the entity being iterated over.
      KineticSpeciesPtrVect* m_base;
   };

   /**
    * @brief Const iterator over the species in a KineticSpeciesPtrVect that
    *        are distributed to the processor.
    */
   class ConstIterator
   {
   public:
      /**
       * @brief Copy constructor.
       *
       * @param[in] a_other The thing to be copied.
       */
      ConstIterator(
         const ConstIterator& a_other);

      /**
       * @brief Assignment operator.
       *
       * @param[in] a_rhs Right hand side of assignment.
       */
      ConstIterator&
      operator = (
         const ConstIterator& a_rhs)
      {
         // If this and a_rhs are different, copy their internals.
         if (this != &a_rhs) {
            m_local_indicies = a_rhs.m_local_indicies;
            m_current = a_rhs.m_current;
            m_from_end = a_rhs.m_from_end;
            m_base = a_rhs.m_base;
         }
         return *this;
      }

      /**
       * @brief Reset iteration from its start.
       */
      void
      reset()
      {
         // Reset the current iteration index to the proper end of
         // m_local_indices.
         if (m_from_end) {
            m_current = m_local_indicies.end();
         }
         else {
            m_current = m_local_indicies.begin();
         }
      }

      /**
       * @brief Pre-increment operator.
       */
      ConstIterator&
      operator ++ ()
      {
         ++m_current;
         return *this;
      }

      /**
       * @brief Post-increment operator.
       */
      ConstIterator
      operator ++ (int)
      {
         ConstIterator tmp(*this);
         ++m_current;
         return tmp;
      }

      /**
       * @brief Dereference operator.
       */
      const KineticSpeciesPtr&
      operator * () const
      {
         return (*m_base)[*m_current];
      }

      /**
       * @brief Equality operator.
       *
       * @param[in] a_other ConstIterator equality is being checked against.
       */
      bool
      operator == (
         const ConstIterator& a_other) const
      {
         bool same_bases = m_base == a_other.m_base;
         bool this_at_end = m_current == m_local_indicies.end();
         bool other_at_end =
            a_other.m_current == a_other.m_local_indicies.end();
         bool both_at_end = this_at_end && other_at_end;
         bool neither_at_end = !this_at_end && !other_at_end;
         return (same_bases &&
            (both_at_end || (neither_at_end && (*(*this) == *a_other))));
      }

      /**
       * @brief Inequality operator.
       *
       * @param[in] a_other ConstIterator inequality is being checked against.
       */
      bool
      operator != (
         const ConstIterator& a_other) const
      {
         bool same_bases = m_base == a_other.m_base;
         bool this_at_end = m_current == m_local_indicies.end();
         bool other_at_end =
            a_other.m_current == a_other.m_local_indicies.end();
         bool both_at_end = this_at_end && other_at_end;
         bool neither_at_end = !this_at_end && !other_at_end;
         return (!same_bases ||
            (!both_at_end && (!neither_at_end || (*(*this) != *a_other))));
      }

   private:
      friend class KineticSpeciesPtrVect;

      // Unimplemented default constructor.
      ConstIterator();

      // Constructor only accessed by KineticSpeciesPtrVect.
      ConstIterator(
         const KineticSpeciesPtrVect& a_in,
         bool a_from_end);

      // The indices of the KineticSpecies local to this processor.
      list<int> m_local_indicies;

      // The fundamental iterator.
      list<int>::iterator m_current;

      // If true, iterate from the end.
      bool m_from_end;

      // Pointer to the entity being iterated over.
      const KineticSpeciesPtrVect* m_base;
   };

   /**
    * @brief Returns an iterator pointing to the beginning.
    */
   Iterator
   begin_locals()
   {
      return Iterator(*this, false);
   }

   /**
    * @brief Returns a const-iterator pointing to the beginning.
    */
   ConstIterator
   begin_locals() const
   {
      return ConstIterator(*this, false);
   }

   /**
    * @brief Returns an iterator pointing to the end.
    */
   Iterator
   end_locals()
   {
      return Iterator(*this, true);
   }

   /**
    * @brief Returns a const-iterator pointing to the end.
    */
   ConstIterator
   end_locals() const
   {
      return ConstIterator(*this, true);
   }

   /**
    * @brief Search the species to find the maximal bounding box.
    */
   void
   FindSpeciesVelocityBoundingBox(
     double* vlo,
     double* vhi) const;

private:
   friend Iterator;
   friend ConstIterator;

   // Unimplemented assignment operator.
   KineticSpeciesPtrVect&
   operator = (
      const KineticSpeciesPtrVect& rhs);
};

typedef KineticSpeciesPtrVect KSPV;
  
} // end namespace Loki

#endif
