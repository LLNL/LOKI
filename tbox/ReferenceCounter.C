/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "ReferenceCounter.H"
#include "../Loki_Defines.H"

#include <cstdlib>
#include <cstdio>

namespace Loki {
namespace tbox {

ReferenceCounter * ReferenceCounter::s_free_list = 0;
bool ReferenceCounter::s_is_finalized = false;


ReferenceCounter::ReferenceCounter()
   : d_references(1),
     d_next(0)
{
}

ReferenceCounter::~ReferenceCounter()
{
   if ((d_next) && (--d_next->d_references == 0)) {
      delete d_next;
   }
}

void *ReferenceCounter::operator new (
   size_t bytes)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
   /* Since low level class; tbox Utilities may not function here */
   if (ReferenceCounter::isFinalized()) {
      LOKI_ABORT("Attempting to new a ReferenceCounter after finalizing");
   }
#endif

   if (s_free_list) {
      ReferenceCounter* node = s_free_list;
      s_free_list = s_free_list->d_next;
      return node;
   } else {
      return ::operator new (
                bytes);
   }
}

void ReferenceCounter::operator delete (
   void* what)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
   /* Since low level class; tbox Utilities may not function here */
   if (ReferenceCounter::isFinalized()) {
      LOKI_ABORT("Attempting to delete a ReferenceCounter after finalizing");
   }
#endif

   ReferenceCounter* node = (ReferenceCounter *)what;
   node->d_next = s_free_list;
   s_free_list = node;
}

void ReferenceCounter::finalizeCallback()
{
   while (s_free_list) {
      void* byebye = s_free_list;
      s_free_list = s_free_list->d_next
      ;
      ::operator delete (
         byebye);
   }

   s_is_finalized = true;
}

}
}
