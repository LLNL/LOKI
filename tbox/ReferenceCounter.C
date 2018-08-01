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
#include "ReferenceCounter.H"
#include "Overture.h"

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
      OV_ABORT("Attempting to new a ReferenceCounter after finalizing");
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
      OV_ABORT("Attempting to delete a ReferenceCounter after finalizing");
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
