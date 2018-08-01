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
