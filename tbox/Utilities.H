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
#include <sstream>
#include "Overture.h"

/**
 * Description:   Useful dimension related macros.
 */
#define TBOX_DIM_ASSERT_CHECK_DIM(dim) \
   TBOX_DIM_ASSERT(                    \
   (dim).isValid()                     \
   )

#define TBOX_DIM_ASSERT_CHECK_DIM_ALLOW_UNINITIALIZED(dim) \
   TBOX_DIM_ASSERT(                                        \
   !(dim).isInitialized || (dim).isValid()                 \
   )

#define TBOX_DIM_ASSERT_CHECK_ARGS1(arg1) \
   TBOX_DIM_ASSERT(                       \
   (arg1).getDim().isValid()              \
   )

#define TBOX_DIM_ASSERT_CHECK_ARGS2(arg1, \
                                    arg2) \
   TBOX_DIM_ASSERT(                       \
   (arg1).getDim().isValid() &&           \
   ((arg1).getDim() == (arg2).getDim())   \
   )

#define TBOX_DIM_ASSERT_CHECK_ARGS3(arg1,  \
                                    arg2,  \
                                    arg3)  \
   TBOX_DIM_ASSERT(                        \
   (arg1).getDim().isValid() &&            \
   ((arg1).getDim() == (arg2).getDim()) && \
   ((arg1).getDim() == (arg3).getDim())    \
   )

#define TBOX_DIM_ASSERT_CHECK_ARGS4(arg1,  \
                                    arg2,  \
                                    arg3,  \
                                    arg4)  \
   TBOX_DIM_ASSERT(                        \
   (arg1).getDim().isValid() &&            \
   ((arg1).getDim() == (arg2).getDim()) && \
   ((arg1).getDim() == (arg3).getDim()) && \
   ((arg1).getDim() == (arg4).getDim())    \
   )

#define TBOX_DIM_ASSERT_CHECK_ARGS5(arg1,  \
                                    arg2,  \
                                    arg3,  \
                                    arg4,  \
                                    arg5)  \
   TBOX_DIM_ASSERT(                        \
   (arg1).getDim().isValid() &&            \
   ((arg1).getDim() == (arg2).getDim()) && \
   ((arg1).getDim() == (arg3).getDim()) && \
   ((arg1).getDim() == (arg4).getDim()) && \
   ((arg1).getDim() == (arg5).getDim())    \
   )

#define TBOX_DIM_ASSERT_CHECK_ARGS6(arg1,  \
                                    arg2,  \
                                    arg3,  \
                                    arg4,  \
                                    arg5,  \
                                    arg6)  \
   TBOX_DIM_ASSERT(                        \
   (arg1).getDim().isValid() &&            \
   ((arg1).getDim() == (arg2).getDim()) && \
   ((arg1).getDim() == (arg3).getDim()) && \
   ((arg1).getDim() == (arg4).getDim()) && \
   ((arg1).getDim() == (arg5).getDim()) && \
   ((arg1).getDim() == (arg6).getDim())    \
   )

#define TBOX_DIM_ASSERT_CHECK_ARGS7(arg1,  \
                                    arg2,  \
                                    arg3,  \
                                    arg4,  \
                                    arg5,  \
                                    arg6,  \
                                    arg7)  \
   TBOX_DIM_ASSERT(                        \
   (arg1).getDim().isValid() &&            \
   ((arg1).getDim() == (arg2).getDim()) && \
   ((arg1).getDim() == (arg3).getDim()) && \
   ((arg1).getDim() == (arg4).getDim()) && \
   ((arg1).getDim() == (arg5).getDim()) && \
   ((arg1).getDim() == (arg6).getDim()) && \
   ((arg1).getDim() == (arg7).getDim())    \
   )

#define TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(dim,  \
                                        arg1) \
   TBOX_DIM_ASSERT(                           \
   (dim).isValid() &&                         \
   ((dim) == (arg1).getDim())                 \
   )                                          \

#define TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(dim,  \
                                        arg1, \
                                        arg2) \
   TBOX_DIM_ASSERT(                           \
   (dim).isValid() &&                         \
   ((dim) == (arg1).getDim()) &&              \
   ((dim) == (arg2).getDim())                 \
   )

#define TBOX_DIM_ASSERT_CHECK_DIM_ARGS3(dim,  \
                                        arg1, \
                                        arg2, \
                                        arg3) \
   TBOX_DIM_ASSERT(                           \
   (dim).isValid() &&                         \
   ((dim) == (arg1).getDim()) &&              \
   ((dim) == (arg2).getDim()) &&              \
   ((dim) == (arg3).getDim())                 \
   )

#define TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dim,  \
                                        arg1, \
                                        arg2, \
                                        arg3, \
                                        arg4) \
   TBOX_DIM_ASSERT(                           \
   (dim).isValid() &&                         \
   ((dim) == (arg1).getDim()) &&              \
   ((dim) == (arg2).getDim()) &&              \
   ((dim) == (arg3).getDim()) &&              \
   ((dim) == (arg4).getDim())                 \
   )

#define TBOX_DIM_ASSERT_CHECK_DIM_ARGS5(dim,  \
                                        arg1, \
                                        arg2, \
                                        arg3, \
                                        arg4, \
                                        arg5) \
   TBOX_DIM_ASSERT(                           \
   (dim).isValid() &&                         \
   ((dim) == (arg1).getDim()) &&              \
   ((dim) == (arg2).getDim()) &&              \
   ((dim) == (arg3).getDim()) &&              \
   ((dim) == (arg4).getDim()) &&              \
   ((dim) == (arg5).getDim())                 \
   )

#define TBOX_DIM_ASSERT_CHECK_DIM_ARGS6(dim,  \
                                        arg1, \
                                        arg2, \
                                        arg3, \
                                        arg4, \
                                        arg5, \
                                        arg6) \
   TBOX_DIM_ASSERT(                           \
   (dim).isValid() &&                         \
   ((dim) == (arg1).getDim()) &&              \
   ((dim) == (arg2).getDim()) &&              \
   ((dim) == (arg3).getDim()) &&              \
   ((dim) == (arg4).getDim()) &&              \
   ((dim) == (arg5).getDim()) &&              \
   ((dim) == (arg6).getDim())                 \
   )

#define TBOX_DIM_ASSERT_CHECK_DIM_ARGS7(dim,  \
                                        arg1, \
                                        arg2, \
                                        arg3, \
                                        arg4, \
                                        arg5, \
                                        arg6, \
                                        arg7) \
   TBOX_DIM_ASSERT(                           \
   (dim).isValid() &&                         \
   ((dim) == (arg1).getDim()) &&              \
   ((dim) == (arg2).getDim()) &&              \
   ((dim) == (arg3).getDim()) &&              \
   ((dim) == (arg4).getDim()) &&              \
   ((dim) == (arg5).getDim()) &&              \
   ((dim) == (arg6).getDim()) &&              \
   ((dim) == (arg7).getDim())                 \
   )

#define TBOX_DIM_ASSERT_CHECK_DIM_ARGS8(dim,  \
                                        arg1, \
                                        arg2, \
                                        arg3, \
                                        arg4, \
                                        arg5, \
                                        arg6, \
                                        arg7, \
                                        arg8) \
   TBOX_DIM_ASSERT(                           \
   (dim).isValid() &&                         \
   ((dim) == (arg1).getDim()) &&              \
   ((dim) == (arg2).getDim()) &&              \
   ((dim) == (arg3).getDim()) &&              \
   ((dim) == (arg4).getDim()) &&              \
   ((dim) == (arg5).getDim()) &&              \
   ((dim) == (arg6).getDim()) &&              \
   ((dim) == (arg7).getDim()) &&              \
   ((dim) == (arg8).getDim())                 \
   )

#ifdef DEBUG

#define TBOX_DIM_ASSERT(EXP)                                             \
   do {                                                                  \
      if (!(EXP)) {                                                      \
         std::ostringstream tboxos;                                      \
         tboxos << "Failed dimension assertion: " << # EXP << std::ends; \
         OV_ABORT(tboxos.str());                                         \
      }                                                                  \
   } while (0)

#else

/*
 * No dimensional assertion checking
 */
#define TBOX_DIM_ASSERT(EXP)

#endif

#ifdef DEBUG_CHECK_ASSERTIONS

#define TBOX_ASSERT(EXP)                                       \
   do {                                                        \
      if (!(EXP)) {                                            \
         std::ostringstream tboxos;                            \
         tboxos << "Failed assertion: " << # EXP << std::ends; \
         OV_ABORT(tboxos.str());                               \
      }                                                        \
   } while (0)
#else

/*
 * No assertion checking
 */
#define TBOX_ASSERT(EXP)

#endif
