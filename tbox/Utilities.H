/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include <sstream>
#include "../Loki_Defines.H"

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
         LOKI_ABORT(tboxos.str());                                         \
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
         LOKI_ABORT(tboxos.str());                               \
      }                                                        \
   } while (0)
#else

/*
 * No assertion checking
 */
#define TBOX_ASSERT(EXP)

#endif
