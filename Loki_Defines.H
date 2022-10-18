/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _LOKI_DEFINES_H
#define _LOKI_DEFINES_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

/*!
 * A null use of a variable, use to avoid GNU compiler
 * warnings about unused variables.
 */
#define NULL_USE(variable)                               \
   do {                                                  \
      if (0) { char* temp = (char *)&variable; temp++; } \
   } while (0)

/**
 * Abort the program when an unrecoverable error is detected with information
 * about the file, function, and line where the error occurred.
 */
#define LOKI_ABORT(message) \
{ \
   printf("Error occured in file %s, function %s, line %d.\n", __FILE__, __func__, __LINE__); \
   cout << message << endl; \
   printf("Loki::abort: I am now going to purposely abort so that you can get a traceback from a debugger\n"); \
   abort(); \
}

#endif
