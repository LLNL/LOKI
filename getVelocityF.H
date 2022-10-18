/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _GETVELOCITY_FORT_H_
#define _GETVELOCITY_FORT_H_

#define FORT_GET_VELOCITY getvelocity_

/**
 * Prototypes of fortran function to compute relativistic velocity.
 */
extern "C" {
   void
   FORT_GET_VELOCITY(
      const double &px,
      const double &py,
      const double &mass,
      const double &c,
      double & vx,
      double & vy);

};

#endif
