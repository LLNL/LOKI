/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _GETMOMENTUM_FORT_H_
#define _GETMOMENTUM_FORT_H_

#define FORT_GET_MOMENTUM getmomentum_

/**
 * Prototypes of fortran function to compute relativistic momentum.
 */
extern "C" {
   void
   FORT_GET_MOMENTUM(
      const double &vx,
      const double &vy,
      const double &mass,
      const double &c,
      double & px,
      double & py);

};

#endif
