/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _SHAPED_RAMPED_COSINE_POTENTIAL_DRIVER_FORT_H_
#define _SHAPED_RAMPED_COSINE_POTENTIAL_DRIVER_FORT_H_

#define evaluateShapedRampedPotentialDriver evaluateshapedrampedpotentialdriver_

/**
 * Prototypes for fortran functions called by
 * ShapedRampedCosineElectricPotentialDriver.
 */
extern "C" {
   void
   evaluateShapedRampedPotentialDriver(
      double &phi,
      const int &nd1a,
      const int &nd1b,
      const int &nd2a,
      const int &nd2b,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const double &t,
      const double &e_ext_param);
}

#endif
