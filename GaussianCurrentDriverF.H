/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _GAUSSIAN_CURRENT_DRIVER_FORT_H_
#define _GAUSSIAN_CURRENT_DRIVER_FORT_H_

#define evaluateGaussianCurrentDriver evaluategaussiancurrentdriver_

/**
 * Prototypes of fortran functions called by GaussianCurrentDriver.
 */
extern "C" {
   void
   evaluateGaussianCurrentDriver(
      double &antenna_source,
      const int &nd1a,
      const int &nd1b,
      const int &nd2a,
      const int &nd2b,
      const int &n1a,
      const int &n1b,
      const int &n2a,
      const int &n2b,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const double &t,
      const int &shape_lo,
      const int &shape_hi,
      const double &e_ext_param,
      const double &shaping,
      const double &phase,
      const int &antenna_loc_idx);
}

#endif
