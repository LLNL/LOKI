/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _SIMPLE_VEL_IC_FORT_H_
#define _SIMPLE_VEL_IC_FORT_H_

#define SET_SIMPLE_VELIC  setsimplevelic_

/**
 * Prototypes for fortran functions called by SimpleVELIC.
 */
extern "C" {
   void
   SET_SIMPLE_VELIC(
      double &vz,
      const int &nd1a,
      const int &nd1b,
      const int &nd2a,
      const int &nd2b,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const int& num_waves,
      const double& amp,
      const double& kx,
      const double& ky,
      const double& phi);
}

#endif
