/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _ELECTRON_TZ_SOURCE_FORT_H_
#define _ELECTRON_TZ_SOURCE_FORT_H_

#define FORT_SET_ELECTRON_TRIG_TZ_SOURCE            setelectrontrigtzsource_
#define FORT_COMPUTE_ELECTRON_TRIG_TZ_SOURCE_ERROR  computeelectrontrigtzsourceerror_

/**
 * Prototypes for fortran functions called by twilight zone sources.
 */
extern "C" {
   void
   FORT_SET_ELECTRON_TRIG_TZ_SOURCE(
      double &f,
      const int& nd1a,
      const int& nd1b,
      const int& nd2a,
      const int& nd2b,
      const int& nd3a,
      const int& nd3b,
      const int& nd4a,
      const int& nd4b,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const double& time,
      const double& velocities,
      const double& dparams);

   void
   FORT_COMPUTE_ELECTRON_TRIG_TZ_SOURCE_ERROR(
      double &error,
      const double& soln,
      const int& nd1a,
      const int& nd1b,
      const int& nd2a,
      const int& nd2b,
      const int& nd3a,
      const int& nd3b,
      const int& nd4a,
      const int& nd4b,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const double& time,
      const double& velocities,
      const double& dparams);
}

#endif
