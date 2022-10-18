/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _UNIDIRECTIONAL_CURRENT_DRIVER_FORT_H_
#define _UNIDIRECTIONAL_CURRENT_DRIVER_FORT_H_

#define evaluateUnidirectionalCurrentDriver evaluateunidirectionalcurrentdriver_

/**
 * Prototypes of fortran functions called by UnidirectionalCurrentDriver.
 */
extern "C" {
   void
   evaluateUnidirectionalCurrentDriver(
      double& antenna_source,
      int& status,
      const double& omega_eff2,
      const int& nd1a,
      const int& nd1b,
      const int& nd2a,
      const int& nd2b,
      const int& n1a,
      const int& n1b,
      const int& n2a,
      const int& n2b,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const double& t,
      const double& e_ext_param,
      const double& R);
}

#endif
