/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _SIMPLE_EM_IC_FORT_H_
#define _SIMPLE_EM_IC_FORT_H_

#define SET_SIMPLE_EMIC  setsimpleemic_

/**
 * Prototypes for fortran functions called by SimpleEMIC.
 */
extern "C" {
   void
   SET_SIMPLE_EMIC(
      const int &nd1a,
      const int &nd1b,
      const int &nd2a,
      const int &nd2b,
      double &u,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const int& num_waves,
      const double& xamp,
      const double& yamp,
      const double& zamp,
      const double& kx,
      const double& ky,
      const double& phi,
      const int& iparams);
}

#endif
