/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _SHAPED_RAMPED_COSINE_DRIVER_FORT_H_
#define _SHAPED_RAMPED_COSINE_DRIVER_FORT_H_

#define evaluateShapedRampedDriver evaluateshapedrampeddriver_
#define evaluateShapedRampedDriverEnvelope evaluateshapedrampeddriverenvelope_

/**
 * Prototypes of fortran functions called by ShapedRampedCosineDriver.
 */
extern "C" {
   void
   evaluateShapedRampedDriver(
      double &em_vars,
      double &ext_efield,
      const int &nd1a,
      const int &nd1b,
      const int &nd2a,
      const int &nd2b,
      const int &num_em_vars,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const int &sums_into,
      const double &t,
      const double &e_ext_param,
      const double &phase,
      const int &shape_type);

   void
   evaluateShapedRampedDriverEnvelope(
      double &envel,
      const double &t,
      const double &t0,
      const double &t_rampup,
      const double &t_hold,
      const double &t_rampdown,
      const double &E_0);
}

#endif
