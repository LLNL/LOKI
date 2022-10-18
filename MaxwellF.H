/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _MAXWELL_FORT_H_
#define _MAXWELL_FORT_H_

#define FORT_ZERO_GHOST_2D                   zeroghost2d_
#define FORT_XPBY_2D                         xpby2d_
#define FORT_MAXWELL_EVAL_RHS                maxwellevalrhs_
#define FORT_MAXWELL_ADD_ANTENNA_SOURCE      maxwelladdantennasource_
#define FORT_MAXWELL_EVAL_VZ_RHS             maxwellevalvzrhs_
#define FORT_SET_MAXWELL_EM_BCS              maxwellsetembcs_
#define FORT_SET_MAXWELL_VZ_BCS              maxwellsetvzbcs_

/**
 * Prototypes for fortran functions called by Maxwell.
 */
extern "C" {
   void
   FORT_ZERO_GHOST_2D(
      double &u,
      const int& n1a,
      const int& n1b,
      const int& n2a,
      const int& n2b,
      const int& nd1a,
      const int& nd1b,
      const int& nd2a,
      const int& nd2b,
      const int& dim);

   void
   FORT_XPBY_2D(
      double &x,
      const double& y,
      const double& b,
      const int& nd1a,
      const int& nd1b,
      const int& nd2a,
      const int& nd2b,
      const int& n1a,
      const int& n1b,
      const int& n2a,
      const int& n2b,
      const int& dim);

   void
   FORT_MAXWELL_EVAL_RHS(
      const int& md1a,
      const int& md1b,
      const int& md2a,
      const int& md2b,
      const int& m1a,
      const int& m1b,
      const int& m2a,
      const int& m2b,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const double& light_speed,
      const double& avWeak,
      const double& avStrong,
      const int& solution_order,
      const double& supergrid_lo,
      const double& supergrid_hi,
      const double& EMvars,
      const double& Jx,
      const double& Jy,
      const double& Jz,
      double& dEMvars);

   void
   FORT_MAXWELL_ADD_ANTENNA_SOURCE(
      const int& md1a,
      const int& md1b,
      const int& md2a,
      const int& md2b,
      const int& m1a,
      const int& m1b,
      const int& m2a,
      const int& m2b,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const double& antenna_source,
      double& dEMvars);

   void
   FORT_MAXWELL_EVAL_VZ_RHS(
      const int& md1a,
      const int& md1b,
      const int& md2a,
      const int& md2b,
      const int& m1a,
      const int& m1b,
      const int& m2a,
      const int& m2b,
      const double& charge_per_mass,
      const double& EMvars,
      double& dvz);

   void
   FORT_SET_MAXWELL_EM_BCS(
      const int& md1a,
      const int& md1b,
      const int& md2a,
      const int& md2b,
      const int& m1a,
      const int& m1b,
      const int& m2a,
      const int& m2b,
      double& EMvars,
      const int& nx,
      const int& ny,
      const int& xPeriodic,
      const int& yPeriodic,
      const int& solution_order,
      const double& light_speed);

   void
   FORT_SET_MAXWELL_VZ_BCS(
      const int& md1a,
      const int& md1b,
      const int& md2a,
      const int& md2b,
      const int& m1a,
      const int& m1b,
      const int& m2a,
      const int& m2b,
      double& vz,
      const int& nx,
      const int& ny,
      const int& xPeriodic,
      const int& yPeriodic,
      const int& solution_order);
};

#endif
