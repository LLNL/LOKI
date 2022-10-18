/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _COLLISION_OPERATOR_FORT_H_
#define _COLLISION_OPERATOR_FORT_H_

#define FORT_COMPUTE_DIST_FUNC_AVG computedistfuncavg_
#define FORT_COMPUTE_VTH_LOCAL     computevthlocal_
#define FORT_COMPUTE_VTH_GLOBAL    computevthglobal_

/**
 * Prototypes of fortran functions called by CollisionOperator.
 */
extern "C" {
   void
   FORT_COMPUTE_DIST_FUNC_AVG(
      const int& nd1a,
      const int& nd1b,
      const int& nd2a,
      const int& nd2b,
      const int& nd3a,
      const int& nd3b,
      const int& nd4a,
      const int& nd4b,
      const int& n1a,
      const int& n1b,
      const int& n2a,
      const int& n2b,
      const int& n3a,
      const int& n3b,
      const int& n4a,
      const int& n4b,
      const int& vel_space_size,
      const int& config_space_size,
      const double& u,
      double& uxybar);

   void
   FORT_COMPUTE_VTH_LOCAL(
      const int& nd1a,
      const int& nd1b,
      const int& nd2a,
      const int& nd2b,
      const int& nd3a,
      const int& nd3b,
      const int& nd4a,
      const int& nd4b,
      const int& n1a,
      const int& n1b,
      const int& n2a,
      const int& n2b,
      const int& n3a,
      const int& n3b,
      const int& n4a,
      const int& n4b,
      const double& vflowx,
      const double& vflowy,
      const double& velocities,
      const int& comm,
      const double& u,
      double& vth);

   void
   FORT_COMPUTE_VTH_GLOBAL(
      const int& nd1a,
      const int& nd1b,
      const int& nd2a,
      const int& nd2b,
      const int& nd3a,
      const int& nd3b,
      const int& nd4a,
      const int& nd4b,
      const int& n1a,
      const int& n1b,
      const int& n2a,
      const int& n2b,
      const int& n3a,
      const int& n3b,
      const int& n4a,
      const int& n4b,
      const int& vel_space_size,
      const double& vflowx,
      const double& vflowy,
      const double& velocities,
      const double& uxybar,
      const int& is_config_space_head_node,
      const int& comm,
      double& vth);
};

#endif
