/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _ORIGINAL_ROSENBLUTH_COLLISION_OPERATOR_FORT_H_
#define _ORIGINAL_ROSENBLUTH_COLLISION_OPERATOR_FORT_H_

#define FORT_APPEND_ORIGINAL_ROSENBLUTH_COLLISION appendoriginalrosenbluthcollision_
#define FORT_COMPUTE_ORIGINAL_ROSENBLUTH_DIFFUSION_TENSOR computeoriginalrosenbluthdiffusiontensor_
#define FORT_COMPUTE_ORIGINAL_ROSENBLUTH_TEMPS computeoriginalrosenbluthtemps_
#define FORT_COMPUTE_ORIGINAL_ROSENBLUTH_DENOMS computeoriginalrosenbluthdenoms_
#define FORT_COMPUTE_ORIGINAL_ROSENBLUTH_NUMERATORS computeoriginalrosenbluthnumerators_
#define FORT_APPEND_ORIGINAL_ROSENBLUTH_BACK_REACTION appendoriginalrosenbluthbackreaction_

/**
 * Prototypes of fortran functions called by RosenbluthCollisionOperator.
 */
extern "C" {

   void
   FORT_APPEND_ORIGINAL_ROSENBLUTH_COLLISION(
      double &rhs,
      const double &u,
      const double &velocities,
      const double &d,
      const int &nd1lo,
      const int &nd1hi,
      const int &nd2lo,
      const int &nd2hi,
      const int &nd3lo,
      const int &nd3hi,
      const int &nd4lo,
      const int &nd4hi,
      const int &n1lo,
      const int &n1hi,
      const int &n2lo,
      const int &n2hi,
      const int &n3lo,
      const int &n3hi,
      const int &n4lo,
      const int &n4hi,
      const double& dx,
      const double& dparams,
      const int& iparams);

   void
   FORT_COMPUTE_ORIGINAL_ROSENBLUTH_DIFFUSION_TENSOR(
      double &d,
      const int &nd1a,
      const int &nd1b,
      const int &nd2a,
      const int &nd2b,
      const int &nd3a,
      const int &nd3b,
      const int &nd4a,
      const int &nd4b,
      const int &n1lo,
      const int &n1hi,
      const int &n2lo,
      const int &n2hi,
      const int &n3lo,
      const int &n3hi,
      const int &n4lo,
      const int &n4hi,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const double& velocities,
      const double& range_lo,
      const double& range_hi,
      const double& dparams,
      const int& iparams);

   void
   FORT_COMPUTE_ORIGINAL_ROSENBLUTH_TEMPS(
      double &fM,
      double &uMomx,
      double &uMomy,
      double &uKE,
      const int &nd1lo,
      const int &nd1hi,
      const int &nd2lo,
      const int &nd2hi,
      const int &nd3lo,
      const int &nd3hi,
      const int &nd4lo,
      const int &nd4hi,
      const double& velocities,
      const double& dparams);

   void
   FORT_COMPUTE_ORIGINAL_ROSENBLUTH_DENOMS(
      double &rMomx,
      double &rMomy,
      double &rKE,
      const double &cMomx,
      const double &cMomy,
      const double &cKE,
      const int &nd1lo,
      const int &nd1hi,
      const int &nd2lo,
      const int &nd2hi,
      const int &nd3lo,
      const int &nd3hi,
      const int &nd4lo,
      const int &nd4hi,
      const int &n1lo,
      const int &n1hi,
      const int &n2lo,
      const int &n2hi,
      const int &n3lo,
      const int &n3hi,
      const int &n4lo,
      const int &n4hi,
      const double& velocities,
      const double& dparams);

   void
   FORT_COMPUTE_ORIGINAL_ROSENBLUTH_NUMERATORS(
      double &rMomx,
      double &rMomy,
      double &rKE,
      const double &cMomx,
      const double &cMomy,
      const double &cKE,
      const double &fM,
      const double &u,
      const int &nd1lo,
      const int &nd1hi,
      const int &nd2lo,
      const int &nd2hi,
      const int &nd3lo,
      const int &nd3hi,
      const int &nd4lo,
      const int &nd4hi,
      const int &n1lo,
      const int &n1hi,
      const int &n2lo,
      const int &n2hi,
      const int &n3lo,
      const int &n3hi,
      const int &n4lo,
      const int &n4hi);

   void
   FORT_APPEND_ORIGINAL_ROSENBLUTH_BACK_REACTION(
      double& rhs,
      const double& cMomx,
      const double& cMomy,
      const double& cKE,
      const double& IMomxN,
      const double& IMomyN,
      const double& IKEN,
      const double& IMomxD,
      const double& IMomyD,
      const double& IKED,
      const int &nd1lo,
      const int &nd1hi,
      const int &nd2lo,
      const int &nd2hi,
      const int &nd3lo,
      const int &nd3hi,
      const int &nd4lo,
      const int &nd4hi,
      const int &n1lo,
      const int &n1hi,
      const int &n2lo,
      const int &n2hi,
      const int &n3lo,
      const int &n3hi,
      const int &n4lo,
      const int &n4hi,
      const double& dparams);
};

#endif
