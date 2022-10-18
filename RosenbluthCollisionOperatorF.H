/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _ROSENBLUTH_COLLISION_OPERATOR_FORT_H_
#define _ROSENBLUTH_COLLISION_OPERATOR_FORT_H_

#define FORT_APPEND_ROSENBLUTH_COLLISION appendrosenbluthcollision_
#define FORT_GET_MAX_RELAM getmaxrelam_
#define FORT_COMPUTE_ROSENBLUTH_TEMPS computerosenbluthtemps_
#define FORT_COMPUTE_ROSENBLUTH_DENOMS computerosenbluthdenoms_
#define FORT_COMPUTE_ROSENBLUTH_NUMERATORS computerosenbluthnumerators_
#define FORT_APPEND_ROSENBLUTH_BACK_REACTION appendrosenbluthbackreaction_

#define FORT_COMPUTE_ROSENBLUTH_SPECIES_MOMENTS computerosenbluthspeciesmoments_
#define FORT_COMPUTE_ROSENBLUTH_DIFFUSION_TENSOR computerosenbluthdiffusiontensor_
#define FORT_COMPUTE_ROSENBLUTH_SPECIES_KEC computerosenbluthspecieskec_
#define FORT_GET_COMMUNICATED_SPECIES_MOMENTS getcommunicatedspeciesmoments_
#define FORT_GET_COMMUNICATED_BR_SPECIES_MOMENTS getcommunicatedbrspeciesmoments_
#define FORT_SET_COMMUNICATED_SPECIES_MOMENTS setcommunicatedspeciesmoments_
#define FORT_SET_COMMUNICATED_SPECIES_BR_MOMENTS setcommunicatedspeciesbrmoments_
#define FORT_COMPUTE_ROSENBLUTH_SPECIES_REDUCED_FIELDS computerosenbluthspeciesreducedfields_
#define FORT_COMPUTE_ROSENBLUTH_SPECIES_VTHERMAL computerosenbluthspeciesvthermal_

#define FORT_ADD_TO_DOUBLEARRAY_4D addtodoublearray4d_

#define FORT_COMPUTE_ROSENBLUTH_ENTROPY_CHANGE computerosenbluthentropychange_

/**
 * Prototypes of fortran functions called by RosenbluthCollisionOperator.
 */
extern "C" {

   void
   FORT_APPEND_ROSENBLUTH_COLLISION(
      double &rhs,
      const double &u,
      const double &velocities,
      const double &IVx,
      const double &IVy,
      const double &d,
      const double &Vth,
      const double &N,
      const double &massR,
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
      const double &dx,
      const double& dparams,
      const int& iparams);

   void
   FORT_GET_MAX_RELAM(
      double& reLamMax,
      const double& reLam,
      const int& n1lo,
      const int& n1hi,
      const int& n2lo,
      const int& n2hi,
      const int& n3lo,
      const int& n3hi,
      const int& n4lo,
      const int& n4hi,
      const int& iparams);

   void
   FORT_COMPUTE_ROSENBLUTH_DIFFUSION_TENSOR(
      double &d,
      double &reLam,
      const double &mass_ratio,
      const double &rVx0,
      const double &rVy0,
      const double &IvTh,
      const double &N,
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
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const double &velocities,
      const int& BR,
      const int& compute_re_lam,
      const double& range_lo,
      const double& range_hi,
      const double& dparams,
      const int& iparams);
  
   void
   FORT_ADD_TO_DOUBLEARRAY_4D(
      double &a,
      const double &b,
      const double &ca,
      const double &cb,
      const int &nd1lo,
      const int &nd1hi,
      const int &nd2lo,
      const int &nd2hi,
      const int &nd3lo,
      const int &nd3hi,
      const int &nd4lo,
      const int &nd4hi);

   void
   FORT_COMPUTE_ROSENBLUTH_ENTROPY_CHANGE(
      double &rEnt,
      const double &u,
      const double &dudt,
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
   FORT_COMPUTE_ROSENBLUTH_TEMPS(
      double &uMomx,
      double &uMomy,
      double &uKE,
      const double &IVx,
      const double &IVy,
      const double &Vth,
      const double &n,
      const int &nd1lo,
      const int &nd1hi,
      const int &nd2lo,
      const int &nd2hi,
      const int &nd3lo,
      const int &nd3hi,
      const int &nd4lo,
      const int &nd4hi,
      const double& velocities);

   void
   FORT_COMPUTE_ROSENBLUTH_DENOMS(
      double &rMomx,
      double &rMomy,
      double &rKE,
      const double &cMomx,
      const double &cMomy,
      const double &cKE,
      const double &IVx,
      const double &IVy,
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
      const double& velocities);

   void
   FORT_COMPUTE_ROSENBLUTH_NUMERATORS(
      double &rMomx,
      double &rMomy,
      double &rKE,
      const double &IVx,
      const double &IVy,
      const double &dfa,
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
      const double& velocities);

  void
  FORT_COMPUTE_ROSENBLUTH_SPECIES_MOMENTS(
      double &rN,
      double &rGammax,
      double &rGammay,
      double &rKE,
      const double &u,
      const int &diagnostics,
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
      const double& velocities);

  void
  FORT_COMPUTE_ROSENBLUTH_SPECIES_KEC(
      double &rKEC,
      const double &rVx0,
      const double &rVy0,
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
      const int &n4hi,
      const double& velocities);

  void
  FORT_COMPUTE_ROSENBLUTH_SPECIES_REDUCED_FIELDS(
      double &vx,
      double &vy,
      const double &n,
      const double &Gammax,
      const double &Gammay,
      const int &nd1lo,
      const int &nd1hi,
      const int &nd2lo,
      const int &nd2hi,
      const int &nd3lo,
      const int &nd3hi,
      const int &nd4lo,
      const int &nd4hi);
  
  void
  FORT_GET_COMMUNICATED_SPECIES_MOMENTS(
      double &data,
      const double &n,
      const double &vx,
      const double &vy,
      const double &vthsq,
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
      const int &nf);

  void
  FORT_GET_COMMUNICATED_BR_SPECIES_MOMENTS(
      double &data,
      const double &IMomxN,
      const double &IMomyN,
      const double &IKEN,
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
      const int &nf,
      const int &other_species_loc);

  void
  FORT_SET_COMMUNICATED_SPECIES_MOMENTS(
      double &n,
      double &vx,
      double &vy,
      double &vthsq,
      const double &data,
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
      const int &nf);
  
  void
  FORT_SET_COMMUNICATED_SPECIES_BR_MOMENTS(
      double &IMomxN,
      double &IMomyN,
      double &IKEN,
      const double &data,
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
      const int &nf);
  
  void
  FORT_COMPUTE_ROSENBLUTH_SPECIES_VTHERMAL(
      double &vthsq,
      const double &KEC,
      const double &n,
      const int &nd1lo,
      const int &nd1hi,
      const int &nd2lo,
      const int &nd2hi,
      const int &nd3lo,
      const int &nd3hi,
      const int &nd4lo,
      const int &nd4hi);
  
   void
   FORT_APPEND_ROSENBLUTH_BACK_REACTION(
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
      const double& massR);
};

#endif
