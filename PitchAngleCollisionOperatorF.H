/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _PITCH_ANGLE_COLLISION_OPERATOR_FORT_H_
#define _PITCH_ANGLE_COLLISION_OPERATOR_FORT_H_

#define FORT_APPEND_PITCH_ANGLE_COLLISION appendpitchanglecollision_
#define FORT_COMPUTE_PITCH_ANGLE_SPECIES_MOMENTS computepitchanglespeciesmoments_
#define FORT_COMPUTE_PITCH_ANGLE_SPECIES_REDUCED_FIELDS computepitchanglespeciesreducedfields_
#define FORT_COMPUTE_PITCH_ANGLE_SPECIES_KEC computepitchanglespecieskec_
#define FORT_COMPUTE_PITCH_ANGLE_SPECIES_VTHERMAL computepitchanglespeciesvthermal_

/**
 * Prototypes of fortran functions called by PitchAngleCollisionOperator.
 */
extern "C" {

   void
   FORT_APPEND_PITCH_ANGLE_COLLISION(
      double &rhs,
      const double &u,
      const double& a_velocities,
      const double &IVx,
      const double &IVy,
      const double &Vth,
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
      const double &range_lo,
      const double &range_hi,
      const double& dparams,
      const int& iparams);

  void
  FORT_COMPUTE_PITCH_ANGLE_SPECIES_MOMENTS(
      double &rN,
      double &rGammax,
      double &rGammay,
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
  FORT_COMPUTE_PITCH_ANGLE_SPECIES_REDUCED_FIELDS(
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
  FORT_COMPUTE_PITCH_ANGLE_SPECIES_KEC(
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
  FORT_COMPUTE_PITCH_ANGLE_SPECIES_VTHERMAL(
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

};

#endif
