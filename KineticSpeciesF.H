/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _KINETIC_SPECIES_FORT_H_
#define _KINETIC_SPECIES_FORT_H_

#define FORT_XPBY_4D                          xpby4d_
#define FORT_SET_ADVECTION_BCS_4D             setadvectionbcs4d_
#define FORT_SET_ADVECTION_BCS_4D_JB          setadvectionbcs4djb_
#define FORT_SET_ACCELERATION_BCS_4D          setaccelerationbcs4d_
#define FORT_SET_ACCELERATION_BCS_4D_JB       setaccelerationbcs4djb_
#define FORT_SET_PHASE_SPACE_VEL_4D           setphasespacevel4d_
#define FORT_SET_PHASE_SPACE_VEL_MAXWELL_4D   setphasespacevelmaxwell4d_
#define FORT_ZERO_GHOST_4D                    zero_ghost4d_
#define FORT_COMPUTE_ADVECTION_FLUXES_4D      computeadvectionfluxes4d_
#define FORT_COMPUTE_ADVECTION_DERIVATIVES_4D computeadvectionderivatives4d_
#define FORT_COMPUTE_ACCELERATION_DERIVATIVES_4D computeaccelerationderivatives4d_
#define FORT_COMPUTE_ACCELERATION_FLUXES_4D   computeaccelerationfluxes4d_
#define FORT_ACCUM_FLUX_DIV_4D                accumfluxdiv4d_
#define FORT_COMPUTE_CURRENTS_4D              computecurrents_
#define FORT_COMPUTE_KE                       computeke_
#define FORT_COMPUTE_KE_MAXWELL               computekemaxwell_
#define FORT_COMPUTE_KE_E_DOT                 computekeedot_
#define FORT_COMPUTE_MOM_4D                   computemom4d_
#define FORT_INTEGRATE_2D                     integrate2d_
#define FORT_COMPUTE_MOM_2D                   computemom2d_
#define FORT_COMPUTE_KE_FLUX                  computekeflux_
#define FORT_COMPUTE_KE_VEL_SPACE_FLUX        computekevelspaceflux_
#define FORT_APPEND_KROOK                     appendkrook_

/**
 * Prototypes of fortran functions called by KineticSpecies.
 */
extern "C" {
   void
   FORT_XPBY_4D(
      double &x,
      const double &y,
      const double &b,
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
   FORT_SET_ADVECTION_BCS_4D(
      double &u,
      const int &ng1lo,
      const int &ng1hi,
      const int &ng2lo,
      const int &ng2hi,
      const int &ng3lo,
      const int &ng3hi,
      const int &ng4lo,
      const int &ng4hi,
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
      const int &solution_order,
      const double &vel1,
      const double &vel2,
      const int &xPeriodic,
      const int &yPeriodic,
      const int64_t &ic);

   void
   FORT_SET_ACCELERATION_BCS_4D(
      double &u,
      const int &ng1lo,
      const int &ng1hi,
      const int &ng2lo,
      const int &ng2hi,
      const int &ng3lo,
      const int &ng3hi,
      const int &ng4lo,
      const int &ng4hi,
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
      const int &solution_order,
      const double &vel3,
      const double &vel4,
      const int64_t &ic);

   void
   FORT_SET_ADVECTION_BCS_4D_JB(
      double &u,
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
      const int &solution_order,
      const int &nGhost,
      const double &vel1,
      const double &vel2,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const int &xPeriodic,
      const int &yPeriodic,
      const int64_t &ic);

   void
   FORT_SET_ACCELERATION_BCS_4D_JB(
      double &u,
      const int &ng1lo,
      const int &ng1hi,
      const int &ng2lo,
      const int &ng2hi,
      const int &ng3lo,
      const int &ng3hi,
      const int &ng4lo,
      const int &ng4hi,
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
      const int &solution_order,
      const int &nGhost,
      const double &vel3,
      const double &vel4,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const int64_t &ic);

   void
   FORT_SET_PHASE_SPACE_VEL_4D(
      double &vel3,
      double &vel4,
      const int &nv1lo,
      const int &nv1hi,
      const int &nv2lo,
      const int &nv2hi,
      const int &nv3lo,
      const int &nv3hi,
      const int &nv4lo,
      const int &nv4hi,
      const int &ni1lo,
      const int &ni1hi,
      const int &ni2lo,
      const int &ni2hi,
      const int &ni3lo,
      const int &ni3hi,
      const int &ni4lo,
      const int &ni4hi,
      const double &vxface_velocities,
      const double &vyface_velocities,
      const double &charge_per_mass, //IEO
      const double &bz_const, //IEO
      const double &accel,
      const int &na1lo,
      const int &na1hi,
      const int &na2lo,
      const int &na2hi,
      double &axmax,
      double &aymax);

   void
   FORT_SET_PHASE_SPACE_VEL_MAXWELL_4D(
      double& vel3,
      double& vel4,
      const int& nv1a,
      const int& nv1b,
      const int& nv2a,
      const int& nv2b,
      const int& nv3a,
      const int& nv3b,
      const int& nv4a,
      const int& nv4b,
      const int& ni1a,
      const int& ni1b,
      const int& ni2a,
      const int& ni2b,
      const int& ni3a,
      const int& ni3b,
      const int& ni4a,
      const int& ni4b,
      const double& vxface_velocities,
      const double& vyface_velocities,
      const double& normalization,
      const double& bz_const, //IEO
      const double& em_vars,
      const double& vz,
      double& axmax,
      double& aymax);

   void
   FORT_ZERO_GHOST_4D(
      double &u,
      const int& n1a,
      const int& n1b,
      const int& n2a,
      const int& n2b,
      const int& n3a,
      const int& n3b,
      const int& n4a,
      const int& n4b,
      const int& nd1a,
      const int& nd1b,
      const int& nd2a,
      const int& nd2b,
      const int& nd3a,
      const int& nd3b,
      const int& nd4a,
      const int& nd4b);

   void
   FORT_COMPUTE_ADVECTION_FLUXES_4D(
      double &flux1,
      double &flux2,
      const int &nd1lo,
      const int &nd1hi,
      const int &nd2lo,
      const int &nd2hi,
      const int &nd3lo,
      const int &nd3hi,
      const int &nd4lo,
      const int &nd4hi,
      const double &vel1,
      const double &vel2,
      double &uface1,
      double &uface2,
      const double &u,
      const double &dx,
      const int &solution_order);

   void
   FORT_COMPUTE_ADVECTION_DERIVATIVES_4D(
      double &rhs,
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
      const double &vel1,
      const double &vel2,
      const double &dx,
      const int &solution_order);
  
   void
   FORT_COMPUTE_ACCELERATION_DERIVATIVES_4D(
      double &rhs,
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
      const double &vel3,
      const double &vel4,
      const double &dx,
      const int &solution_order);

   void
   FORT_COMPUTE_ACCELERATION_FLUXES_4D(
      double &flux3,
      double &flux4,
      const int &nd1lo,
      const int &nd1hi,
      const int &nd2lo,
      const int &nd2hi,
      const int &nd3lo,
      const int &nd3hi,
      const int &nd4lo,
      const int &nd4hi,
      const double &vel3,
      const double &vel4,
      double &uface3,
      double &uface4,
      const double &u,
      const double &dx,
      const int &solution_order);

   void
   FORT_ACCUM_FLUX_DIV_4D(
      double &rhs,
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
      const double &flux1,
      const double &flux2,
      const double &flux3,
      const double &flux4,
      const double& dx);

   void
   FORT_COMPUTE_CURRENTS_4D(
      const int& nd1lo,
      const int& nd1hi,
      const int& nd2lo,
      const int& nd2hi,
      const int& nd3lo,
      const int& nd3hi,
      const int& nd4lo,
      const int& nd4hi,
      const int& n1lo,
      const int& n1hi,
      const int& n2lo,
      const int& n2hi,
      const int& n3lo,
      const int& n3hi,
      const int& n4lo,
      const int& n4hi,
      const double& velocities,
      const double& u,
      const double& vz,
      double& JxLocal,
      double& JyLocal,
      double& JzLocal);

   void
   FORT_COMPUTE_KE(
      const int& nd1lo,
      const int& nd1hi,
      const int& nd2lo,
      const int& nd2hi,
      const int& nd3lo,
      const int& nd3hi,
      const int& nd4lo,
      const int& nd4hi,
      const int& n1lo,
      const int& n1hi,
      const int& n2lo,
      const int& n2hi,
      const int& n3lo,
      const int& n3hi,
      const int& n4lo,
      const int& n4hi,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const double& u,
      const double& mass,
      const double& velocities,
      double& ke,
      double& ke_x,
      double& ke_y,
      double& px,
      double& py);

   void
   FORT_COMPUTE_KE_MAXWELL(
      const int& nd1lo,
      const int& nd1hi,
      const int& nd2lo,
      const int& nd2hi,
      const int& nd3lo,
      const int& nd3hi,
      const int& nd4lo,
      const int& nd4hi,
      const int& n1lo,
      const int& n1hi,
      const int& n2lo,
      const int& n2hi,
      const int& n3lo,
      const int& n3hi,
      const int& n4lo,
      const int& n4hi,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const double& u,
      const double& mass,
      const double& velocities,
      const double& vz,
      double& ke,
      double& ke_x,
      double& ke_y,
      double& px,
      double& py);

   void
   FORT_COMPUTE_KE_E_DOT(
      const int& nd1lo,
      const int& nd1hi,
      const int& nd2lo,
      const int& nd2hi,
      const int& nd3lo,
      const int& nd3hi,
      const int& nd4lo,
      const int& nd4hi,
      const int& n1lo,
      const int& n1hi,
      const int& n2lo,
      const int& n2hi,
      const int& n3lo,
      const int& n3hi,
      const int& n4lo,
      const int& n4hi,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const double& u,
      const double& normalization,
      const double& velocities,
      const double& Ex,
      double& ke_e_dot);

   void
   FORT_COMPUTE_MOM_4D(
      const int& nd1lo,
      const int& nd1hi,
      const int& nd2lo,
      const int& nd2hi,
      const int& nd3lo,
      const int& nd3hi,
      const int& nd4lo,
      const int& nd4hi,
      const int& n1lo,
      const int& n1hi,
      const int& n2lo,
      const int& n2hi,
      const int& n3lo,
      const int& n3hi,
      const int& n4lo,
      const int& n4hi,
      const double& velocities,
      const double& u,
      double& momx,
      double& momy,
      double& ke,
      double& ent);

   void
   FORT_INTEGRATE_2D(
      const int& nd1lo,
      const int& nd1hi,
      const int& nd2lo,
      const int& nd2hi,
      const int& n1lo,
      const int& n1hi,
      const int& n2lo,
      const int& n2hi,
      const double& dx,
      const double& ke_2d,
      double& ke);

   void
   FORT_COMPUTE_MOM_2D(
      const int& nd1lo,
      const int& nd1hi,
      const int& nd2lo,
      const int& nd2hi,
      const int& n1lo,
      const int& n1hi,
      const int& n2lo,
      const int& n2hi,
      const double& dx,
      const double& momx_2d,
      const double& momy_2d,
      const double& ke_2d,
      const double& ent_2d,
      double& momx,
      double& momy,
      double& ke,
      double& ent);

   void
   FORT_COMPUTE_KE_FLUX(
      const int& nd1lo,
      const int& nd1hi,
      const int& nd2lo,
      const int& nd2hi,
      const int& nd3lo,
      const int& nd3hi,
      const int& nd4lo,
      const int& nd4hi,
      const int& n1lo,
      const int& n1hi,
      const int& n2lo,
      const int& n2hi,
      const int& n3lo,
      const int& n3hi,
      const int& n4lo,
      const int& n4hi,
      const int& ng1lo,
      const int& ng1hi,
      const int& ng2lo,
      const int& ng2hi,
      const int& ng3lo,
      const int& ng3hi,
      const int& ng4lo,
      const int& ng4hi,
      const double& dx,
      const double& face_flux1,
      const double& face_flux2,
      const double& face_flux3,
      const double& face_flux4,
      const double& velocities,
      const double& vxface_velocities,
      const double& vyface_velocities,
      const int& dir,
      const int& side,
      const double& mass,
      double& ke_flux);

   void
   FORT_COMPUTE_KE_VEL_SPACE_FLUX(
      const int& nd1lo,
      const int& nd1hi,
      const int& nd2lo,
      const int& nd2hi,
      const int& nd3lo,
      const int& nd3hi,
      const int& nd4lo,
      const int& nd4hi,
      const int& n1lo,
      const int& n1hi,
      const int& n2lo,
      const int& n2hi,
      const int& n3lo,
      const int& n3hi,
      const int& n4lo,
      const int& n4hi,
      const int& ng1lo,
      const int& ng1hi,
      const int& ng2lo,
      const int& ng2hi,
      const int& ng3lo,
      const int& ng3hi,
      const int& ng4lo,
      const int& ng4hi,
      const double& dx,
      const double& face_flux3,
      const double& face_flux4,
      double& ke_flux,
      const double& mass,
      const double& vxface_velocities,
      const double& vyface_velocities,
      const int& side,
      const int& dir);

   void
   FORT_APPEND_KROOK(
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
      const double& a_dt,
      const int64_t& ic,
      const double& nu,
      const double& u,
      double& rhs);

};

#endif
