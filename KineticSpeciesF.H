/*************************************************************************
 *
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * Written by Jeffrey Banks banksj3@rpi.edu (Rensselaer Polytechnic Institute,
 * Amos Eaton 301, 110 8th St., Troy, NY 12180); Jeffrey Hittinger
 * hittinger1@llnl.gov, William Arrighi arrighi2@llnl.gov, Richard Berger
 * berger5@llnl.gov, Thomas Chapman chapman29@llnl.gov (LLNL, P.O Box 808,
 * Livermore, CA 94551); Stephan Brunner stephan.brunner@epfl.ch (Ecole
 * Polytechnique Federale de Lausanne, EPFL SB SPC-TH, PPB 312, Station 13,
 * CH-1015 Lausanne, Switzerland).
 * CODE-744849
 *
 * All rights reserved.
 *
 * This file is part of Loki.  For details, see.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
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
#define FORT_COMPUTE_KE_4D                    computeke4d_
#define FORT_COMPUTE_KE_2D                    computeke2d_
#define FORT_COMPUTE_KE_FLUX_4D               computekeflux4d_
#define FORT_COMPUTE_KE_FLUX_2D               computekeflux2d_
#define FORT_COMPUTE_KE_VEL_SPACE_FLUX        computekevelspaceflux_
#define FORT_APPEND_KROOK                     appendkrook_

/**
 * Prototypes of fortran functions called by KineticSpecies.\
 */
extern "C" {
   void
   FORT_XPBY_4D(
      real &x,
      const int &nx1lo,
      const int &nx1hi,
      const int &nx2lo,
      const int &nx2hi,
      const int &nx3lo,
      const int &nx3hi,
      const int &nx4lo,
      const int &nx4hi,
      const real &y,
      const int &ny1lo,
      const int &ny1hi,
      const int &ny2lo,
      const int &ny2hi,
      const int &ny3lo,
      const int &ny3hi,
      const int &ny4lo,
      const int &ny4hi,
      const int &n1lo,
      const int &n1hi,
      const int &n2lo,
      const int &n2hi,
      const int &n3lo,
      const int &n3hi,
      const int &n4lo,
      const int &n4hi,
      const real &b);

   void
   FORT_SET_ADVECTION_BCS_4D(
      real &u,
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
      const real &vel1,
      const real &vel2,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const real &ic_phase_param,
      const int &xPeriodic,
      const int &yPeriodic);

   void
   FORT_SET_ACCELERATION_BCS_4D(
      real &u,
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
      const real &vel3,
      const real &vel4,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const real &ic_phase_param);

   void
   FORT_SET_ADVECTION_BCS_4D_JB(
      real &u,
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
      const real &vel1,
      const real &vel2,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const real &ic_phase_param,
      const int &xPeriodic,
      const int &yPeriodic);

   void
   FORT_SET_ACCELERATION_BCS_4D_JB(
      real &u,
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
      const real &vel3,
      const real &vel4,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const real &ic_phase_param);

   void
   FORT_SET_PHASE_SPACE_VEL_4D(
      real &vel3,
      real &vel4,
      const int &nv1lo,
      const int &nv1hi,
      const int &nv2lo,
      const int &nv2hi,
      const int &nv3lo,
      const int &nv3hi,
      const int &nv4lo,
      const int &nv4hi,
      const real &accel,
      const int &na1lo,
      const int &na1hi,
      const int &na2lo,
      const int &na2hi,
      const int &na3lo,
      const int &na3hi,
      const int &na4lo,
      const int &na4hi,
      real &axmax,
      real &aymax);

   void
   FORT_SET_PHASE_SPACE_VEL_MAXWELL_4D(
      real& vel3,
      real& vel4,
      const int& nv1a,
      const int& nv1b,
      const int& nv2a,
      const int& nv2b,
      const int& nv3a,
      const int& nv3b,
      const int& nv4a,
      const int& nv4b,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const real& charge_per_mass,
      const real& em_vars,
      const real& vz,
      real& axmax,
      real& aymax);

   void
   FORT_ZERO_GHOST_4D(
      real &u,
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
      real &flux1,
      real &flux2,
      const int &nd1lo,
      const int &nd1hi,
      const int &nd2lo,
      const int &nd2hi,
      const int &nd3lo,
      const int &nd3hi,
      const int &nd4lo,
      const int &nd4hi,
      const real &vel1,
      const real &vel2,
      real &uface1,
      real &uface2,
      const real &u,
      const real &dx,
      const int &solution_order);

   void
   FORT_COMPUTE_ADVECTION_DERIVATIVES_4D(
      real &rhs,
      const real &u,
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
      const real &vel1,
      const real &vel2,
      const real &dx,
      const int &solution_order);

   void
   FORT_COMPUTE_ACCELERATION_DERIVATIVES_4D(
      real &rhs,
      const real &u,
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
      const real &vel3,
      const real &vel4,
      const real &dx,
      const int &solution_order);

   void
   FORT_COMPUTE_ACCELERATION_FLUXES_4D(
      real &flux3,
      real &flux4,
      const int &nd1lo,
      const int &nd1hi,
      const int &nd2lo,
      const int &nd2hi,
      const int &nd3lo,
      const int &nd3hi,
      const int &nd4lo,
      const int &nd4hi,
      const real &vel3,
      const real &vel4,
      real &uface3,
      real &uface4,
      const real &u,
      const real &dx,
      const int &solution_order);

   void
   FORT_ACCUM_FLUX_DIV_4D(
      real &rhs,
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
      const real &flux1,
      const real &flux2,
      const real &flux3,
      const real &flux4,
      PROBLEMDOMAIN_TO_FORT_PROTO);

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
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const real& u,
      const real& vz,
      real& JxLocal,
      real& JyLocal,
      real& JzLocal);

   void
   FORT_COMPUTE_KE_4D(
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
      const real& u,
      real& ke);

   void
   FORT_COMPUTE_KE_2D(
      const int& nd1lo,
      const int& nd1hi,
      const int& nd2lo,
      const int& nd2hi,
      const int& n1lo,
      const int& n1hi,
      const int& n2lo,
      const int& n2hi,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const real& ke_2d,
      real& ke);

   void
   FORT_COMPUTE_KE_FLUX_4D(
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
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const real& face_flux1,
      const real& face_flux2,
      const real& face_vel1,
      const real& face_vel2,
      real& ke_flux,
      const int& dir);

   void
   FORT_COMPUTE_KE_FLUX_2D(
      const int& nd1lo,
      const int& nd1hi,
      const int& nd2lo,
      const int& nd2hi,
      const int& n1lo,
      const int& n1hi,
      const int& n2lo,
      const int& n2hi,
      const int& ng1lo,
      const int& ng1hi,
      const int& ng2lo,
      const int& ng2hi,
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const real& ke_flux_2d,
      real& ke_flux,
      const int& dir,
      const int& side);

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
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const real& face_flux3,
      const real& face_flux4,
      real& ke_flux,
      const real& mass,
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
      PROBLEMDOMAIN_TO_FORT_PROTO,
      const real& xlo_krook,
      const real& xhi_krook,
      const real& krookPower,
      const real& krookCoeff,
      const real& u,
      real& rhs);

};

#endif
