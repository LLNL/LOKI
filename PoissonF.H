/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _POISSON_F_H_
#define _POISSON_F_H_

#define neutralizeCharge4D            neutralizecharge4d_
#define computeEFieldFromPotential    computeefieldfrompotential_

/**
 * Prototypes of fortran functions called by Poisson.
 */
extern "C"
{
   void
   neutralizeCharge4D(
      const int &md1a,
      const int &md1b,
      const int &md2a,
      const int &md2b,
      const int &n1a,
      const int &n1b,
      const int &n2a,
      const int &n2b,
      double &rhs,
      int &comm);

   void
   computeEFieldFromPotential(
      const int &nd1a,
      const int &nd1b,
      const int &nd2a,
      const int &nd2b,
      const int &n1a,
      const int &n1b,
      const int &n2a,
      const int &n2b,
      const int &solution_order,
      const int &em_vars_dim,
      const double &dx,
      double &em_vars,
      double &phi);
};

#endif
