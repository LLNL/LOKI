/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "LokiPoissonSolveSuperLU.H"
#include "Directions.H"
#include "Loki_Utilities.H"

namespace Loki {

const string LokiPoissonSolveSuperLU::s_CLASS_NAME("superlu");


bool
LokiPoissonSolveSuperLU::isType(
   const string& a_name)
{
   if (a_name.compare(s_CLASS_NAME) == 0) {
      return true;
   }
   return false;
}


LokiPoissonSolveSuperLU::LokiPoissonSolveSuperLU(
   const ProblemDomain& a_domain,
   int a_order,
   MPI_Comm a_comm)
   : m_nx(a_domain.numberOfCells(X1)),
     m_ny(a_domain.numberOfCells(X2))
{
   // Set up general problem parameters.
   if (a_order == 4) {
      m_nghosts = 2;
   }
   else {
      m_nghosts = 3;
   }
   m_nx_w_bdy = m_nx+2*m_nghosts;
   m_ny_w_bdy = m_ny+2*m_nghosts;
   m_nx_tot = m_nx_w_bdy + 1;
   m_ny_tot = m_ny_w_bdy + 1;

   // Set up PETSC_COMM_WORLD to be the Poisson object's communicator.  Then
   // initialize Petsc.
   int argc = 0;
   if (!PetscInitializeCalled) {
      PETSC_COMM_WORLD = a_comm;
      PetscInitialize(&argc, NULL, NULL, NULL);
   }

   // Create the linear system matrix, the rhs, and solution vectors.
   buildAndFillSystemMat(a_domain, a_order);

   // Create the linear system context.
   createLinearSystemContext();
}


LokiPoissonSolveSuperLU::~LokiPoissonSolveSuperLU()
{
   // Destroy all PETSc objects and finalize PETSc.
   KSPDestroy(&m_ksp);
   VecDestroy(&m_x);
   VecDestroy(&m_b);
   MatDestroy(&m_A);
   PetscFinalize();
}


void
LokiPoissonSolveSuperLU::solve(
   ParallelArray& a_solution,
   const ParallelArray& a_rhs)
{
   // Dump the interior of a_rhs into m_b.
   const double* a_rhs_vals = a_rhs.getData();
   int offset = m_nghosts*(m_nx_w_bdy+1);
   for (int j = m_nghosts; j < m_ny+m_nghosts; ++j) {
      for (int i = m_nghosts; i < m_nx+m_nghosts; ++i) {
         VecSetValue(m_b,
            j*(m_nx_tot)+i,
            a_rhs_vals[offset+i-m_nghosts]*m_scale,
            INSERT_VALUES);
      }
      offset += m_nx_w_bdy;
   }

   // Solve system.
   KSPSolve(m_ksp, m_b, m_x);

   // Dump m_x into a_solution.
   int elements[m_nx_w_bdy];
   double* a_solution_vals = a_solution.getData();
   offset = 0;
   for (int j = 0; j < m_ny_w_bdy; ++j) {
      for (int i = 0; i < m_nx_w_bdy; ++i) {
         elements[i] = j*m_nx_tot+i;
      }
      VecGetValues(m_x, m_nx_w_bdy, elements, a_solution_vals+offset);
      offset += m_nx_w_bdy;
   }
   elements[0] = m_nx_tot*m_ny_tot-1;
   double alpha;
   VecGetValues(m_x, 1, elements, &alpha);
}


void
LokiPoissonSolveSuperLU::printParameters() const
{
   // Write the solver's parameters.  For this solver there are not parameters,
   // just the type of solver.
   Loki_Utilities::printF("  Using SuperLU Poisson solver\n" );
}


void
LokiPoissonSolveSuperLU::buildAndFillSystemMat(
   const ProblemDomain& a_domain,
   int a_order)
{
   // Total number of equations.
   int neqns = m_nx_tot*m_ny_tot;

   // Create and structure the linear sytem matrix.
   MatCreate(PETSC_COMM_WORLD, &m_A);
   MatSetSizes(m_A, neqns, neqns, neqns, neqns);
   MatSetFromOptions(m_A);
   int preallocation[neqns];
   int row = 0;
   for (int i = 0; i < m_nghosts*m_nx_tot; ++i) {
      preallocation[row++] = 2;
   }
   for (int j = m_nghosts; j < m_nghosts+m_ny; ++j) {
      for (int i = 0; i < m_nghosts; ++i) {
         preallocation[row++] = 2;
      }
      for (int i = m_nghosts; i < m_nghosts+m_nx; ++i) {
         preallocation[row++] = 2*a_order+2;
      }
      for (int i = m_nghosts+m_nx; i < m_nx_tot; ++i) {
         preallocation[row++] = 2;
      }
   }
   for (int i = 0; i < m_nghosts*m_nx_tot+m_nx_w_bdy; ++i) {
      preallocation[row++] = 2;
   }
   preallocation[row] = m_nx*m_ny+1;
   int tmp;
   MatMPIAIJSetPreallocation(m_A, tmp, preallocation, tmp, preallocation);
   MatSeqAIJSetPreallocation(m_A, tmp, preallocation);

   // Create the RHS and solution vectors.
   VecCreate(PETSC_COMM_WORLD, &m_x);
   VecSetSizes(m_x, neqns, neqns);
   VecSetFromOptions(m_x);
   VecDuplicate(m_x, &m_b);
   VecZeroEntries(m_b);

   // Fill the linear system matrix.
   double d2x = a_domain.dx(X1)*a_domain.dx(X1);
   double d2y = a_domain.dx(X2)*a_domain.dx(X2);
   int eqn = 0;
   if (a_order == 4) {
      // 4th order Laplacian stencil.
      double xfactor = 1.0/(12.0*d2x);
      double yfactor = 1.0/(12.0*d2y);
      m_scale = 1.0/(30.0*(xfactor+yfactor));
      double stencilLower[2][2];
      double stencilUpper[2][2];
      stencilLower[X1][0] = -1.0*xfactor*m_scale;
      stencilLower[X1][1] = 16.0*xfactor*m_scale;
      stencilLower[X2][0] = -1.0*yfactor*m_scale;
      stencilLower[X2][1] = 16.0*yfactor*m_scale;
      stencilUpper[X1][0] = 16.0*xfactor*m_scale;
      stencilUpper[X1][1] = -1.0*xfactor*m_scale;
      stencilUpper[X2][0] = 16.0*yfactor*m_scale;
      stencilUpper[X2][1] = -1.0*yfactor*m_scale;
      double stencilCenter = -30.0*(xfactor+yfactor)*m_scale;

      // Array of columns to fill with stencil.
      int colsLower[2];
      int colsUpper[2];
      int colsCenter;

      // Ghosts at y<ya.
      for (int j = 0; j < m_nghosts; ++j) {
         for (int i = 0; i < m_nghosts; ++i) {
            MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
            MatSetValue(m_A, eqn, eqn+m_nx_tot*m_ny+m_nx, -1.0, INSERT_VALUES);
            ++eqn;
         }
         for (int i = m_nghosts; i < m_nghosts+m_nx; ++i) {
            MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
            MatSetValue(m_A, eqn, eqn+m_nx_tot*m_ny, -1.0, INSERT_VALUES);
            ++eqn;
         }
         for (int i = m_nghosts+m_nx; i < m_nx_tot; ++i) {
            MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
            MatSetValue(m_A, eqn, eqn+m_nx_tot*m_ny-m_nx, -1.0, INSERT_VALUES);
            ++eqn;
         }
      }
      // Lines between y=ya and y=yb.
      for (int j = m_nghosts; j < m_nghosts+m_ny; ++j) {
         // Ghosts at x<xa.
         for (int i = 0; i < m_nghosts; ++i) {
            MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
            MatSetValue(m_A, eqn, eqn+m_nx, -1.0, INSERT_VALUES);
            ++eqn;
         }
         // Interior, between y=ya and y=yb and between x=xa and x=xb.
         for (int i = m_nghosts; i < m_nghosts+m_nx; ++i) {
            // d2x
            if (i < 2*m_nghosts) {
               if (i == m_nghosts) {
                  colsLower[0] = eqn-2+m_nx;
                  colsLower[1] = eqn-1+m_nx;
               }
               else if (i == m_nghosts+1) {
                  colsLower[0] = eqn-2+m_nx;
                  colsLower[1] = eqn-1;
               }
               colsUpper[0] = eqn+1;
               colsUpper[1] = eqn+2;
            }
            else if (i > m_nx-1) {
               if (i == m_nx+1) {
                  colsUpper[0] = eqn+1-m_nx;
                  colsUpper[1] = eqn+2-m_nx;
               }
               else if (i == m_nx) {
                  colsUpper[0] = eqn+1;
                  colsUpper[1] = eqn+2-m_nx;
               }
               colsLower[0] = eqn-2;
               colsLower[1] = eqn-1;
            }
            else {
               colsLower[0] = eqn-2;
               colsLower[1] = eqn-1;
               colsUpper[0] = eqn+1;
               colsUpper[1] = eqn+2;
            }
            MatSetValues(m_A,
               1,
               &eqn,
               2,
               colsLower,
               &stencilLower[X1][0],
               INSERT_VALUES);
            MatSetValues(m_A,
               1,
               &eqn,
               2,
               colsUpper,
               &stencilUpper[X1][0],
               INSERT_VALUES);
            // d2y
            if (j < 2*m_nghosts) {
               if (j == m_nghosts) {
                  colsLower[0] = eqn+m_nx_tot*(m_ny-2);
                  colsLower[1] = eqn+m_nx_tot*(m_ny-1);
               }
               else if (j == m_nghosts+1) {
                  colsLower[0] = eqn+m_nx_tot*(m_ny-2);
                  colsLower[1] = eqn-m_nx_tot;
               }
               colsUpper[0] = eqn+m_nx_tot;
               colsUpper[1] = eqn+2*m_nx_tot;
            }
            else if (j > m_ny-1) {
               if (j == m_ny+1) {
                  colsUpper[0] = eqn-m_nx_tot*(m_ny-1);
                  colsUpper[1] = eqn-m_nx_tot*(m_ny-2);
               }
               else if (j == m_ny) {
                  colsUpper[0] = eqn+m_nx_tot;
                  colsUpper[1] = eqn-m_nx_tot*(m_ny-2);
               }
               colsLower[0] = eqn-2*m_nx_tot;
               colsLower[1] = eqn-m_nx_tot;
            }
            else {
               colsLower[0] = eqn-2*m_nx_tot;
               colsLower[1] = eqn-m_nx_tot;
               colsUpper[0] = eqn+m_nx_tot;
               colsUpper[1] = eqn+2*m_nx_tot;
            }
            MatSetValues(m_A,
               1,
               &eqn,
               2,
               colsLower,
               &stencilLower[X2][0],
               INSERT_VALUES);
            MatSetValues(m_A,
               1,
               &eqn,
               2,
               colsUpper,
               &stencilUpper[X2][0],
               INSERT_VALUES);
            // dx2 & dy2
            colsCenter = eqn;
            MatSetValue(m_A, eqn, colsCenter, stencilCenter, INSERT_VALUES);
            // Extra equation.
            MatSetValue(m_A, eqn, neqns-1, 1.0, INSERT_VALUES);
            ++eqn;
         }
         // Ghosts at x>xb.
         for (int i = m_nghosts+m_nx; i < m_nx_tot; ++i) {
            MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
            MatSetValue(m_A, eqn, eqn-m_nx, -1.0, INSERT_VALUES);
            ++eqn;
         }
      }
      // Ghosts at y>yb.
      for (int j = m_nghosts+m_ny; j < m_ny_w_bdy; ++j) {
         for (int i = 0; i < m_nghosts; ++i) {
            MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
            MatSetValue(m_A, eqn, eqn-m_nx_tot*m_ny+m_nx, -1.0, INSERT_VALUES);
            ++eqn;
         }
         for (int i = m_nghosts; i < m_nghosts+m_nx; ++i) {
            MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
            MatSetValue(m_A, eqn, eqn-m_nx_tot*m_ny, -1.0, INSERT_VALUES);
            ++eqn;
         }
         for (int i = m_nghosts+m_nx; i < m_nx_tot; ++i) {
            MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
            MatSetValue(m_A, eqn, eqn-m_nx_tot*m_ny-m_nx, -1.0, INSERT_VALUES);
            ++eqn;
         }
      }
      // Extra equation.
      for (int i = 0; i < m_nghosts; ++i) {
         MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
         MatSetValue(m_A, eqn, eqn-m_nx_tot*m_ny+m_nx, -1.0, INSERT_VALUES);
         ++eqn;
      }
      for (int i = m_nghosts; i < m_nghosts+m_nx; ++i) {
         MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
         MatSetValue(m_A, eqn, eqn-m_nx_tot*m_ny, -1.0, INSERT_VALUES);
         ++eqn;
      }
      for (int i = m_nghosts+m_nx; i < m_nx_w_bdy; ++i) {
         MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
         MatSetValue(m_A, eqn, eqn-m_nx_tot*m_ny-m_nx, -1.0, INSERT_VALUES);
         ++eqn;
      }
      int loc = m_nghosts*m_nx_tot+m_nghosts;
      for (int j = 0; j < m_ny; ++j) {
         for (int i = 0; i < m_nx; ++i) {
            if (loc+i == eqn-m_nx_tot*m_ny-m_nx) {
               MatSetValue(m_A, eqn, loc+i, -1.0, INSERT_VALUES);
            }
            else {
               MatSetValue(m_A, eqn, loc+i, 1.0, INSERT_VALUES);
            }
         }
         loc += m_nx_tot;
      }
      MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
   }
   else {
      // 6th order Laplacian stencil.
      double xfactor = 1.0/(180.0*d2x);
      double yfactor = 1.0/(180.0*d2y);
      m_scale = 1.0/(490.0*(xfactor+yfactor));
      double stencilLower[2][3];
      double stencilUpper[2][3];
      stencilLower[X1][0] =   2.0*xfactor*m_scale;
      stencilLower[X1][1] = -27.0*xfactor*m_scale;
      stencilLower[X1][2] = 270.0*xfactor*m_scale;
      stencilLower[X2][0] =   2.0*yfactor*m_scale;
      stencilLower[X2][1] = -27.0*yfactor*m_scale;
      stencilLower[X2][2] = 270.0*yfactor*m_scale;
      stencilUpper[X1][0] = 270.0*xfactor*m_scale;
      stencilUpper[X1][1] = -27.0*xfactor*m_scale;
      stencilUpper[X1][2] =   2.0*xfactor*m_scale;
      stencilUpper[X2][0] = 270.0*yfactor*m_scale;
      stencilUpper[X2][1] = -27.0*yfactor*m_scale;
      stencilUpper[X2][2] =   2.0*yfactor*m_scale;
      double stencilCenter = -490.0*(xfactor+yfactor)*m_scale;

      // Array of columns to fill with stencil.
      int colsLower[3];
      int colsUpper[3];
      int colsCenter;

      // Ghosts at y<ya.
      for (int j = 0; j < m_nghosts; ++j) {
         for (int i = 0; i < m_nghosts; ++i) {
            MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
            MatSetValue(m_A, eqn, eqn+m_nx_tot*m_ny+m_nx, -1.0, INSERT_VALUES);
            ++eqn;
         }
         for (int i = m_nghosts; i < m_nghosts+m_nx; ++i) {
            MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
            MatSetValue(m_A, eqn, eqn+m_nx_tot*m_ny, -1.0, INSERT_VALUES);
            ++eqn;
         }
         for (int i = m_nghosts+m_nx; i < m_nx_tot; ++i) {
            MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
            MatSetValue(m_A, eqn, eqn+m_nx_tot*m_ny-m_nx, -1.0, INSERT_VALUES);
            ++eqn;
         }
      }
      // Lines between y=ya and y=yb.
      for (int j = m_nghosts; j < m_ny+m_nghosts; ++j) {
         // Ghosts at x<xa.
         for (int i = 0; i < m_nghosts; ++i) {
            MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
            MatSetValue(m_A, eqn, eqn+m_nx, -1.0, INSERT_VALUES);
            ++eqn;
         }
         // Interior, between y=ya and y=yb and between x=xa and x=xb.
         for (int i = m_nghosts; i < m_nx+m_nghosts; ++i) {
            // d2x
            if (i < 2*m_nghosts) {
               if (i == m_nghosts) {
                  colsLower[0] = eqn-3+m_nx;
                  colsLower[1] = eqn-2+m_nx;
                  colsLower[2] = eqn-1+m_nx;
               }
               else if (i == m_nghosts+1) {
                  colsLower[0] = eqn-3+m_nx;
                  colsLower[1] = eqn-2+m_nx;
                  colsLower[2] = eqn-1;
               }
               else if (i == m_nghosts+2) {
                  colsLower[0] = eqn-3+m_nx;
                  colsLower[1] = eqn-2;
                  colsLower[2] = eqn-1;
               }
               colsUpper[0] = eqn+1;
               colsUpper[1] = eqn+2;
               colsUpper[2] = eqn+3;
            }
            else if (i > m_nx-1) {
               if (i == m_nx+2) {
                  colsUpper[0] = eqn+1-m_nx;
                  colsUpper[1] = eqn+2-m_nx;
                  colsUpper[2] = eqn+3-m_nx;
               }
               else if (i == m_nx+1) {
                  colsUpper[0] = eqn+1;
                  colsUpper[1] = eqn+2-m_nx;
                  colsUpper[2] = eqn+3-m_nx;
               }
               else if (i == m_nx) {
                  colsUpper[0] = eqn+1;
                  colsUpper[1] = eqn+2;
                  colsUpper[2] = eqn+3-m_nx;
               }
               colsLower[0] = eqn-3;
               colsLower[1] = eqn-2;
               colsLower[2] = eqn-1;
            }
            else {
               colsLower[0] = eqn-3;
               colsLower[1] = eqn-2;
               colsLower[2] = eqn-1;
               colsUpper[0] = eqn+1;
               colsUpper[1] = eqn+2;
               colsUpper[2] = eqn+3;
            }
            MatSetValues(m_A,
               1,
               &eqn,
               3,
               colsLower,
               &stencilLower[X1][0],
               INSERT_VALUES);
            MatSetValues(m_A,
               1,
               &eqn,
               3,
               colsUpper,
               &stencilUpper[X1][0],
               INSERT_VALUES);
            // d2y
            if (j < 2*m_nghosts) {
               if (j == m_nghosts) {
                  colsLower[0] = eqn+m_nx_tot*(m_ny-3);
                  colsLower[1] = eqn+m_nx_tot*(m_ny-2);
                  colsLower[2] = eqn+m_nx_tot*(m_ny-1);
               }
               else if (j == m_nghosts+1) {
                  colsLower[0] = eqn+m_nx_tot*(m_ny-3);
                  colsLower[1] = eqn+m_nx_tot*(m_ny-2);
                  colsLower[2] = eqn-m_nx_tot;
               }
               else if (j == m_nghosts+2) {
                  colsLower[0] = eqn+m_nx_tot*(m_ny-3);
                  colsLower[1] = eqn-m_nx_tot*2;
                  colsLower[2] = eqn-m_nx_tot;
               }
               colsUpper[0] = eqn+m_nx_tot;
               colsUpper[1] = eqn+2*m_nx_tot;
               colsUpper[2] = eqn+3*m_nx_tot;
            }
            else if (j > m_ny-1) {
               if (j == m_ny+2) {
                  colsUpper[0] = eqn-m_nx_tot*(m_ny-1);
                  colsUpper[1] = eqn-m_nx_tot*(m_ny-2);
                  colsUpper[2] = eqn-m_nx_tot*(m_ny-3);
               }
               else if (j == m_ny+1) {
                  colsUpper[0] = eqn+m_nx_tot;
                  colsUpper[1] = eqn-m_nx_tot*(m_ny-2);
                  colsUpper[2] = eqn-m_nx_tot*(m_ny-3);
               }
               else if (j == m_ny) {
                  colsUpper[0] = eqn+m_nx_tot;
                  colsUpper[1] = eqn+m_nx_tot*2;
                  colsUpper[2] = eqn-m_nx_tot*(m_ny-3);
               }
               colsLower[0] = eqn-3*m_nx_tot;
               colsLower[1] = eqn-2*m_nx_tot;
               colsLower[2] = eqn-m_nx_tot;
            }
            else {
               colsLower[0] = eqn-3*m_nx_tot;
               colsLower[1] = eqn-2*m_nx_tot;
               colsLower[2] = eqn-m_nx_tot;
               colsUpper[0] = eqn+m_nx_tot;
               colsUpper[1] = eqn+2*m_nx_tot;
               colsUpper[2] = eqn+3*m_nx_tot;
            }
            MatSetValues(m_A,
               1,
               &eqn,
               3,
               colsLower,
               &stencilLower[X2][0],
               INSERT_VALUES);
            MatSetValues(m_A,
               1,
               &eqn,
               3,
               colsUpper,
               &stencilUpper[X2][0],
               INSERT_VALUES);
            // dx2 & dy2
            colsCenter = eqn;
            MatSetValue(m_A, eqn, colsCenter, stencilCenter, INSERT_VALUES);
            // Extra equation.
            MatSetValue(m_A, eqn, neqns-1, 1.0, INSERT_VALUES);
            ++eqn;
         }
         // Ghosts at x>xb.
         for (int i = m_nghosts+m_nx; i < m_nx_tot; ++i) {
            MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
            MatSetValue(m_A, eqn, eqn-m_nx, -1.0, INSERT_VALUES);
            ++eqn;
         }
      }
      // Ghosts at y>yb.
      for (int j = m_nghosts+m_ny; j < m_ny_w_bdy; ++j) {
         for (int i = 0; i < m_nghosts; ++i) {
            MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
            MatSetValue(m_A, eqn, eqn-m_nx_tot*m_ny+m_nx, -1.0, INSERT_VALUES);
            ++eqn;
         }
         for (int i = m_nghosts; i < m_nghosts+m_nx; ++i) {
            MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
            MatSetValue(m_A, eqn, eqn-m_nx_tot*m_ny, -1.0, INSERT_VALUES);
            ++eqn;
         }
         for (int i = m_nghosts+m_nx; i < m_nx_tot; ++i) {
            MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
            MatSetValue(m_A, eqn, eqn-m_nx_tot*m_ny-m_nx, -1.0, INSERT_VALUES);
            ++eqn;
         }
      }
      // Extra equation.
      for (int i = 0; i < m_nghosts; ++i) {
         MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
         MatSetValue(m_A, eqn, eqn-m_nx_tot*m_ny+m_nx, -1.0, INSERT_VALUES);
         ++eqn;
      }
      for (int i = m_nghosts; i < m_nghosts+m_nx; ++i) {
         MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
         MatSetValue(m_A, eqn, eqn-m_nx_tot*m_ny, -1.0, INSERT_VALUES);
         ++eqn;
      }
      for (int i = m_nghosts+m_nx; i < m_nx_w_bdy; ++i) {
         MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
         MatSetValue(m_A, eqn, eqn-m_nx_tot*m_ny-m_nx, -1.0, INSERT_VALUES);
         ++eqn;
      }
      int loc = m_nghosts*m_nx_tot+m_nghosts;
      for (int j = 0; j < m_ny; ++j) {
         for (int i = 0; i < m_nx; ++i) {
            if (loc+i == eqn-m_nx_tot*m_ny-m_nx) {
               MatSetValue(m_A, eqn, loc+i, -1.0, INSERT_VALUES);
            }
            else {
               MatSetValue(m_A, eqn, loc+i, 1.0, INSERT_VALUES);
            }
         }
         loc += m_nx_tot;
      }
      MatSetValue(m_A, eqn, eqn, 1.0, INSERT_VALUES);
   }

   // Assemble the linear system matrix.
   MatAssemblyBegin(m_A, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(m_A, MAT_FINAL_ASSEMBLY);
}

  
void
LokiPoissonSolveSuperLU::createLinearSystemContext()
{
   // Create linear solver context.
   KSPCreate(PETSC_COMM_WORLD, &m_ksp);
   KSPSetOperators(m_ksp, m_A, m_A, DIFFERENT_NONZERO_PATTERN);

   // Specify SuperLU solver.
   KSPSetType(m_ksp, KSPPREONLY);
   PC  pc;
   KSPGetPC(m_ksp, &pc);
   PCSetType(pc, PCLU);
   PCFactorSetMatSolverPackage(pc, MATSOLVERSUPERLU);
   PCFactorSetUpMatSolverPackage(pc);
   KSPSetFromOptions(m_ksp);
}

} // end namespace Loki
