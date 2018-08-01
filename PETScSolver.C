#ifdef OVERTURE_USE_PETSC

// Loki private version of PETScSolver.C supplied by Overture.
// The version supplied by Overture does not properly support superlu or
// superlu_dist via PETSc.  This file addresses those and a few other
// errors.  However, it would be wise to completely rewrite this as
// it is exceptionally difficult to understand and most likely has other
// issues lurking.

// **************************************************************************
// *********** PETScSolver: Oges interface to parallel PETsc ****************
// ***********         PARALLEL VERSION                      ****************
// **************************************************************************

#include "PETScSolver.h"
#include "display.h"
#include "SparseRep.h"

// =============================================================================
//  \brief Here is the function that Overture::finish() calls to shutdown PETSc
// =============================================================================
static void 
finalizePETSc()
{
  #ifdef OVERTURE_USE_PETSC
  if( Oges::OGES_COMM_WORLD!=MPI_COMM_NULL )
  {
    // const int myid=max(0,Communication_Manager::My_Process_Number);
    // printf("--PETScSolver:  call finalizePETSc, myid=%i\n",myid);
    int ierr = PetscFinalize(); 
  }
  #endif
}

// ========================================================================================
/// \brief Call this function (before using Oges) if you want to use the PETSc solvers
// =======================================================================================
void
initPETSc()
{
  Oges::petscIsAvailable=true;                            // set to true if PETSc is available
  Oges::createPETSc=PETScSolver::newPETScSolver;  // pointer to a function that can "new" a PETSc instance
  Overture::shutDownPETSc = &finalizePETSc;               // set the function that will shut down PETSc
}


// =======================================================================================
/// \brief This function is called by buildEquationSolvers to create a PETScEquationSolver
// =======================================================================================
EquationSolver* PETScSolver::newPETScSolver(Oges &oges)
{
  return new PETScSolver(oges);
}




int PETScSolver::debug=0;
int PETScSolver::instancesOfPETSc=0;

static int NUM_KSP_CONVERGED_ITS_ERRORS=0;  // counts number of these errors

// static int numberOfKsp=0;
// static int numberOfMats=0;
// static int numberOfVects=0;

// From help:
// SuperLU_Dist Options -------------------------------------------------
//   -mat_superlu_dist_r <1>: Number rows in processor partition (None)
//   -mat_superlu_dist_c <1>: Number columns in processor partition (None)
//   -mat_superlu_dist_matinput <0>: Matrix input mode (0: GLOBAL; 1: DISTRIBUTED) (None)
//   -mat_superlu_dist_equil: <true> Equilibrate matrix (None)
//   -mat_superlu_dist_rowperm <LargeDiag> (one of) LargeDiag NATURAL
//   -mat_superlu_dist_colperm <MMD_AT_PLUS_A> (one of) MMD_AT_PLUS_A NATURAL MMD_ATA COLAMD
//   -mat_superlu_dist_replacetinypivot: <true> Replace tiny pivots (None)
//   -mat_superlu_dist_iterrefine: <false> Use iterative refinement (None)
//   -mat_superlu_dist_statprint: <true> Print factorization information (None)

#define FOR_3D(i1,i2,i3,I1,I2,I3) \
    int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  \
    int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); \
    for(i3=I3Base; i3<=I3Bound; i3++) \
    for(i2=I2Base; i2<=I2Bound; i2++) \
    for(i1=I1Base; i1<=I1Bound; i1++)

PETScSolver::
PETScSolver(Oges & oges_) : EquationSolver(oges_) 
{
  initialized=false;
  reInitialize=false;
  
  pnab=NULL;
  pnoffset=NULL;
  numberOfProcessors=Communication_Manager::Number_Of_Processors;
  numberOfGridPoints=0;
  numberOfGridPointsThisProcessor=0;
  relativeTol=1.e-4;
  numberOfComponents=1;
  pCoeff=NULL;
  problemIsSingular=notSingular;
  nullVector=NULL;
  diagonalScale=NULL; // for scaling the equations

  turnOnPETScMemoryTracing=false;
  
  solverMethod              =-1;
  preconditioner            =-1;   // initialize to an invalid value so we assign it later
  matrixOrdering            =-1;
  numberOfIncompleteLULevels=-1;
  gmresRestartLength        =-1;
  // parameters.rescaleRowNorms
  // old: useDiagonalScaling=false; // true;
  adjustPeriodicCoefficients=true;
  
  processorIsActive=true;       // set to true if this processor is in PETSc communicator
  allProcessorsAreActive=true;  // allProcessorsAreActive==true if all processors are in MPI_COMM_WORLD   
}

PETScSolver::
~PETScSolver()
{
  destroy();  // we need to destroy before closing down PETSc

  finalizePETSc(); 
  
  delete [] pnab;
  delete [] pnoffset;
  delete [] nullVector;
  delete diagonalScale;
}

int PETScSolver::
destroy()
{

  if( !processorIsActive )
    return 0;

  int ierr;
  /*
    Free work space.  All PETSc objects should be destroyed when they
    are no longer needed.
  */
//   printF("VectDestroy x, vect object %i.\n",numberOfVects);
//   numberOfVects--;
  ierr = VecDestroy(&x);   CHKERRQ(ierr);

//   printF("VectDestroy b, vect object %i.\n",numberOfVects);
//   numberOfVects--;
  ierr = VecDestroy(&b);   CHKERRQ(ierr);  

//   printF("MatDestroy object %i.\n",numberOfMats);
//   numberOfMats--;
  ierr = MatDestroy(&A);   CHKERRQ(ierr);
  assert( A==NULL );

  // printF("KSPDestroy object %i.\n",numberOfKsp);
  // numberOfKsp--;

  ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
  if( problemIsSingular!=notSingular && nullVector!=NULL )
  {
//     printF("VectDestroy nullVector[0], vect object %i.\n",numberOfVects);
//     numberOfVects--;
    ierr = VecDestroy(&nullVector[0]);  CHKERRQ(ierr);
  }
  delete diagonalScale; diagonalScale=NULL;
  
}


static char _p_ov_help[]="PETSc being used by the Overture `Oges' equation solver.";
int PETScSolver::
initializePETSc()
// ===================================================================================
//   /Description: Initialize PETSc (if necessary)
// ===================================================================================
{
  if( !initialized )
  { 
    // *new* way 10026
    Overture::incrementReferenceCountForPETSc();

    // We are another instance of PETSc, increment the global counter:
    instancesOfPETSc++;
  
    const int myid=max(0,Communication_Manager::My_Process_Number);

    if( !PetscInitializeCalled ) 
    {
      int ierr;
      if( Oges::debug & 2 && Communication_Manager::My_Process_Number<=0)
         printf("PETScSolver::initializePETSc, oges.argc=%i...\n",oges.argc);
      // int argc=0;
      // char **args=NULL;
       
      // tell petsc to use the communicator we have already initialized:
      // 2.2.1 PetscSetCommWorld(MPI_COMM_WORLD); 

      // ----- Use the OGES_WORLD_COMM parallel communicator -----  *wdh* 2015/06/30  
      // PETSC_COMM_WORLD = MPI_COMM_WORLD;

      PETSC_COMM_WORLD = Oges::OGES_COMM_WORLD;   // Here is the PETSC WORLD 

      MPI_Comm & OGES_COMM = parameters.getCommunicator(); // Here is the sub-communicator for this instance of OGES

      if( OGES_COMM != MPI_COMM_NULL )    
      {
	processorIsActive=true;  // this process IS in the PETSc sub-communicator
	if( PETSC_COMM_WORLD == MPI_COMM_NULL )
	{
	  printf("--PETScSolver-- initializePETSc:ERROR: myid=%i is in OGES_COMM but not in OGES_COMM_WORLD\n",
		 myid);
	  OV_ABORT("ERROR");
	}
	
	if( true && Oges::debug & 4 )
	{
	  int size=-1, rank=-1; 
	  MPI_Comm_size(OGES_COMM,&size);
	  MPI_Comm_rank(OGES_COMM,&rank);
      
	  printf(" --PETScSolver-- initializePETSc: OGES_COMM: myid=%i size=%i, rank=%i, processorIsActive=%i.\n",myid,size,rank,(int)processorIsActive);
	}
      }
      else
      {
	processorIsActive=false; // this process is NOT in the PETSc sub-communicator
	if( true && Oges::debug & 4 )
	{
      	  printf(" --PETScSolver--initializePETSc:  PETSC_COMM_WORLD==NULL, myid=%i, processorIsActive=%i.\n",myid,(int)processorIsActive);
	}
      }
      

      if( PETSC_COMM_WORLD!=MPI_COMM_NULL )
      {
        // -- Initialize PETSc on any processor that lives in PETSC_COMM_WORLD --
	if( true && Oges::debug & 4 )
	{
	  int size=-1, rank=-1; 
	  MPI_Comm_size(PETSC_COMM_WORLD,&size);
	  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
      
	  printf(" --PETScSolver-- initializePETSc: PETSC_COMM_WORLD: myid=%i size=%i, rank=%i.\n",myid,size,rank);
	}
	
	PetscInitialize(&oges.argc,&oges.argv,(char *)0,_p_ov_help);

	if( turnOnPETScMemoryTracing )
	{
	  // Activate logging of PETSC's malloc call
	  // 2.2.1 ierr = PetscTrLog();CHKERRQ( ierr );
	  ierr = PetscMallocDumpLog(stdout); CHKERRQ( ierr );
	}
      }
      
      // allProcessorsAreActive==true if PETSC_COMM_WORLD == MPI_COMM_WORLD
      allProcessorsAreActive = ParallelUtility::getMinValue( (int)processorIsActive );
      
    }
  }
  
  return 0;
}

int PETScSolver::
finalizePETSc()
// =======================================================================
//  This function will close PETSc if it the last instance.
// =======================================================================
{
  if( initialized )
  {

    // *new* way 10026
    Overture::decrementReferenceCountForPETSc();

    instancesOfPETSc--;  // decrement the global count

    if( false &&   // *wdh* 10026 -- PETSc is now shut down by Overture::finish()
        instancesOfPETSc==0 && PetscInitializeCalled )
    {
      if( true || Oges::debug & 2 )
        printF("PETScSolver::INFO: PetscFinalize CALLED -- This should only happen at the end of the program.\n");

      int ierr;
      ierr = PetscFinalize();CHKERRQ(ierr);
    }

    initialized=false;  // do this to prevent errors from multiple calls to this function
  }
  return 0;
}

// =======================================================================================
/// \brief return the name of the PETSc solver (KSP method, PC type etc.) 
// =======================================================================================
aString PETScSolver::
getSolverName() const
{
  // extract a name of the solver -- here we get the name from petsc:
  aString name="";
  const int maxLen=100;
  char buff[maxLen+1];
  PetscBool     flg;
  PetscOptionsGetString(PETSC_NULL,"-ksp_type",buff,maxLen,&flg);
  if( flg )
    name = name + "ksp[" + buff + "]";
  else
    name = "ksp[unknown]";
  
  PetscOptionsGetString(PETSC_NULL,"-pc_type",buff,maxLen,&flg);
  aString pcType=buff;
  name = name + " pc[" + pcType;

  if( pcType=="hypre" )
  {
    PetscOptionsGetString(PETSC_NULL,"-pc_hypre_type",buff,maxLen,&flg);
    aString hypreType=buff;
    // if( hypreType=="boomeramg" ) hypreType="AMG";
    name = name + "-" + hypreType;
  }
  name = name + "]";
  

  if( pcType!="hypre" )
  {
    PetscOptionsGetString(PETSC_NULL,"-sub_ksp_type",buff,maxLen,&flg);
    name = name + " sub-ksp[" + buff + "]";

    PetscOptionsGetString(PETSC_NULL,"-sub_pc_type",buff,maxLen,&flg);
    aString subPCType=buff;
    if( subPCType=="ilu" )
    {
//        PetscOptionsGetString(PETSC_NULL,"-sub_pc_ilu_levels",buff,maxLen,&flg);
      PetscOptionsGetString(PETSC_NULL,"-sub_pc_factor_levels",buff,maxLen,&flg);
      if( flg )
	name = name + "-" + subPCType + "(" + buff + ")";
    }
    else
    {
      name = name + " sub-pc[" + subPCType + "]";
    }
      
  }

  // printF("--PETSc-- solver: %s\n",(const char*)name);

  return name;
}



int  PETScSolver::
setProblemIsSingular( SingularProblemEnum singularOption /* =specifyConstantNullVector */ )
// ===============================================================================
// /Description:
//    There are currently 3 ways to solve a singular system.
// 
// 0=no, 1=use special PETSc option, 2=give PETSc a null vector, 3=add an extra equation
// ===============================================================================
{
  problemIsSingular=singularOption; 
  return 0;
}


//\begin{>>EquationSolverInclude.tex}{\subsection{printStatistics}} 
int PETScSolver::
printStatistics(FILE *file /* = stdout */ ) const
//===================================================================================
// /Description:
//   Output any relevant statistics
//\end{>>EquationSolverInclude.tex}
//===================================================================================
{
  return 0;
}


real PETScSolver::
getCurrentMemoryUsage()
// =====================================================
// Return the current memory usage in Mb
// =====================================================
{
  PetscLogDouble mem;
  int ierr = PetscMemoryGetCurrentUsage(&mem);
  return real(mem/(1024*1024));
}

// static int numberOfProcessors; 
// static int *pnab,*pnoffset; 
  
#define nab(side,axis,p,grid) pnab[(side)+2*( (axis) + 3*( (p) + numberOfProcessors*( (grid) ) ) )]
#define ndab(axis,p,grid) (nab(1,axis,p,grid)-nab(0,axis,p,grid)+1)
#define noffset(p,grid) pnoffset[(p)+numberOfProcessors*(grid)]


int PETScSolver::
getGlobalIndex( int n, int *iv, int grid, int p ) const
// ===============================================================================
// /Description: Return the global index given the point, grid and processor
// ===============================================================================
{
  return  n + numberOfComponents*(
         (iv[0]-nab(0,axis1,p,grid))+ndab(0,p,grid)*(
          iv[1]-nab(0,axis2,p,grid) +ndab(1,p,grid)*(
	  iv[2]-nab(0,axis3,p,grid))) + noffset(p,grid) );
}

int PETScSolver::
getGlobalIndex( int n, int *iv, int grid, realArray & ug ) const
// ===============================================================================
// /Description: Return the global index given the point, grid and grid function
// ===============================================================================
{
  int p= ug.Array_Descriptor.findProcNum( iv );  // processor number
  return  n + numberOfComponents*(
         (iv[0]-nab(0,axis1,p,grid))+ndab(0,p,grid)*(
          iv[1]-nab(0,axis2,p,grid) +ndab(1,p,grid)*(
	  iv[2]-nab(0,axis3,p,grid))) + noffset(p,grid) ); 
}

// =====================================================================================
// \brief Return the maximum residual.
// =====================================================================================
real PETScSolver::
getMaximumResidual()
{
  // This is not the max norm!
  double rnorm=0.;
  if( processorIsActive )
    KSPGetResidualNorm(ksp,&rnorm);

  // all processors need to know the residual: 
  if( !allProcessorsAreActive )
    rnorm = ParallelUtility::getMaxValue(rnorm);
  
  return rnorm;
}

// =====================================================================================
// \brief Return the number of iterations used in the last solve.
// =====================================================================================
int PETScSolver::
getNumberOfIterations() const
{
  if( !allProcessorsAreActive )
  {
    // If not all processors are active we need to make sure that we return the correct
    // answer on non-active processors.
    oges.numberOfIterations=ParallelUtility::getMaxValue(oges.numberOfIterations);
  }
  return oges.numberOfIterations;
  
} 





int PETScSolver::
buildGlobalIndexing(CompositeGrid & cg, realCompositeGridFunction & uu )
// ============================================================================
//  Build the arrays needed for defining the global indexing scheme
//
//   nab(side,axis,p,grid) : bounds of grid on processor p.
// ============================================================================
{
 if( !processorIsActive )
    return 0;

  const int myid=max(0,Communication_Manager::My_Process_Number);

  int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2]; 
  int jv[3], &j1=jv[0], &j2=jv[1], &j3=jv[2]; 

  const int nabDimension=2*3*numberOfProcessors*cg.numberOfComponentGrids();
  delete [] pnab;  // *wdh* 091128 -- for refactor
  pnab = new int[nabDimension];  
  delete [] pnoffset; // *wdh* 091128 -- for refactor
  pnoffset = new int [numberOfProcessors*cg.numberOfComponentGrids()]; 


  numberOfGridPoints=0;
  numberOfGridPointsThisProcessor=0;

  int grid;
  for( grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {
    realArray & ug= uu[grid];
      
    int nd1a = ug.getBase(0), nd1b=ug.getBound(0);
    int nd2a = ug.getBase(1), nd2b=ug.getBound(1);
    int nd3a = ug.getBase(2), nd3b=ug.getBound(2);

    IndexBox uBox;
    for( int p=0; p<numberOfProcessors; p++ )
    {
      CopyArray::getLocalArrayBox( p, ug, uBox ); // find the array bounds on proc. p (no ghost)
      for( int axis=0; axis<3; axis++ )
      {
	nab(0,axis,p,grid)=uBox.base(axis);
	nab(1,axis,p,grid)=uBox.bound(axis);
      }
      if( debug & 1 )
	printF("PETScSolver::NEW: grid=%i p=%i local bounds=[%i,%i][%i,%i][%i,%i] \n",grid,p,
	       nab(0,axis1,p,grid),nab(1,axis1,p,grid),
	       nab(0,axis2,p,grid),nab(1,axis2,p,grid),
	       nab(0,axis3,p,grid),nab(1,axis3,p,grid));
    }
    
    numberOfGridPointsThisProcessor+=ndab(0,myid,grid)*ndab(1,myid,grid)*ndab(2,myid,grid);
    numberOfGridPoints+= (nd1b-nd1a+1)*(nd2b-nd2a+1)*(nd3b-nd3a+1);   
  } // end for grid


  int offset=0;
  for( int p=0; p<numberOfProcessors; p++ )
  {
    for( grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      noffset(p,grid)=offset;
      offset+=ndab(0,p,grid)*ndab(1,p,grid)*ndab(2,p,grid);  // add number of pts on this local array
    }
  }
  
  if( debug & 4 )
  {
    printf(" >>>>>>>>>>> myid=%i numberOfGridPointsThisProcessor=%i, numberOfGridPoints=%i<<<<<<<<<<<<<<<<<<<\n",
	   myid,numberOfGridPointsThisProcessor,numberOfGridPoints);
    fflush(0);
  }
  
  return 0;
}


// ============================================================================
/// \brief Print a description of the solver and options used.
/// \note: NOTE: the corresponding function in class Oges will first print some info,
///    before this routine is called.
// ============================================================================
int PETScSolver::
printSolverDescription( const aString & label, FILE *file /* = stdout */ ) const
{
  PetscBool flg;
  aString name;
  const int maxLen=100;
  char buff[maxLen+1];
  // PetscOptionsGetString(PETSC_NULL,"-ksp_type",buff,maxLen,&flg);
  PetscOptionsGetString(PETSC_NULL,"-pc_type",buff,maxLen,&flg);
  aString pcType=buff;
  if( pcType=="hypre" )
  {
    PetscOptionsGetString(PETSC_NULL,"-pc_hypre_type",buff,maxLen,&flg);
    name=buff;
    if( name=="boomeramg" )
    {
      aString cycle="?",threshold="?",coarsenType="?";
      PetscOptionsGetString(PETSC_NULL,"-pc_hypre_boomeramg_cycle_type",buff,maxLen,&flg); if( flg) cycle=buff;
      PetscOptionsGetString(PETSC_NULL,"-pc_hypre_boomeramg_strong_threshold",buff,maxLen,&flg); if( flg) threshold=buff;
      PetscOptionsGetString(PETSC_NULL,"-pc_hypre_boomeramg_coarsen_type",buff,maxLen,&flg); if( flg)coarsenType=buff;
      
      fPrintF(file," Hypre: AMG (boomeramg) cycle=%s, strong-threshold=%s, coarsen-type=%s \n",
              (const char*)cycle,
              (const char*)threshold,
	      (const char*)coarsenType);
    }

  }
  
  // Results fron -help (.petscrc)
  // HYPRE BoomerAMG Options
  // -pc_hypre_boomeramg_cycle_type <V> (choose one of) V W (None)
  // -pc_hypre_boomeramg_max_levels <25>: Number of levels (of grids) allowed (None)
  // -pc_hypre_boomeramg_max_iter <1>: Maximum iterations used PER hypre call (None)
  // -pc_hypre_boomeramg_tol <0>: Convergence tolerance PER hypre call (0.0 = use a fixed number of iterations) (None)
  // -pc_hypre_boomeramg_truncfactor <0>: Truncation factor for interpolation (0=no truncation) (None)
  // -pc_hypre_boomeramg_P_max <0>: Max elements per row for interpolation operator (0=unlimited) (None)
  // -pc_hypre_boomeramg_agg_nl <0>: Number of levels of aggressive coarsening (None)
  // -pc_hypre_boomeramg_agg_num_paths <1>: Number of paths for aggressive coarsening (None)
  // -pc_hypre_boomeramg_strong_threshold <0.25>: Threshold for being strongly connected (None)
  // -pc_hypre_boomeramg_max_row_sum <0.9>: Maximum row sum (None)
  // -pc_hypre_boomeramg_grid_sweeps_all <1>: Number of sweeps for the up and down grid levels (None)
  // -pc_hypre_boomeramg_grid_sweeps_down <1>: Number of sweeps for the down cycles (None)
  // -pc_hypre_boomeramg_grid_sweeps_up <1>: Number of sweeps for the up cycles (None)
  // -pc_hypre_boomeramg_grid_sweeps_coarse <1>: Number of sweeps for the coarse level (None)
  // -pc_hypre_boomeramg_relax_type_all <symmetric-SOR/Jacobi> (choose one of) Jacobi sequential-Gauss-Seidel seqboundary-Gauss-Seidel SOR/Jacobi backward-SOR/Jacobi  symmetric-SOR/Jacobi  l1scaled-SOR/Jacobi Gaussian-elimination      CG Chebyshev FCF-Jacobi l1scaled-Jacobi (None)
  // -pc_hypre_boomeramg_relax_type_down <symmetric-SOR/Jacobi> (choose one of) Jacobi sequential-Gauss-Seidel seqboundary-Gauss-Seidel SOR/Jacobi backward-SOR/Jacobi  symmetric-SOR/Jacobi  l1scaled-SOR/Jacobi Gaussian-elimination      CG Chebyshev FCF-Jacobi l1scaled-Jacobi (None)
  // -pc_hypre_boomeramg_relax_type_up <symmetric-SOR/Jacobi> (choose one of) Jacobi sequential-Gauss-Seidel seqboundary-Gauss-Seidel SOR/Jacobi backward-SOR/Jacobi  symmetric-SOR/Jacobi  l1scaled-SOR/Jacobi Gaussian-elimination      CG Chebyshev FCF-Jacobi l1scaled-Jacobi (None)
  // -pc_hypre_boomeramg_relax_type_coarse <Gaussian-elimination> (choose one of) Jacobi sequential-Gauss-Seidel seqboundary-Gauss-Seidel SOR/Jacobi backward-SOR/Jacobi  symmetric-SOR/Jacobi  l1scaled-SOR/Jacobi Gaussian-elimination      CG Chebyshev FCF-Jacobi l1scaled-Jacobi (None)
  // -pc_hypre_boomeramg_relax_weight_all <1>: Relaxation weight for all levels (0 = hypre estimates, -k = determined with k CG steps) (None)
  // -pc_hypre_boomeramg_relax_weight_level <1>: Set the relaxation weight for a particular level (weight,level) (None)
  // -pc_hypre_boomeramg_outer_relax_weight_all <1>: Outer relaxation weight for all levels (-k = determined with k CG steps) (None)
  // -pc_hypre_boomeramg_outer_relax_weight_level <1>: Set the outer relaxation weight for a particular level (weight,level) (None)
  // -pc_hypre_boomeramg_no_CF: <FALSE> Do not use CF-relaxation (None)
  // -pc_hypre_boomeramg_measure_type <local> (choose one of) local global (None)
  // -pc_hypre_boomeramg_coarsen_type <Falgout> (choose one of) CLJP Ruge-Stueben  modifiedRuge-Stueben   Falgout  PMIS  HMIS (None)
  // -pc_hypre_boomeramg_interp_type <classical> (choose one of) classical   direct multipass multipass-wts ext+i ext+i-cc standard standard-wts   FF FF1 (None)
  // -pc_hypre_boomeramg_print_statistics: Print statistics (None)
  // -pc_hypre_boomeramg_print_debug: Print debug information (None)
  // -pc_hypre_boomeramg_nodal_coarsen: <FALSE> HYPRE_BoomerAMGSetNodal() (None)
  // -pc_hypre_boomeramg_nodal_relaxation: <FALSE> Nodal relaxation via Schwarz (None)

  return 0;
}

int PETScSolver::
buildMatrix( realCompositeGridFunction & coeff, realCompositeGridFunction & uu )
// ===========================================================================================
// /Description:
//    Build the PETSc matrix A from a coefficient grid function
// 
// ==========================================================================================
{
  if( !processorIsActive )
    return 0;

  debug=Oges::debug;
  
  MPI_Comm & OGES_COMM = parameters.getCommunicator();

  CompositeGrid & cg = *coeff.getCompositeGrid();
  const int myid=Communication_Manager::My_Process_Number;
  real time=getCPU();
  
  if( debug & 1 )
    printF("PETScSolver:: build matrix... Oges::debug=%i\n",Oges::debug);

  if( false ) 
    printF("*********** PETScSolver:: build matrix : parameters.rescaleRowNorms=%i ****************\n",parameters.rescaleRowNorms);
  

  if( debug & 4 )
  {
    coeff.display("PETScSolver::buildMatrix: Here is coeff on input","%4.1f ");
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      display(coeff[grid].sparse->classify,sPrintF("classify on grid=%i",grid),"%3i ");
    }
  }
  

  buildGlobalIndexing( cg,uu );

  pCoeff=&coeff;  // save this -- we need the classify array in solve

  int grid0=0;
  numberOfComponents=coeff[grid0].sparse->numberOfComponents;
  const int stencilSize=coeff[grid0].sparse->stencilSize;
  const int stencilDim=stencilSize*numberOfComponents;  // this is how many coeff's per equation

  // printF("PETScSolver:: build matrix: stencilSize=%i\n",stencilSize);

  const int fullStencilDimension=coeff[0].getLength(0);
  
  
  int grid;
//   int ngp=0;
//   for( grid=0; grid<cg.numberOfComponentGrids(); grid++ )
//   {  
//     realArray & coeffg= coeff[grid];
//     int nd1a = coeffg.getBase(1), nd1b=coeffg.getBound(1);
//     int nd2a = coeffg.getBase(2), nd2b=coeffg.getBound(2);
//     int nd3a = coeffg.getBase(3), nd3b=coeffg.getBound(3);

//     ngp+= (nd1b-nd1a+1)*(nd2b-nd2a+1)*(nd3b-nd3a+1);
//   }
//   numberOfGridPoints=ngp;
  
  if( Oges::debug & 2 )
    printF("PETScSolver:INFO: Total number of grid points is numberOfGridPoints=%i\n",numberOfGridPoints);
  

  numberOfUnknowns=numberOfGridPoints*numberOfComponents;
  numberOfUnknownsThisProcessor=numberOfGridPointsThisProcessor*numberOfComponents;
  

  int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2]; 
  int jv[3], &j1=jv[0], &j2=jv[1], &j3=jv[2]; 

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
         Compute the matrix and right-hand-side vector that define
         the linear system, Ax = b.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /* 
     Create parallel matrix, specifying only its global dimensions.
     When using MatCreate(), the matrix format can be specified at
     runtime. Also, the parallel partitioning of the matrix is
     determined by PETSc at runtime.

     Performance tuning note:  For problems of substantial size,
     preallocation of matrix memory is crucial for attaining good 
     performance.  Since preallocation is not possible via the generic
     matrix creation routine MatCreate(), we recommend for practical 
     problems instead to use the creation routine for a particular matrix
     format, e.g.,
         MatCreateMPIAIJ() - parallel AIJ (compressed sparse row)
         MatCreateMPIBAIJ() - parallel block AIJ
     See the matrix chapter of the users manual for details.
  */
  //  ierr = MatCreate(OGES_COMM,PETSC_DECIDE,PETSC_DECIDE,numberOfUnknowns,numberOfUnknowns,&A);CHKERRQ(ierr);
  if( parameters.externalSolver!=OgesParameters::defaultExternalSolver )
  {
    // Is this needed???
    // 2.2.1 ierr = MatCreate(OGES_COMM,numberOfUnknownsThisProcessor,numberOfUnknownsThisProcessor,
    //                  numberOfUnknowns,numberOfUnknowns,&A);CHKERRQ(ierr);
    ierr = MatCreate(OGES_COMM,&A);CHKERRQ(ierr);
    ierr = MatSetSizes(A,numberOfUnknownsThisProcessor,numberOfUnknownsThisProcessor,numberOfUnknowns,numberOfUnknowns);CHKERRQ(ierr);
  }
  else
  {
    PetscInt blockSize=parameters.blockSize;
    PetscBool optionWasSet;
    PetscOptionsGetInt(PETSC_NULL,"-mat_block_size",&blockSize,&optionWasSet);
    if( optionWasSet )
    {
      if( Oges::debug & 2 )
	printF("PETScSolver:Using -mat_block_size option: block size = %i\n"
	       "                    This option over-rides the value in OgesParameters\n",blockSize);
    }
    else
    {
      blockSize=parameters.blockSize;
      if( Oges::debug & 2 )
	printF("PETScSolver:Using blockSize=%i from OgesParameters\n",blockSize);
    }

    // ****** fix these:
    int d_nz = fullStencilDimension;  // expected number of non-zero entries on this processor ("diagonal block")
    int o_nz = 2;  // expected number of non-zero entries off processor
    int *d_nnz=PETSC_NULL;
    int *o_nnz=PETSC_NULL;
   
    // --- We should first determine the actual number of non-zeros in each row to be more efficient ----

    if( blockSize==1 )
    {
      // 	numberOfMats++;
      // 	printF("MatCreateMPIAIJ: create object %i.\n",numberOfMats);

      if( true )
      {
	// d_nz = fullStencilDimension;  // expected number of non-zero entries on this processor ("diagonal block")
	// o_nz = fullStencilDimension;
	ierr = MatCreateAIJ(OGES_COMM,numberOfUnknownsThisProcessor,numberOfUnknownsThisProcessor,
			    numberOfUnknowns,numberOfUnknowns,
			    d_nz,d_nnz,o_nz,o_nnz,
			    &A); CHKERRQ(ierr);
        // This next line is needed to avoid error when malloc'ing an additional entry that was more
        // than the estimated number
        MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
      }
      else
      {
	// v 2.3.2
	// ierr = MatCreateMPIAIJ(OGES_COMM,numberOfUnknownsThisProcessor,numberOfUnknownsThisProcessor,
	// 		       numberOfUnknowns,numberOfUnknowns,
	// 		       d_nz,d_nnz,o_nz,o_nnz,
	ierr = MatCreateAIJ(OGES_COMM,numberOfUnknownsThisProcessor,numberOfUnknownsThisProcessor,
			    numberOfUnknowns,numberOfUnknowns,
			    d_nz,d_nnz,o_nz,o_nnz,
			    &A); CHKERRQ(ierr);
      }

      // // PETSc documentation recommends doing this instead **FINISH ME**
      // ierr = MatCreate(OGES_COMM,&A);  CHKERRQ(ierr);
      // ierr = MatSetType(A,MATMPIAIJ);  CHKERRQ(ierr);
      // ierr = MatMPIAIJSetPreallocation(A,d_nz,d_nnz,o_nz,o_nnz);  CHKERRQ(ierr);
      
    }
    else
    {
      // for block matrices:
      // bs=block size
      // nz=number of non-zero blocks per block row (if same per row)
      // nnz[] =number of non-zero blocks per block row (different for each row)

      printF("\n *** PETScSolver: build a block matrix BAIJ with blockSize=%i **** \n",blockSize);

      // v 2.3.2
      // ierr = MatCreateMPIBAIJ(PETSC_COMM_SELF,blockSize,numberOfUnknownsThisProcessor,numberOfUnknownsThisProcessor,
      // 			     numberOfUnknowns,numberOfUnknowns,d_nz,d_nnz,o_nz,o_nnz,&A); CHKERRQ(ierr);

      ierr = MatCreateBAIJ(PETSC_COMM_SELF,blockSize,numberOfUnknownsThisProcessor,numberOfUnknownsThisProcessor,
			   numberOfUnknowns,numberOfUnknowns,d_nz,d_nnz,o_nz,o_nnz,&A); CHKERRQ(ierr);

    }
    
  }
  
  if( parameters.parallelExternalSolver==OgesParameters::superlu_dist ||
      parameters.externalSolver==OgesParameters::superlu )
  {
    // v2.3.2 ierr = MatSetType(A,MATSUPERLU_DIST);CHKERRQ(ierr);
  }
  else if( parameters.parallelExternalSolver!=OgesParameters::defaultExternalSolver &&
          parameters.parallelExternalSolver!=OgesParameters::hypre )
  {
    printf("*** PETScSolver:ERROR: un-expected externalSolver=%i ***\n",parameters.externalSolver);
  }
  
  
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  if( oges.getCompatibilityConstraint() &&
      problemIsSingular==notSingular )
  {
    problemIsSingular=addExtraEquation;
    if( Oges::debug & 1 )
      printF("PETScSolver:buildMatrix:INFO: problem is singular -- using addExtraEquation option\n");
  }
  

  int extraEquation=-1;
  if( problemIsSingular==addExtraEquation )
  {
    // find a spot to put the extra equation
    int grid=cg.numberOfComponentGrids()-1;
    const IntegerArray & d = cg[grid].dimension();

    i1=d(1,0), i2=d(1,1), i3=d(1,2); 
    int p= uu[grid].Array_Descriptor.findProcNum( iv );  // processor number
    int n=numberOfComponents-1;
    extraEquation=getGlobalIndex( n, iv, grid, p );  // get the global index

  }
  


  /* 
     Currently, all PETSc parallel matrix formats are partitioned by
     contiguous chunks of rows across the processors.  Determine which
     rows of the matrix are locally owned. 
  */
  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);

  /* 
     Set matrix elements for the 2-D, five-point stencil in parallel.
      - Each processor needs to insert only elements that it owns
        locally (but any non-local elements will be sent to the
        appropriate processor during matrix assembly). 
      - Always specify global rows and columns of matrix entries.

     Note: this uses the less common natural ordering that orders first
     all the unknowns for x = h then for x = 2h etc; Hence you see J = I +- n
     instead of J = I +- m as you might expect. The more standard ordering
     would first do all variables for y = h, then y = 2h etc.

   */
  real coeffScale=1.;  // fix this 
  const real eps= coeffScale*REAL_EPSILON*100.;  // cutoff tolerance for keeping coefficients

  if( parameters.rescaleRowNorms )
  {
    Range all;
    delete diagonalScale;  // *wdh* 091128 
    diagonalScale = new realCompositeGridFunction(cg,all,all,all,numberOfComponents);
  }
  
  for( grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {  
    const IntegerArray & gid = cg[grid].gridIndexRange();
    const IntegerArray & ir = cg[grid].indexRange();

    realArray & coeffg= coeff[grid];
    realSerialArray coeffLocal; getLocalArrayWithGhostBoundaries(coeffg,coeffLocal); 
    SparseRepForMGF & sparseRep = *coeff[grid].sparse;
    

    realArray & ug= uu[grid];
    realSerialArray uLocal; getLocalArrayWithGhostBoundaries(ug,uLocal);

    intSerialArray mask; getLocalArrayWithGhostBoundaries(cg[grid].mask(),mask);

    IntegerDistributedArray & equationNumberX = coeff[grid].sparse->equationNumber;
    intSerialArray equationNumber; getLocalArrayWithGhostBoundaries(equationNumberX,equationNumber);

    intSerialArray classify; getLocalArrayWithGhostBoundaries(coeff[grid].sparse->classify,classify);
    realSerialArray ds;
    if( parameters.rescaleRowNorms )
    {
      getLocalArrayWithGhostBoundaries((*diagonalScale)[grid],ds);
    }
    

    int n1a = uLocal.getBase(0) +ug.getGhostBoundaryWidth(0); 
    int n1b = uLocal.getBound(0)-ug.getGhostBoundaryWidth(0);
    int n2a = uLocal.getBase(1) +ug.getGhostBoundaryWidth(1); 
    int n2b = uLocal.getBound(1)-ug.getGhostBoundaryWidth(1);
    int n3a = uLocal.getBase(2) +ug.getGhostBoundaryWidth(2); 
    int n3b = uLocal.getBound(2)-ug.getGhostBoundaryWidth(2);

    if( debug & 4 )
      printf("PETScSolver:setMatrix: myid=%i n1a,n1b,n2a,n2b,n3a,n3b=[%i,%i][%i,%i][%i,%i]\n",
             myid,n1a,n1b,n2a,n2b,n3a,n3b);
    

#define arrayDims(grid,side,axis) cg[grid].dimension(side,axis)
#define arraySize(grid,axis) (arrayDims(grid,1,axis)-arrayDims(grid,0,axis)+1)

#define FOR3N(i1,i2,i3,n,n1a,n1b,n2a,n2b,n3a,n3b)\
    for( i3=n3a; i3<=n3b; i3++ )\
    for( i2=n2a; i2<=n2b; i2++ )\
    for( i1=n1a; i1<=n1b; i1++ )\
    for( n=0; n<numberOfComponents; n++ )

    int n,n1,eqn;     // component number

#define EQUATIONNUMBER(m,n,i1,i2,i3) equationNumber(m+stencilDim*(n),i1,i2,i3)
#define COEFF(m,n,i1,i2,i3) coeffLocal(m+stencilDim*(n),i1,i2,i3)

    const real eps=REAL_MIN*1000.;

    FOR3N(i1,i2,i3,n,n1a,n1b,n2a,n2b,n3a,n3b)
    {
      if( classify(i1,i2,i3,n)==-1 ) continue; // interp. pt, these are done later --- watch out for periodic

      // 	int p= ug.Array_Descriptor.findProcNum( iv );  // processor number
      // 	assert( p==myid );

      int ig=getGlobalIndex( n, iv, grid, myid );  // get the global index

      if( classify(i1,i2,i3,n)==0 )
      {
	// unused
	int jg=ig;
	v=1.;      // eqn at unused pts is the identity
        if( ig!=extraEquation )
	{
	  ierr = MatSetValues(A,1,&ig,1,&jg,&v,INSERT_VALUES);CHKERRQ(ierr);
	}
      }
      else
      {
	real dScale;
	if( parameters.rescaleRowNorms )
	{
          dScale=0.;
	  for( int m=0; m<stencilDim; m++ )
	  {
	    dScale=max(dScale,fabs(COEFF(m,n,i1,i2,i3)));
	  }
          if( dScale>eps )
            dScale=1./dScale;
	  else
            dScale=1.;  // all coefficients are small, just scale by 1.
          
          ds(i1,i2,i3,n)=dScale;  // save inverse of the scale factor
	}
	else
	{
	  dScale=1.;
	}

	if( i1<ir(0,0) || i1>ir(1,0) || 
            i2<ir(0,1) || i2>ir(1,1) || 
            i3<ir(0,2) || i3>ir(1,2) ) // *wdh* 051118
	{
	  // ***** fill-in equations at ghost points or periodic points *********

	  // int jg=ig;
	  // v=1.;      // eqn at ghost points is the identity for now
	  // ierr = MatSetValues(A,1,&ig,1,&jg,&v,INSERT_VALUES);CHKERRQ(ierr);
	  // printf("p=%i: i1,i2=%i,%i   ig,jg=%i,%i  set to identity.\n",myid,i1,i2,ig,jg);
	  
	  if( debug & 8 ) 
	    printf(" Boundary Eqn: grid=%i classify=%i (i1,i2,i3)=(%i,%i,%i) mask=%i:\n",
		   grid,classify(i1,i2,i3,n),i1,i2,i3,mask(i1,i2,i3)); 

	

	  for( int m=0; m<stencilDim; m++ )
	  {
	    if( COEFF(m,n,i1,i2,i3)!=0. )
	    {
	      // eqn = equationNumber(m,i1,i2,i3) -1;
	      eqn = EQUATIONNUMBER(m,n,i1,i2,i3) -1;

	      sparseRep.equationToIndex(eqn,n1,j1,j2,j3); 

		
	      int p= ug.Array_Descriptor.findProcNum( jv );  // processor number
	      int jg=getGlobalIndex( n1, jv, grid, p );  // get the global index


	      v=COEFF(m,n,i1,i2,i3)*dScale; 
	      if( fabs(v)<eps ) continue;
	    
	      ierr = MatSetValues(A,1,&ig,1,&jg,&v,INSERT_VALUES);CHKERRQ(ierr);

	      if( debug & 8 ) 
		printf("      m=%i n=%i j=(%i,%i,%i) value=%6.2f [eqn=%i,ig=%i,jg=%i],\n",
		       m,n,j1,j2,j3,COEFF(m,n,i1,i2,i3),EQUATIONNUMBER(m,n,i1,i2,i3) -1,ig,jg);

	    }
	  }


	}
	else 
	{
	  for( int m=0; m<stencilDim; m++ )
	  {
	    // compute the offset for this stencil pt  --------------- fix this ---------------
	    // int m1 = (m % 3) -1;
	    // int m2 = int(m/3) -1;
	  
	    if( COEFF(m,n,i1,i2,i3)!=0. )
	    {

	      // int n1=0;
	      // j1= i1+m1;
	      // j2= i2+m2;
	      // j3= i3;
		
	      eqn = EQUATIONNUMBER(m,n,i1,i2,i3) -1;

	      sparseRep.equationToIndex(eqn,n1,j1,j2,j3); 
	      int p= ug.Array_Descriptor.findProcNum( jv );  // processor number

              if( adjustPeriodicCoefficients && p==myid && classify(j1,j2,j3,n1)==-2 )
	      {
                // Here we explicitly change the matrix to use the periodic image that is inside the grid

                // printf(" grid=%i, i=(%i,%i,%i), j=(%i,%i,%i) --> periodic pt",grid,i1,i2,i3,j1,j2,j3);
		
		// This is a periodic point
                // eqn = EQUATIONNUMBER(1,n1,j1,j2,j3) -1;  // *wdh* 061014  bug found for systems
		// The periodic point must be located at m=1: 
                eqn = EQUATIONNUMBER(1,n1,j1,j2,j3) -1;
		sparseRep.equationToIndex(eqn,n1,j1,j2,j3); 
		p= ug.Array_Descriptor.findProcNum( jv );  // processor number

                // printf(" j=(%i,%i,%i)\n",j1,j2,j3);
		
	      }
	      
	      int jg=getGlobalIndex( n1, jv, grid, p );  // get the global index

	      v=COEFF(m,n,i1,i2,i3)*dScale;

	      if( fabs(v)<eps ) continue;
	      ierr = MatSetValues(A,1,&ig,1,&jg,&v,INSERT_VALUES);CHKERRQ(ierr);


	      // printf("p=%i: i1,i2=%i,%i   ig,jg=%i,%i  m=%i value=%6.4f\n",myid,i1,i2,ig,jg,m,v);

	    }
	  }
	  if( problemIsSingular==addExtraEquation )
	  {
	    // add entries to this row and to the extra equation
	    v=1.;
	    ierr = MatSetValues(A,1,&ig,1,&extraEquation,&v,INSERT_VALUES);CHKERRQ(ierr);
	    ierr = MatSetValues(A,1,&extraEquation,1,&ig,&v,INSERT_VALUES);CHKERRQ(ierr);
	  }

	}
      }
    }
    
  }  // end for grid

  fillInterpolationCoefficients(uu); 

  
  /* 
     Assemble matrix, using the 2-step process:
       MatAssemblyBegin(), MatAssemblyEnd()
     Computations can be done while messages are in transition
     by placing code between these two statements.
  */
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  if( false && myid==0 )
  {
    for( int ig=0; ig<numberOfUnknowns; ig++ )
    {
      ierr = MatGetValues(A,1,&ig,1,&ig,&v);CHKERRQ(ierr);
      if( fabs(v)<1.e-10 )
      {
	printf("WARNING: matrix diagonal a(%i,%i)=%e is very small\n",ig,ig,v);
      }
    }
  }

  if( debug & 4 )
  {
    printF("============= PETScSolver: Here is the sparse matrix ==============\n");
    ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }

  if( false )
  { // save the matrix in a file that can me read in matlab by typing 'petscMatrix'
    printF("Saving matrix in matlab format file=`petscMatrix.m'\n");
    PetscViewer viewer;
    PetscViewerASCIIOpen(OGES_COMM,"petscMatrix.m",&viewer);
    PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
    ierr = MatView(A,viewer);CHKERRQ(ierr);

  }
  
  /* 
     Create parallel vectors.
      - We form 1 vector from scratch and then duplicate as needed.
      - When using VecCreate(), VecSetSizes and VecSetFromOptions()
        in this example, we specify only the
        vector's global dimension; the parallel partitioning is determined
        at runtime. 
      - When solving a linear system, the vectors and matrices MUST
        be partitioned accordingly.  PETSc automatically generates
        appropriately partitioned matrices and vectors when MatCreate()
        and VecCreate() are used with the same communicator.  
      - The user can alternatively specify the local vector and matrix
        dimensions when more sophisticated partitioning is needed
        (replacing the PETSC_DECIDE argument in the VecSetSizes() statement
        below).
  */
//   numberOfVects++;
//   printF("VecCreate b, vect object %i.\n",numberOfVects);

  ierr = VecCreate(OGES_COMM,&b); CHKERRQ(ierr);
  //  ierr = VecSetSizes(b,PETSC_DECIDE,numberOfUnknowns);CHKERRQ(ierr);
  ierr = VecSetSizes(b,numberOfUnknownsThisProcessor,numberOfUnknowns); CHKERRQ(ierr);
  ierr = VecSetFromOptions(b); CHKERRQ(ierr);

//   numberOfVects++;
//   printF("VecDup x, vect object %i.\n",numberOfVects);
  ierr = VecDuplicate(b,&x); CHKERRQ(ierr);

  
  if( problemIsSingular==specifyNullVector )
  {
    // Create the null vector

//   cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
//   cout << " 1=interior, 2=bndry, 3=ghost1, 4=ghost2, -1=interp, -2=periodic, -3=extrap, 0=unused\n";
//   cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

    assert( nullVector==NULL );
    nullVector = new Vec [1];

//     numberOfVects++;
//     printF("VecDup nullVector[0], vect object %i.\n",numberOfVects);
    ierr = VecDuplicate(b,&nullVector[0]);CHKERRQ(ierr);
    
    PetscScalar *nv;
    VecGetArray(nullVector[0],&nv);  // get the local array from Petsc
    ierr = VecGetOwnershipRange(nullVector[0],&Istart,&Iend);CHKERRQ(ierr);

    for( grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {  
      MappedGrid & mg = cg[grid];
      intSerialArray classify; getLocalArrayWithGhostBoundaries(coeff[grid].sparse->classify,classify);
      const int *pir = &mg.indexRange(0,0);
      #define IR(side,axis) pir[side+2*(axis)]


      realArray & ug= uu[grid];
      realSerialArray uLocal; getLocalArrayWithGhostBoundaries(ug,uLocal); 
      int n1a = uLocal.getBase(0) +ug.getGhostBoundaryWidth(0); 
      int n1b = uLocal.getBound(0)-ug.getGhostBoundaryWidth(0);
      int n2a = uLocal.getBase(1) +ug.getGhostBoundaryWidth(1); 
      int n2b = uLocal.getBound(1)-ug.getGhostBoundaryWidth(1);
      int n3a = uLocal.getBase(2) +ug.getGhostBoundaryWidth(2); 
      int n3b = uLocal.getBound(2)-ug.getGhostBoundaryWidth(2);
      
      i1=n1a, i2=n2a, i3=n3a;
      int n=0;
      int ig=getGlobalIndex( n,iv, grid, myid );  // get the global index for the first point

      FOR3N(i1,i2,i3,n,n1a,n1b,n2a,n2b,n3a,n3b)
      {
	assert( ig>= Istart && ig<=Iend );
	  
	if( classify(i1,i2,i3,n)==0 ) 
//             i1<IR(0,0) || i1>IR(1,0) || i2<IR(0,1) || i2>IR(1,1) || i3<IR(0,2) || i3>IR(1,2) )
	{
	  nv[ig-Istart]=0.; // Null vector is zero at unused pts
	}
	else
	{
	  nv[ig-Istart]=1.; 
          if( debug & 2 ) printf(" Null vector: grid=%i, (i1,i2,i3)=(%i,%i,%i) ig=%i\n",grid,i1,i2,i3,ig);
	  
	}
	ig++;
      }
    }
    #undef IR
    VecRestoreArray(nullVector[0],&nv);
    if( debug & 8 )
    {
      if( myid==0 ) printf("PETScSolver::Null vector\n");
      ierr = VecView(nullVector[0],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    }
  
    if( debug & 16 )
    {
      ierr = MatMult(A,nullVector[0],b);CHKERRQ(ierr);
      if( myid==0 ) printf("PETScSolver::b = A*v[0]\n");
      ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    }
    
  }

  if( parameters.solveForTranspose )
  {
    if( debug & 16 )
    {
      printF("============= PETScSolver: Here is A BEFORE TRANSPOSE ==============\n");
      ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    }

    // Transpose the matrix
    printF("--- PETScSolver: TRANSPOSE the matrix ...\n");
    if( parameters.rescaleRowNorms )
    {
      printF("--- PETScSolver: ERROR: solve for transpose but rescaleRowNorms=true\n"
             "    This is very likely an error. To turn off row scaling use: \n"
             " solver.set(OgesParameters::THErescaleRowNorms,false);\n" );
      OV_ABORT("error");
    }
    
    real cpu0=getCPU();
    
    // PetscErrorCode  MatTranspose(Mat mat,MatReuse reuse,Mat *B)
    // -- transpose in place: 
    // typedef enum {MAT_INITIAL_MATRIX,MAT_REUSE_MATRIX} MatReuse;
    // MatTranspose( A, MAT_INITIAL_MATRIX, &A );
 #if (PETSC_VERSION_MAJOR==3)
    MatTranspose( A, MAT_INITIAL_MATRIX, &A );
#else
    MatTranspose( A, &A );
#endif

    if( debug & 16 )
    {
      printF("============= PETScSolver: Here is A AFTER TRANSPOSE ==============\n");
      ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    }

    printF("--- PETScSolver: .... done TRANSPOSE. cpu=%8.2e(s) ----\n",getCPU()-cpu0);

  }
  



  if( Oges::debug & 2 )
  {
    time=ParallelUtility::getMaxValue(getCPU()-time,-1,OGES_COMM);
    printF("PETScSolver:: ... done build matrix, cpu=%8.2e\n",time);
  }
  
  
  return 0;

}


//==============================================================================================
/// \brief Add options to PETSc's list of options
//==============================================================================================
int PETScSolver::
buildSolver()
{

  if( !processorIsActive )
    return 0;

  MPI_Comm & OGES_COMM = parameters.getCommunicator();

  ListOfShowFileParameters & petscOptions = oges.parameters.petscOptions; 
  std::list<ShowFileParameter>::iterator iter; 
  for(iter = petscOptions.begin(); iter!=petscOptions.end(); iter++ )
  {
    ShowFileParameter & param = *iter;
    aString name; ShowFileParameter::ParameterType type; int ivalue; real rvalue; aString stringValue;
    param.get( name, type, ivalue, rvalue, stringValue );

    if( type==ShowFileParameter::stringParameter )
    {
      if( debug & 1 )
	printF("PETScSolver::buildSolver: INFO: adding option=[%s] value=[%s]\n",(const char*)name,(const char*)stringValue);
      PetscOptionsSetValue(name,stringValue);

      if( name=="-pc_type" )
      {
        
      }
      

    }
//     else if( type==ShowFileParameter::realParameter )
//     {
//       // textStrings[nt]=sPrintF("%e",rvalue);
//     }
//     else if( type==ShowFileParameter::intParameter )
//     {
//       // textStrings[nt]=sPrintF("%i",ivalue);
//     }
    else
    {
      Overture::abort("error"); 
    }
  }


  /* 
     Create linear solver context
  */

  // numberOfKsp++;
  // printF("KSPCreate object %i.\n",numberOfKsp);
  
  ierr = KSPCreate(OGES_COMM,&ksp);CHKERRQ(ierr);

  /* 
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
  ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);

  /* 
     Set linear solver defaults for this problem (optional).
     - By extracting the KSP and PC contexts from the KSP context,
       we can then directly call any KSP and PC routines to set
       various options.
     - The following two statements are optional; all of these
       parameters could alternatively be specified at runtime via
       KSPSetFromOptions().  All of these defaults can be
       overridden at runtime, as indicated below.
  */

//   ierr = KSPSetTolerances(ksp,relativeTol/numberOfUnknowns,1.e-50,PETSC_DEFAULT,
//                           PETSC_DEFAULT);CHKERRQ(ierr);

  /* 
    Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
  */

  // PetscOptionsSetValue("-mg_coarse_ksp_max_it","6");
//   PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold",".5");
//   PetscOptionsSetValue("-pc_hypre_boomeramg_max_levels","20");
// -pc_hypre_boomeramg_coarsen_type <Falgout> (one of) CLJP Ruge-Stueben  modifiedRuge-Stueben   Falgout
//   PetscOptionsSetValue("-pc_hypre_boomeramg_coarsen_type","Falgout");
  // PetscOptionsSetValue("-pc_hypre_boomeramg_coarsen_type","modifiedRuge-Stueben");
//   PetscOptionsSetValue("-pc_hypre_boomeramg_coarsen_type","CLJP");
//   PetscOptionsSetValue("-pc_hypre_boomeramg_coarsen_type","Ruge-Stueben");
//  PetscOptionsSetValue("-pc_hypre_boomeramg_max_row_sum","1.");  // turn off
   
  // PetscOptionsSetValue("-pc_hypre_boomeramg_grid_sweeps","4, 4, 4, 4");

//   PetscOptionsSetValue("-pc_hypre_boomeramg_sweep_all","true");
  
  // set the parameters such as the preconditioner and krylov method
  setPetscParameters();


  /* 
    Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
  */
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  // ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  // ierr = KSPGetPC(ksp, &pc);    CHKERRQ(ierr);

  

/* ---

  KSPType krylovSpaceMethod;
 //  krylovSpaceMethod=KSPGMRES;
  krylovSpaceMethod=KSPBCGS;  // biConjugateGradientStabilized

//  ierr = KSPSetType(ksp, krylovSpaceMethod); CHKERRQ(ierr);

   ierr = KSPGetPC(ksp, &pc);    CHKERRQ(ierr);
   PCType  petscPreconditioner; 
   petscPreconditioner=PCILU;
//   ierr = PCSetType(pc,      petscPreconditioner);  CHKERRQ(ierr);

   MatOrderingType  matOrdering; 
   matOrdering=MATORDERING_RCM;  // reverseCuthillMcKeeOrdering:
   if( petscPreconditioner==PCILU )
   {
//     ierr = PCILUSetMatOrdering(pc, matOrdering); CHKERRQ(ierr);
   }
  
   ----- */

   if( parameters.parallelExternalSolver==OgesParameters::superlu_dist ||
       parameters.externalSolver==OgesParameters::superlu) // is this the right way to do this?
   {
     ierr = KSPSetType(ksp,KSPPREONLY); CHKERRQ(ierr);
     ierr = PCSetType(pc,PCLU);  CHKERRQ(ierr);

   }
  
  
   if( problemIsSingular==specifyConstantNullVector || problemIsSingular==specifyNullVector )
   {
     // Specify that null vector
     if( problemIsSingular==specifyNullVector )
     {
       int numberOfNullVectors=1;
       MatNullSpaceCreate (OGES_COMM,PETSC_FALSE,numberOfNullVectors,nullVector,&nsp); 
     }
     else if( problemIsSingular==specifyConstantNullVector )
     {
       MatNullSpaceCreate (OGES_COMM,PETSC_TRUE,0,NULL,&nsp); 
     }
     
     KSPSetNullSpace(ksp,nsp); 

     PetscReal damping=1.e-10;
     // 2.2.1 PCLUSetDamping(pc,damping);
#if (PETSC_VERSION_MAJOR==3)
     PCFactorSetShiftType(pc,MAT_SHIFT_NONZERO);
     PCFactorSetShiftAmount(pc,damping);
#else
     PCFactorSetShiftNonzero(pc,damping);
#endif
   }

//    // The initial guess is non-zero unless we use a direct solver
//    if( krylovSpaceMethod!=KSPPREONLY ) // petscPreconditioner!=PCLU )
//    {
//      ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE); CHKERRQ(ierr); 
//    }
//    else
//    {
//      ierr = KSPSetInitialGuessNonzero(ksp,PETSC_FALSE); CHKERRQ(ierr); 
//    }
   
   return 0;
  
}

static PCType getPetscPreconditioner( OgesParameters::PreconditionerEnum preconditioner )
// =========================================================================================
// Convert an OgesParameters preconditioner type into a Petsc type.
// =========================================================================================
{
  PCType petscPreconditioner;
  switch ( preconditioner )
  {
  case OgesParameters::noPreconditioner:
    petscPreconditioner=PCNONE ;
    break;
  case OgesParameters::jacobiPreconditioner:
    petscPreconditioner=PCJACOBI;
    break;
  case OgesParameters::sorPreconditioner:
    petscPreconditioner=PCSOR;
    break;
  case OgesParameters::luPreconditioner:
    petscPreconditioner=PCLU;
    break;
  case OgesParameters::shellPreconditioner:
    petscPreconditioner=PCSHELL;
    break;
  case OgesParameters::blockJacobiPreconditioner:
    petscPreconditioner=PCBJACOBI;
    break;
  case OgesParameters::multigridPreconditioner:
    petscPreconditioner=PCMG;
    break;
  case OgesParameters::eisenstatPreconditioner:
    petscPreconditioner=PCEISENSTAT;
    break;
  case OgesParameters::incompleteCholeskyPreconditioner:
    petscPreconditioner=PCICC;
    break;
  case OgesParameters::incompleteLUPreconditioner:
    petscPreconditioner=PCILU;
    break;
  case OgesParameters::additiveSchwarzPreconditioner:
    petscPreconditioner=PCASM;
    break;
  case OgesParameters::kspPreconditioner:
    petscPreconditioner=PCKSP;
    break;
  case OgesParameters::compositePreconditioner:
    petscPreconditioner=PCCOMPOSITE;
    break;
  case OgesParameters::redundantPreconditioner:
    petscPreconditioner=PCREDUNDANT;
    break;
  case OgesParameters::hyprePreconditioner:
    petscPreconditioner=PCHYPRE;
    break;
  default:
    printf("****WARNING**** Unknown preconditionner for PETSc\n");
    petscPreconditioner=PCILU;
  }
  return petscPreconditioner;
}

static KSPType getPetscKrylovSpaceMethod( OgesParameters::SolverMethodEnum solverMethod )
// =========================================================================================
// Convert an OgesParameters krylov space method into a Petsc type.
// =========================================================================================
{
    
  KSPType krylovSpaceMethod;

  switch( solverMethod )
  {
  case OgesParameters::richardson:
    krylovSpaceMethod=KSPRICHARDSON;
    break;
  case OgesParameters::chebychev:
    krylovSpaceMethod=KSPCHEBYSHEV;
    break;
  case OgesParameters::conjugateGradient:
    krylovSpaceMethod=KSPCG;
    break;
  case OgesParameters::gmres:
    krylovSpaceMethod=KSPGMRES;
    break;
  case OgesParameters::biConjugateGradientStabilized: // bcgs
    krylovSpaceMethod=KSPBCGS;
    break;
  case OgesParameters::tcqmr:
    krylovSpaceMethod=KSPTCQMR;
    break;
  case OgesParameters::tfqmr:
    krylovSpaceMethod=KSPTFQMR;
    break;
  case OgesParameters::conjugateResidual:
    krylovSpaceMethod=KSPCR;
    break;
  case OgesParameters::leastSquares:
    krylovSpaceMethod=KSPLSQR;
    break;
  case OgesParameters::preonly:
    krylovSpaceMethod=KSPPREONLY;
    break;
  case OgesParameters::qcg:
    krylovSpaceMethod=KSPQCG;
    break;
  case OgesParameters::biConjugateGradient:
    krylovSpaceMethod=KSPBICG;
    break;
  case OgesParameters::biConjugateGradientSquared:   // do this ***
    printf("PETScSolver: WARNING: no biConjugateGradientSquared, using conjugateGradientSquared\n");
    krylovSpaceMethod=KSPCGS;
    break;
  case OgesParameters::conjugateGradientSquared:
    krylovSpaceMethod=KSPCGS;
    break;
  default:
    printf("****getPetscKrylovSpaceMethod:WARNING: Unknown Krylov space method for PETSc\n");
    krylovSpaceMethod=KSPGMRES;
  }
  
  return krylovSpaceMethod;
}

#undef __FUNC__
#define __FUNC__ "PETScSolver::setPetscParameters"
int PETScSolver::
setPetscParameters() 
  //....Solver options (from PETSc, see the documentation, Chap. 4)
  // KSP Options:
  //   KSPRICHARDSON   = Richardson iter.
  //   KSPCHEBYSHEV    = Chebyshev iter.
  //   KSPCG           = Conj. Gradients (only for SPD systems)
  //   KSPGMRES        = Restarted GMRes
  //   KSPTCQMR        = Transp. Free QMR (T.F.Chan version); SLOW!
  //   KSPBCGS         = BiConj. Gradient method
  //   KSPCGS          = Conj. Gradient Squared
  //   KSPTFQMR        = Transp. Free QMR (Freund's version): Fast, _default_
  //   KSPCR           = Conj. Resid.
  //   KSPPREONLY      = Use only the preconditioner
{
  const int myid=Communication_Manager::My_Process_Number;

  bool parametersHaveChanged=false;  // set to true if any parameters have changed
  PetscBool optionWasSet;

  ierr = KSPGetPC(ksp, &pc);    CHKERRQ(ierr);

  // Choose the KSP type:
//  OgesParameters::SolverMethodEnum newSolverMethod = numberOfProcessors>1 ? parameters.parallelSolverMethod :
//                                                                            parameters.solverMethod;
  // *wdh* 061018 : always do this:
  OgesParameters::SolverMethodEnum newSolverMethod = parameters.parallelSolverMethod;
  PetscOptionsHasName(PETSC_NULL,"-ksp_type",&optionWasSet);  // check if user has set
  if( !optionWasSet && newSolverMethod )
  {
    parametersHaveChanged=true;
    KSPType krylovSpaceMethod = getPetscKrylovSpaceMethod(newSolverMethod);
    
    if( false && parameters.solveForTranspose ) // this didn't seem to work
    {
      krylovSpaceMethod = KSPPREONLY; 
      if (Oges::debug & 2) 
        printF("Solution of transpose requires krylovSpaceMethod = KSPPREONLY\n");
    }
    if( Oges::debug & 2 ) 
      printF(" ********** PETScSolver:INFO set krylov space method =%i (p=%i) ***********\n",newSolverMethod,myid);

    ierr = KSPSetType(ksp, krylovSpaceMethod); CHKERRQ(ierr);
    // solverMethod=parameters.solverMethod;
    
  }
  else
  {
    // *wdh* 070130 -- should we set parameters.parallelSolverMethod to match -ksp_type ?? -------------------------
  }
  

  // *** Here we set the preconditioner (or parallelPreconditioner for parallel jobs)
  OgesParameters::PreconditionerEnum newPreconditioner = numberOfProcessors>1 ? parameters.parallelPreconditioner :
                                                                                parameters.preconditioner;
  // PetscOptionsHasName(PETSC_NULL,"-pc_type",&optionWasSet);  // check if user has set
  
  // In parallel we always need to set the preconditioner so we can later set the sub_ksp etc. for block methods
  if( newPreconditioner!=preconditioner )
  {
    parametersHaveChanged=true;

    PCType petscPreconditioner = getPetscPreconditioner(newPreconditioner);

    const int maxLen=100;
    char buff[maxLen];
    PetscOptionsGetString(PETSC_NULL,"-pc_type",buff,maxLen,&optionWasSet);
    if( optionWasSet )
    {
      aString pcType=buff;

        // #define PCNONE            "none"
        // #define PCJACOBI          "jacobi"
        // #define PCSOR             "sor"
        // #define PCLU              "lu"
        // #define PCSHELL           "shell"
        // #define PCBJACOBI         "bjacobi"
        // #define PCMG              "mg"
        // #define PCEISENSTAT       "eisenstat"
        // #define PCILU             "ilu"
        // #define PCICC             "icc"
        // #define PCASM             "asm"
        // #define PCKSP             "ksp"
        // #define PCCOMPOSITE       "composite"
        // #define PCREDUNDANT       "redundant"
        // #define PCSPAI            "spai"
        // #define PCNN              "nn"
        // #define PCCHOLESKY        "cholesky"
        // #define PCGAMG            "gamg"
        // #define PCPBJACOBI        "pbjacobi"
        // #define PCMAT             "mat"
        // #define PCHYPRE           "hypre"
        // #define PCFIELDSPLIT      "fieldsplit"
        // #define PCTFS             "tfs"
        // #define PCML              "ml"
        // #define PCGALERKIN        "galerkin"

      petscPreconditioner=PCASM;
      bool found=true;
      if( pcType=="none" ) petscPreconditioner=PCNONE;
      else if( pcType=="jacobi" ) petscPreconditioner=PCJACOBI;
      else if( pcType=="sor" ) petscPreconditioner=PCSOR;
      else if( pcType=="lu" ) petscPreconditioner=PCLU;
      else if( pcType=="shell" ) petscPreconditioner=PCSHELL;
      else if( pcType=="bjacobi" ) petscPreconditioner=PCBJACOBI;
      else if( pcType=="mg" ) petscPreconditioner=PCMG;
      else if( pcType=="eisenstat" ) petscPreconditioner=PCEISENSTAT;
      else if( pcType=="ilu" ) petscPreconditioner=PCILU;
      else if( pcType=="icc" ) petscPreconditioner=PCICC;
      else if( pcType=="asm" ) petscPreconditioner=PCASM;
      else if( pcType=="ksp" ) petscPreconditioner=PCKSP;
      else if( pcType=="composite" ) petscPreconditioner=PCCOMPOSITE;
      else if( pcType=="redundant" ) petscPreconditioner=PCREDUNDANT;
      else if( pcType=="spai" ) petscPreconditioner=PCSPAI;
      else if( pcType=="nn" ) petscPreconditioner=PCNN;
      else if( pcType=="cholesky" ) petscPreconditioner=PCCHOLESKY;
      else if( pcType=="gamg" ) petscPreconditioner=PCGAMG;
      else if( pcType=="pbjacobi" ) petscPreconditioner=PCPBJACOBI;
      else if( pcType=="mat" ) petscPreconditioner=PCMAT;
      else if( pcType=="hypre" ) petscPreconditioner=PCHYPRE;
      else if( pcType=="fieldsplit" ) petscPreconditioner=PCFIELDSPLIT;
      else if( pcType=="tfs" ) petscPreconditioner=PCTFS;
      else if( pcType=="ml" ) petscPreconditioner=PCML;
      else if( pcType=="galerkin" ) petscPreconditioner=PCGALERKIN;
      else
      {
        found=false;
	printF("PETScSolver::setPetscParameters:ERROR: unknown result from -pc_type found in the "
	       "options database = [%s]\n",(const char*)pcType);
      }

//       printf("PETScSolver::setPetscParameters: pcType= [%s]\n",(const char*)pcType);
      if( found )
        printf("*** -pc_type found in the options database: [%s]\n",petscPreconditioner);

    }

    // In parallel we can use block-jacobi, additive-Schwatrz PC or Hypre

    if( numberOfProcessors==1 ||
        !strcmp(petscPreconditioner, PCBJACOBI) || 
        !strcmp(petscPreconditioner, PCASM) ||
        !strcmp(petscPreconditioner, PCLU) ||
        !strcmp(petscPreconditioner, PCHYPRE))
    {
      if( Oges::debug & 2 ) 
	printF(" ********** PETScSolver: set preconditioner ***********\n");
      ierr = PCSetType(pc,      petscPreconditioner);  CHKERRQ(ierr);
    }
    else if( numberOfProcessors>1 )
    {
      printF("PETScSolver:ERROR: cannot set the parallel preconditioner to be %i\n",
	     newPreconditioner);
    }
    
  }
  
  MatOrderingType  matOrdering;           //  == ORDER_RCM default;  
  PetscOptionsHasName(PETSC_NULL,"-pc_factor_mat_ordering_type",&optionWasSet);  // check if user has set
  if( !optionWasSet &&
      parameters.matrixOrdering!=matrixOrdering && 
      parameters.preconditioner==OgesParameters::incompleteLUPreconditioner)
  {
    parametersHaveChanged=true;
    switch( parameters.matrixOrdering )
    {
    case OgesParameters::naturalOrdering:
      matOrdering=(char*)MATORDERINGNATURAL;
      break;
    case OgesParameters::nestedDisectionOrdering:
      matOrdering=(char*)MATORDERINGND;
      break;
    case OgesParameters::oneWayDisectionOrdering:
      matOrdering=(char*)MATORDERING1WD;
      break;
    case OgesParameters::reverseCuthillMcKeeOrdering:
      matOrdering=(char*)MATORDERINGRCM;
      break;
    case OgesParameters::quotientMinimumDegreeOrdering:
      matOrdering=(char*)MATORDERINGQMD;
      break;
    case OgesParameters::rowlengthOrdering:
      matOrdering=(char*)MATORDERINGROWLENGTH;
      break;
    default:
      printf("****WARNING**** Unknown matrix ordering PETSc\n");
      matOrdering=(char*)MATORDERINGNATURAL;
    }
    if( numberOfProcessors==1 )
    {
      if( Oges::debug & 2 ) 
	printF(" ********** PETScSolver: set matrix ordering ***********\n");

      // 2.2.1 ierr = PCILUSetMatOrdering(pc, matOrdering); CHKERRQ(ierr);
#if (PETSC_VERSION_MAJOR==3)
      ierr = PCFactorSetMatOrderingType(pc, matOrdering); CHKERRQ(ierr);
#else
      ierr = PCFactorSetMatOrdering(pc, matOrdering); CHKERRQ(ierr);
#endif

    }
  }
  
  PetscOptionsHasName(PETSC_NULL,"-ksp_gmres_restart",&optionWasSet);  // check if user has set
  if( !optionWasSet &&
      parameters.parallelSolverMethod==OgesParameters::gmres &&
      parameters.gmresRestartLength!=gmresRestartLength )
  {
    ierr = KSPGMRESSetRestart(ksp, parameters.gmresRestartLength); CHKERRQ(ierr);
    gmresRestartLength=parameters.gmresRestartLength;
  }

  PetscOptionsHasName(PETSC_NULL,"-pc_ilu_levels",&optionWasSet);  // check if user has set
  if( !optionWasSet &&
      parameters.numberOfIncompleteLULevels!=numberOfIncompleteLULevels &&
      parameters.preconditioner==OgesParameters::incompleteLUPreconditioner )
  {
    parametersHaveChanged=true;
    if( numberOfProcessors==1 )
    {
      if( Oges::debug & 2 ) 
	printF(" ********** PETScSolver: set ilu levels ***********\n");

      // 2.2.1 ierr = PCILUSetLevels(pc, parameters.numberOfIncompleteLULevels);  CHKERRQ(ierr);
      // 2.2.1 ierr = PCILUSetFill(pc,   parameters.incompleteLUExpectedFill);    CHKERRQ(ierr);
      ierr = PCFactorSetLevels(pc, parameters.numberOfIncompleteLULevels);  CHKERRQ(ierr);
      ierr = PCFactorSetFill(pc,   parameters.incompleteLUExpectedFill);    CHKERRQ(ierr);

      numberOfIncompleteLULevels=parameters.numberOfIncompleteLULevels;
    }
  }

//  PCASM:    Set the overlap, using the default PETSc decomposition via
//      PCASMSetOverlap(pc,overlap);
//      Could instead use the option -pc_asm_overlap <ovl> 

  if( true && numberOfProcessors>=1 && ( true || parametersHaveChanged)  )
  {
    // ***this causes memory problems when viewing the ksp for some reason ***

    // Parallel preconditioners are:
    // petscPreconditioner==PCBJACOBI ||
    // petscPreconditioner==PCASM )

    // we can set the block preconditioner:

    // This code comes from ex7.c and ex8.c (ASM)
    PetscBool isbjacobi=PETSC_FALSE, isasm=PETSC_FALSE;
    // v2.3.2 ierr = PetscTypeCompare((PetscObject)pc,PCBJACOBI,&isbjacobi);CHKERRQ(ierr);
    ierr = PetscObjectTypeCompare((PetscObject)pc,PCBJACOBI,&isbjacobi);CHKERRQ(ierr);
    if( Oges::debug & 1 && isbjacobi )
      printF(" PETScSolver::setPetscParameters:INFO: Using Block Jacobi\n");
    if( !isbjacobi ) 
    {
      // ierr = PetscTypeCompare((PetscObject)pc,PCASM,&isasm);CHKERRQ(ierr);
      ierr = PetscObjectTypeCompare((PetscObject)pc,PCASM,&isasm);CHKERRQ(ierr);
      if( Oges::debug & 1 && isasm )
        printF(" PETScSolver::setPetscParameters:INFO: Using the Additive Schwarz Method\n"); 
    }
    
    if( isbjacobi || isasm ) 
    {
      ierr = KSPSetUp(ksp);CHKERRQ(ierr);
  
      KSP *subksp;     /* array of local KSP contexts on this processor */
      // Extract the array of KSP contexts for the local blocks
      PetscInt nlocal,first; 
     
      if( isbjacobi )
      {
	ierr = PCBJacobiGetSubKSP(pc,&nlocal,&first,&subksp);CHKERRQ(ierr);
      }
      else
      {
	ierr = PCASMGetSubKSP(pc,&nlocal,&first,&subksp);CHKERRQ(ierr);
      }
      if( Oges::debug & 1 )
	printF(" PETScSolver::setPetscParameters:INFO: number of local blocks=%i\n",nlocal); 

      // Loop over the local blocks, setting various KSP options for each block.  
      PC subpc;        /* PC context for subdomain */
      for (i=0; i<nlocal; i++) 
      {
	ierr = KSPGetPC(subksp[i],&subpc);CHKERRQ(ierr);
        PetscOptionsHasName(PETSC_NULL,"-sub_pc_type",&optionWasSet);  // check if user has set
        if( !optionWasSet && parameters.preconditioner!=preconditioner )
	{
	  if( Oges::debug & 1 )
	    printF(" PETScSolver::setPetscParameters:INFO: setting block preconditioner for block %i (p=%i) to %i\n",
		   i,myid,parameters.preconditioner);

          PCType petscPreconditioner = getPetscPreconditioner(parameters.preconditioner);
	  
	  ierr = PCSetType(subpc,petscPreconditioner);CHKERRQ(ierr);
	}

        PetscOptionsHasName(PETSC_NULL,"-sub_ksp_type",&optionWasSet);  // check if user has set
	if( !optionWasSet && parameters.solverMethod!=solverMethod )
	{
	  if( Oges::debug & 1 )
	    printF(" PETScSolver::setPetscParameters:INFO: setting KSP type for block %i (p=%i) to %i\n",
		   i,myid,parameters.solverMethod);
          
	  KSPType krylovSpaceMethod = getPetscKrylovSpaceMethod(parameters.solverMethod);
	  ierr = KSPSetType(subksp[i],krylovSpaceMethod);CHKERRQ(ierr);
	}

	PetscOptionsHasName(PETSC_NULL,"-sub_pc_factor_mat_ordering_type",&optionWasSet);  // check if user has set
        if( !optionWasSet &&
	    parameters.matrixOrdering!=matrixOrdering && 
	    parameters.preconditioner==OgesParameters::incompleteLUPreconditioner)
	{
	  // 2.2.1 ierr = PCILUSetMatOrdering(subpc, matOrdering); CHKERRQ(ierr);
#if (PETSC_VERSION_MAJOR==3)
	  ierr = PCFactorSetMatOrderingType(subpc, matOrdering); CHKERRQ(ierr);
#else
	  ierr = PCFactorSetMatOrdering(subpc, matOrdering); CHKERRQ(ierr);
#endif
	}
  
        PetscOptionsHasName(PETSC_NULL,"-sub_pc_factor_levels",&optionWasSet);  // check if user has set
	if( !optionWasSet &&
            parameters.numberOfIncompleteLULevels!=numberOfIncompleteLULevels &&
	    parameters.preconditioner==OgesParameters::incompleteLUPreconditioner )
	{
	  ierr = PCFactorSetLevels(subpc, parameters.numberOfIncompleteLULevels);  CHKERRQ(ierr);
	}
        PetscOptionsHasName(PETSC_NULL,"-sub_pc_factor_fill",&optionWasSet);  // check if user has set
	if( !optionWasSet &&
            parameters.numberOfIncompleteLULevels!=numberOfIncompleteLULevels &&
	    parameters.preconditioner==OgesParameters::incompleteLUPreconditioner )
	{
	  ierr = PCFactorSetFill(subpc,   parameters.incompleteLUExpectedFill);    CHKERRQ(ierr);
	}
	// ierr = KSPSetTolerances(subksp[i],1.e-6,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

// 	if( parameters.solverMethod!=OgesParameters::preonly ) // petscPreconditioner!=PCLU )
// 	{
// 	  ierr = KSPSetInitialGuessNonzero(subksp[i],PETSC_TRUE); CHKERRQ(ierr); 
// 	}
// 	else
// 	{
// 	  ierr = KSPSetInitialGuessNonzero(subksp[i],PETSC_FALSE); CHKERRQ(ierr); 
// 	}

      }
    }
  }

  solverMethod=parameters.solverMethod;
  matrixOrdering=parameters.matrixOrdering;
  preconditioner=parameters.preconditioner;
  numberOfIncompleteLULevels=parameters.numberOfIncompleteLULevels;


}
#undef __FUNC__
#define __FUNC__ "PETScSolver::setPetscParameters"
int PETScSolver::
setPetscRunTimeParameters() 
{
  if( !processorIsActive )
    return 0;

  const int myid=Communication_Manager::My_Process_Number;

  bool parametersHaveChanged=false;  // set to true if any parameters have changed
  
//   if( true ) return 0; // ************************************************************************ 051112 * temp
  
  // rtol : reduce residual by this factor.
  // atol : absolute tolerance
  // dtol : divergence detector
  //    convergence:  | r_k |_2 < max( rtol*| r_0 |_2, atol )

  double rtol=parameters.relativeTolerance>0. ? parameters.relativeTolerance : REAL_EPSILON*1000.;
  double atol=parameters.absoluteTolerance>0. ? parameters.absoluteTolerance : 
              max( real(numberOfEquations),500.)*REAL_EPSILON;
  double dtol=parameters.maximumAllowableIncreaseInResidual;
  int maxits = parameters.maximumNumberOfIterations > 0 ?  parameters.maximumNumberOfIterations : 900;
  
//   if( parameters.solveForTranspose )
//     dtol=DBL_MAX; // ** assume we are solving for the left null vector

  if( Oges::debug & 4 ) 
     printF(" PETScSolver: rtol=%e, atol=%e, dtol=%e\n",rtol,atol,dtol);
  
  ierr = KSPSetTolerances(ksp, rtol, atol, dtol, maxits); CHKERRQ(ierr);
  
  // DO NOT assign the sub KSP's here -- let the user set through "define petscOption -sub_ksp_rtol 1.e-5" etc.
  if( false && numberOfProcessors>1 )
  {

    ierr = KSPGetPC(ksp, &pc);    CHKERRQ(ierr);

    PetscBool isbjacobi;
    // ierr = PetscTypeCompare((PetscObject)pc,PCBJACOBI,&isbjacobi);CHKERRQ(ierr);
    ierr = PetscObjectTypeCompare((PetscObject)pc,PCBJACOBI,&isbjacobi);CHKERRQ(ierr);

    if( Oges::debug & 4 )
      printF(" PETScSolver::setPetscRunTimeParameters:INFO: isbjacobi=%i\n",isbjacobi); 

    if (isbjacobi) 
    {
      ierr = KSPSetUp(ksp);CHKERRQ(ierr);
  
      KSP *subksp;     /* array of local KSP contexts on this processor */
      // Extract the array of KSP contexts for the local blocks
      PetscInt nlocal,first; 
      ierr = PCBJacobiGetSubKSP(pc,&nlocal,&first,&subksp);CHKERRQ(ierr);
      for (i=0; i<nlocal; i++) 
      {
        ierr = KSPSetTolerances(subksp[i], rtol, atol, dtol, maxits); CHKERRQ(ierr);
      }
      

    }
  }
  


}
static bool firstSolve=true;

#undef __FUNC__
#define __FUNC__ "PETScSolver::sizeOf"
real PETScSolver::
sizeOf( FILE *file /* =NULL */  )
// return number of bytes allocated 
{

  FILE *outputFile = file==NULL ? stdout : file;
  
  real size=0.;  
  PetscLogDouble mem=0;
  if( initialized )
  {

    // 2.2.1  PetscLogDouble space=0, fragments=0, maximumBytes=0, mem=0;
    if( turnOnPETScMemoryTracing )
    {
      // 2.2.1 ierr = PetscTrSpace( &space, &fragments, &maximumBytes);  CHKERRQ(ierr);
      // 2.2.1 PetscGetResidentSetSize(&mem); //  maximum memory used
      PetscMemoryGetMaximumUsage(&mem);       //  maximum memory used
    
      size=mem;
    }
    else
    {
      if( processorIsActive )
      {
	MatInfo matInfo;
	ierr=MatGetInfo(A,MAT_GLOBAL_SUM,&matInfo); CHKERRQ(ierr);
	size=matInfo.memory;
      }
      
    }
  }
  
  if( Oges::debug & 2 )
  {
    // 2.2.1  fprintf(outputFile,
    // 2.2.1          ">> PETSC: %e Kbytes in use.  (maximum used=%e Kbytes, total memory use=%e Kbytes)\n"
    // 2.2.1  	    "  matrix=%e Kbytes, fragments=%e\n",
    // 2.2.1  	    space*.001,maximumBytes*.001,mem*.001,matInfo.memory*.001,fragments);
    fprintf(outputFile,">> PETSC: maximum memory used=%e (Kbytes)\n",mem/1000);

  }

  if( Oges::debug & 4 )
  {
    // 2.2.1 PetscTrDump(stdout);
    PetscMallocDump(stdout);  // Dumps the allocated memory blocks to a file
  }

  return size;
}

int PETScSolver::
solve( realCompositeGridFunction & uu, realCompositeGridFunction & f )
// ===================================================================================================
//  /Description:
//      Solve the equations.
// 
// ===================================================================================================
{
  
  CompositeGrid & cg = *uu.getCompositeGrid();
  const int myid=Communication_Manager::My_Process_Number;



  // new way 090707 -- we need to fix this so tcm3 works when solving first a dirichlet then a neumann problem
  bool shouldUpdateMatrix=oges.refactor || !oges.initialized || oges.shouldBeInitialized; // *wdh* added 090706 
  if( initialized && shouldUpdateMatrix )
  {
    // -- rebuild the matrix after it has already been built ---
    // do this for now: 
    if( debug & 2 )
      printF("PETScSolver::solve: rebuilding an existing matrix...\n");

    // We probably don't need to destroy everything -- could keep vectors ??
    destroy(); // *wdh* 091128 

    // ierr = MatZeroEntries(A);CHKERRQ(ierr);
      
    // *wdh* 091128 -- these must be reset ---
    solverMethod=-1;
    preconditioner=-1;
    matrixOrdering=-1;
    numberOfIncompleteLULevels=-1;
    gmresRestartLength=-1;


  }

  if( !initialized || reInitialize )
  {
    initializePETSc();

    initialized=true;  // This should remain true from now on (for PETSc instances)
    reInitialize=false;
    shouldUpdateMatrix=true;
  }

  // printf(" --PETSCsSolver: solve: myid=%i, processorIsActive=%i \n",myid,(int)processorIsActive);
  

  if( !processorIsActive )
  { // --- return here if this processor is not involved ---
    return 0;
  }

  MPI_Comm & OGES_COMM = parameters.getCommunicator();

  if( shouldUpdateMatrix )  
  {

    buildMatrix(oges.coeff,uu);

    if( debug & 2 )
      printF("PETScSolver::solve: done buildMatrix.\n");

    buildSolver(); 

    if( debug & 2 )
      printF("PETScSolver::solve: done buildSolver.\n");

    oges.refactor=false;
    oges.initialized=true;
    oges.shouldBeInitialized=false;
  }
  

  // set any run time parameters that have changed (e.g. tolerances)
  setPetscRunTimeParameters();
  
  // The initial guess is non-zero unless we use a direct solver
//   if( (numberOfProcessors<=1 && parameters.solverMethod!=OgesParameters::preonly) ||
//       (numberOfProcessors>1  &&  parameters.parallelSolverMethod!=OgesParameters::preonly ) )
  int procCt;
  if (oges.parameters.OGES_COMM == MPI_COMM_NULL) {
    procCt = 0;
  }
  else {
    MPI_Comm_size(oges.parameters.OGES_COMM, &procCt);
  }
  if( (procCt > 1 && parameters.parallelSolverMethod != OgesParameters::preonly) ||
      (procCt == 1 && parameters.solverMethod != OgesParameters::preonly) )
  {
    ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE); CHKERRQ(ierr); 
  }
  else
  {
    ierr = KSPSetInitialGuessNonzero(ksp,PETSC_FALSE); CHKERRQ(ierr); 
  }


  int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2]; 
  int grid;

  int gride=cg.numberOfComponentGrids()-1;  // for the extra eqution


  // optimized version for filling in the RHS
  PetscScalar *bv, *xv;
  VecGetArray(x,&xv);  // get the local array from Petsc
  VecGetArray(b,&bv);  // get the local array from Petsc

  if( parameters.rescaleRowNorms )
  {
    assert( diagonalScale!=NULL ); 
  }
  ierr = VecGetOwnershipRange(b,&Istart,&Iend);CHKERRQ(ierr);
  
  assert( pCoeff!=NULL );
  realCompositeGridFunction & coeff = *pCoeff;

  for( grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {  
    MappedGrid & mg = cg[grid];
    // mg.update( MappedGrid::THEcenter ); // *wdh* 2011/08/20 - this is not needed
    intSerialArray mask; getLocalArrayWithGhostBoundaries(cg[grid].mask(),mask);
    intSerialArray classify; getLocalArrayWithGhostBoundaries(coeff[grid].sparse->classify,classify);

    realArray & ug= uu[grid];
    realSerialArray uLocal; getLocalArrayWithGhostBoundaries(ug,uLocal); 
    int n1a = uLocal.getBase(0) +ug.getGhostBoundaryWidth(0); 
    int n1b = uLocal.getBound(0)-ug.getGhostBoundaryWidth(0);
    int n2a = uLocal.getBase(1) +ug.getGhostBoundaryWidth(1); 
    int n2b = uLocal.getBound(1)-ug.getGhostBoundaryWidth(1);
    int n3a = uLocal.getBase(2) +ug.getGhostBoundaryWidth(2); 
    int n3b = uLocal.getBound(2)-ug.getGhostBoundaryWidth(2);
      
    // realSerialArray xLocal; getLocalArrayWithGhostBoundaries(mg.center(),xLocal); 
    realSerialArray fLocal; getLocalArrayWithGhostBoundaries(f[grid],fLocal); 
    realSerialArray ds;
    if( parameters.rescaleRowNorms )
    {
      getLocalArrayWithGhostBoundaries((*diagonalScale)[grid],ds);
    }

    const IntegerArray & gid = cg[grid].gridIndexRange(); 

    i1=n1a, i2=n2a, i3=n3a;
    int n=0;
    int ig=getGlobalIndex( n,iv, grid, myid );  // get the global index for the first point

    const int nb=fLocal.getBase(3);
    FOR3N(i1,i2,i3,n,n1a,n1b,n2a,n2b,n3a,n3b)
    {
      // int ig=getGlobalIndex( iv, grid, myid );  // get the global index 
      // assert( ig>= Istart && ig<=Iend );
      if( ig<Istart || ig>Iend )
      {
        int igg=getGlobalIndex( n,iv, grid, myid );
	printf("PETScSolver:solve:ERROR filling in the rhs and initial conditions\n"
	       " ig=%i is not in the range [Istart,Iend]=[%i,%i]\n"
               " myid=%i, grid=%i, getGlobalIndex( n,iv, grid, myid )=%i (should equal ig)\n" 
               " (i1,i2,i3,n)=(%i,%i,%i,%i) [n1a,n1b][n2a,n2b][n3a,n3b]=[%i,%i][%i,%i][%i,%i]\n",
	       ig,Istart,Iend, myid,grid, igg, i1,i2,i3,n, n1a,n1b,n2a,n2b,n3a,n3b);

// 	i1=n1a, i2=n2a, i3=n3a;
// 	int ig0=getGlobalIndex( n,iv, grid, myid ); 
//         printf(" intial value: (i1,i2,i3)=(%i,%i,%i) getGlobalIndex( n,iv, grid, myid )=%i\n",
// 	       i1,i2,i3,ig0);

// 	int p=myid;
//         printf("getGlobalIndex: nab(0,axis1,p,grid)=%i, ndab(0,p,grid)=%i\n"
//                "                nab(0,axis2,p,grid)=%i, ndab(1,p,grid)=%i\n" 
//                "                nab(0,axis3,p,grid)=%i \n" 
//                "  noffset(p,grid)=%i\n",
// 	       nab(0,axis1,p,grid),ndab(0,p,grid),nab(0,axis2,p,grid),ndab(1,p,grid),nab(0,axis3,p,grid),
//                noffset(p,grid));

	Overture::abort("fatal error");
      }
	  
      if( classify(i1,i2,i3,n)<=0 ) // mask(i1,i2,i3)<0 )
      {
        // -1=interp, -2=periodic, -3=extrap, 0=unused
	bv[ig-Istart]=0.; // VecSetValues( b, 1, &ig, &v, INSERT_VALUES );
        xv[ig-Istart]=uLocal(i1,i2,i3,n+nb);  
      }
      else
      {
        if( parameters.rescaleRowNorms )
	  bv[ig-Istart]=fLocal(i1,i2,i3,n+nb)*ds(i1,i2,i3,n); // VecSetValues( b, 1, &ig, &v, INSERT_VALUES );
	else
	  bv[ig-Istart]=fLocal(i1,i2,i3,n+nb); // VecSetValues( b, 1, &ig, &v, INSERT_VALUES );
	
        xv[ig-Istart]=uLocal(i1,i2,i3,n+nb);
      }

      ig++;
    } // end FOR3N

    if( grid==gride && problemIsSingular==addExtraEquation )
    {
      // fill in the RHS for the extra equation
      int grid=cg.numberOfComponentGrids()-1;
      const IntegerArray & d = cg[grid].dimension();
      i1=d(1,0), i2=d(1,1), i3=d(1,2); 
      int p= uu[grid].Array_Descriptor.findProcNum( iv );  // processor number
      if( myid==p )
      {
	int n=numberOfComponents-1;
	int extraEquation=getGlobalIndex( n, iv, grid, p );  // get the global index
	assert( extraEquation>=Istart && extraEquation<=Iend );
        // NOTE: no diagonal scaling on this equation
        bv[extraEquation-Istart]=fLocal(i1,i2,i3,n+nb);
        // printf("PETScSolver: Assign RHS for constraint : value=%g\n",bv[extraEquation-Istart]);
      }
    }
  }


  VecRestoreArray(b,&bv);
  VecRestoreArray(x,&xv);

  // from ex29.c
  if( false && problemIsSingular!=notSingular ) // Make the RHS compatible with the null space
  {
    ierr = KSPGetNullSpace(ksp,&nsp);CHKERRQ(ierr);
    ierr = MatNullSpaceRemove(nsp,b,PETSC_NULL);CHKERRQ(ierr);
  }

  if( debug & 8 )
  {
    printF("PETScSolver::b (RHS before solve)\n");
    ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }

  // How do we solve for the transpose?
  // PetscErrorCode  MatSolveTranspose(Mat mat,Vec b,Vec x)
  // PetscErrorCode  MatTranspose(Mat mat,MatReuse reuse,Mat *B)

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  real time=getCPU();

  if( debug & 2 )
    printF("PETScSolver::solve: before KSPSolve...\n");

  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

  if( debug & 2 )
    printF("PETScSolver::solve: after KSPSolve...\n");

  time=getCPU()-time;
  time=ParallelUtility::getMaxValue(time,-1,OGES_COMM);

  KSPConvergedReason reason;
  ierr = KSPGetConvergedReason(ksp,&reason);
  if( reason<0 && !( reason==-3 && NUM_KSP_CONVERGED_ITS_ERRORS>10 ) )
  {
    // typedef enum {/* converged */
    //               KSP_CONVERGED_RTOL               =  2,
    //               KSP_CONVERGED_ATOL               =  3,
    //               KSP_CONVERGED_ITS                =  4,
    //               KSP_CONVERGED_STCG_NEG_CURVE     =  5,
    //               KSP_CONVERGED_STCG_CONSTRAINED   =  6,
    //               KSP_CONVERGED_STEP_LENGTH        =  7,
    //               /* diverged */
    //               KSP_DIVERGED_NULL                = -2,
    //               KSP_DIVERGED_ITS                 = -3,
    //               KSP_DIVERGED_DTOL                = -4,
    //               KSP_DIVERGED_BREAKDOWN           = -5,
    //               KSP_DIVERGED_BREAKDOWN_BICG      = -6,
    //               KSP_DIVERGED_NONSYMMETRIC        = -7,
    //               KSP_DIVERGED_INDEFINITE_PC       = -8,
    //               KSP_DIVERGED_NAN                 = -9,
    //               KSP_DIVERGED_INDEFINITE_MAT      = -10,
    //  
    //               KSP_CONVERGED_ITERATING          =  0} KSPConvergedReason;

    printF("PETScEquationSolver:ERROR: Solution diverged! reason=%i : \n",(int)reason);
    printF("     KSP_DIVERGED_NULL                = -2,\n"
           "     KSP_DIVERGED_ITS                 = -3,\n"
           "     KSP_DIVERGED_DTOL                = -4,\n"
           "     KSP_DIVERGED_BREAKDOWN           = -5,\n"
           "     KSP_DIVERGED_BREAKDOWN_BICG      = -6,\n"
           "     KSP_DIVERGED_NONSYMMETRIC        = -7,\n"
           "     KSP_DIVERGED_INDEFINITE_PC       = -8,\n"
           "     KSP_DIVERGED_NAN                  = -9,\n"
           "     KSP_DIVERGED_INDEFINITE_MAT      = -10\n");
    printF("NOTE: to see more information turn on the '-info' PETSc option (e.g. in your .petscrc)\n");
    printF("NOTE 2: to avoid the divergence error '-4' you can set the Oges option 'maximum allowable increase in the residual' \n");

    if( reason==-3 )
      NUM_KSP_CONVERGED_ITS_ERRORS++;
    if( NUM_KSP_CONVERGED_ITS_ERRORS>=10 )
    {
      printF("\n ******** PETScEquationSolver: TOO MANY KSP_DIVERGED_ITS errors. I will not print anymore. *********\n");
    }
    

    // OV_ABORT("error");

  }

  PetscInt its;
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  oges.numberOfIterations=its;   // tell Oges how many iterations we used.
  
  if( Oges::debug & 1 || debug & 2 )
  {
    printF("PETScSolver::Time for solve=%8.2e (its=%i) \n",time,its);
  }
  
  if( debug & 8 )
  {
    printF("PETScSolver::x (solution after solve)\n");
    ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }




    //  === copy the solution from Petsc to Overture ====

  VecGetArray(x,&xv);  // get the local array from Petsc

  ierr = VecGetOwnershipRange(x,&Istart,&Iend);CHKERRQ(ierr);

  for( grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {  
    realArray & ug= uu[grid];
    realSerialArray uLocal; getLocalArrayWithGhostBoundaries(ug,uLocal); 

    int n1a = uLocal.getBase(0) +ug.getGhostBoundaryWidth(0), 
        n1b = uLocal.getBound(0)-ug.getGhostBoundaryWidth(0);
    int n2a = uLocal.getBase(1) +ug.getGhostBoundaryWidth(1), 
        n2b = uLocal.getBound(1)-ug.getGhostBoundaryWidth(1);
    int n3a = uLocal.getBase(2) +ug.getGhostBoundaryWidth(2), 
        n3b = uLocal.getBound(2)-ug.getGhostBoundaryWidth(2);
      
    if( false )
      printf("PETScSolver:solver myid=%i local array bounds = [%i,%i][%i,%i][%i,%i]\n",myid,n1a,n1b,n2a,n2b,n3a,n3b);

    i1=n1a, i2=n2a, i3=n3a;
    int n=0;
    int ig=getGlobalIndex( n, iv, grid, myid );  // get the global index for the first point

    const int nb=uLocal.getBase(3);
    FOR3N(i1,i2,i3,n,n1a,n1b,n2a,n2b,n3a,n3b)
    {

      // ******** NOTE we can probably just increment ig by 1 if we start correctly
      // int ig=getGlobalIndex( iv, grid, myid );  // get the global index

      if( ig>=Istart && ig<=Iend )
      {
	if( false ) printf(" myid=%i: i1,i2=%i,%i, ig=%i xv[ig]=%6.4f\n",myid,i1,i2,ig,xv[ig-Istart]);
	uLocal(i1,i2,i3,n+nb)=xv[ig-Istart];
      }
      else
      {
        int p=myid;
        printf("ERROR: myid=%i, i1,i2=%i,%i, ig=%i Istart,Iend=[%i,%i]\n", myid,i1,i2,ig,Istart,Iend);
      }
      ig++;
    }
    
  }
  for( grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {  
    uu[grid].updateGhostBoundaries();
    if( debug & 8 )
      display(uu[grid],sPrintF("Solution: uu[%i]",grid),"%6.3f ");
  }
  
  VecRestoreArray(x,&xv);

  if( debug & 4 )
  {
    // compute the maximum residual
    // printF(" ***PETScSolver:compute max residual...\n");

    Vec res;
//     numberOfVects++;
//     printF("VecDup res, vect object %i.\n",numberOfVects);
    ierr = VecDuplicate(b,&res);CHKERRQ(ierr);

    ierr = MatMult(A,x,res);CHKERRQ(ierr);   // res = A*x

//     if( myid==0 ) printf("PETScSolver::res=A*x\n");
//     ierr = VecView(res,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

//     if( myid==0 ) printf("PETScSolver::b\n");
//     ierr = VecView(res,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    // Computes y = x + alpha y.
    PetscScalar alpha; alpha=-1.;
    // 2.2.1 VecAYPX(&alpha,b,res);  // res = b - res
    VecAYPX(res,alpha,b);  // res = b - res

//     if( myid==0 ) printf("PETScSolver::res=b-A*x\n");
//     ierr = VecView(res,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    
//     if( myid==0 ) printf("PETScSolver::b\n");
//     ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    
 
    PetscReal maxNorm;
    VecNorm(res, NORM_INFINITY, &maxNorm);

    ierr = VecDestroy(&res);CHKERRQ(ierr);    

    printF(" ***PETScSolver:solve: max residual = %8.2e (relativeTol=%8.2e)\n",maxNorm,
              parameters.relativeTolerance);

  }
  


  if( firstSolve && Oges::debug & 1  )
  {
    aString name = getSolverName();
    printF("--PETSc-- solver: %s\n",(const char*)name);
  }


  if( firstSolve && debug & 4 )
  { // trouble when this is turned on ??
    firstSolve=false;
    printF("PETScSolver::ksp\n");
    ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  return 0;
}


#define initExplicitInterp EXTERN_C_NAME(initexplicitinterp)
extern "C"
{
  void initExplicitInterp(const int&ndc1,const int&ndc2,const int&ndc3,const int&ndci,
          const int&ipar,real&coeff,const real&ci,real&pr,real&ps,real&pt,
          const real&gridSpacing,const int&indexStart,
	  const int&variableInterpolationWidth,const int&interpoleeLocation,const int&interpoleeGrid);
}

int PETScSolver::
fillInterpolationCoefficients(realCompositeGridFunction & uu)
//===================================================================
// /Description:
//   Fill the matrix with the interpolation coefficients
// (this routine started from Interpolant::initializeExplicitInterpolation
//===================================================================
{
  CompositeGrid & cg = *uu.getCompositeGrid();
  const int myid=Communication_Manager::My_Process_Number;

  if( Oges::debug & 4 )
  {
    printF("PETScSolver::fillInterpolationCoefficients: max(cg.numberOfInterpolationPoints)=%i\n",max(cg.numberOfInterpolationPoints));
  }
  

  if( cg.numberOfBaseGrids() ==0 || max(cg.numberOfInterpolationPoints) <= 0 )
    return 0;

  real time0=getCPU();

  const int numberOfDimensions=cg.numberOfDimensions();

  // for now we use only one width per grid
  int axis,grid;
  IntegerArray width(3,cg.numberOfComponentGrids()); width=1;
  Range NG(0,cg.numberOfComponentGrids()-1);
  for( grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  for( axis=axis1; axis<numberOfDimensions; axis++ ) 
    width(axis,grid)=max(width(axis,grid),max(cg.interpolationWidth(axis,grid,NG)));

  const int maxWidth=max(width);

  IntegerArray indexStart(3,cg.numberOfComponentGrids());
  RealArray gridSpacing(3,cg.numberOfComponentGrids());
  for( grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {
    MappedGrid & mg = cg[grid];
    for( axis=0; axis<3; axis++ )
    {
      indexStart(axis,grid)=mg.gridIndexRange(0,axis);
      gridSpacing(axis,grid)=mg.gridSpacing(axis);
    }
      
  }
    
  int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2]; 
  int ivd[3], &i1d=ivd[0], &i2d=ivd[1], &i3d=ivd[2];  // for donor position

  realSerialArray coeff;
  
  for( grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {
    MappedGrid & mg = cg[grid];
    const realArray & ug = uu[grid];

    int ni=cg.numberOfInterpolationPoints(grid);
    if( ni==0 ) continue;
	
    intArray & interpoleeLocation = cg.interpoleeLocation[grid];
    intSerialArray il; 
    
    // *wdh* 091129 -- fix to use local interpolation arrays if they are there.
    if( ( grid<cg.numberOfBaseGrids() && 
	  cg->localInterpolationDataState==CompositeGridData::localInterpolationDataForAMR ) || 
	cg->localInterpolationDataState==CompositeGridData::noLocalInterpolationData )
    {

      // use the interpolation data in the parallel arrays
      getLocalArrayWithGhostBoundaries(cg.interpoleeLocation[grid],il);
    }
    else
    {
      // use the interpolation data in the serial arrays (for now these are refinement grids)
      il.reference(cg->interpoleeLocationLocal[grid]);
    }

    int n1a = il.getBase(0) +interpoleeLocation.getGhostBoundaryWidth(0), 
        n1b = il.getBound(0)-interpoleeLocation.getGhostBoundaryWidth(0);

    ni = n1b-n1a+1;  // number of interpolation points on this processor
    
    if( ni==0 ) continue;
    

    Range R(n1a,n1b);

    coeff.redim(R,width(axis1,grid),width(axis2,grid),width(axis3,grid));

    intSerialArray ip; 
    intSerialArray interpoleeGrid; 
    intSerialArray viw; 
    realSerialArray ci; 

    if( ( grid<cg.numberOfBaseGrids() && 
	  cg->localInterpolationDataState==CompositeGridData::localInterpolationDataForAMR ) || 
	cg->localInterpolationDataState==CompositeGridData::noLocalInterpolationData )
    {

      // use the interpolation data in the parallel arrays
      getLocalArrayWithGhostBoundaries(cg.interpolationPoint[grid],ip);
      getLocalArrayWithGhostBoundaries(cg.interpoleeGrid[grid],interpoleeGrid);
      getLocalArrayWithGhostBoundaries(cg.variableInterpolationWidth[grid],viw);
      getLocalArrayWithGhostBoundaries(cg.interpolationCoordinates[grid],ci);

    }
    else
    {
      // use the interpolation data in the serial arrays (for now these are refinement grids)
      // printf("PETScSolver::USE LOCAL INTERP ARRAYS\n");
      
      ip.reference( cg->interpolationPointLocal[grid]);
      il.reference( cg->interpoleeLocationLocal[grid]);
      interpoleeGrid.reference( cg->interpoleeGridLocal[grid]);
      viw.reference( cg->variableInterpolationWidthLocal[grid]);
      ci.reference(cg->interpolationCoordinatesLocal[grid]);
    }
      
    int debug = 0; // Oges::debug;
    if( debug & 4 )
    {
      printf("PETScSolver::fillInterp: myid=%i: grid=%i [n1a,n1b]=[%i,%i] numberOfComponents=%i interp-width=%i\n",myid,grid,n1a,n1b,
	     numberOfComponents,width(axis1,grid));
      // for( int i=n1a; i<=n1b; i++ )
      // {
      // 	printf(" grid=%i i=%i ip=(%i,%i) il=(%i,%i) viw=%i\n",grid,i,ip(i,0),ip(i,1),il(i,0),il(i,1),viw(i));
      // }
    }

    int ipar[7]={numberOfDimensions,
		 grid,
		 ni,
		 mg.isCellCentered(0),
		 (int)Interpolant::precomputeAllCoefficients,
		 maxWidth,
		 0}; //  This last position is saved for a return value of useVariableWidthInterpolation

    realSerialArray & cc = coeff;

    // ************ warning -- watch out for indexing of local arrays -- base and bound : n1a,..,n1b	
    //   use interpoleeGrid(n1a,0)
    RealArray pr(R),ps(R),pt(R);
    initExplicitInterp(cc.getLength(0),cc.getLength(1),cc.getLength(2),il.getLength(0),
		       ipar[0], 
		       *cc.getDataPointer(),
		       ci(n1a,0),pr(n1a),ps(n1a),pt(n1a),gridSpacing(0,0),indexStart(0,0),
		       viw(n1a),il(n1a,0),interpoleeGrid(n1a,0));
	
    int useVariableWidthInterpolation=ipar[6]; 

    // printf("PETScSolver::useVariableWidthInterpolation=%i\n",useVariableWidthInterpolation);

    // cc.display("coeff after initExplicitInterp");

    const real epsForInterpCoeff=REAL_EPSILON*100.; // neglect interpolation coeff smaller than this

    assert( width(0,grid)==width(1,grid) );

    i3=0; i3d=0;
    int n=0;
    for( int i=n1a; i<=n1b; i++ )
    {
      i1=ip(i,0);
      i2=ip(i,1);
      if( numberOfDimensions==3 ) i3=ip(i,2);
      
      const int gridi = interpoleeGrid(i);
      const int iw = useVariableWidthInterpolation ? viw(i) : width(0,grid);  // *wdh* 100113 -- added support for VIW
      const int iw3 = numberOfDimensions==2 ? 1 : iw;

      int p= ug.Array_Descriptor.findProcNum( iv );  // processor number
      for( int n=0; n<numberOfComponents; n++ )
      {

        // fill in value -1 for grid,(i1,i2,i3)
	const int ig=getGlobalIndex( n, iv, grid, p );  // get the global index
	v=-1.; 
	ierr = MatSetValues(A,1,&ig,1,&ig,&v,INSERT_VALUES);CHKERRQ(ierr);

	// if( debug & 1 ) printf("interp: grid=%i p=%i i=%i ip=(%i,%i,%i) n=%i coeff[m1,m2]:\n",grid,p,i,i1,i2,i3,n);


	for( int m3=0; m3<iw3; m3++ )
	{
	  if( numberOfDimensions==3 ) i3d=il(i,2)+m3;
	  for( int m2=0; m2<iw; m2++ )
	  {
	    i2d=il(i,1)+m2;
	    for( int m1=0; m1<iw; m1++ )
	    {
	      i1d=il(i,0)+m1;
	
	      int p= uu[gridi].Array_Descriptor.findProcNum( ivd );  // processor number
	      int jg=getGlobalIndex( n, ivd, gridi, p );  // get the global index

	      // fill in value  grid, (i1,i2,i3) coeff(i,m1,m2,m3)
	      v=coeff(i,m1,m2,m3);

	      if( fabs(v)>epsForInterpCoeff )
	      {
		// if( debug & 1 ) printf(" (ig,jg)=(%i,%i) [m1=%i,m2=%i,m3=%i](gridi=%i,p=%i)=%4.2f \n",ig,jg,m1,m2,m3,gridi,p,v);
		ierr = MatSetValues(A,1,&ig,1,&jg,&v,INSERT_VALUES);CHKERRQ(ierr);

	      }
	    
	    }
	  }
	}
      }
    }

  }
  real time=getCPU()-time0;
  time0=getCPU();
  if( Mapping::debug & 1 ) printF("*** time for new init explicit = %8.2e\n",time);

    
  return 0;

}
  
int PETScSolver::
setExtraEquationValues( realCompositeGridFunction & f, real *value )
//==================================================================================
// /Description:
//   Assign values to the right-hand-side for the extra equations
//
// /f (input/output) : fill in rhs values here
// /values[i] (input) : values for each extra equation, i=0,1,2,...,
//
// /Return values: 0=success
//==================================================================================
{
  const int myid=max(0,Communication_Manager::My_Process_Number);

//  assert( problemIsSingular==addExtraEquation );
//  assert( oges.numberOfExtraEquations==1 );  // we only do this case so far

  CompositeGrid & cg = *f.getCompositeGrid();
  // find a spot to put the extra equation
  int grid=cg.numberOfComponentGrids()-1;
  const IntegerArray & d = cg[grid].dimension();

  int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2]; 
  i1=d(1,0), i2=d(1,1), i3=d(1,2); 
  int n=numberOfComponents-1;

  int p= f[grid].Array_Descriptor.findProcNum( iv );  // processor number

  const realSerialArray & fLocal = f[grid].getLocalArray();
  if( p==myid )
  {
    fLocal(i1,i2,i3,n)=value[0];
  }
  if( Oges::debug & 4 )
    printF("PETScSolver::setExtraEquationValues: f[%i](%i,%i,%i,n=%i)= %14.10e \n",grid,i1,i2,i3,n,value[0]);

  return 0;
}

int PETScSolver::
getExtraEquationValues( const realCompositeGridFunction & u, real *value, const int maxNumberToReturn /* =1 */ )
//==================================================================================
/// \brief Return solution values from the extra equations
///
/// \param u (input) : grid function holding the solution.
/// \param value[i] (output) : values for each extra equation, i=0,1,2,...,
/// \param maxNumberToReturn (input) : max number of values to return.
//
//==================================================================================
{
  const int myid=max(0,Communication_Manager::My_Process_Number);

//  assert( problemIsSingular==addExtraEquation );

  const int numExtra=min(maxNumberToReturn,oges.numberOfExtraEquations);
  
  assert( numExtra==1 );  // we only do this case so far

  const CompositeGrid & cg = *u.getCompositeGrid();
  // find a spot to put the extra equation
  int grid=cg.numberOfComponentGrids()-1;
  const IntegerArray & d = cg[grid].dimension();

  int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2]; 
  i1=d(1,0), i2=d(1,1), i3=d(1,2); 
  int n=numberOfComponents-1;

  int p= u[grid].Array_Descriptor.findProcNum( iv );  // processor number

  const realSerialArray & uLocal = u[grid].getLocalArray();
  if( p==myid )
  {
    value[0]=uLocal(i1,i2,i3,n);
  }

  broadCast(value[0],p); // broadcast value from processor p

  if( Oges::debug & 4 )
    printF("PETScSolver::getExtraEquationValues: value= u[%i](%i,%i,%i,n=%i)= %14.10e \n",
           grid,i1,i2,i3,n,value[0]);

  
  return 0;
}

int PETScSolver::
evaluateExtraEquation( const realCompositeGridFunction & u, real & value, int extraEquation /* =0 */ )
//==================================================================================
// /Description:
//    Evaluate the dot product of the coefficients of an extra equation times u 
//
// /u (input) : grid function to dot with the extra equation
// /value (output) : the dot product
// /extraEquation (input) : the number of the extra equation (0,1,...,numberOfExtraEquations-1)
// 
// /Return values: 0=success
// /Author: wdh
//==================================================================================
{
  real sumOfExtraEquationCoefficients=0.;
  return evaluateExtraEquation(u,value,sumOfExtraEquationCoefficients,extraEquation);
}

int PETScSolver::
evaluateExtraEquation( const realCompositeGridFunction & u, real & value, real & sumOfExtraEquationCoefficients,
                       int extraEquation /* =0 */ )
//==================================================================================
// /Description:
//    Evaluate the dot product of the coefficients of an extra equation times u 
//  Also return the sum of the coefficients of the extra equation (i.e. the dot product with the "1" vector)
//
// /u (input) : grid function to dot with the extra equation
// /value (output) : the dot product
// /extraEquation (input) : the number of the extra equation (0,1,...,numberOfExtraEquations-1)
// 
// /Return values: 0=success
// /Author: wdh
//==================================================================================
{
  const int myid=max(0,Communication_Manager::My_Process_Number);
  assert( extraEquation==0 );

  if( oges.getCompatibilityConstraint() &&
      problemIsSingular==notSingular )
  {
    problemIsSingular=addExtraEquation;
    if( Oges::debug & 1 )
      printF("PETScSolver:evaluateExtraEquation:INFO: problem is singular -- using addExtraEquation option\n");
  }

  const CompositeGrid & cg = *u.getCompositeGrid();

  Index I1,I2,I3;
  value=0.;
  sumOfExtraEquationCoefficients=0.;
  if( problemIsSingular==addExtraEquation )
  {
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      const MappedGrid & mg = cg[grid];

      // *** fix this for parallel ***
      const realSerialArray & uLocal = u[grid].getLocalArray();
      const intSerialArray & mask = mg.mask().getLocalArray();
      const intSerialArray & classify = oges.coeff[grid].sparse->classify.getLocalArray();

      getIndex(mg.indexRange(),I1,I2,I3);  // do not include ghost or periodic
      int includeGhost=0;  // do NOT include ghost pts
      bool ok = ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3,includeGhost);
      if( !ok ) continue;      

      // This needs to match how the constraint is filled in function buildMatrix
      int i1,i2,i3,n=0;
      FOR_3D(i1,i2,i3,I1,I2,I3)
      {
	if( classify(i1,i2,i3,n)>0 )  // include interior and boundary 
	{
	  value+=uLocal(i1,i2,i3);
	  sumOfExtraEquationCoefficients+=1;
	}
      }
      
    }
    // printf("PETScSolver::evaluateExtraEquation: myid=%i value=%g\n",myid,value);

    // note: get sum over all processors since all processors are assumed to have a need to know 
    value=ParallelUtility::getSum(value);
    sumOfExtraEquationCoefficients=ParallelUtility::getSum(sumOfExtraEquationCoefficients);
    

    // printF("PETScSolver::evaluateExtraEquation: total-value=%14.9g\n",value);
  }
  else
  {
    printF("PETScSolver::evaluateExtraEquation:ERROR: not implemented for problemIsSingular=%i\n",problemIsSingular);
    return 1;
  }
  return 0;
}


#endif


