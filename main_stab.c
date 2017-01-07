static char help[] = "Performs stability analysis.\n";

#define PI    3.141592653589793
#define EXP1  2.718281828459045

#include "petscsnes.h"
#include "slepceps.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include "main.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    /* ====================================         MPI_INIT       ==================================== */
    SlepcInitialize(&argc,&argv,(char*)0,help);
    int ierr;
    int nproc,rank;
    ierr=MPI_Comm_size(MPI_COMM_WORLD,&nproc);CHKERRQ(ierr);
    ierr=MPI_Comm_rank(MPI_COMM_WORLD,&rank);CHKERRQ(ierr);
    
    int ierror,ia;
    clock_t time_start, time_end; time_start=clock();
    /* ====================================      PROBLEM SETUP     ==================================== */
    int jobnumber = 100;
    int ndim ; set_ndim(&ndim)  ;
    int nddim; set_nddim(&nddim);

    /* ------------------------------------          MESH          -------------------------------------*/
    int mref_base= 4;
    int nref     = 0;
    int mref     = mref_base+nref;
    int porder   = 2;
    int nelem_global[ndim];
    double *knotVector_global[ndim];
    if (rank==0) { mesh_uniform(nelem_global,knotVector_global,mref,porder,ndim); }
    ierror=MPI_Bcast(nelem_global,ndim,MPI_INT,0,MPI_COMM_WORLD);assert(ierror==MPI_SUCCESS);
    // ----  MPI: PARTITION
    int npart[ndim],ipart[ndim];
    int nelem[ndim],nknot[ndim],ielem_displ[ndim];
    double *knotVector[ndim];
    set_mpi_part(npart,ipart,nelem,ielem_displ,nknot,knotVector,nelem_global,knotVector_global,ndim,porder,nproc,rank);
    if (rank==0) { for (int idim=0; idim<ndim; idim++) { free(knotVector_global[idim]); } } // allocated in "mesh_uniform.c"
    int nbasis[ndim]; for (int idim=0; idim<ndim; idim++) { nbasis[idim] = ( ipart[idim]==npart[idim]-1 ? nelem[idim]+porder : nelem[idim]); }

    /* ------------------------------------           BVP          -------------------------------------*/
    int ndof ; set_ndof(&ndof)  ;
    //                    r ,  curv,   e0 , e345 , le,
    double par_mat[]={  0.25,  0.50,  585.,  180.,  0.1 }; assert_par_mat(par_mat);
    //              ux,   uy,   uz, ux_n, uy_n, uz_n,       //                     ______________________
    int bc_type[]={  0,    1,    1,    1,    1,    1,       // Face 0: YZ         /                     /|
                     0,    1,    1,    1,    1,    1,       // Face 1: YZ        /          5          / |
                     1,    0,    1,    1,    1,    1,       // Face 2: ZX       /_____________________/  |
                     1,    0,    1,    1,    1,    1,       // Face 3: ZX       |                     |  |
                     1,    1,    0,    1,    1,    1,       // Face 4: XY     0 |                     | 1|    Z
                     1,    1,    0,    1,    1,    1  };    // Face 5: XY       |          2          |  /    |  Y
                                                            //                  |                     | /     | /        0: Dirichlet
                                                            //                  | ____________________|/      |/____X    1: Neumann
    // Grad u average:
    double par_dirichlet[]={0.01,0.01,0.0,
                            0.02,0.02,0.0,
                            0.03,0.03,0.0 };
    double par_neumann[]  ={0.00};
    // ---- MPI: 
    int nboun = (int)(ipart[0]==0)+(int)(ipart[0]==npart[0]-1)+(int)(ipart[1]==0)+(int)(ipart[1]==npart[1]-1)+(int)(ipart[2]==0)+(int)(ipart[2]==npart[2]-1);
    int boun[nboun];
    ia=0;
    if (ipart[0]==0) { boun[ia++]=0; } if (ipart[0]==npart[0]-1) { boun[ia++]=1; }
    if (ipart[1]==0) { boun[ia++]=2; } if (ipart[1]==npart[1]-1) { boun[ia++]=3; }
    if (ipart[2]==0) { boun[ia++]=4; } if (ipart[2]==npart[2]-1) { boun[ia++]=5; }

    /* ------------------------------------          IGA           -------------------------------------*/
    int            nquad=4;
    double         cquad[nquad];
    double         wquad[nquad];
    legendre_handle(cquad,wquad,nquad,0.,1.);
    // ---- MPI : COMMUNICATION
    int nselfIdx,nsendPart,nsendIdx,nrecvPart,nrecvIdx,*selfIdx,*sendPart,*sendIdx,*sendPtr,*recvPart,*recvIdx,*recvPtr;
    set_mpi_comm(&nselfIdx,&selfIdx,&nsendPart,&sendPart,&nsendIdx,&sendIdx,&sendPtr,&nrecvPart,&recvPart,&nrecvIdx,&recvIdx,&recvPtr,ndof,porder,rank,npart,ipart,nelem,nbasis);
    // ---- universalIdx (for partition independant data storage)
    int nDof_global = ndof*((nelem_global[0]+porder)*(nelem_global[1]+porder)*(nelem_global[2]+porder));
    int nDof        = ndof*(nbasis[0]*nbasis[1]*nbasis[2]);
    int *nDof_array;          if (rank==0) { nDof_array      =(int*) malloc(nproc*sizeof(int)); }
    int *iDof_displ_array;    if (rank==0) { iDof_displ_array=(int*) malloc(nproc*sizeof(int)); }
    int *universalIdx;        if (rank==0) { universalIdx    =(int*) malloc(nDof_global*sizeof(int)); };
    set_mpi_universalIdx(universalIdx,nDof_array,iDof_displ_array,nDof,nselfIdx,nbasis,ielem_displ,nelem_global,ndof,porder,nproc,rank);
    double *U_universal     ; if (rank==0) { U_universal     =(double*) malloc(nDof_global*sizeof(double)); }
    double *U_universal_temp; if (rank==0) { U_universal_temp=(double*) malloc(nDof_global*sizeof(double)); }
    // ---- globalIdx (maps local idx to global idx)
    int *globalIdx      = (int*) malloc((nselfIdx+nrecvIdx)*sizeof(int));
    set_mpi_globalIdx(globalIdx,iDof_displ_array,nselfIdx,selfIdx,nsendPart,sendPart,nsendIdx,sendIdx,sendPtr,nrecvPart,recvPart,nrecvIdx,recvIdx,recvPtr,nbasis,ndof,nproc,rank);

    /* ------------------------------------          DATA          -------------------------------------*/
    double *U_;
    FILE                *fptr0,*fptr1;
    char                fname0[128],fname1[128];
    int                 lSize;
    char par_lE_char[32];


    /* ====================================    ASSEMBLE / SOLVE    ==================================== */

    /* ------------------------------------           APP          -------------------------------------*/
    App app;

    app.ndim         = ndim;
    app.nddim        = nddim;
    // ---- BVP
    app.ndof         = ndof;
    app.par_mat      = par_mat;
    app.bc_type      = bc_type;
    app.nboun        = nboun;
    app.boun         = boun;
    app.par_dirichlet= par_dirichlet;
    app.par_neumann  = par_neumann;
    // ---- MESH
    app.porder       = porder;
    app.nelem        = nelem;
    app.nbasis       = nbasis;
    app.nknot        = nknot;
    app.knotVector   = knotVector;
    // ---- IGA
    app.nquad        = nquad;
    app.cquad        = cquad;
    app.wquad        = wquad;
    app.nselfIdx     = nselfIdx;
    app.selfIdx      = selfIdx;
    app.nsendPart    = nsendPart;
    app.sendPart     = sendPart;
    app.nsendIdx     = nsendIdx;
    app.sendIdx      = sendIdx;
    app.sendPtr      = sendPtr;
    app.nrecvPart    = nrecvPart;
    app.recvPart     = recvPart;
    app.nrecvIdx     = nrecvIdx;
    app.recvIdx      = recvIdx;
    app.recvPtr      = recvPtr;
    app.globalIdx    = globalIdx;

    /* ------------------------------------          PETSc         -------------------------------------*/
    Vec            U;
    Mat            Tangent;
    // ---- Vec
    ierr = VecCreateMPI(MPI_COMM_WORLD,nDof,nDof_global,&U);CHKERRQ(ierr);
    // ---- Mat
    int *d_nz=(int*) malloc(nDof*sizeof(int)),*o_nz=(int*) malloc(nDof*sizeof(int));
    set_nz(d_nz,o_nz,nbasis,npart,ipart,ndof,porder);
    ierr = MatCreateAIJ(MPI_COMM_WORLD,nDof,nDof,nDof_global,nDof_global,0,d_nz,0,o_nz,&Tangent);CHKERRQ(ierr);
    free(d_nz);free(o_nz);
    ierr = MatSetOption(Tangent,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
    ierr = MatSetOption(Tangent,MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE);CHKERRQ(ierr);
    ierr = MatSetOption(Tangent,MAT_SUBSET_OFF_PROC_ENTRIES,PETSC_TRUE);CHKERRQ(ierr);
    ierr = MatSetOption(Tangent,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);
    ierr = MatSetOption(Tangent,MAT_SYMMETRY_ETERNAL,PETSC_TRUE);CHKERRQ(ierr);
    
    /* ====================================        STABILITY       ==================================== */
    if (rank==0) { printf("Stability anlysis:\n\n"); }

    /* ------------------------------------          SLEPc         -------------------------------------*/
    EPS                eps;
    EPSType            eps_type;
    PetscReal          eps_error,eps_vreal,eps_vimag;
    PetscInt           eps_its;
    EPSConvergedReason eps_reason;
    int canfgets=0;

    /* ------------------------------------       EIGENSOLVER      -------------------------------------*/
    for (int imode=0; imode<1; imode++)
    {
        sprintf(fname0,"data/%03d/PI_cube_%d_%d_%d_%08.8f_%08.8f_%08.8f_%08.8f_lE.txt",jobnumber,mref,porder,imode,par_mat[0],par_mat[1],par_mat[3],par_mat[2]/par_mat[3]);
    	if (file_exist(fname0,rank)==1)
        {
            if (rank==0) { fptr0=fopen(fname0,"r"); } 
            while (1)
            {
                if (rank==0) { canfgets = ( (fgets(par_lE_char,32,fptr0) != NULL) ? 1:0 ); } ierror=MPI_Bcast(&canfgets,1,MPI_INT,0,MPI_COMM_WORLD);assert(ierror==MPI_SUCCESS);
                if (canfgets==1)
                {
                    if (rank == 0) { par_mat[4]=atof(par_lE_char); } ierror=MPI_Bcast(par_mat+4,1,MPI_DOUBLE,0,MPI_COMM_WORLD);assert(ierror==MPI_SUCCESS);
                    sprintf(fname0,"data/%03d/smallestEigval_%d_%d_%d_%08.8f_%08.8f_%08.8f_%08.8f_%08.8f.bin",jobnumber,mref,porder,imode,par_mat[0],par_mat[1],par_mat[3],par_mat[2]/par_mat[3],par_mat[4]);
                    sprintf(fname1,"data/%03d/U_cube_%d_%d_%d_%08.8f_%08.8f_%08.8f_%08.8f_%08.8f.bin",jobnumber,mref,porder,imode,par_mat[0],par_mat[1],par_mat[3],par_mat[2]/par_mat[3],par_mat[4]);
                    if ( (file_exist(fname0,rank)==0) && (file_exist(fname1,rank)==1) )
                    {
                        if (rank==0)
                        {
                            fptr1=fopen(fname1,"rb");lSize=fread(U_universal,sizeof(double),nDof_global,fptr1);assert(lSize==nDof_global);fclose(fptr1);
                            for (int i=0; i<nDof_global; i++) { U_universal_temp[i]=U_universal[universalIdx[i]]; }
                        }
                        ierr = VecGetArray(U,&U_);CHKERRQ(ierr);
                        ierror=MPI_Scatterv(U_universal_temp,nDof_array,iDof_displ_array,MPI_DOUBLE,U_,nselfIdx,MPI_DOUBLE,0,MPI_COMM_WORLD);assert(ierror==MPI_SUCCESS);
                        ierr = VecRestoreArray(U,&U_);CHKERRQ(ierr);

                        compute_Tangent(NULL,U,Tangent,Tangent,&app);
                        ierr = EPSCreate(MPI_COMM_WORLD,&eps);CHKERRQ(ierr);
                        ierr = EPSSetOperators(eps,Tangent,NULL);CHKERRQ(ierr);
                        ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
                        ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);CHKERRQ(ierr);

                        ierr = EPSSetConvergenceTest(eps,EPS_CONV_ABS);CHKERRQ(ierr);
                        ierr = EPSSetTolerances(eps,1.e-6,2000000);CHKERRQ(ierr);
                        ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
                        ierr = EPSSetUp(eps);CHKERRQ(ierr);
                        ierr = EPSSolve(eps);CHKERRQ(ierr);
                        // ---- output eigs.
                        if (rank==0) { printf("mode %2d,  ",imode); }
                        ierr = EPSGetType(eps,&eps_type);CHKERRQ(ierr);
                        ierr = PetscPrintf(MPI_COMM_WORLD,"%s,  ",eps_type);CHKERRQ(ierr);
                        ierr = EPSGetConvergedReason(eps,&eps_reason);CHKERRQ(ierr);
                        ierr = PetscPrintf(MPI_COMM_WORLD,"conv.reason=%D,  ",eps_reason);CHKERRQ(ierr);
                        ierr = EPSGetIterationNumber(eps,&eps_its);CHKERRQ(ierr);
                        ierr = PetscPrintf(MPI_COMM_WORLD,"num.its=%10d,  ",eps_its);CHKERRQ(ierr);
                        if (eps_reason>0)
                        {
                            ierr = EPSGetEigenpair(eps,0,&eps_vreal,&eps_vimag,NULL,NULL);CHKERRQ(ierr);
                            ierr = EPSComputeError(eps,0,EPS_ERROR_RELATIVE,&eps_error);CHKERRQ(ierr);
                            ierr = PetscPrintf(MPI_COMM_WORLD,"l = %16.16f,  smallest eigv = %16.16f%+16.16f j,  error = %g",par_mat[4],(double)eps_vreal,(double)eps_vimag,(double)eps_error);CHKERRQ(ierr);
                            if (rank==0)
                            {
                                sprintf(fname0,"data/%03d/smallestEigval_%d_%d_%d_%08.8f_%08.8f_%08.8f_%08.8f_%08.8f.bin",jobnumber,mref,porder,imode,par_mat[0],par_mat[1],par_mat[3],par_mat[2]/par_mat[3],par_mat[4]);
                                fptr1=fopen(fname0,"wb");fwrite(&eps_vreal,sizeof(double),1,fptr1);fclose(fptr1);
                            }
                        }
                        if (rank==0) { printf("\n"); }
                        ierr = EPSDestroy(&eps);CHKERRQ(ierr);
                    }
                }
                else { break; }
            } // while
            if (rank==0) { fclose(fptr0); }
        } // if file_exist()
    } // for imode

    /* ====================================         CLEANUP         ==================================== */
    for (int idim=0; idim<ndim; idim++) { free(knotVector[idim]); } // allocated in "set_mpi_part.c"
    free(selfIdx);                                                  // allocated in "set_mpi_comm.c"
    free(sendPart);                                                 //
    free(sendIdx);                                                  //
    free(sendPtr);                                                  //
    free(recvPart);                                                 //
    free(recvIdx);                                                  //
    free(recvPtr);                                                  //
    free(globalIdx);
    if (rank==0)
    {
        free(nDof_array);
        free(iDof_displ_array);
        free(universalIdx);
        free(U_universal);
        free(U_universal_temp);
    }
    ierr = VecDestroy(&U);CHKERRQ(ierr);
    ierr = MatDestroy(&Tangent);CHKERRQ(ierr);

    time_end=clock(); if (rank==0) { printf("Time consumed: %f[min].\n", ((float)(time_end-time_start))/CLOCKS_PER_SEC/60.0 ); }

    /* ====================================        MPI_FINAL       ==================================== */
    SlepcFinalize();

    return 0;
}
