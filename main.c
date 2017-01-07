static char help[] = "BVP of gradient elasticity.\n";

#define PI    3.141592653589793
#define EXP1  2.718281828459045

#include "petscsnes.h"
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
    PetscInitialize(&argc,&argv,(char*)0,help);
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
    int mref_base= 3;
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
    Vec            Residual;
    Mat            Tangent;
    SNES           snes;        
    KSP            ksp;         
    PC             pc;
    SNESLineSearch ls;
    int            snes_niter,ksp_niter;
    SNESConvergedReason snes_reason;
    // ---- Vec
    ierr = VecCreateMPI(MPI_COMM_WORLD,nDof,nDof_global,&Residual);CHKERRQ(ierr);
    ierr = VecDuplicate(Residual,&U);CHKERRQ(ierr);
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
    // ---- SNES
    ierr = SNESCreate(MPI_COMM_WORLD,&snes);CHKERRQ(ierr);
    ierr = SNESSetFunction(snes,Residual,compute_Residual,&app);CHKERRQ(ierr);
    ierr = SNESSetJacobian(snes,Tangent,Tangent,compute_Tangent,&app);CHKERRQ(ierr);
    ierr = SNESSetType(snes,SNESNEWTONLS);CHKERRQ(ierr);
    ierr = SNESGetLineSearch(snes,&ls);CHKERRQ(ierr);
    ierr = SNESLineSearchSetType(ls,SNESLINESEARCHBT);CHKERRQ(ierr); //SNESLINESEARCHBT, SNESLINESEARCHL2, SNESLINESEARCHCP, SNESLINESEARCHBASIC, SNESLINESEARCHSHELL
    ierr = SNESLineSearchSetOrder(ls,SNES_LINESEARCH_ORDER_CUBIC);CHKERRQ(ierr); //LINEAR, QUADRATIC, CUBIC
    // ---- KSP/PC
    ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr); 
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(ksp,1);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,1e8,1e8);CHKERRQ(ierr);
    ierr = KSPSetType(ksp,KSPMINRES);CHKERRQ(ierr); // The operator and the preconditioner must be symmetric and the preconditioner must be positive definite for this method.
    ierr = KSPSetPCSide(ksp,PC_LEFT);CHKERRQ(ierr); //
    ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
    ierr = PCJacobiSetUseAbs(pc,1);CHKERRQ(ierr);
    /*    ierr = KSPSetType(ksp,KSPGMRES);CHKERRQ(ierr);
    ierr = KSPGMRESSetRestart(ksp,10000);CHKERRQ(ierr); // eats memory. 2e4 doesn't work on edison; srun error.
    ierr = KSPGMRESSetOrthogonalization(ksp,KSPGMRESModifiedGramSchmidtOrthogonalization);CHKERRQ(ierr);//KSPGMRESClassicalGramSchmidtOrthogonalization
    ierr = KSPSetPCSide(ksp,PC_RIGHT);CHKERRQ(ierr);//PC_RIGHT, PC_LEFT, PC_SYMMETRIC
    ierr = PCSetType(pc,PCASM);CHKERRQ(ierr);
    ierr = PCASMSetType(pc,PC_ASM_RESTRICT);CHKERRQ(ierr);//PC_ASM_BASIC, PC_ASM_INTERPOLATE, PC_ASM_RESTRICT, PC_ASM_NONE*/
    // ---- overwrite with runtime option
    ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
    
    /* ------------------------------------ PATH-FOLLOW/REFINE ETC.-------------------------------------*/
    int problem_type=0;
    int snes_max_it;
    int ipar_step,ipar_step_min,ipar_step_max;
    switch (problem_type)                             // 2^13=       8192
    {                                                 // 2^14=     1,6384
        case 0: // arbitrary init. guess              // 2^15=     3,2768
            ipar_step     = 1;                        // 2^16=     6,5536
            ipar_step_min =-1;                        // 2^17=    13,1072
            ipar_step_max =-1;                        // 2^18=    26,2144
            snes_max_it   = 20;                       // 2^19=    52,4288
            break;                                    // 2^20=   104,8576
        case 1: // path-following                     // 2^21=   209,7152
            ipar_step     = 1000000+0*int_pow(2,13);            // 2^22=   419,4304
            ipar_step_min = 1000000+0*int_pow(2,13);            // 2^23=   838,8608
            ipar_step_max = int_pow(2,24);            // 2^24=  1677,7216
            snes_max_it   = 36;                          
            break;                                   
        case 2: // refinement
            ipar_step     = 1048576;
            ipar_step_min =-1;
            ipar_step_max =-1;
            snes_max_it   = 1;
            break;
    }
    int    ipar_begin =10000000;
    int    ipar_end   =10000000;
    int    dir        = (ipar_begin < ipar_end ? 1:-1);
    int    converge_count=0,diverge_count=0;
    FILE   *fptr0;
    char   fname0[128];
    double *U_,Energy,Energy_global;
    // ---- 
    for (int imode=0; imode<1; imode++) {
    // ---- increment ipar
    for (int ipar=ipar_begin; dir*(ipar_end-ipar)>=0; )
    {
        par_mat[4]=1.e-8*ipar; assert_par_mat(par_mat);
        
        sprintf(fname0,"data/%03d/U_cube_%d_%d_%d_%08.8f_%08.8f_%08.8f_%08.8f_%08.8f.bin",jobnumber,mref_base,porder,imode,par_mat[0],par_mat[1],par_mat[3],par_mat[2]/par_mat[3],par_mat[4]);
        if ( (problem_type != 2) || ( (problem_type == 2) && ((fptr0=fopen(fname0,"rb")) != NULL) ) )
        {
            if (rank==0) { printf("r=%04.4f, beta=%04.4f, e3=%04.4f, er=%04.4f, lE=%08.8f, ipar=%d, ipar_step=%d\n",par_mat[0],par_mat[1],par_mat[3],par_mat[2]/par_mat[3],par_mat[4],ipar,ipar_step); }
            // ---- load/set initial guess.
            ierr = VecGetArray(U,&U_);CHKERRQ(ierr);
            if (rank==0)
            {
                switch (problem_type) {
                    case 0: for (int i=0; i<nDof_global; i++) { U_universal[i]=0.0+3.0e-2*sin(5.0*i); } break;
                    case 1: sprintf(fname0,"data/%03d/U_cube_%d_%d_temp%d.bin",jobnumber,mref,porder,imode);
                            fptr0=fopen(fname0,"rb"); ierror=fread(U_universal,sizeof(double),nDof_global,fptr0); assert(ierror==nDof_global); fclose(fptr0); break;
                    case 2: insert_knot_uniform(U_universal,fptr0,nelem_global,nref,porder,ndim,ndof); break;
                    default: exit(0); }
                for (int i=0; i<nDof_global; i++) { U_universal_temp[i]=U_universal[universalIdx[i]]; }
            }
            if (problem_type==2) { fclose(fptr0); }
            ierror=MPI_Scatterv(U_universal_temp,nDof_array,iDof_displ_array,MPI_DOUBLE,U_,nselfIdx,MPI_DOUBLE,0,MPI_COMM_WORLD);assert(ierror==MPI_SUCCESS);
            ierr = VecRestoreArray(U,&U_);CHKERRQ(ierr);
            // ---- subiterate: to save U at every snes_max_it nonlin. iterations.
            for (int isubiter=0; isubiter<50; isubiter++)
            {
                ierr = SNESSetTolerances(snes,1.e-12,1.e-12,1.e-20,snes_max_it,1e9);CHKERRQ(ierr);
                // ---- solve nonlinear system.
                ierr = SNESSolve(snes,NULL,U);CHKERRQ(ierr);
                ierr = KSPGetIterationNumber(ksp,&ksp_niter);CHKERRQ(ierr);
                ierr = SNESGetIterationNumber(snes,&snes_niter);CHKERRQ(ierr);
                ierr = SNESGetConvergedReason(snes,&snes_reason);CHKERRQ(ierr);
                time_end=clock();
                if (rank==0) { printf("Time=%f[min], reason = %s, #nonlin.ite.= %d, #lin.ite.= %d.\n",((float)(time_end-time_start))/CLOCKS_PER_SEC/60.0,SNESConvergedReasons[snes_reason],snes_niter,ksp_niter); }
                if ( (problem_type != 2) || (snes_reason != SNES_DIVERGED_MAX_IT) ) { break; }
                // ---- write solution
                ierr = VecGetArray(U,&U_);CHKERRQ(ierr);
                ierror=MPI_Gatherv(U_,nselfIdx,MPI_DOUBLE,U_universal_temp,nDof_array,iDof_displ_array,MPI_DOUBLE,0,MPI_COMM_WORLD);assert(ierror==MPI_SUCCESS);
                if (rank==0)
                {
                    for (int i=0; i<nDof_global; i++) { U_universal[universalIdx[i]]=U_universal_temp[i]; }
                    sprintf(fname0,"data/%03d/U_cube_%d_%d_%d_%08.8f_%08.8f_%08.8f_%08.8f_%08.8f.bin",jobnumber,mref,porder,imode,par_mat[0],par_mat[1],par_mat[3],par_mat[2]/par_mat[3],par_mat[4]);
                    fptr0=fopen(fname0,"wb");fwrite(U_universal,sizeof(double),nDof_global,fptr0);fclose(fptr0);
                }
                ierr = VecRestoreArray(U,&U_);CHKERRQ(ierr);
            }
            ierr = VecGetArray(U,&U_);CHKERRQ(ierr);
            // ---- write solution
            ierror=MPI_Gatherv(U_,nselfIdx,MPI_DOUBLE,U_universal_temp,nDof_array,iDof_displ_array,MPI_DOUBLE,0,MPI_COMM_WORLD);assert(ierror==MPI_SUCCESS);
            if (rank==0 && snes_reason > 0)
            {
                for (int i=0; i<nDof_global; i++) { U_universal[universalIdx[i]]=U_universal_temp[i]; }
                sprintf(fname0,"data/%03d/U_cube_%d_%d_%d_%08.8f_%08.8f_%08.8f_%08.8f_%08.8f.bin",jobnumber,mref,porder,imode,par_mat[0],par_mat[1],par_mat[3],par_mat[2]/par_mat[3],par_mat[4]);
                fptr0=fopen(fname0,"wb");fwrite(U_universal,sizeof(double),nDof_global,fptr0);fclose(fptr0);
                if (problem_type==1){
                sprintf(fname0,"data/%03d/U_cube_%d_%d_temp%d.bin",jobnumber,mref,porder,imode);
                fptr0=fopen(fname0,"wb");fwrite(U_universal,sizeof(double),nDof_global,fptr0);fclose(fptr0);
                }
            }
            // ---- compute/write Energy
            compute_Energy(&Energy,U_,&app);
            ierror=MPI_Reduce(&Energy,&Energy_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);assert(ierror==MPI_SUCCESS);
            if (rank==0 && snes_reason > 0)
            {
                printf("%e\n",Energy_global);
                sprintf(fname0,"data/%03d/PI_cube_%d_%d_%d_%08.8f_%08.8f_%08.8f_%08.8f_%08.8f.bin",jobnumber,mref,porder,imode,par_mat[0],par_mat[1],par_mat[3],par_mat[2]/par_mat[3],par_mat[4]);
                fptr0=fopen(fname0,"wb");fwrite(&Energy_global,sizeof(double),1,fptr0);fclose(fptr0);
            }
            ierr = VecRestoreArray(U,&U_);CHKERRQ(ierr);
            // ---- control step-size for path-following.
            if (problem_type==1)
            {
                if (snes_reason > 0) { converge_count++; diverge_count=0; }
                else            { converge_count=0; diverge_count++; }

                if ( (converge_count==16) && (ipar_step < ipar_step_max) )
                {
                    converge_count=0;
                    ipar_step*=2;
                }
                else if (diverge_count == 1)
                {
                    if (ipar_step > ipar_step_min)
                    {
                        diverge_count=0;
                        ipar-=1*(dir*ipar_step);
                        ipar_step/=2;
                    }
                }
                else if (diverge_count ==   4) { ipar_step*=2; snes_max_it=30; }
                else if (diverge_count ==   8) { ipar_step*=2; snes_max_it=30; }
                else if (diverge_count ==  16) { ipar_step*=2; snes_max_it=30; }
                else if (diverge_count ==  32) { ipar_step*=2; snes_max_it=30; }
                else if (diverge_count ==  64) { ipar_step*=2; snes_max_it=30; }
                else if (diverge_count == 128) { ipar_step*=2; snes_max_it=30; }
                else if (diverge_count == 1000) { break; }
            }
        } // if ( (problem_type != 2) || ( (problem_type == 2) && ((fptr0=fopen(fname0,"rb")) != NULL) ) )
        // ---- increment ipar
        if ( dir*(ipar_end-ipar)>=ipar_step )  { ipar+=(dir*ipar_step); }
        else if ( dir*(ipar_end-ipar)>1.e-15 ) { ipar=ipar_end; }
        else                                   { break; }
    } // for ipar
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
    ierr = VecDestroy(&Residual);CHKERRQ(ierr);
    ierr = MatDestroy(&Tangent);CHKERRQ(ierr);
    ierr = SNESDestroy(&snes);CHKERRQ(ierr);

    time_end=clock(); if (rank==0) { printf("Time consumed: %f[min].\n", ((float)(time_end-time_start))/CLOCKS_PER_SEC/60.0 ); }

    /* ====================================        MPI_FINAL       ==================================== */
    PetscFinalize();

    return 0;
}
