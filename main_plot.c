#define PI    3.141592653589793
#define EXP1  2.718281828459045

#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include "main_plot.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    /* ====================================         MPI_INIT       ==================================== */
    MPI_Init(&argc,&argv);
    int ierror;
    int nproc,rank;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    int ia;
    clock_t time_start, time_end; time_start=clock();
    /* ====================================      PROBLEM SETUP     ==================================== */
    int jobnumber = 150;
    int ndim ; set_ndim(&ndim)  ;
    int nddim; set_nddim(&nddim);

    /* ------------------------------------          MESH          -------------------------------------*/
    int mref_base= 7;
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
    double par_mat[]={  0.25,  0.50,  585.,  180.,  0.  }; assert_par_mat(par_mat);
    //              ux,   uy,   uz, ux_n, uy_n, uz_n,       //                     ______________________
    int bc_type[]={  0,    0,    0,    0,    0,    0,       // Face 0: YZ         /                     /|
                     0,    1,    1,    0,    1,    1,       // Face 1: YZ        /          5          / |
                     1,    1,    1,    1,    1,    1,       // Face 2: ZX       /_____________________/  |
                     1,    1,    1,    1,    1,    1,       // Face 3: ZX       |                     |  |
                     1,    1,    1,    1,    1,    1,       // Face 4: XY     0 |                     | 1|    Z
                     1,    1,    1,    1,    1,    1  };    // Face 5: XY       |          2          |  /    |  Y
                                                            //                  |                     | /     | /        0: Dirichlet
                                                            //                  | ____________________|/      |/____X    1: Neumann
    // Grad u average:
    double par_dirichlet[]={0.0,0.0,0.0,
                            0.0,0.0,0.0,
                            0.0,0.0,0.0 };
    double par_neumann[]  ={0.01};
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

    /* ====================================       POSTPROCESS      ==================================== */
    double *U_;
    double  u;
    FILE   *fptr0,*fptr1;
    char   fname0[128],fname1[128];
    int    lSize;
    char   par_lE_char[32];
    
    /* ------------------------------------          MATHGL        -------------------------------------*/
    mgl_def_font("STIX","/opt/mathgl/mathgl-2.3/fonts/");//adventor, bonum, chorus, cursor, heros, heroscn, pagella, schola, STIX, termes.
    HMGL    gr;
    HMDT    xplot;
    HMDT    yplot;

    U_=(double*) malloc(nselfIdx*sizeof(double));
    for (int imode=21; imode<35; imode++)
    {
        for (int ipar=6250000; ipar<6250000+1; ipar+=1250000)
        {
            par_mat[4]=1.e-8*ipar; 
            sprintf(fname0,"data/%03d/U_cube_%d_%d_%d_%08.8f_%08.8f_%08.8f_%08.8f_%08.8f.bin",jobnumber,mref,porder,imode,par_mat[0],par_mat[1],par_mat[3],par_mat[2]/par_mat[3],par_mat[4]);
            if (file_exist(fname0,rank)==1)
            {
                // ---- e1 e2
                if (rank==0)
                {
                    fptr0=fopen(fname0,"rb");
                    lSize=fread(U_universal,sizeof(double),nDof_global,fptr0); assert(lSize==nDof_global);
                    fclose(fptr0);
                    for (int i=0; i<nDof_global; i++) { U_universal_temp[i]=U_universal[universalIdx[i]]; }
                }
                ierror=MPI_Scatterv(U_universal_temp,nDof_array,iDof_displ_array,MPI_DOUBLE,U_,nselfIdx,MPI_DOUBLE,0,MPI_COMM_WORLD);assert(ierror==MPI_SUCCESS);
                // ---- dens. plot e1
                //gr=mgl_create_graph(1200,1200);
                //mgl_set_font_size(gr,3); // 3 (paper) or 4 (slide)
                //plot_solution(gr,U_,&app,4,nproc,rank);
                //sprintf(fname0,"plot/%03d/e1_cube_%d_%d_%d_%08.8f_%08.8f_%08.8f_%08.8f_%08.8f.png",jobnumber,mref,porder,imode,par_mat[0],par_mat[1],par_mat[3],par_mat[2]/par_mat[3],par_mat[4]);
                //if (rank==0) { mgl_write_frame(gr,fname0,""); }
                //mgl_delete_graph(gr);
                // ---- dens. plot e2 
                //gr=mgl_create_graph(1200,1200);
                //mgl_set_font_size(gr,3); // 3 (paper) or 4 (slide)
                //plot_solution(gr,U_,&app,5,nproc,rank);
                //sprintf(fname0,"plot/%03d/e2_cube_%d_%d_%d_%08.8f_%08.8f_%08.8f_%08.8f_%08.8f.png",jobnumber,mref,porder,imode,par_mat[0],par_mat[1],par_mat[3],par_mat[2]/par_mat[3],par_mat[4]);
                //if (rank==0) { mgl_write_frame(gr,fname0,""); }
                //mgl_delete_graph(gr);
                // ---- dens. plot phase
                gr=mgl_create_graph(1200,1200);
                mgl_set_font_size(gr,3); // 3 (paper) or 4 (slide)
                plot_solution(gr,U_,&app,6,nproc,rank);
                sprintf(fname0,"plot/%03d/phase_cube_%d_%d_%d_%08.8f_%08.8f_%08.8f_%08.8f_%08.8f_xy.png",jobnumber,mref,porder,imode,par_mat[0],par_mat[1],par_mat[3],par_mat[2]/par_mat[3],par_mat[4]);
                if (rank==0) { mgl_write_frame(gr,fname0,""); }
                mgl_delete_graph(gr);
            }
        }
    }
    // ----
    free(U_);

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

    time_end=clock(); if (rank==0) { printf("Time consumed: %f[min].\n", ((float)(time_end-time_start))/CLOCKS_PER_SEC/60.0 ); }

    /* ====================================        MPI_FINAL       ==================================== */
    MPI_Finalize();

    return 0;
}
