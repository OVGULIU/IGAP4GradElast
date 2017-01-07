/******** definition ********/
void plot_solution(HMGL gr, double *U_self, void *app_, int plotID, int nproc, int rank)
{
    //mgl_set_num_thr(1);
    /* ================ input parameters ================ */
    App *app=(App*)app_;
    int ndim             = app->ndim;
    int nddim            = app->nddim;
    int ndof             = app->ndof;
    double *par_mat      = app->par_mat;
    //int *bc_type         = app->bc_type;
    //double *par_dirichlet= app->par_dirichlet;
    //double *par_neumann  = app->par_neumann;
    // ---- MESH
    int porder           = app->porder;
    int *nelem           = app->nelem;
    int *nbasis          = app->nbasis;
    int *nknot           = app->nknot;
    double **knotVector  = app->knotVector;
    // ---- IGA
    //int nquad            = app->nquad;
    //double *cquad        = app->cquad;
    //double *wquad        = app->wquad;
    // ---- MPI
    int nselfIdx         = app->nselfIdx;
    int *selfIdx         = app->selfIdx;
    int nsendPart        = app->nsendPart;
    int *sendPart        = app->sendPart;
    int nsendIdx         = app->nsendIdx;
    int *sendIdx         = app->sendIdx;
    int *sendPtr         = app->sendPtr;
    int nrecvPart        = app->nrecvPart;
    int *recvPart        = app->recvPart;
    int nrecvIdx         = app->nrecvIdx;
    int *recvIdx         = app->recvIdx;
    int *recvPtr         = app->recvPtr;
    //int *globalIdx       = app->globalIdx;
    //int nboun            = app->nboun;
    //int *boun            = app->boun;

    /* ================ assemble U and Ui ================ */
    double *U           = (double*) malloc((nselfIdx+nrecvIdx)*sizeof(double));
    for (int i=0; i<nselfIdx; i++) { U[selfIdx[i]]=U_self[i]; }
    MPI_Status status;
    int ierror;
    double *U_send=(double*) malloc(nsendIdx*sizeof(double));
    double *U_recv=(double*) malloc(nrecvIdx*sizeof(double));
    for (int i=0; i<nsendIdx; i++) { U_send[i]=U[sendIdx[i]]; }
    for (int isendPart=0; isendPart<nsendPart; isendPart++) { ierror=MPI_Send(U_send+sendPtr[isendPart],sendPtr[isendPart+1]-sendPtr[isendPart],MPI_DOUBLE,sendPart[isendPart],1,MPI_COMM_WORLD);assert(ierror==MPI_SUCCESS); }
    for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) { ierror=MPI_Recv(U_recv+recvPtr[irecvPart],recvPtr[irecvPart+1]-recvPtr[irecvPart],MPI_DOUBLE,recvPart[irecvPart],1,MPI_COMM_WORLD,&status);assert(ierror==MPI_SUCCESS); }
    for (int i=0; i<nrecvIdx; i++) { U[recvIdx[i]]=U_recv[i]; }
    free(U_send);
    free(U_recv);

    /* ================ set up gr ================ */
    char titleStr[36];
    HMGL gr_temp;
    gr_temp = mgl_create_graph(1200,1200);
    HMDT cTickVal=mgl_create_data();
    HMDT xTickVal=mgl_create_data();
    HMDT yTickVal=mgl_create_data();
    
    HMDT zroArray,oneArray;
    HMDT tArray;
    HMDT eArray;
    
    HMDT xArray,yArray,zArray;
    HMDT cArray,Psi_c;
    double rc,xc,yc,zc,x0,y0,Lcube,a0,a1,Lx,Ly,Lz,Ld,Lbase,theta,drel,del,mrgn,x,y,z,*twinBoun,U11,U22,theta_p,theta_m;
    int ntwinBoun;
    double emax;
    
    switch (plotID)
    {
    case 0:
        break;
    case 1:
        break;
    case 2:
        break;
    case 3:
        break;
    case 4:
    case 5:
    case 6:
        //sprintf(titleStr,"t=%8.6f, \\Pi=",t);
        //mgl_title(gr,titleStr,"k:rC",FS);
        if (0) // standard view
        {
            mgl_set_origin(gr,0,0,0);
            mgl_relplot(gr,-0.5,1.5,-0.5,1.5); // Warning: relplot changes fontsize, too. -1.4 would results in a diff. size
            mgl_rotate(gr,60,30,0);
            mgl_aspect(gr,1,1,1);
            mgl_set_ranges(gr,0.5-0.64,0.5+0.64,0.5-0.64,0.5+0.64,0.5-0.64,0.5+0.64);
            //mgl_set_ticks(gr,'x',0.1,4,0.0);
            //mgl_set_ticks(gr,'y',0.1,4,0.0);
            //mgl_set_ticks(gr,'z',0.1,4,0.0);
            //mgl_axis(gr,"xyzU","","");
        }
        else // top/side view
        {
            mgl_set_origin(gr,0,0,0);
            mgl_relplot(gr,-0.45,1.45,-0.45,1.45);
            mgl_rotate(gr,0,0,0);
            mgl_aspect(gr,1,1,1);
            mgl_set_ranges(gr,-0.075,1.075,-0.075,1.075,-0.075,1.075);
            mgl_set_ticks(gr,'x',999.0,0,0.0);
            mgl_set_ticks(gr,'y',999.0,0,0.0);
            mgl_set_ticks(gr,'z',999.0,0,0.0);
            //mgl_axis(gr,"xyzU~","","");
        }
        mgl_set_range_val(gr,'c',-0.25,0.25);
        if (plotID==6) { mgl_set_range_val(gr,'c',0,3); }
        mgl_set_alpha(gr,0);
        
        if (rank < nproc) { plot_volume(gr,U,ndim,nddim,ndof,par_mat,porder,nelem,nbasis,nknot,knotVector,plotID,1.0); }
        
        break;
    default:
        assert(plotID<999);
    } // switch (plotID)
    
    // ---- combine and output graphics
    if (rank!=0)
    { mgl_mpi_send(gr,0); }
    else
    {
        for (int iproc=1; iproc<nproc; iproc++)
        {
            mgl_mpi_recv(gr_temp,iproc);
            mgl_combine_gr(gr,gr_temp);
        }
    }
    // ---- finalize
    mgl_delete_graph(gr_temp);
    mgl_delete_data(xTickVal);
    mgl_delete_data(yTickVal);
    mgl_delete_data(cTickVal);
    free(U);
}

/*
 * ---- plot_time_energy.c & plot_convergence.c
 * LineWidth = '2' : paper
 *             '3' : slides
 * 
 * ---- plot_solution.c
 * MarkerSize= '1' : paper
 *             '1' : slides
 * 
 * ---- main.c
 * Fontsize = 3 (x1.4) : paper
 *            4 (x1.4) : slide
 * 
*/



