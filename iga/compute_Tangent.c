#undef __FUNCT__
#define __FUNCT__ "compute_Tangent"

PetscErrorCode compute_Tangent(SNES snes, Vec Ui_, Mat Tangent_, Mat Pmat_, void *app_)
{
    int ierr;

    /* ================ input parameters ================ */
    App *app=(App*)app_;
    // ---- PDE
    int ndim             = app->ndim;
    int nddim            = app->nddim;
    int ndof             = app->ndof;
    double *par_mat      = app->par_mat;
    int *bc_type         = app->bc_type;
    //double *par_dirichlet= app->par_dirichlet;
    //double *par_neumann  = app->par_neumann;
    int nboun            = app->nboun;
    int *boun            = app->boun;
    // ---- MESH
    int porder           = app->porder;
    int *nelem           = app->nelem;
    int *nbasis          = app->nbasis;
    int *nknot           = app->nknot;
    double **knotVector  = app->knotVector;
    // ---- IGA
    int nquad            = app->nquad;
    double *cquad        = app->cquad;
    double *wquad        = app->wquad;
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
    int *globalIdx       = app->globalIdx;
    /* ================ assemble Ui ================ */
    double *Ui           = (double*) malloc((nselfIdx+nrecvIdx)*sizeof(double));
    const double *Ui_self;
    ierr = VecGetArrayRead(Ui_,&Ui_self);CHKERRQ(ierr);
    for (int i=0; i<nselfIdx; i++)
    {
        Ui[selfIdx[i]]=Ui_self[i];
    }
    ierr = VecRestoreArrayRead(Ui_,&Ui_self);CHKERRQ(ierr);
    MPI_Status status;
    int ierror;
    double *U_send=(double*) malloc(nsendIdx*sizeof(double));
    double *U_recv=(double*) malloc(nrecvIdx*sizeof(double));
    for (int i=0; i<nsendIdx; i++) { U_send[i]=Ui[sendIdx[i]]; }
    for (int isendPart=0; isendPart<nsendPart; isendPart++) 
    { ierror=MPI_Send(U_send+sendPtr[isendPart],sendPtr[isendPart+1]-sendPtr[isendPart],MPI_DOUBLE,sendPart[isendPart],32,MPI_COMM_WORLD); assert(ierror==MPI_SUCCESS);}
    for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) 
    { ierror=MPI_Recv(U_recv+recvPtr[irecvPart],recvPtr[irecvPart+1]-recvPtr[irecvPart],MPI_DOUBLE,recvPart[irecvPart],32,MPI_COMM_WORLD,&status); assert(ierror==MPI_SUCCESS);}
    for (int i=0; i<nrecvIdx; i++) { Ui[recvIdx[i]]=U_recv[i]; }
    free(U_send);
    free(U_recv);
    /* ================ assemble Tangent_ ================ */
    // ---- Mesh
    int ibpe,lbpe;
    int nbpe = int_pow(porder+1,3);
    int ibasis;
    double X0[ndim],X1[ndim];
    // ---- PDE
    double *tangent = (double*) malloc(int_pow(ndof*nddim,2)*sizeof(double));
    double ui[ndof*nddim];
    // ---- quadrature
    double xq[ndim],xi;
    double weight;
    // ---- IGA
    double N[nbpe*nddim];
    ierr = MatZeroEntries(Pmat_);CHKERRQ(ierr);
    double val[int_pow(ndof*nbpe,2)];
    int         idx[ndof*nbpe];
    double      temp;
    int         ia,row,col;
    int nbasis_y_active=nelem[1]+porder;
    int nbasis_z_active=nelem[2]+porder;
    /* ================================================================   VOLUME INTEGRAL   ================================================================ */
    // ---- loop over elements
    for (int ielem_x=0; ielem_x<nelem[0]; ielem_x++) {
    for (int ielem_y=0; ielem_y<nelem[1]; ielem_y++) {
    for (int ielem_z=0; ielem_z<nelem[2]; ielem_z++) {
    for (int i=0; i<int_pow(ndof*nbpe,2); i++) { val[i]=0.0; }
    // ---- loop over quadrature points
    for (int iquad_x=0; iquad_x<nquad; iquad_x++) {
    for (int iquad_y=0; iquad_y<nquad; iquad_y++) {
    for (int iquad_z=0; iquad_z<nquad; iquad_z++) {
        // ---- evaluate quadrature coord.
        X0[0]=knotVector[0][ielem_x+porder]; X1[0]=knotVector[0][ielem_x+porder+1];
        X0[1]=knotVector[1][ielem_y+porder]; X1[1]=knotVector[1][ielem_y+porder+1];
        X0[2]=knotVector[2][ielem_z+porder]; X1[2]=knotVector[2][ielem_z+porder+1];
        xi=cquad[iquad_x];xq[0]=(-xi+1.0)*X0[0]+xi*X1[0];
        xi=cquad[iquad_y];xq[1]=(-xi+1.0)*X0[1]+xi*X1[1];
        xi=cquad[iquad_z];xq[2]=(-xi+1.0)*X0[2]+xi*X1[2];
        // ---- evaluate N
        ia=0;
        for (int ibpe_x=0; ibpe_x<porder+1; ibpe_x++){
        for (int ibpe_y=0; ibpe_y<porder+1; ibpe_y++){
        for (int ibpe_z=0; ibpe_z<porder+1; ibpe_z++){
            N[ia++]=evalN(knotVector[0],nknot[0],0,ielem_x+ibpe_x,porder,xq[0])
                   *evalN(knotVector[1],nknot[1],0,ielem_y+ibpe_y,porder,xq[1])
                   *evalN(knotVector[2],nknot[2],0,ielem_z+ibpe_z,porder,xq[2]);
            for (int JJ=0; JJ<ndim; JJ++)
            {
                N[ia++]=evalN(knotVector[0],nknot[0],(int)(JJ==0),ielem_x+ibpe_x,porder,xq[0])
                       *evalN(knotVector[1],nknot[1],(int)(JJ==1),ielem_y+ibpe_y,porder,xq[1])
                       *evalN(knotVector[2],nknot[2],(int)(JJ==2),ielem_z+ibpe_z,porder,xq[2]);
            }
            for (int JJ=0; JJ<ndim; JJ++)
            {
                for (int KK=JJ; KK<ndim; KK++)
                {
                    N[ia++]=evalN(knotVector[0],nknot[0],(int)(JJ==0)+(int)(KK==0),ielem_x+ibpe_x,porder,xq[0])
                           *evalN(knotVector[1],nknot[1],(int)(JJ==1)+(int)(KK==1),ielem_y+ibpe_y,porder,xq[1])
                           *evalN(knotVector[2],nknot[2],(int)(JJ==2)+(int)(KK==2),ielem_z+ibpe_z,porder,xq[2]);
                }
            }
        }}}
        // ---- evaluate ui
        for (ia=0; ia<nddim*ndof; ia++) { ui[ia]=0.0; }
        for (int ibpe_x=0; ibpe_x<porder+1; ibpe_x++){
        for (int ibpe_y=0; ibpe_y<porder+1; ibpe_y++){
        for (int ibpe_z=0; ibpe_z<porder+1; ibpe_z++){
            ibpe=(porder+1)*(porder+1)*ibpe_x+(porder+1)*ibpe_y+ibpe_z;
            ibasis=(nbasis_y_active*nbasis_z_active)*(ielem_x+ibpe_x)+nbasis_z_active*(ielem_y+ibpe_y)+(ielem_z+ibpe_z);
            for (int idof=0; idof<ndof; idof++)
            {
                ia=ndof*ibasis+idof;
                for (int iddim=0; iddim<nddim; iddim++)
                {
                    ui[nddim*idof+iddim]+=N[nddim*ibpe+iddim]*Ui[ia];
                }
            }
        }}}
        // ---- evaluate tangent
        eval_tangent(tangent,ui,par_mat);
        // ---- add to val
        weight=wquad[iquad_x]*wquad[iquad_y]*wquad[iquad_z]*(X1[0]-X0[0])*(X1[1]-X0[1])*(X1[2]-X0[2]);
        for (int ibpe_x=0; ibpe_x<porder+1; ibpe_x++){
        for (int ibpe_y=0; ibpe_y<porder+1; ibpe_y++){
        for (int ibpe_z=0; ibpe_z<porder+1; ibpe_z++){
            ibpe=(porder+1)*(porder+1)*ibpe_x+(porder+1)*ibpe_y+ibpe_z;
            for (int idof=0; idof<ndof; idof++)
            {
                row=ndof*ibpe+idof;
                for (int lbpe_x=0; lbpe_x<porder+1; lbpe_x++){
                for (int lbpe_y=0; lbpe_y<porder+1; lbpe_y++){
                for (int lbpe_z=0; lbpe_z<porder+1; lbpe_z++){
                    lbpe=(porder+1)*(porder+1)*lbpe_x+(porder+1)*lbpe_y+lbpe_z;
                    for (int ldof=0; ldof<ndof; ldof++)
                    {
                        col=ndof*lbpe+ldof;
                        temp=0.0;
                        for (int iddim=0; iddim<nddim; iddim++)
                        {
                            for (int lddim=0; lddim<nddim; lddim++)
                            {
                                temp+=N[nddim*ibpe+iddim]*tangent[(ndof*nddim)*(nddim*idof+iddim)+(nddim*ldof+lddim)]*N[nddim*lbpe+lddim];
                            } // lddim
                        } // iddim
                        val[(ndof*nbpe)*row+col]+=weight*temp;
                    } // ldof
                }}} // lbpe
            } // idof
        }}} // ibpe
    }}} // iquad
    // ---- set idx
    for (int ibpe_x=0; ibpe_x<porder+1; ibpe_x++){
    for (int ibpe_y=0; ibpe_y<porder+1; ibpe_y++){
    for (int ibpe_z=0; ibpe_z<porder+1; ibpe_z++){
        ibpe=(porder+1)*(porder+1)*ibpe_x+(porder+1)*ibpe_y+ibpe_z;
        ibasis=(nbasis_y_active*nbasis_z_active)*(ielem_x+ibpe_x)+nbasis_z_active*(ielem_y+ibpe_y)+(ielem_z+ibpe_z);
        for (int idof=0; idof<ndof; idof++)
        {
            idx[ndof*ibpe+idof]=globalIdx[ndof*ibasis+idof];
        }
    }}}
    // ---- add to Pmat_
    ierr = MatSetValues(Pmat_,ndof*nbpe,idx,ndof*nbpe,idx,val,ADD_VALUES);CHKERRQ(ierr);
    }}} // ielem
    // ---- 
    ierr = MatAssemblyBegin(Pmat_,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Pmat_,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    
    /* ================================================================ BOUNDARY CONDITIONS ================================================================ */
    int nbasis_x_,nbasis_y_,nbasis_z_,ibasis_z_0,ibasis_z_1;
    int ax_,ay_,az_;
    /* ----------------------------------------------------------------    Dirichlet BCs    ---------------------------------------------------------------- */
    int nrow_dirichlet=0;
    int itemp;
    for (int iboun=0; iboun<nboun; iboun++)
    {
        if      (boun[iboun]/2==0) { itemp=(2*ndof); for (int i=0; i<(2*ndof); i++) { itemp-=bc_type[boun[iboun]*(2*ndof)+i]; } nrow_dirichlet += itemp*nbasis[1]*nbasis[2]; }
        else if (boun[iboun]/2==1) { itemp=(2*ndof); for (int i=0; i<(2*ndof); i++) { itemp-=bc_type[boun[iboun]*(2*ndof)+i]; } nrow_dirichlet += itemp*nbasis[2]*nbasis[0]; }
        else                       { itemp=(2*ndof); for (int i=0; i<(2*ndof); i++) { itemp-=bc_type[boun[iboun]*(2*ndof)+i]; } nrow_dirichlet += itemp*nbasis[0]*nbasis[1]; }
    }
    int *row_dirichlet = (int*) malloc(nrow_dirichlet*sizeof(int));
    ia=0;
    for (int iboun=0; iboun<nboun; iboun++)
    {
        switch (boun[iboun]/2)
        {
        case 0:
            nbasis_x_=nbasis[1]; ax_=nbasis_z_active;
            nbasis_y_=nbasis[2]; ay_=1;
            nbasis_z_=nbasis[0]; az_=nbasis_y_active*nbasis_z_active;
            break;
        case 1:
            nbasis_x_=nbasis[0]; ax_=nbasis_y_active*nbasis_z_active;
            nbasis_y_=nbasis[2]; ay_=1;
            nbasis_z_=nbasis[1]; az_=nbasis_z_active;
            break;
        case 2:
            nbasis_x_=nbasis[0]; ax_=nbasis_y_active*nbasis_z_active;
            nbasis_y_=nbasis[1]; ay_=nbasis_z_active;
            nbasis_z_=nbasis[2]; az_=1;
            break;
        default:
            printf("error: \n");
            exit(0);
        }
        switch (boun[iboun]%2)
        {
        case 0:
            ibasis_z_0=0;
            ibasis_z_1=1;
            break;
        case 1:
            ibasis_z_0=nbasis_z_-1;
            ibasis_z_1=nbasis_z_-2;
            break;
        default:
            printf("error: \n");
            exit(0);
        }
        // standard
        for (int idof=0; idof<ndof; idof++) {
            if (bc_type[boun[iboun]*(2*ndof)+0*ndof+idof]==0) {
                for (int ibasis_x_=0; ibasis_x_<nbasis_x_; ibasis_x_++) {
                for (int ibasis_y_=0; ibasis_y_<nbasis_y_; ibasis_y_++) {
                    ibasis=ax_*ibasis_x_+ay_*ibasis_y_+az_*ibasis_z_0;
                    row_dirichlet[ia++]=globalIdx[ndof*ibasis+idof];
        }}}}
        // high-order
        for (int idof=0; idof<ndof; idof++) {
            if (bc_type[boun[iboun]*(2*ndof)+1*ndof+idof]==0) {
                for (int ibasis_x_=0; ibasis_x_<nbasis_x_; ibasis_x_++) {
                for (int ibasis_y_=0; ibasis_y_<nbasis_y_; ibasis_y_++) {
                    ibasis=ax_*ibasis_x_+ay_*ibasis_y_+az_*ibasis_z_1;
                    row_dirichlet[ia++]=globalIdx[ndof*ibasis+idof];
        }}}}
    } // for iboun
    assert(ia==nrow_dirichlet);

    ierr = MatZeroRowsColumns(Pmat_,nrow_dirichlet,row_dirichlet,1.0,PETSC_NULL,PETSC_NULL); // symmetrized //CHKERRQ(ierr);
    //ierr = MatZeroRows(Pmat_,nrow_dirichlet,row_dirichlet,1.0,PETSC_NULL,PETSC_NULL); // unsymmetrized //CHKERRQ(ierr);
    free(row_dirichlet);
    
    /* ---------------------------------------------------------------- FINAL ASSEMBLY ---------------------------------------------------------------- */
    if (Tangent_ != Pmat_)
    {
        ierr = MatAssemblyBegin(Tangent_,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(Tangent_,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }

    free(Ui);
    free(tangent);
    return 0;
}

