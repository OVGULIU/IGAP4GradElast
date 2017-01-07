void plot_volume(
HMGL gr,
double *U0,
int ndim,
int nddim,
int ndof,
double *par_mat,
int porder,
int nelem[],
int nbasis[],
int nknot[],
double *knotVector[],
int vflag,
double uscale)
{
    double machine_tol=1.e-15;
    
    //mgl_set_num_thr(1);
    // ---- BVP
    double *u  = (double*) malloc(nddim*ndof*sizeof(double));
    //double cmin_temp=1;
    // ---- Mesh
    int ibpe;
    int nbpe = (porder+1)*(porder+1)*(porder+1);                   // number of active basis per element
    int ibasis;                                                    // local basis id
    double X0,Y0,Z0,X1,Y1,Z1;
    // ---- IGA
    double *N = (double*) malloc(nddim*nbpe*sizeof(double));
    int         ia;
    
    //int nbasis_x_active=nelem[0]+porder;
    int nbasis_y_active=nelem[1]+porder;
    int nbasis_z_active=nelem[2]+porder;
    
    // ---- Plot
    double X,Y,Z;
    double c,ux,uy,uz;
    double e1,e2,vwell;

    double *avalPtr;
    double *cvalPtr;

    switch (vflag)
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
        avalPtr=&e1;
        cvalPtr=&c;
        break;
    case 5:
        avalPtr=&e2;
        cvalPtr=&c;
        break;
    case 6:
        avalPtr=&vwell;
        cvalPtr=&c;
        break;
    case 9:
        avalPtr=&e1;
        cvalPtr=&e2;
        break;
    default:
        exit(0);
    }
    int nppe=4;
    HMDT cplot=mgl_create_data(); mgl_data_create(cplot,nppe*nelem[0]+1,nppe*nelem[1]+1,nppe*nelem[2]+1);
    HMDT aplot=mgl_create_data(); mgl_data_create(aplot,nppe*nelem[0]+1,nppe*nelem[1]+1,nppe*nelem[2]+1);
    HMDT xplot=mgl_create_data(); mgl_data_create(xplot,nppe*nelem[0]+1,nppe*nelem[1]+1,nppe*nelem[2]+1);
    HMDT yplot=mgl_create_data(); mgl_data_create(yplot,nppe*nelem[0]+1,nppe*nelem[1]+1,nppe*nelem[2]+1);
    HMDT zplot=mgl_create_data(); mgl_data_create(zplot,nppe*nelem[0]+1,nppe*nelem[1]+1,nppe*nelem[2]+1);
    // ---- make contour data
    for (int ielem_x=0; ielem_x<nelem[0]; ielem_x++) {
    for (int ielem_y=0; ielem_y<nelem[1]; ielem_y++) {
    for (int ielem_z=0; ielem_z<nelem[2]; ielem_z++) {
    X0=knotVector[0][ielem_x+porder]; X1=knotVector[0][ielem_x+porder+1];
    Y0=knotVector[1][ielem_y+porder]; Y1=knotVector[1][ielem_y+porder+1];
    Z0=knotVector[2][ielem_z+porder]; Z1=knotVector[2][ielem_z+porder+1];
    for (int ippe_x=0; ippe_x<nppe+(int)(ielem_x==nelem[0]-1); ippe_x++) {
    for (int ippe_y=0; ippe_y<nppe+(int)(ielem_y==nelem[1]-1); ippe_y++) {
    for (int ippe_z=0; ippe_z<nppe+(int)(ielem_z==nelem[2]-1); ippe_z++) {
        X=X0+(X1-X0)*ippe_x/nppe-machine_tol*(int)(ippe_x==nppe);
        Y=Y0+(Y1-Y0)*ippe_y/nppe-machine_tol*(int)(ippe_y==nppe);
        Z=Z0+(Z1-Z0)*ippe_z/nppe-machine_tol*(int)(ippe_z==nppe);
        // ---- evaluate N
        ia=0;
        for (int ibpe_x=0; ibpe_x<porder+1; ibpe_x++){
        for (int ibpe_y=0; ibpe_y<porder+1; ibpe_y++){
        for (int ibpe_z=0; ibpe_z<porder+1; ibpe_z++){
            N[ia++]=evalN(knotVector[0],nknot[0],0,ielem_x+ibpe_x,porder,X)
                   *evalN(knotVector[1],nknot[1],0,ielem_y+ibpe_y,porder,Y)
                   *evalN(knotVector[2],nknot[2],0,ielem_z+ibpe_z,porder,Z);
            for (int JJ=0; JJ<ndim; JJ++)
            {
                N[ia++]=evalN(knotVector[0],nknot[0],(int)(JJ==0),ielem_x+ibpe_x,porder,X)
                       *evalN(knotVector[1],nknot[1],(int)(JJ==1),ielem_y+ibpe_y,porder,Y)
                       *evalN(knotVector[2],nknot[2],(int)(JJ==2),ielem_z+ibpe_z,porder,Z);
            }
            for (int JJ=0; JJ<ndim; JJ++)
            {
                for (int KK=JJ; KK<ndim; KK++)
                {
                    N[ia++]=evalN(knotVector[0],nknot[0],(int)(JJ==0)+(int)(KK==0),ielem_x+ibpe_x,porder,X)
                           *evalN(knotVector[1],nknot[1],(int)(JJ==1)+(int)(KK==1),ielem_y+ibpe_y,porder,Y)
                           *evalN(knotVector[2],nknot[2],(int)(JJ==2)+(int)(KK==2),ielem_z+ibpe_z,porder,Z);
                }
            }
        }}}
        // ---- evaluate u
        for (ia=0; ia<nddim*ndof; ia++) { u[ia]=0.0; }
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
                    u[nddim*idof+iddim]+=N[nddim*ibpe+iddim]*U0[ia];
                }
            }
        }}}
        // ---- evaluate e's
        eval_e1(&e1,u,par_mat);
        eval_e2(&e2,u,par_mat);
        eval_energy_3well(&vwell,u,par_mat);
        // ----
        if      ((vwell<-0.5) && (e2<0)) {vwell=3.0;}
        else if ((vwell<-0.5) && (e1>0)) {vwell=1.0;}
        else if ((vwell<-0.5) && (e1<0)) {vwell=2.0;}
        else {vwell=0.0;}
        // ----
        c =-999;
        ux=u[nddim*0];
        uy=u[nddim*1];
        uz=u[nddim*2];
        mgl_data_set_value(cplot,*cvalPtr,nppe*ielem_x+ippe_x,nppe*ielem_y+ippe_y,nppe*ielem_z+ippe_z); // *cvalPtr 0 <-> 1 : transparent <-> opaque (if used for transparency)
        mgl_data_set_value(aplot,*avalPtr,nppe*ielem_x+ippe_x,nppe*ielem_y+ippe_y,nppe*ielem_z+ippe_z);
        mgl_data_set_value(xplot,X+uscale*ux,nppe*ielem_x+ippe_x,nppe*ielem_y+ippe_y,nppe*ielem_z+ippe_z);
        mgl_data_set_value(yplot,Y+uscale*uy,nppe*ielem_x+ippe_x,nppe*ielem_y+ippe_y,nppe*ielem_z+ippe_z);
        mgl_data_set_value(zplot,Z+uscale*uz,nppe*ielem_x+ippe_x,nppe*ielem_y+ippe_y,nppe*ielem_z+ippe_z);
    }}}  // ippe
    }}} // ielem

    //
    int    naval;
    double *aval;
    int    space=nelem[0]*nppe;//upto nelem*nppe;
    HMDT e1plot  =mgl_create_data(); mgl_data_create(e1plot,1,1,1);
    HMDT e2plot  =mgl_create_data(); mgl_data_create(e2plot,1,1,1);
    switch (vflag)
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
        naval=nppe*nelem[0]/space+1;
        aval=(double*) malloc(naval*sizeof(double));
        for (int iaval=0; iaval<naval; iaval++) { aval[iaval]=space*iaval; }
        for (int iaval=0; iaval<naval; iaval++) { mgl_dens3_xyz(gr,xplot,yplot,zplot,aplot,"x#",aval[iaval],"meshnum 15"); }
        free(aval);
        naval=nppe*nelem[1]/space+1;
        aval=(double*) malloc(naval*sizeof(double));
        for (int iaval=0; iaval<naval; iaval++) { aval[iaval]=space*iaval; }
        for (int iaval=0; iaval<naval; iaval++) { mgl_dens3_xyz(gr,xplot,yplot,zplot,aplot,"#",aval[iaval],"meshnum 15"); }
        free(aval);
        naval=nppe*nelem[2]/space+1;
        aval=(double*) malloc(naval*sizeof(double));
        for (int iaval=0; iaval<naval; iaval++) { aval[iaval]=space*iaval; }
        for (int iaval=0; iaval<naval; iaval++) { mgl_dens3_xyz(gr,xplot,yplot,zplot,aplot,"z#",aval[iaval],"meshnum 15"); }
        free(aval);
        break;
    case 6:
        naval=nppe*nelem[0]/space+1;
        aval=(double*) malloc(naval*sizeof(double));
        for (int iaval=0; iaval<naval; iaval++) { aval[iaval]=space*iaval; }
        for (int iaval=0; iaval<naval; iaval++) { mgl_dens3_xyz(gr,xplot,yplot,zplot,aplot,"{x202020}{xFFA500}{xFFA500}{x3CB371}{x3CB371}{xD2691E}|x#",aval[iaval],"meshnum 15"); }
        free(aval);
        naval=nppe*nelem[1]/space+1;
        aval=(double*) malloc(naval*sizeof(double));
        for (int iaval=0; iaval<naval; iaval++) { aval[iaval]=space*iaval; }
        for (int iaval=0; iaval<naval; iaval++) { mgl_dens3_xyz(gr,xplot,yplot,zplot,aplot,"{x202020}{xFFA500}{xFFA500}{x3CB371}{x3CB371}{xD2691E}|#",aval[iaval],"meshnum 15"); }
        free(aval);
        naval=nppe*nelem[2]/space+1;
        aval=(double*) malloc(naval*sizeof(double));
        for (int iaval=0; iaval<naval; iaval++) { aval[iaval]=space*iaval; }
        for (int iaval=0; iaval<naval; iaval++) { mgl_dens3_xyz(gr,xplot,yplot,zplot,aplot,"{x202020}{xFFA500}{xFFA500}{x3CB371}{x3CB371}{xD2691E}|z#",aval[iaval],"meshnum 15"); }
        free(aval);
        // x3CB371 or x20B2AA 
        break;
    case 9:
        for (int i=0; i<nppe*nelem[0]+1; i++) {
        for (int j=0; j<nppe*nelem[1]+1; j++) {
        for (int k=0; k<nppe*nelem[2]+1; k++) {
            if ( (i%4==0) && (j%4==0) && (k%4==0) )
            {
                mgl_data_set_value(e1plot,mgl_data_get_value(aplot,i,j,k),0,0,0);
                mgl_data_set_value(e2plot,mgl_data_get_value(cplot,i,j,k),0,0,0);
                mgl_plot_xy(gr,e1plot,e2plot," H.","size 0.05");
            }
        }}}
        break;
    default:
        exit(0);
    }

    // ---- finalize
    free(u);
    free(N);
    //free(aval);
    mgl_delete_data(cplot);
    mgl_delete_data(aplot);
    mgl_delete_data(xplot);
    mgl_delete_data(yplot);
    mgl_delete_data(zplot);
    mgl_delete_data(e1plot);
    mgl_delete_data(e2plot);

}
