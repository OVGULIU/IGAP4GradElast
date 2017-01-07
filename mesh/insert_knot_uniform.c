void insert_knot_uniform(double *U_universal, FILE *filePtr, const int nelem[], int nref, int porder, int ndim, int ndof)
{    
    int lSize;

    double vx[4],vy[4],vz[4];
    double nelem_temp[ndim];
    int nbasis_temp[ndim];
    int nbasis_next[ndim];
    int nDof_base;
    double *U_universal_base,*Ui,*Ui_temp;

    if (nref==0)
    {
        lSize=fread(U_universal,sizeof(double),ndof*(nelem[0]+porder)*(nelem[1]+porder)*(nelem[2]+porder),filePtr);
        assert(lSize==ndof*(nelem[0]+porder)*(nelem[1]+porder)*(nelem[2]+porder));
    }
    else if (nref>0)
    {
        assert(porder==2); // currently only supports porder==2
        // load base solution
        nDof_base=ndof*(nelem[0]/int_pow(2,nref)+porder)*(nelem[1]/int_pow(2,nref)+porder)*(nelem[2]/int_pow(2,nref)+porder);
        U_universal_base =(double*) malloc(nDof_base*sizeof(double));
        lSize=fread(U_universal_base,sizeof(double),nDof_base,filePtr); assert(lSize==nDof_base);
        //
        Ui=(double*) malloc((nelem[0]+porder)*(nelem[1]+porder)*(nelem[2]+porder)*sizeof(double));
        for (int idof=0; idof<ndof; idof++)
        {
            for (int idim=0; idim<ndim; idim++) { nelem_temp[idim]=nelem[idim]/int_pow(2,nref); }
            Ui_temp=(double*) malloc((nelem_temp[0]+porder)*(nelem_temp[1]+porder)*(nelem_temp[2]+porder)*sizeof(double));
            for (int ibasis=0; ibasis<(nelem_temp[0]+porder)*(nelem_temp[1]+porder)*(nelem_temp[2]+porder); ibasis++) { Ui_temp[ibasis]=U_universal_base[ndof*ibasis+idof]; }
            for (int iref=0; iref<nref; iref++)
            {
                for (int idim=0; idim<ndim; idim++) { nbasis_temp[idim]=  nelem_temp[idim]+porder; }
                for (int idim=0; idim<ndim; idim++) { nbasis_next[idim]=2*nelem_temp[idim]+porder; }
                for (int i=0; i<nbasis_next[0]*nbasis_next[1]*nbasis_next[2]; i++) { Ui[i]=0; }

                for (int ix=0; ix<nbasis_temp[0]; ix++) {
                for (int iy=0; iy<nbasis_temp[1]; iy++) {
                for (int iz=0; iz<nbasis_temp[2]; iz++) {
                if (ix==0)                     { vx[0]=0    ; vx[1]=0    ; vx[2]=1    ; vx[3]=1./2.; }
                else if (ix==1)                { vx[0]=0    ; vx[1]=1./2.; vx[2]=3./4.; vx[3]=1./4.; }
                else if (ix<nelem_temp[0])     { vx[0]=1./4.; vx[1]=3./4.; vx[2]=3./4.; vx[3]=1./4.; }
                else if (ix==nelem_temp[0])    { vx[0]=1./4.; vx[1]=3./4.; vx[2]=1./2.; vx[3]=0    ; }
                else if (ix==nelem_temp[0]+1)  { vx[0]=1./2.; vx[1]=1    ; vx[2]=0    ; vx[3]=0    ; }
                if (iy==0)                     { vy[0]=0    ; vy[1]=0    ; vy[2]=1    ; vy[3]=1./2.; }
                else if (iy==1)                { vy[0]=0    ; vy[1]=1./2.; vy[2]=3./4.; vy[3]=1./4.; }
                else if (iy<nelem_temp[1])     { vy[0]=1./4.; vy[1]=3./4.; vy[2]=3./4.; vy[3]=1./4.; }
                else if (iy==nelem_temp[1])    { vy[0]=1./4.; vy[1]=3./4.; vy[2]=1./2.; vy[3]=0    ; }
                else if (iy==nelem_temp[1]+1)  { vy[0]=1./2.; vy[1]=1    ; vy[2]=0    ; vy[3]=0    ; }
                if (iz==0)                     { vz[0]=0    ; vz[1]=0    ; vz[2]=1    ; vz[3]=1./2.; }
                else if (iz==1)                { vz[0]=0    ; vz[1]=1./2.; vz[2]=3./4.; vz[3]=1./4.; }
                else if (iz<nelem_temp[2])     { vz[0]=1./4.; vz[1]=3./4.; vz[2]=3./4.; vz[3]=1./4.; }
                else if (iz==nelem_temp[2])    { vz[0]=1./4.; vz[1]=3./4.; vz[2]=1./2.; vz[3]=0    ; }
                else if (iz==nelem_temp[2]+1)  { vz[0]=1./2.; vz[1]=1    ; vz[2]=0    ; vz[3]=0    ; }
                for (int jx=0; jx<4; jx++) {
                for (int jy=0; jy<4; jy++) {
                for (int jz=0; jz<4; jz++) {
                if (vx[jx]*vy[jy]*vz[jz] != 0) { Ui[(nbasis_next[1]*nbasis_next[2])*(2*ix-2+jx)+nbasis_next[2]*(2*iy-2+jy)+(2*iz-2+jz)]+=Ui_temp[(nbasis_temp[1]*nbasis_temp[2])*ix+nbasis_temp[2]*iy+iz]*vx[jx]*vy[jy]*vz[jz]; }
                }}}
                }}}
                free(Ui_temp);
                Ui_temp=(double*) malloc(nbasis_next[0]*nbasis_next[1]*nbasis_next[2]*sizeof(double));
                for (int i=0; i<nbasis_next[0]*nbasis_next[1]*nbasis_next[2]; i++) { Ui_temp[i]=Ui[i]; }
                for (int idim=0; idim<ndim; idim++) { nelem_temp[idim]*=2; }
            }
            free(Ui_temp);
            for (int ibasis=0; ibasis<(nelem[0]+porder)*(nelem[1]+porder)*(nelem[2]+porder); ibasis++) { U_universal[ndof*ibasis+idof]=Ui[ibasis]; }
        } // for ndof
        free(Ui);
        free(U_universal_base);
    }
}
