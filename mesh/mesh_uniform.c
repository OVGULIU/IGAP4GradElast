void mesh_uniform(int nelem_global[], double *knotVector_global[], int mref, int porder, int ndim)
{

    double L[]={1.0,1.0,1.0};

    for (int idim=0; idim<ndim; idim++) { nelem_global[idim]=int_pow(2,mref); }
    for (int idim=0; idim<ndim; idim++) { knotVector_global[idim]=(double*) malloc((nelem_global[idim]+1+2*porder)*sizeof(double)); }

    double l[ndim]; for (int idim=0; idim<ndim; idim++) { l[idim]=L[idim]/nelem_global[idim]; }
    int iknot;
    for (int idim=0; idim<ndim; idim++)
    {
        iknot=0;
        while (iknot < porder+1)                      { knotVector_global[idim][iknot++] = 0.; }
        while (iknot < nelem_global[idim]+porder)     { knotVector_global[idim][iknot]   = (iknot-porder)*l[idim]+0*sin(iknot)*0.25*l[idim]; iknot++; }
        while (iknot < nelem_global[idim]+1+2*porder) { knotVector_global[idim][iknot++] = L[idim]; }
    }
}
