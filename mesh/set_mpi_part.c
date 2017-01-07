void set_mpi_part(int npart[], int ipart[], int nelem[], int ielem_displ[], int nknot[], double *knotVector[], int nelem_global[], double *knotVector_global[], int ndim, int porder, int nproc, int rank)
{
    int ierror;

    double nepp_approx;
    nepp_approx=pow((double)(nelem_global[0]*nelem_global[1]*nelem_global[2])/(double)nproc,1./3.);
    for (int i=round((double)nelem_global[0]/nepp_approx); i>0; i--)  { if (nproc%i ==0) { npart[0]=i; break;} } 
    nepp_approx=pow((double)(nelem_global[1]*nelem_global[2])/(double)(nproc/npart[0]),1./2.);
    for (int i=round((double)nelem_global[1]/nepp_approx); i>0; i--)  { if (nproc/npart[0]%i ==0) { npart[1]=i; break;} }
    npart[2]=nproc/npart[0]/npart[1];
    if (rank==0) { printf("#proc = %d, Domain partition: %d x %d x %d.\n",nproc,npart[0],npart[1],npart[2]); }
    
    ipart[0]=rank/(npart[1]*npart[2]);
    ipart[1]=(rank%(npart[1]*npart[2]))/npart[2];
    ipart[2]=rank%npart[2];

    int *nelem_w_array;       if (rank==0) { nelem_w_array=(int*) malloc(nproc*sizeof(int)); }
    int *nknot_w_array;       if (rank==0) { nknot_w_array=(int*) malloc(nproc*sizeof(int)); }
    int *ielem_w_displ_array; if (rank==0) { ielem_w_displ_array=(int*) malloc(nproc*sizeof(int)); }
    int ipart_w;
    int *nelem_w_1d;
    int *ielem_w_displ_1d;
    for (int idim=0; idim<ndim; idim++)
    {
        if (rank==0)
        {
            nelem_w_1d      =(int*) malloc(npart[idim]*sizeof(int));
            ielem_w_displ_1d=(int*) malloc(npart[idim]*sizeof(int));
            for (int i=0; i<npart[idim]; i++)                    { nelem_w_1d[i]=nelem_global[idim]/npart[idim]; }
            for (int i=0; i<nelem_global[idim]%npart[idim]; i++) { nelem_w_1d[i]++; }
            //nelem_w_1d[npart[idim]-1]-=porder;
            //for (int i=0; i<porder; i++)                         { nelem_w_1d[((npart[idim]-1)-i)%npart[idim]]++; }
            ielem_w_displ_1d[0]=0;
            for (int i=1; i<npart[idim]; i++) { ielem_w_displ_1d[i]=ielem_w_displ_1d[i-1]+nelem_w_1d[i-1]; }
            for (int iproc=0; iproc<nproc; iproc++)
            {
                switch (idim) {
                case 0: ipart_w=iproc/(npart[1]*npart[2]); break;
                case 1: ipart_w=(iproc%(npart[1]*npart[2]))/npart[2]; break;
                case 2: ipart_w=iproc%npart[2]; break;
                }
                nelem_w_array[iproc]=nelem_w_1d[ipart_w];
                nknot_w_array[iproc]=nelem_w_1d[ipart_w]+1+2*porder;
                ielem_w_displ_array[iproc]=ielem_w_displ_1d[ipart_w];
            }
            free(nelem_w_1d);
            free(ielem_w_displ_1d);
        }
        ierror=MPI_Scatter(nelem_w_array,1,MPI_INT,nelem+idim,1,MPI_INT,0,MPI_COMM_WORLD);
        ierror=MPI_Scatter(nknot_w_array,1,MPI_INT,nknot+idim,1,MPI_INT,0,MPI_COMM_WORLD);
        ierror=MPI_Scatter(ielem_w_displ_array,1,MPI_INT,ielem_displ+idim,1,MPI_INT,0,MPI_COMM_WORLD);
        knotVector[idim]=(double*) malloc(nknot[idim]*sizeof(double));
        ierror=MPI_Scatterv(knotVector_global[idim],nknot_w_array,ielem_w_displ_array,MPI_DOUBLE,knotVector[idim],nknot[idim],MPI_DOUBLE,0,MPI_COMM_WORLD);assert(ierror==MPI_SUCCESS);
    }
    if (rank==0)
    {
        free(nelem_w_array);
        free(nknot_w_array);
        free(ielem_w_displ_array);
    }
}
