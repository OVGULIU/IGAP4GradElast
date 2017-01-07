void set_mpi_universalIdx(int *universalIdx, int *nDof_array, int *iDof_displ_array, int nDof, int nselfIdx, int nbasis[], int ielem_displ[], int nelem_global[], int ndof, int porder, int nproc, int rank)
{
    int ierror;
    ierror=MPI_Gather(&nDof,1,MPI_INT,nDof_array,1,MPI_INT,0,MPI_COMM_WORLD);assert(ierror==MPI_SUCCESS);
    if (rank==0) { iDof_displ_array[0]=0; for (int iproc=1; iproc<nproc; iproc++) { iDof_displ_array[iproc]=iDof_displ_array[iproc-1]+nDof_array[iproc-1]; } }
    
    int *universalIdx_self=(int*) malloc(nselfIdx*sizeof(int));
    int ia=0, ibasis;
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++) {
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++) {
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++) {
        ibasis=(nelem_global[1]+porder)*(nelem_global[2]+porder)*(ielem_displ[0]+ibasis_x)+(nelem_global[2]+porder)*(ielem_displ[1]+ibasis_y)+(ielem_displ[2]+ibasis_z);
        for (int idof=0; idof<ndof; idof++) { universalIdx_self[ia++]=ndof*ibasis+idof; }
    }}}
    ierror=MPI_Gatherv(universalIdx_self,nselfIdx,MPI_INT,universalIdx,nDof_array,iDof_displ_array,MPI_INT,0,MPI_COMM_WORLD);assert(ierror==MPI_SUCCESS);
    free(universalIdx_self);
}
