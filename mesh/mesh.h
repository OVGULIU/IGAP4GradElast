#if !defined (MESH)
#define MESH

void mesh_uniform(int nelem_global[], double *knotVector_global[], int mref, int porder, int ndim);

void insert_knot_uniform(double *U_universal, FILE *filePtr, const int nelem[], int nref, int porder, int ndim, int ndof);


// ---- MPI
void set_mpi_part(int npart[], int ipart[], int nelem[], int ielem_displ[], int nknot[], double *knotVector[], int nelem_global[], double *knotVector_global[], int ndim, int porder, int nproc, int rank);

#endif
