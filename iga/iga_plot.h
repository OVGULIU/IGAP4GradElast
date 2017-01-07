#if !defined (IGA)
#define IGA

//PetscErrorCode compute_Residual(SNES snes, Vec Ui_, Vec Residual_, void *app_);
//PetscErrorCode compute_Tangent(SNES snes, Vec Ui_, Mat Tangent_, Mat Pmat_, void *app_);

double evalN(double *kV, int nk, int k, int i, int p, double xi);
void compute_Energy(double *Energy, double *U0_self, void *app_);
void appl_dirichlet(double *dirichlet, double par[], int face_id, int bc_order, int idof, double *knotVector_[], int nknot_[], int nbasis_[], int ndim, int porder);

// ---- MPI
void set_mpi_comm(int *nselfIdx_,int **selfIdx_,int *nsendPart_,int **sendPart_,int *nsendIdx_,int **sendIdx_,int **sendPtr_,int *nrecvPart_,int **recvPart_,int *nrecvIdx_,int **recvIdx_,int **recvPtr_,int ndof,int porder,int rank,int npart[],int ipart[],int nelem[],int nbasis[]);
void set_mpi_universalIdx(int *universalIdx, int *nDof_array, int *iDof_displ_array, int nDof, int nselfIdx, int nbasis[], int ielem_displ[], int nelem_global[], int ndof, int porder, int nproc, int rank);
void set_mpi_globalIdx(int *globalIdx, int *iDof_displ_array, int nselfIdx, int *selfIdx, int nsendPart, int *sendPart, int nsendIdx, int *sendIdx, int *sendPtr, int nrecvPart, int *recvPart, int nrecvIdx, int *recvIdx, int *recvPtr, int nbasis[], int ndof, int nproc, int rank);
void set_nz(int *d_nz, int *o_nz, int nbasis[], int npart[], int ipart[], int ndof, int porder);
// ---- QUADRATURE
#include "quadrature.h"

#endif
