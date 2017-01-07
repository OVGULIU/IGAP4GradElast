#if !defined (PLOT)
#define PLOT

#include <mgl2/mgl_cf.h>
#include <mgl2/mpi.h>

//void plot_N(double *kV, int nk, int n, int p, double Lend);

void plot_solution(HMGL gr, double *U_self, void *app_, int plotID, int nproc, int rank);
void plot_volume(HMGL gr,double *U_,int ndim,int nddim,int ndof,double *par,int porder,int nelem[],int nbasis[],int nknot[],double *knotVector[],int vflag,double uscale);

#endif
