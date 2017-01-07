#if !defined (MAINDOTH)
#define MAINDOTH

/* ---------------- declaration (structure) ---------------- */
//typedef int Size_d;

typedef struct{
    int ndim;
    int nddim;
    // ---- MESH
    int porder;
    int *nelem;
    int *nbasis;
    int *nknot;
    double **knotVector;
    // ---- PDE
    int ndof;
    double *par_mat;
    int *bc_type;
    double *par_dirichlet;
    double *par_neumann;
    int nboun;
    int *boun;
    // ---- IGA
    //Vec U0;
    int nquad;
    double *cquad;
    double *wquad;
    int nselfIdx;
    int *selfIdx;
    int nsendPart;
    int *sendPart;
    int nsendIdx;
    int *sendIdx;
    int *sendPtr;
    int nrecvPart;
    int *recvPart;
    int nrecvIdx;
    int *recvIdx;
    int *recvPtr;
    int *globalIdx;
} App;

/* ---------------- declarations ---------------- */
#include "util.h"
#include "mathutil.h"
#include "mesh/mesh.h"
#include "bvp/bvp.h"
#include "iga/iga_plot.h"
#include "plot/plot.h"

/* ---------------- definitions ---------------- */
#include "util.c"
#include "mathutil.c"
#include "mesh/mesh.c"
#include "bvp/bvp.c"
#include "iga/iga_plot.c"
#include "plot/plot.c"

#endif
