#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "iga.h"

#include "evalN.c"
#include "compute_Residual.c"
#include "compute_Tangent.c"
#include "compute_Energy.c"
#include "appl_dirichlet.c"

// ---- MPI
#include "set_mpi_comm.c"
#include "set_mpi_universalIdx.c"
#include "set_mpi_globalIdx.c"
#include "set_nz.c"
// ---- QUADRATURE
#include "quadrature.c"
