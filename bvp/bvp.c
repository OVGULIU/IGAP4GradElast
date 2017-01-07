#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "bvp.h"

#include "set_ndim.c"
#include "set_nddim.c"
#include "set_ndof.c"

#include "eval_neumann.c"
#include "eval_residual.c"
#include "eval_tangent.c"
#include "eval_energy.c"
#include "eval_e1.c"
#include "eval_e2.c"

#include "assert_par_mat.c"

#include "eval_energy_3well.c"

