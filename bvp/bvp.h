#if !defined (BVP)
#define BVP

void set_ndim(int *ndim);
void set_nddim(int *nddim);
void set_ndof(int *ndof);

void eval_residual(double *residual, double *u, double *par);
void eval_tangent(double *tangent, double *u, double *par);

void eval_neumann(double *traction, double par_neumann[], int face_id, int bc_order, int idof, double xq[]);

void eval_energy(double *energy, double *u, double *par);
void eval_e1(double *e1, double *u, double *par);
void eval_e2(double *e2, double *u, double *par);

void assert_par_mat(double *par);


void eval_energy_3well(double *vwell, double *u, double *par);

#endif
