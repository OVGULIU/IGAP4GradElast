#if !defined (MATHUTIL)
#define MATHUTIL

inline int factorial(int n);
inline int int_pow(int base, int pow);
inline double double_abs(double a);

void matinv_dot_vec(double *mat, double *vec, int n);
void matinv(double *mat, int n);

#endif
