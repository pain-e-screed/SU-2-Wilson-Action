#include "includes.h"

extern gsl_matrix_complex * gamma_matrix[4];
extern gsl_matrix_complex * identity_gamma_minus[4];
extern gsl_matrix_complex * identity_gamma_plus[4];

void WilsonDirac(gsl_vector_complex * input, gsl_vector_complex * output, void * ctxt);
outer_product(gsl_matrix_complex * w,gsl_matrix_complex * u,gsl_matrix_complex * v);
void Lanczos( void (*MV) (gsl_vector_complex *,gsl_vector_complex *, void * ctxt) ,const int k  , void * data  );