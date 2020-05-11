#include "includes.h"

extern gsl_matrix_complex * gamma_matrix[4];
extern gsl_matrix_complex * identity_gamma_minus[4];
extern gsl_matrix_complex * identity_gamma_plus[4];

void WilsonDirac(gsl_vector_complex * input, gsl_vector_complex * output, void * ctxt);
hermitian_conj(gsl_matrix_complex *out, gsl_matrix_complex *in);
outer_product(gsl_matrix_complex * w,gsl_matrix_complex * u,gsl_matrix_complex * v);
void compute_gamma_matrices();
void compute_gamma_identities();
void Lanczos( void (*MV) (gsl_vector_complex *,gsl_vector_complex *, void * ctxt),gsl_vector_complex * q_0 , const int k, void * data    );
void vectorScaleSub(gsl_vector_complex * v,gsl_vector_complex * u, gsl_complex c);
void grahamSchmidtInsertion(gsl_vector_complex * r,gsl_matrix_complex * Q, int k);
