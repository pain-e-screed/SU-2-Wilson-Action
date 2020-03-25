#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_cblas.h>
#include<gsl/gsl_blas.h>
#include "includes.h"


#ifndef _WRAPPERS_
#define _WRAPPERS_

gsl_matrix_complex * makeIdentity(int n);
gsl_vector * getEigenvalues(gsl_matrix_complex * A);
void print_matrix_complex(gsl_matrix_complex * A);
gsl_complex  matrix_complex_trace(gsl_matrix_complex * A);
int normalize(gsl_vector_complex * n,gsl_vector_complex * u);
gsl_complex inner_product(gsl_vector_complex * u, gsl_vector_complex * v);
int gramSchmidtStep( gsl_matrix_complex * Q, int m);

#endif
