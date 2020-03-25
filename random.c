#include "includes.h"

unsigned long time_seed()
{
  long accum = 0;
  clock_t time1 = clock();
  struct tm * unit;
  unit = localtime(&time1);

  accum =31* unit->tm_mon + unit->tm_mday;
  accum = 24  * accum  + unit->tm_hour;
  accum = 60 * accum + unit->tm_min;
  accum = 60 * accum + unit->tm_sec;
  return accum;
}

//Generate random double from (0, a]
double gen_rand(double a, gsl_rng * r)
{
  return a * gsl_rng_uniform(r);
}

void gen_su2_matrix(gsl_matrix_complex *B, gsl_rng * r)
{
  double a0,a1,a2,a3,a;
  gsl_complex a11, a12, a21, a22;
  gsl_matrix_complex *A = gsl_matrix_complex_alloc(2,2);

  a0 =  gen_rand(2.0,r) - 1.0;
  a1 =  gen_rand(2.0,r) - 1.0;
  a2 =  gen_rand(2.0,r) - 1.0;
  a3 =  gen_rand(2.0,r) - 1.0;

  a = 1.0/(sqrt(a0*a0 + a1*a1 + a2*a2 + a3*a3));
  a0*=a;
  a1*=a;
  a2*=a;
  a3*=a;

  a11 = gsl_complex_rect( a0, a3);
  a22 = gsl_complex_rect( a0,-a3);
  a12 = gsl_complex_rect( a2, a1);
  a21 = gsl_complex_rect( -a2,a1);

  //Insert those values in A
  gsl_matrix_complex_set(A,0,0,a11);
  gsl_matrix_complex_set(A,0,1,a12);
  gsl_matrix_complex_set(A,1,0,a21);
  gsl_matrix_complex_set(A,1,1,a22);

  gsl_matrix_complex_memcpy(B,A);
  gsl_matrix_complex_free(A);
}

gsl_matrix_complex * gen_rotation_matrix(gsl_rng * r,double step)
{
  gsl_matrix_complex * temp = gsl_matrix_complex_alloc(2,2);
  gsl_complex a11, a12, a21, a22;
  double a0, a1, a2, a3;
  double a,b;
  double epsilon;

  epsilon = gen_rand(1.0, r) * step;
  a0 = sqrt(1.0 - epsilon*epsilon);
  a1 =  gen_rand(2.0,r) - 1.0;
  a2 =  gen_rand(2.0,r) - 1.0;
  a3 =  gen_rand(2.0,r) - 1.0;
  b = a1*a1 + a2*a2 + a3*a3;
  a = epsilon*sqrt(1.0/b);

  a11 = gsl_complex_rect(a0,a*a1);
  a12 = gsl_complex_rect(-a*a2,a*a3);
  a21 = gsl_complex_rect(a*a2,a*a3);
  a22 = gsl_complex_rect(a0,-a*a1);

  gsl_matrix_complex_set(temp,0,0,a11);
  gsl_matrix_complex_set(temp,0,1,a12);
  gsl_matrix_complex_set(temp,1,0,a21);
  gsl_matrix_complex_set(temp,1,1,a22);

  return temp;

}
//This function rotates an input matrix R. It returns the original matrix in B
void rand_rotation(gsl_rng * r, gsl_matrix_complex * R, gsl_matrix_complex * B,double step)
{

  gsl_complex alpha = gsl_complex_rect(1.0,0.0);
  gsl_complex beta = gsl_complex_rect(0.0,0.0);
  gsl_matrix_complex * A = gen_rotation_matrix(r,step);
  gsl_matrix_complex * C = gsl_matrix_complex_alloc(2,2);

  gsl_matrix_complex_memcpy(B,R);
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,  alpha, R, A, beta, C);
  gsl_matrix_complex_memcpy(R,C);
  gsl_matrix_complex_free(A);
  gsl_matrix_complex_free(C);
}
