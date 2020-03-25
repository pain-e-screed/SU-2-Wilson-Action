#ifndef _RANDOM_
#define _RANDOM_


unsigned long time_seed();
double gen_rand(double a, gsl_rng * r);
void gen_su2_matrix(gsl_matrix_complex *B, gsl_rng * r);
gsl_matrix_complex * gen_rotation_matrix(gsl_rng * r,double step);
void rand_rotation(gsl_rng * r, gsl_matrix_complex * R, gsl_matrix_complex * B,double step);


#endif
