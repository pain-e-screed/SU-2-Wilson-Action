#include "includes.h"

gsl_matrix_complex * makeIdentity(int n)
{
    gsl_matrix_complex * I = gsl_matrix_complex_calloc(n,n);
    gsl_complex unit1  ;
    unit1 = gsl_complex_rect(1.0,0.0);
    int i;
    for(i=0;i<n;i++)
      gsl_matrix_complex_set(I,i,i, unit1 );

    return I;
}


gsl_vector * getEigenvalues(gsl_matrix_complex * A)
{
  gsl_eigen_genherm_workspace * work = gsl_eigen_genherm_alloc(A->size1);
  gsl_matrix_complex * I = makeIdentity(A->size1);
  gsl_matrix_complex * temp = gsl_matrix_complex_calloc(A->size1,A->size2);
  gsl_vector * eval = gsl_vector_calloc(A->size1);
  gsl_matrix_complex_memcpy(temp, A);
  gsl_eigen_genherm(A,I,eval,work);
  printf("The eigenvalues are:\n");
  for(int i=0;i<A->size1;i++)
    printf("%g\n",gsl_vector_get(eval,i));
  gsl_eigen_genherm_free(work);
  gsl_matrix_complex_free(I);
  gsl_matrix_complex_free(temp);
  return eval;
}


void print_matrix_complex(gsl_matrix_complex * A)
{
  int i,j;
  int m,n;
  m = A->size1;
  n = A->size2;
  gsl_complex temp;

  for(i = 0;i<m;i++)
  {
    for(j=0;j<n;j++)
    {
      temp = gsl_matrix_complex_get(A,i,j);
      printf("%g+i%g      ", GSL_REAL(temp), GSL_IMAG(temp)   );
    }
    printf("\n");
  }
  printf("\n");
}

gsl_complex  matrix_complex_trace(gsl_matrix_complex * A)
{
  int i;
  gsl_complex temp;
  gsl_complex total = gsl_complex_rect(0.0,0.0);
  for(i=0;i<A->size1;i++)
  {
    temp = gsl_matrix_complex_get(A,i,i);
    total = gsl_complex_add(total,temp);
  }
  return total;
}


int normalize(gsl_vector_complex * n,gsl_vector_complex * u)
{
  double temp;
  gsl_vector_complex * a = gsl_vector_complex_calloc(u->size);
  gsl_vector_complex_memcpy(a,u);
  temp = gsl_blas_dznrm2(u);
  temp = sqrt(temp);
  if( gsl_vector_complex_scale( a, gsl_complex_rect(1.0/temp,0.0) ) == GSL_SUCCESS)
  {
    gsl_vector_complex_memcpy(n,a);
    gsl_vector_complex_free(a);
  }
  return 0;
}

gsl_complex inner_product(gsl_vector_complex * u, gsl_vector_complex * v)
{
  gsl_complex  temp;
  if(GSL_SUCCESS == gsl_blas_zdotc(u,v,&temp))
    return temp;
}


//orthogonalizes the mth column vectors against the m-1 vectors before it
//it then normalizes the mth column as well
int gramSchmidtStep(gsl_matrix_complex * Q, int m)
{
  int i;
  //If it's the first vector in the series just normalize it
  if(m==1)
  {
    gsl_vector_complex * v = gsl_vector_complex_alloc(Q->size1);
    gsl_vector_complex_view temp = gsl_matrix_complex_column(Q, m-1);
    normalize(v,&temp.vector);
    gsl_matrix_complex_set_col(Q,m-1,v);
    gsl_vector_complex_free(v);
    return m;
  }
  else
  {
    gsl_vector_complex* q = gsl_vector_complex_calloc(Q->size1);
    for( i = 0; i<m-1;i++)
    {
      gsl_vector_complex_view temp = gsl_matrix_complex_column(Q, i);
      if(gsl_vector_complex_sub(q, &temp.vector)== GSL_SUCCESS);
    }
    normalize(q,q);
    gsl_matrix_complex_set_col(Q,m-1,q);
    gsl_vector_complex_free(q);
    return m;
  }
}
