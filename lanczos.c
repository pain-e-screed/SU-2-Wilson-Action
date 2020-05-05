
#include "includes.h"
void WilsonDirac( gsl_vector_complex * output,gsl_vector_complex * input, void * ctxt)
{
  lattice * L = (lattice * ) ctxt;
  gsl_matrix_complex * gamma_matrix[4];
  gsl_vector_complex_view output_vector_view;

  long unsigned site_index, vector_site_index, temp_site_index, temp_vector_site_index;
  //spin_index should either be a 0 or 1 for down/up respectively
  //color_index runs from 0 to 3, indicating t to z
  unsigned int spin_index, color_index, sub_index;

  if(!gamma_matrix[0])
    compute_gamma_matrices(&gamma_matrix[0]);

  if(!output)
    output = gsl_vector_complex_alloc(8*N4);

  if(!input)
  {
    printf("Input vector has not been allocated.\n");
    exit(0);
  }

  if(input->size != output->size )
  {
    printf("Input vector and output vector have mismatched lengths\n");
    exit(0);
  }

  gsl_vector_complex_set_all(output,gsl_complex_rect(5.0,0.0));
  FORALLSITES(site_index)
  {
    vector_site_index = 8 * site_index;
    output_vector_view = gsl_vector_complex_subvector(output, vector_site_index, 8 );
    for(spin_index=0;spin_index<2;spin_index++)
    {
      FORALLDIR(color_index)
      {
        sub_index = 4 * spin_index + color_index;

      }
    }

  }








}
//w = u cross v
outer_product(gsl_matrix_complex * w,gsl_matrix_complex * u,gsl_matrix_complex * v)
{
  int i,j;
  int xshift, yshift;
  gsl_matrix_complex_view temp_matrix_view;
  gsl_matrix_complex temp = gsl_matrix_complex_alloc(v->size1,v->size2);
  if(!w)
  {
    printf("Return matrix has not been properly allocated\n");
    exit(0);
  }
  if( !(u->size1 * v->size1 == w->size1 && u->size2 * v->size2 == w->size2))
  {
    printf("Return matrix does not have proper dimensions to compute outer product\n");
    exit(0);
  }
  xshift = v->size1;
  yshift = v->size2;

  for(i=0;i<u->size1;i++)
  {
    for(j=0;j<u->size2;j++)
    {
      temp_matrix_view = gsl_matrix_complex_submatrix(w,i*xshift,j*yshift,xshift,yshift);
      gsl_matrix_complex_memcpy(temp, v);
      gsl_matrix_complex_scale(temp, gsl_matrix_complex_get(u,i,j));
      gsl_matrix_complex_memcpy(&temp_matrix_view->matrix, temp);
    }
  }
}

compute_gamma_matrices(gsl_matrix_complex ** gamma)
{
  int i;
  FORALLDIR(i)
    gamma[i] = gsl_matrix_complex_calloc(4,4);
}

void Lanczos( void (*MV) (gsl_vector_complex *,gsl_vector_complex *, void * ctxt) , const int k, void * data    )
{
  //The ctxt argument is a generic pointer so I can re-use Lanczos later for other things
  int j;
  //1. Choose initial vector r = q_0     beta_0 = ||q_0||
  constructT(T, alpha, beta);
  //2. Begin Loop
  for(j=1;j<k;j++)
  {
    //3. Normalize vector q_j = r /beta_{j-1}
    //4. Compute r = Aq_j - q_{j-1} beta_{j-1}
    MV(r,q,data);
    vectorSub(r,qp,betap);
    //5. alpha_j = q^T_j  r
    alpha = innerProduct(q,r);
    //6. r = r - q_j alpha _j
    vectorSub(r,q,alpha);
    //7. Orthogonalize r against Q
    grahamSchmidt(r,Q);
    //8. beta_j = ||r||
    beta = Magnitude(r);
    //9. Compute eigenvalues of T_j and Test for convergence
    appendT(T, alpha, beta);
    QRalgorithm(S,eigenvalues, T,j);
    convergenceTest(S,beta);
  }
  //11. End loop

}
