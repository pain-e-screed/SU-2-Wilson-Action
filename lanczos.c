
#include "includes.h"
void WilsonDirac( gsl_vector_complex * output,gsl_vector_complex * input, void * ctxt)
{
  lattice * L = (lattice * ) ctxt;
  gsl_matrix_complex * gamma_matrix[4];
  gsl_matrix_complex * identity_gamma_minus[4];
  gsl_matrix_complex * identity_gamma_plus[4];
  gsl_matrix_complex * temp_link_conj = gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex * temp_mat1 = gsl_matrix_complex_alloc(8,8);
  gsl_matrix_complex * temp_mat2 = gsl_matrix_complex_alloc(8,8);
  gsl_vector_complex * temp_vec1 = gsl_vector_complex_alloc(8);
  gsl_vector_complex * temp_vec2 = gsl_vector_complex_alloc(8);

  gsl_vector_complex_view output_vector_view;
  gsl_vector_complex_view input_vector_view1;
  gsl_vector_complex_view input_vector_view2;


  long unsigned site_index, vector_site_index, temp_site_index, temp_vector_site_index;
  unsigned int dir;
  //Begin testing input vectors and matrices
  if(!gamma_matrix[0])
    compute_gamma_matrices();
  if(!identity_gamma_plus[0] || !identity_gamma_minus[0])
    compute_gamma_identities();
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
  //End vec/matrix testing
  //If changing the value of m in the WD operator, this is where you'll do it
  gsl_vector_complex_memcpy(output,input)
  gsl_vector_complex_set_all(output,gsl_complex_rect(5.0,0.0));
  FORALLSITES(site_index)
  {
    vector_site_index = 8 * site_index;
    output_vector_view = gsl_vector_complex_subvector(output, vector_site_index, 8 );
    gsl_vector_complex_set_zero(temp_vec1);
    FORALLDIR(dir)
    {
      temp_site_index = hop_index(site_index,1, dir, BACKWARD);
      input_vector_view1 = gsl_vector_complex_subvector(input, temp_site_index, 8 );
      outer_product(temp_mat1,identity_gamma_minus[dir],L->R[site_index]->link[dir]);

      temp_site_index = hop_index(site_index,1, dir, FORWARD);
      input_vector_view2 = gsl_vector_complex_subvector(input, temp_site_index, 8 );
      //Add something to do the Hermitian conjugate to this link
      hermitian_conj(temp_link_conj,L->R[site_index]->link[dir])
      outer_product(temp_mat2,identity_gamma_plus[dir],temp_link_conj);

      gsl_blas_zgemv(CblasNoTrans, gsl_complex_rect(-0.5,0.0), temp_mat1, &input_vector_view1.vector, gsl_complex_rect(1.0,0.0), &output_vector_view.vector);
      gsl_blas_zgemv(CblasNoTrans, gsl_complex_rect(-0.5,0.0), temp_mat2, &input_vector_view2.vector, gsl_complex_rect(1.0,0.0), &output_vector_view.vector);
    }

      gsl_matrix_complex_free(temp_link_conj);
      gsl_matrix_complex_free(temp_mat1);
      gsl_matrix_complex_free(temp_mat2);
      gsl_matrix_complex_free(temp_vec1);
      gsl_matrix_complex_free(temp_vec2);
  }




hermitian_conj(gsl_matrix_complex *out, gsl_matrix_complex *in)
{
  gsl_matrix_complex *temp = gsl_matrix_complex_alloc(in->size1,in->size2);
  int i, j;
  for(i=0;i<in->size1;i++)
  {
    for (j = 0; j < in->size2; j++)
    {
      gsl_matrix_complex_set(temp,i,j,gsl_complex_conjugate(gsl_matrix_complex_get(in,j,i))   );
    }
  }
  gsl_matrix_complex_memcpy(out,temp);
  gsl_matrix_complex_free(temp);
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
  T = constructT(alpha, beta,k);
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
    appendT(T, alpha, beta,j);
    QRalgorithm(S,eigenvalues, T,j);
    convergenceTest(S,beta);
  }
  //11. End loop
}


gsl_matrix * constructT( alpha, beta,k)
{
  gsl_matrix * temp = gsl_matrix_calloc(k,k);
  gsl_matrix_complex_set(temp,0,0,alpha);
  gsl_matrix_complex_set(temp,0,1,beta);
  gsl_matrix_complex_set(temp,1,0,beta);
  return temp;
}


vectorSub(gsl_vector_complex * v,gsl_vector_complex * u, gsl_complex c)
{
  
}


grahamSchmidt(r,Q)
appendT(T, alpha, beta,j)
QRalgorithm(S,eigenvalues, T,j)
convergenceTest(S,beta)
