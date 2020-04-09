
#include "includes.h"
void WilsonDirac( gsl_vector_complex * output,gsl_vector_complex * input, void * ctxt)
{
  lattice * L = ctxt;
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
