#include "includes.h"

int AnnealSite(lattice * L,gsl_rng * r, long unsigned index, int dir, double beta, double epsilon )
{
      double Si, Sf, dS;
      gsl_matrix_complex * B = gsl_matrix_complex_alloc(2,2);
      Si = partial_action(L,index,dir);
      rand_rotation(r, L->R[index]->link[dir], B,epsilon);
      Sf = partial_action(L,index,dir);
      dS = Sf - Si;
      if(gen_rand(1.0, r) < fmin(1.0,exp(-dS*beta) ))
      {
        gsl_matrix_complex_free(B);
        return 1;
      }
      else
      {
        gsl_matrix_complex_memcpy(L->R[index]->link[dir],B);
        gsl_matrix_complex_free(B);
        return 0;
      }
}

double AnnealStepLattice(lattice * L, gsl_rng * r, double beta, double epsilon)
{
  long unsigned i;
  long unsigned nReject,nAccept;
  int j;
  nReject = 0;
  nAccept = 0;

  FORALLSITES(i)
  {
    FORALLDIR(j)
    {
      if(AnnealSite(L,r, i, j,beta, epsilon ))
        nAccept++;
      else
        nReject++;
    }
  }
  return ((double) nAccept )/ ( (double) (nAccept + nReject));
}

void AnnealingSchedule(lattice * L, gsl_rng *r )
{
  double beta;
  double step;
  double AccRej;
  const int iter = 10;
  double B[iter], X[iter],steps[iter];
  int i, j;

  beta = 0.01;
  steps[0] = 0.5;
  for(i=1;i<iter;i++)
    steps[i] = steps[i-1] * 0.8;

  for(j=0;j<iter;j++)
  {
    for(int i=0;i<1000;i++)
    {
      AccRej = AnnealStepLattice( L, r,beta,steps[j]);
      if(i%300 == 1)
      {
        printf("Acceptance = %lf%%\n",100.0*AccRej );
      }
    }
    B[j] = beta;
    X[j] =  CreutzRatio(L,r,1,1);
    beta += 0.05;
    //step *= 0.9;
    printf("beta = %lf, X(1,1) = %lf\n",B[j],X[j]);
    printf("ln(g_0^2):%lf\n",log(4.0/B[j]) );
    printf("g_0^-2:%lf\n",B[j]/4.0 );
    printf("Action:%lf\n\n\n",action(L,B[j]));
  }
}
//Test whether or not the real part of the trace of the link fields is greater than 1 or less than -1
void TestPlaquettes(lattice * L)
{
  long unsigned index;
  int dir1,dir2;
  gsl_matrix_complex * temp = gsl_matrix_complex_alloc(2,2);
  double retr;

  FORALLSITES(index)
    FORALLDIR(dir1)
      FORALLDIRBUT(dir2,dir1)
      {
        plaquette(temp,L,index, dir1, dir2);
        retr = GSL_REAL( matrix_complex_trace(temp) );
        if(  fabs(retr) > 2.0  )
        {
          printf("The plaquette is out of bounds for some reason\n");
          print_matrix_complex(L->R[index]->link[dir1]);
        }
      }

}
