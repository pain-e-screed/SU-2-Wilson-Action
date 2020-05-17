#include "includes.h"

int RunLatticeSimulation(char * file_name)
{
  int return_status=0;
  //rng setup
  gsl_rng * r;
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r,time_seed());
  int r_max = gsl_rng_max(r);
  printf("RNG setup completed.\n");

  //Allocate and Initalize Lattice
  lattice * L = lattice_allocation();
  printf("L occupies %lu bytes\n",sizeof(lattice ));
  printf("Lattice allocated\n The lattice has %d sites and %d link fields.\n",L4,4*L4);
  initialize_lattice(r,L);
  printf("Lattice initalized\n");


  AnnealingSchedule(L,r);
  printf("Lattice annealed.\n");

  latticeWriteToFile(L,file_name);

  lattice_free(L);
  gsl_rng_free(r);
  return return_status;

}

int latticeWriteToFile(lattice * L,char * file_name)
{
  int i;
  long unsigned r;

  FILE * file_pointer;
  file_pointer = fopen(file_name, "wb");
  FORALLSITES(r)
    FORALLDIR(i)
    if(!gsl_matrix_fwrite(file_pointer, L->R[r]->link[i]))
    {
      printf("Matrix fwrite failed at site %ld and link %d\n",r,i );
      return 1;
    }

  return 0;
}


int latticeUnpack(lattice * L,char * file_name)
{
  int i;
  long unsigned r;

  FILE * file_pointer;
  file_pointer = fopen(file_name, "rb");
  FORALLSITES(r)
    FORALLDIR(i)
    if(!gsl_matrix_fread(file_pointer, L->R[r]->link[i]))
    {
      printf("Matrix fread failed at site %ld and link %d\n",r,i );
      return 1;
    }
}

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
  const int iter = 30;
  double B[iter], X[iter],steps[iter];
  int i, j;
  FILE * f_ptr;

  f_ptr = fopen("values.txt","w");
  beta = 0.1;
  steps[0] = 0.5;
  for(i=1;i<iter;i++)
    steps[i] = steps[i-1] * 0.9;

  for(j=0;j<iter;j++)
  {
    for(int i=0;i<500;i++)
    {
      AccRej = AnnealStepLattice( L, r,beta,steps[j]);
    }
    B[j] = beta;
    X[j] =  CreutzRatio(L,r,1,1);
    beta += 0.1;
    //step *= 0.9;
    if(!isnan(X[j]))
      fprintf(f_ptr, "%lf %lf %lf\n",B[j]/4.0,log(4.0/B[j]), X[j]);
    printf("beta = %lf, X(1,1) = %lf\n",B[j],X[j]);
    printf("ln(g_0^2):%lf\n",log(4.0/B[j]) );
    printf("g_0^-2:%lf\n",B[j]/4.0 );
    printf("Action:%lf\n\n\n",action(L,B[j]));
  }
  fclose(f_ptr);
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
