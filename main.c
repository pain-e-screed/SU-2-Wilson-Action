#include "includes.h"

int main()
{

  rng setup
  gsl_rng * r;
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r,time_seed());
  int r_max = gsl_rng_max(r);
  printf("RNG setup completed\n");

  //Allocate and Initalize Lattice
  lattice * L = lattice_allocation();
  printf("L occupies %lu bytes\n",sizeof(lattice ));
  printf("Lattice allocated\n The lattice has %d sites with %d link fields\n",L4,4*L4);
  initialize_lattice(r,L);
  printf("Lattice initalized\n");
  AnnealingSchedule(L,r);
  printf("lattice annealed\n");

  lattice_free(L);
  gsl_rng_free(r);
  return 0;
}
