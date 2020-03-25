#ifndef _SIM_
#define _SIM_

int AnnealSite(lattice * L,gsl_rng * r, long unsigned index, int dir, double beta, double epsilon );
double AnnealStepLattice(lattice * L, gsl_rng * r, double beta, double epsilon);
void AnnealingSchedule(lattice * L, gsl_rng *r );
void TestPlaquettes(lattice * L);

#endif
