#ifndef _SIM_
#define _SIM_

int RunLatticeSimulation(char * file_name);
int latticeWriteToFile(lattice * L,char * file_name);
int latticeUnpack(lattice * L,char * file_name);
int AnnealSite(lattice * L,gsl_rng * r, long unsigned index, int dir, double beta, double epsilon );
double AnnealStepLattice(lattice * L, gsl_rng * r, double beta, double epsilon);
void AnnealingSchedule(lattice * L, gsl_rng *r );
void TestPlaquettes(lattice * L);

#endif
