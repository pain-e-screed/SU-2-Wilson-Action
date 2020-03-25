#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_cblas.h>
#include<gsl/gsl_blas.h>
#include "includes.h"

#ifndef _LATTICE_
#define _LATTICE_

#define TUP 0
#define XUP 1
#define YUP 2
#define ZUP 3
#define NDIR 4

#define FORWARD 1
#define BACKWARD -1
//number of lattice sites along an axis:
#define LN 6
#define L4  LN*LN*LN*LN

#define FORALLDIRBUT(i,j) \
   for (i = TUP; i < NDIR; i++) if (i != j)
#define FORALLDIR(i) for(i = 0;i<NDIR;i++)
#define FORALLSITES(i) for(i=0;i<L4;i++)
#define OUTSIDE(i) (i<0 || i>LN-1)

typedef struct site{
  gsl_matrix_complex * link[4];
}site;

typedef struct lattice{
  site * R[L4];
}lattice;

site * site_allocation();
void site_free(site* s);
lattice * lattice_allocation();
void lattice_free(lattice * L);
void initialize_lattice(gsl_rng *r , lattice * lat1);
int plaquette(gsl_matrix_complex * P,lattice * L,long unsigned index, int mu, int nu);
long unsigned hop_index(long unsigned index,int jump, int dir, int FOB);
double partial_action(lattice * L, long unsigned index, int dir);

double action(lattice * L, double beta);
void WilsonArrow(lattice * L, long unsigned index, int J, int dir, gsl_matrix_complex * output);
double WilsonRectangle(lattice * L,long unsigned index, int I, int J, int mu,int nu);
double WilsonExpectation(lattice * L,gsl_rng *r,  int I, int J);
double average_plaquette(lattice *L,gsl_rng * r);
double CreutzRatio(lattice * L, gsl_rng * r, int I, int J);
#endif
