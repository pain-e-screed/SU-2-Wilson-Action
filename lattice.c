
#include "includes.h"
/*
Below are my definitions for the lattice, site, and loc (de)-allocations.
The hierarchy for these structures is as follows
lattice->site[L4]
                 ->gsl_matrix_complex[4]

The actual call is given by
(lattice) L->R[index]->place
                     ->link[dir]
*/
//==============================================================================
site* site_allocation()
{
  int i;
  site * temp_site = (site * ) malloc(sizeof(site));
  FORALLDIR(i){ temp_site->link[i] = gsl_matrix_complex_alloc(2,2); }
  return temp_site;
}

void site_free(site* s)
{
  int i;
  FORALLDIR(i){ gsl_matrix_complex_free(s->link[i]);}
  free(s);
}

lattice * lattice_allocation()
{
  int i;
  lattice * temp_lattice = (lattice * ) malloc(sizeof(lattice));
  FORALLSITES(i){
    temp_lattice->R[i] = site_allocation();
  }
return temp_lattice;
}

void lattice_free(lattice * L)
{
  long unsigned i;
  FORALLSITES(i){site_free(L->R[i]);}
  free(L);
}
//==============================================================================


/* Now we have initializations of the lattice. I have it where it randomly
generates an SU(2) matrix and copies it into one of the links. It then goes
over every site and repeats the process for all links.
*/
//==============================================================================
void initialize_lattice(gsl_rng *r , lattice * lat1)
{
  gsl_matrix_complex *temp_tup = gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex *temp_xup = gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex *temp_yup = gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex *temp_zup = gsl_matrix_complex_alloc(2,2);

  long unsigned i;
  FORALLSITES(i)
  {
    gen_su2_matrix(temp_tup, r);
    gen_su2_matrix(temp_xup, r);
    gen_su2_matrix(temp_yup, r);
    gen_su2_matrix(temp_zup, r);

    gsl_matrix_complex_memcpy(lat1->R[i]->link[TUP],temp_tup);
    gsl_matrix_complex_memcpy(lat1->R[i]->link[XUP],temp_xup);
    gsl_matrix_complex_memcpy(lat1->R[i]->link[YUP],temp_yup);
    gsl_matrix_complex_memcpy(lat1->R[i]->link[ZUP],temp_zup);
  }
  gsl_matrix_complex_free(temp_tup);
  gsl_matrix_complex_free(temp_xup);
  gsl_matrix_complex_free(temp_yup);
  gsl_matrix_complex_free(temp_zup);
}

/* First, a routine for generating a single plaquette. That is,
  summing along a length 1 Wilson loop at a given site.
*/
//==============================================================================
 int plaquette(gsl_matrix_complex * P,lattice * L,long unsigned index, int mu, int nu)
{
  long unsigned indexmu, indexnu;
  gsl_matrix_complex * u1 = gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex * u2 = gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex * u3 = gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex * u4 = gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex * temp_left = gsl_matrix_complex_calloc(2,2);
  gsl_matrix_complex * temp_right = gsl_matrix_complex_calloc(2,2);
  gsl_matrix_complex * temp_total = gsl_matrix_complex_calloc(2,2);
  gsl_complex alpha, beta;
  if(mu != nu){
    indexmu = hop_index(index,1, mu,FORWARD);
    indexnu = hop_index(index,1, nu,FORWARD);
    alpha = gsl_complex_rect(1.0,0.0);
    beta = gsl_complex_rect(0.0,0.0);
    gsl_matrix_complex_memcpy(u1,L->R[index]->link[mu]);
    gsl_matrix_complex_memcpy(u2,L->R[indexmu]->link[nu]);
    gsl_matrix_complex_memcpy(u3,L->R[indexnu]->link[mu]);
    gsl_matrix_complex_memcpy(u4,L->R[index]->link[nu]);

    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans, alpha, u1, u2, beta, temp_left);
    gsl_blas_zgemm(CblasConjTrans, CblasConjTrans, alpha, u3, u4, beta, temp_right);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, temp_left, temp_right, beta, temp_total);

    gsl_matrix_complex_memcpy(P,temp_total);
    gsl_matrix_complex_free(u1);
    gsl_matrix_complex_free(u2);
    gsl_matrix_complex_free(u3);
    gsl_matrix_complex_free(u4);
    gsl_matrix_complex_free(temp_left);
    gsl_matrix_complex_free(temp_right);
    gsl_matrix_complex_free(temp_total);
    return 0;
  }
  else
  {
    gsl_matrix_complex_set_identity(P);
    gsl_matrix_complex_free(u1);
    gsl_matrix_complex_free(u2);
    gsl_matrix_complex_free(u3);
    gsl_matrix_complex_free(u4);
    gsl_matrix_complex_free(temp_left);
    gsl_matrix_complex_free(temp_right);
    gsl_matrix_complex_free(temp_total);
    return 0;
  }

}

//For an LxLxLxL lattice
//Calculate the index of a lattice site from its (t,x,y,z) coordinates

//If I shift from one site to an adjacent site in a direction, what is the new
//index of the site assuming periodic boundary conditions. if I need to change later, can do it here
long unsigned hop_index(long unsigned index,int jump, int dir, int FOB)
{
  int t,x,y,z;
  long unsigned val = index;

  z = val%LN;
  val = (val - z)/LN;
  y = val%LN;
  val = (val - y)/LN;
  x = val%LN;
  val = (val - x)/LN;
  t = val%LN;

  if(dir == ZUP)
    {
      z+=FOB*jump;
      if(OUTSIDE(z))
        z = (z+LN)%LN;
    }
  else if(dir == YUP)
    {
      y+=FOB*jump;
      if(OUTSIDE(y))
        y = (y+LN)%LN;
    }
  else if(dir == XUP)
    {
      x+=FOB*jump;
      if(OUTSIDE(x))
        x = (x+LN)%LN;
    }
  else
    {
      t+=FOB*jump;
      if(OUTSIDE(t))
        t = (t+LN)%LN;
    }
    val = LN*(LN*(LN*t+x)+y)+z;
  return val;
}

double partial_action(lattice * L, long unsigned index, int dir)
{
  int i;
  gsl_matrix_complex * A = gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex * B = gsl_matrix_complex_alloc(2,2);
  gsl_complex total;
  long unsigned temp_index;

  total = gsl_complex_rect(0.0,0.0);

  FORALLDIRBUT(i,dir)
  {
    temp_index = hop_index(index,1,i,BACKWARD);
    plaquette(A, L,temp_index, dir, i);
    plaquette(B, L,index, dir, i);
    total = gsl_complex_add(matrix_complex_trace(A),total);
    total = gsl_complex_add(matrix_complex_trace(B),total);
  }

  gsl_matrix_complex_free(A);
  gsl_matrix_complex_free(B);
  return -GSL_REAL(total)/2.0;
}

double action(lattice * L, double beta)
{
  long unsigned index;
  int i,j;
  double total = 0.0;
  gsl_matrix_complex * temp = gsl_matrix_complex_alloc(2,2);

  FORALLSITES(index)
  {
    FORALLDIR(i)
    {
      FORALLDIRBUT(j,i)
      {
        plaquette(temp,L,index,i,j);
        total += GSL_REAL( matrix_complex_trace(temp)  );
      }
    }
  }
  gsl_matrix_complex_free(temp);
  return -beta*total/2.0;
}

/* This function starts at a site on the lattice and goes in a particular
direction until it has the right lengths
index is the starting site location. J is how many link fields need to be
multiplied together. Dir is the direction to take. */
void WilsonArrow(lattice * L, long unsigned index, int J, int dir, gsl_matrix_complex * output)
{
  int i;
  unsigned long temp_index = index;
  gsl_matrix_complex * temp1 = gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex * temp2 = gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex * temp3 = gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex * identity = gsl_matrix_complex_alloc(2,2);
  gsl_complex alpha, beta;

  alpha = gsl_complex_rect(1.0,0.0);
  beta = gsl_complex_rect(0.0,0.0);
  gsl_matrix_complex_set_identity(identity);

  gsl_matrix_complex_memcpy(temp1,L->R[temp_index]->link[dir]);
  temp_index = hop_index(temp_index,1, dir,FORWARD);
  for(i=1;i<J;i++)
    {
      gsl_matrix_complex_memcpy(temp2,L->R[temp_index]->link[dir]);
      gsl_blas_zgemm(CblasNoTrans,CblasNoTrans, alpha, temp1,temp2, beta, temp3);
      gsl_matrix_complex_memcpy(temp1,temp3);
      temp_index = hop_index(temp_index,1, dir,FORWARD);
    }
  gsl_matrix_complex_memcpy(output, temp1);
  gsl_matrix_complex_free(temp1);
  gsl_matrix_complex_free(temp2);
  gsl_matrix_complex_free(identity);

}

double WilsonRectangle(lattice * L,long unsigned index, int I, int J, int mu,int nu)
{
  long unsigned index1,index2,index3;
  int i;
  double return_trace;
  gsl_matrix_complex * u1 = gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex * u2 = gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex * u3 = gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex * u4 = gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex * v1 = gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex * v2 = gsl_matrix_complex_alloc(2,2);
  gsl_complex alpha, beta;

  alpha = gsl_complex_rect(1.0,0.0);
  beta = gsl_complex_rect(0.0,0.0);

  index1 = index;
  index2 = hop_index(index,I,mu,FORWARD);
  index3 = hop_index(index,J,nu,FORWARD);
  if(I == 0 || J== 0)
    {gsl_matrix_complex_set_identity(u1); }
  else{
    WilsonArrow(L,index1,I,mu,u1);
    WilsonArrow(L,index2,J,nu,u2);
    WilsonArrow(L,index3,I,mu,u3);
    WilsonArrow(L,index1,J,nu,u4);

    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans, alpha, u1,u2, beta, v1);
    gsl_blas_zgemm(CblasConjTrans,CblasConjTrans, alpha, u3,u4, beta, v2);
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans, alpha, v1,v2, beta, u1);
  }

  return_trace=GSL_REAL(matrix_complex_trace(u1));
  gsl_matrix_complex_free(u1);
  gsl_matrix_complex_free(u2);
  gsl_matrix_complex_free(u3);
  gsl_matrix_complex_free(u4);
  gsl_matrix_complex_free(v1);
  gsl_matrix_complex_free(v2);
  return return_trace/2.0;
}

double WilsonExpectation(lattice * L,gsl_rng *r,  int I, int J)
{
  int i, j;
  int dir1, dir2;
  long unsigned index,N;
  double temp = 0.0;
  N = 0;
  if(I == 0 || J ==0 )
  {
    return 1.0;
  }
  else
  {
    FORALLSITES(index)
    {
      FORALLDIR(dir1)
      {
        FORALLDIRBUT(dir2,dir1)
        {
        temp+=WilsonRectangle(L,index,I,J,dir1,dir2);
        N++;
        }
      }
    }
  }
  return temp / (double) N  ;
}


double average_plaquette(lattice *L,gsl_rng * r)
{
  return WilsonExpectation(L,r,1,1);
}


double CreutzRatio(lattice * L, gsl_rng * r, int I, int J)
{
  double temp;
  temp = WilsonExpectation(L,r,I-1,J-1);
  temp *= WilsonExpectation(L,r,I,J);
  temp /= WilsonExpectation(L,r,I-1,J);
  temp /= WilsonExpectation(L,r,I,J-1);
  if(fabs(temp) > 1.0E-15 )
  {
    return -log(temp);
  }
  else
  {
    printf("Invalid Creutz ratio\n" );
    return 0.0;
  }
}
