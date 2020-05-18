#include "includes.h"

int main()
{

  char file_name[30], overwrite = "y";
  printf("Please enter the filename you wish to use for the lattice\n");
  scanf("%s",file_name);
  if(!access(file_name,F_OK))
  {
    printf("A lattice config file with this name has already been generated." );
    printf("Would you like to overwrite the file? y/n: " );
    scanf(" %c\n",&overwrite );
  }

  if(overwrite = "y")
  {
    if(!RunLatticeSimulation( file_name))
      printf("Lattice written to config file\n" );
  }
  else if(overwrite = "n")
  {
    printf("Skipping simulation and going straight to Lanczos\n" );
  }
  else
  {
    printf("This is not a valid input\n" );
    return 0;
  }

  lattice * L;
  gsl_vector_complex * q = randomComplexVector(8*N4);
  latticeUnpack(L,file_name);
  void (*f_ptr) (gsl_vector_complex *,gsl_vector_complex *, void * ctxt) = &WilsonDirac;

  Lanczos(f_ptr,q , 10, L );

}
