#include "forwardmodel.h"

int main(int argc, char *argv[]) {

  double t, nu;
  
  /* Check arguments... */
  if(argc<3)
    ERRMSG("Give parameters: <t> <nu>");
  
  /* Read arguments... */
  t=atof(argv[1]);
  nu=atof(argv[2]);
  
  /* Compute Planck function... */
  printf("%g\n", planck(t, nu));
  
  return EXIT_SUCCESS;
}
