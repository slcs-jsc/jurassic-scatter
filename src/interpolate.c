#include "jurassic.h"
#include "control.h"
#include "atmosphere.h"

int main(int argc, char *argv[]) {
  
  static atm_t atm_in, atm_pts;
  static ctl_t ctl;
  
  /* Check arguments... */
  if(argc<5) 
    ERRMSG("Give parameters: <ctl> <atm_in> <atm_pts> <atm_out>");
  
  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  
  /* Read atmospheric data... */
  read_atm(NULL, argv[2], &ctl, &atm_in);
  read_atm(NULL, argv[3], &ctl, &atm_pts);
  
  /* Interpolate atmospheric data... */
  intpol_atm(&ctl, &atm_pts, &atm_in);
  
  /* Save interpolated data... */
  write_atm(NULL, argv[4], &ctl, &atm_pts);
  
  return EXIT_SUCCESS;
}
