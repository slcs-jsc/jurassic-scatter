#include "jurassic.h"
#include "control.h"
#include "misc.h"

/* Create observations grid */

int main(int argc, char *argv[]) {
  
  static ctl_t ctl;
  
  static obs_t obs;
  
  double dz, obsz, z, zmax, zmin;
  
  /* Check arguments... */
  if(argc<7)
    ERRMSG("Give parameters: <ctl> <obsz> <zmin> <zmax> <dz> <obs>");
  
  /* Read forward model control parameters... */
  read_ctl(argc, argv, &ctl);
  
  /* Read arguments... */
  obsz=atof(argv[2]);
  zmin=atof(argv[3]);
  zmax=atof(argv[4]);
  dz=atof(argv[5]);
  
  /* Create measurement geometry... */
  for(z=zmin; z<=zmax; z+=dz) {
    obs.obsz[obs.nr]=obsz;
    obs.obslat[obs.nr]=180/M_PI*acos((RE+z)/(RE+obsz));
    obs.vpz[obs.nr]=z;
    if((++obs.nr)>=NRMAX)
      ERRMSG("Too many rays!");
  }
  
  /* Write observation data... */
  write_obs(NULL, argv[6], &ctl, &obs);
  
  return EXIT_SUCCESS;
}
