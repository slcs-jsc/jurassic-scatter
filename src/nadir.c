#include "jurassic.h"
#include "control.h"
#include "misc.h"

int main(int argc, char *argv[]) {
  
  static ctl_t ctl;
  static obs_t obs;
  
  double dlat, lat, latmax, latmin, obsz;
  
  /* Check arguments... */
  if(argc<7)
    ERRMSG("Give parameters: <ctl> <obsz> <latmin> <latmax> <dlat> <obs>");
  
  /* Read forward model control parameters... */
  read_ctl(argc, argv, &ctl);
  
  /* Read arguments... */
  obsz=atof(argv[2]);
  latmin=atof(argv[3]);
  latmax=atof(argv[4]);
  dlat=atof(argv[5]);
  
  /* Create measurement geometry... */
  for(lat=latmin; lat<=latmax; lat+=dlat) {
    obs.obsz[obs.nr]=obsz;
    obs.vplat[obs.nr]=lat;
    if((++obs.nr)>=NRMAX)
      ERRMSG("Too many rays!");
  }
  
  /* Write observation data... */
  write_obs(NULL, argv[6], &ctl, &obs);
  
  return EXIT_SUCCESS;
}
