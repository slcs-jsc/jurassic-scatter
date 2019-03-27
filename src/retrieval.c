#include "retrievalmodel.h"

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(int argc, char *argv[]) {
  
  static ret_t ret;
  
  static ctl_t ctl;
  
  static atm_t atm_i, atm_apr;
  
  static obs_t obs_i, obs_meas;
  
  static aero_t aero_i, aero_apr;

  FILE *dirlist;
  
  int id, ir, nbad;
  
  size_t m;
  
  /* Check arguments... */
  if(argc<3)
    ERRMSG("Give parameters: <ctl> <dirlist>");
  
  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  read_ret(argc, argv, &ctl, &ret);
  
  /* Open directory list... */
  if(!(dirlist=fopen(argv[2], "r")))
    ERRMSG("Cannot open directory list!");
  
  /* Loop over directories... */
  while(fscanf(dirlist, "%s", ret.dir)!=EOF) {
    
    /* Write info... */
    printf("\nRetrieve in directory %s...\n\n", ret.dir);
    
    /* Read atmospheric data... */
    read_atm(ret.dir, "atm_apr.tab", &ctl, &atm_apr);

    /* Read particle data... */
    read_aero(ret.dir, "aero_apr.tab", &ctl, &aero_apr);
    
    /* Read observation data... */
    read_obs(ret.dir, "obs_meas.tab", &ctl, &obs_meas);
    
    /* Retrieval... */
    while(1) {

      /* Run retrieval... */
      optimal_estimation(&ret, &ctl, &obs_meas, &obs_i, &atm_apr, &atm_i, &aero_apr, &aero_i);

      /* Check radiance residuals... */
      nbad=0;
      for(id=0; id<ctl.nd; id++)
	for(ir=0; ir<obs_meas.nr; ir++)
	  if(ret.resmax>0
	     && gsl_finite(obs_i.rad[id][ir])
	     && gsl_finite(obs_meas.rad[id][ir])
	     && fabs(1-obs_i.rad[id][ir]/obs_meas.rad[id][ir])
	     >=ret.resmax/100) {
	    obs_meas.rad[id][ir]=obs_i.rad[id][ir]=GSL_NAN;
	    nbad++;
	  }
      
      /* Redo retrieval? */
      m=obs2y(&ctl, &obs_meas, NULL, NULL, NULL);
      if(nbad>0 && m>0)
	printf("\nFound %d bad measurements. Redo retrieval...\n\n", nbad);
      else
	break;
    }
  }
  
  /* Write info... */
  printf("\nRetrieval done...\n");
  
  return EXIT_SUCCESS;
}

/*****************************************************************************/
