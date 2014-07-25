#include "jurassic.h"
#include "control.h"
#include "atmosphere.h"
#include "forwardmodel.h"
#include "scatter.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Carry out forward model computation in a single directory. */
void call_formod(ctl_t *ctl,
		 const char *wrkdir,
		 const char *obsfile,
		 const char *atmfile,
		 const char *radfile,
		 const char *task,
		 const char *aerofile);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(int argc, char *argv[]) {
  
  static ctl_t ctl;
  
  FILE *in;
  
  char dirlist[LEN], wrkdir[LEN], task[LEN], aerofile[LEN];
  
  /* Check arguments... */
  if(argc<5)
    ERRMSG("Give parameters: <ctl> <obs> <atm> <rad>");
  
  /* Read forward model control parameters... */
  read_ctl(argc, argv, &ctl);
  
  /* Get task... */
  scan_ctl(argc, argv, "TASK", -1, "-", task);
  
  /* Get dirlist... */
  scan_ctl(argc, argv, "DIRLIST", -1, "-", dirlist);

  /* Get aero... */
  scan_ctl(argc, argv, "AEROFILE", -1, "-", aerofile);
  
  /* Single forward calculation... */
  if(dirlist[0]=='-')
    call_formod(&ctl, NULL, argv[2], argv[3], argv[4], task, aerofile);
  
  /* Work on directory list... */
  else {
    
    /* Open directory list... */
    if(!(in=fopen(dirlist, "r")))
      ERRMSG("Cannot open directory list!");
    
    /* Loop over directories... */
    while(fscanf(in, "%s", wrkdir)!=EOF) {
      
      /* Write info... */
      printf("\nWorking directory: %s\n", wrkdir);
      
      /* Call forward model... */
      call_formod(&ctl, wrkdir, argv[2], argv[3], argv[4], task, aerofile);
    }

    /* Close dirlist... */
    fclose(in);
  }
  
  return EXIT_SUCCESS;
}

/*****************************************************************************/

void call_formod(ctl_t *ctl,
		 const char *wrkdir,
		 const char *obsfile,
		 const char *atmfile,
		 const char *radfile,
		 const char *task,
		 const char *aerofile) {
  
  static atm_t atm, atm2;
  
  static obs_t obs;

  static aero_i aeroin;
  
  static aero_t aero;
  
  char filename[LEN];
  
  int ig, ig2, ip, iw;
  
  /* Read observation geometry... */
  read_obs(wrkdir, obsfile, ctl, &obs);
  
  /* Read atmospheric data... */
  read_atm(wrkdir, atmfile, ctl, &atm);
  
  /* Read aerosol and cloud data */
  if(aerofile[0]!='-' && ctl->sca_n>0) {
    read_aero(wrkdir, aerofile, ctl, &aeroin);
    /* Get aerosol/cloud optical properties */
    get_opt_prop(ctl, &aeroin, &aero);
  } 
  else if (aerofile[0]=='-' && ctl->sca_n>0) {
    ERRMSG("Please give aerosol file name or set SCA_N=0 for clear air simulation!");
  }

  /* Call forward model... */
  formod(ctl, &atm, &obs, &aero);
  
  /* Save radiance data... */
  write_obs(wrkdir, radfile, ctl, &obs);
  
  /* Compute contributions... */
  if(task[0]=='c' || task[0]=='C') {
    
    /* Switch off N2 and O2 continuum... */
    ctl->ctm_n2=0;
    ctl->ctm_o2=0;
    
    /* Loop over emitters... */
    for(ig=0; ig<ctl->ng; ig++) {
      
      /* Copy atmospheric data... */
      copy_atm(ctl, &atm2, &atm, 0);
      
      /* Set extinction to zero... */
      for(iw=0; iw<ctl->nw; iw++)
	for(ip=0; ip<atm2.np; ip++)
	  atm2.k[iw][ip]=0;
      
      /* Set volume mixing ratios to zero... */
      for(ig2=0; ig2<ctl->ng; ig2++)
	if(ig2!=ig)
	  for(ip=0; ip<atm2.np; ip++)
	    atm2.q[ig2][ip]=0;
      
      /* Call forward model... */
      formod(ctl, &atm2, &obs, &aero);
      
      /* Save radiance data... */
      sprintf(filename, "%s.%s", radfile, ctl->emitter[ig]);
	write_obs(wrkdir, filename, ctl, &obs);
    }
    
    /* Copy atmospheric data... */
    copy_atm(ctl, &atm2, &atm, 0);
    
    /* Set volume mixing ratios to zero... */
    for(ig=0; ig<ctl->ng; ig++)
      for(ip=0; ip<atm2.np; ip++)
	atm2.q[ig][ip]=0;
    
    /* Call forward model... */
    formod(ctl, &atm2, &obs, &aero);
    
    /* Save radiance data... */
    sprintf(filename, "%s.EXTINCT", radfile);
    write_obs(wrkdir, filename, ctl, &obs);
  }
}
