#include "jurassic.h"
#include "control.h"
#include "misc.h"

/* Create observations grid */

/* Read altitudes file. */
void read_altitudes(const char *dirname,
		    const char *filename,
		    double *altis,
		    int *ii);

int main(int argc, char *argv[]) {
  
  static ctl_t ctl;
  
  static obs_t obs;
  
  double dz, obsz, z, zmax, zmin, altis[NRMAX];
  
  char altfile[LEN];
  
  int nalt=0, i; 

  /* Check arguments... */
  if(argc<6)
    ERRMSG("Give parameters: <ctl> <obsz> <zmin> <zmax> <dz> <obs> \n or  <ctl> <obsz> ALTITUDEFILE <filename> <obs>");
  
  /* Read forward model control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Get altitude file... */
  scan_ctl(argc, argv, "ALTITUDEFILE", -1, "-", altfile);
  
  /* Read altitudes from file. */
  if(altfile[0]!='-') {
    read_altitudes(NULL, altfile, altis, &nalt);

    /* Read observer altitude... */
    obsz=atof(argv[2]);

    for(i=0; i<nalt; i++) {
      obs.obsz[obs.nr]=obsz;
      obs.obslat[obs.nr]=180/M_PI*acos((RE+altis[i])/(RE+obsz));
      obs.vpz[obs.nr]=altis[i];
      if((++obs.nr)>=NRMAX)
	ERRMSG("Too many rays!");
    }

    /* Write observation data... */
    write_obs(NULL, argv[5], &ctl, &obs);
  }
  /* Check arguments... */
  else if(argc<7) {
    ERRMSG("Give parameters: <ctl> <obsz> <zmin> <zmax> <dz> <obs>");
  }
  /* Generate altitudes from top, bottom and dz. */
  else {
    
    /* Read arguments... */
    obsz=atof(argv[2]);
    zmin=atof(argv[3]);
    zmax=atof(argv[4]);
    dz=atof(argv[5]);
    
    /* Create measurement geometry... */
    for(z=zmin; z<=zmax; z+=dz) {
      obs.obsz[obs.nr]=obsz;
      obs.obslat[obs.nr]=180./M_PI*acos((RE+z)/(RE+obsz));
      obs.vpz[obs.nr]=z;
      if((++obs.nr)>=NRMAX)
	ERRMSG("Too many rays!");
    }

    /* Write observation data... */
    write_obs(NULL, argv[6], &ctl, &obs);
  }
  
  return EXIT_SUCCESS;
}

/*****************************************************************************/

void read_altitudes(const char *dirname,
		    const char *filename,
		    double *altis,
		    int *ii){

  FILE *in;
  
  char file[LEN], line[LEN], *tok;
  
  /* Init... */
  *ii=0;
  
  /* Set filename... */
  if(dirname!=NULL)
    sprintf(file, "%s/%s", dirname, filename);
  else
    sprintf(file, "%s", filename);
  
  /* Write info... */
  printf("Read altitude data: %s\n", file);
  
  /* Open file... */
  if(!(in=fopen(file, "r")))
    ERRMSG("Cannot open file!");
  
  /* Read data... */
  while(fgets(line, LEN, in)) {
    /* Read data... */
    TOK(line, tok, "%lg", altis[*ii]);
    /* Increment counter... */
    ++*ii;
  }
  
  /* Close file... */
  fclose(in);
  
  /* Check number of points... */
  if(*ii<1)
    ERRMSG("Could not read any data!");
}
