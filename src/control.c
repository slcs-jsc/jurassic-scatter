#include "control.h"

/* Module to read and handle all input parameters. */

/*****************************************************************************/

void read_ctl(int argc,
	      char *argv[],
	      ctl_t *ctl) {
  
  int id, ig, iw;
  
  /* Write info... */
  printf("\nJURASSIC - Version with Scattering\n"
	 "Executable: %s\nControl parameters: %s\n", argv[0], argv[1]);
  
  /* Emitters... */
  ctl->ng=(int)scan_ctl(argc, argv, "NG", -1, "0", NULL);
  if(ctl->ng<0 || ctl->ng>NGMAX)
    ERRMSG("Set 0 <= NG <= MAX!");
  for(ig=0; ig<ctl->ng; ig++)
    scan_ctl(argc, argv, "EMITTER", ig, "", ctl->emitter[ig]);
  
  /* Radiance channels... */
  ctl->nd=(int)scan_ctl(argc, argv, "ND", -1, "0", NULL);
  if(ctl->nd<0 || ctl->nd>NDMAX)
    ERRMSG("Set 0 <= ND <= MAX!");
  for(id=0; id<ctl->nd; id++)
    ctl->nu[id]=scan_ctl(argc, argv, "NU", id, "", NULL);

  /* Spectral windows... */
  ctl->nw=(int)scan_ctl(argc, argv, "NW", -1, "1", NULL);
  if(ctl->nw<0 || ctl->nw>NWMAX)
    ERRMSG("Set 0 <= NW <= NWMAX!");
  for(id=0; id<ctl->nd; id++)
    ctl->window[id]=(int)scan_ctl(argc, argv, "WINDOW", id, "0", NULL);
  
  /* Emissivity look-up tables... */
  scan_ctl(argc, argv, "TBLBASE", -1, "-", ctl->tblbase);
  
  /* Hydrostatic equilibrium... */
  ctl->hydz=scan_ctl(argc, argv, "HYDZ", -1, "-999", NULL);
  
  /* Continua... */
  ctl->ctm_co2=(int)scan_ctl(argc, argv, "CTM_CO2", -1, "1", NULL);
  ctl->ctm_h2o=(int)scan_ctl(argc, argv, "CTM_H2O", -1, "1", NULL);
  ctl->ctm_n2=(int)scan_ctl(argc, argv, "CTM_N2", -1, "1", NULL);
  ctl->ctm_o2=(int)scan_ctl(argc, argv, "CTM_O2", -1, "1", NULL);
  
  /* Scattering on Aerosol/Clouds ... */
  ctl->sca_n=(int)scan_ctl(argc, argv, "SCA_N", -1, "0", NULL);
  if(ctl->sca_n<0 || ctl->sca_n>SCAMOD)
    ERRMSG("Set 0 <= SCA_NMOD <= MAX!");
  ctl->sca_mult=(int)scan_ctl(argc, argv, "SCA_MULT", -1, "1", NULL);
  scan_ctl(argc, argv, "SCA_EXT", -1, "beta_a", ctl->sca_ext);
  if(ctl->sca_n>0 && ctl->sca_mult==0 && 
     !strcasecmp(ctl->sca_ext, "beta_e") && 
     !strcasecmp(ctl->sca_ext, "beta_a"))
    ERRMSG("Please set extinction to beta_a or beta_e.");

  /* Interpolation of atmospheric data... */
  ctl->ip=(int)scan_ctl(argc, argv, "IP", -1, "1", NULL);
  ctl->cz=scan_ctl(argc, argv, "CZ", -1, "0", NULL);
  ctl->cx=scan_ctl(argc, argv, "CX", -1, "0", NULL);
  
  /* Ray-tracing... */
  ctl->refrac=(int)scan_ctl(argc, argv, "REFRAC", -1, "1", NULL);
  ctl->rayds=scan_ctl(argc, argv, "RAYDS", -1, "10", NULL);
  ctl->raydz=scan_ctl(argc, argv, "RAYDZ", -1, "1", NULL);
  ctl->transs=scan_ctl(argc, argv, "TRANSS", -1, "0.02", NULL);
  
  /* Field of view... */
  scan_ctl(argc, argv, "FOV", -1, "-", ctl->fov);
  
  /* Retrieval interface... */
  ctl->retp_zmin=scan_ctl(argc, argv, "RETP_ZMIN", -1, "-999", NULL);
  ctl->retp_zmax=scan_ctl(argc, argv, "RETP_ZMAX", -1, "-999", NULL);
  ctl->rett_zmin=scan_ctl(argc, argv, "RETT_ZMIN", -1, "-999", NULL);
  ctl->rett_zmax=scan_ctl(argc, argv, "RETT_ZMAX", -1, "-999", NULL);
  for(ig=0; ig<ctl->ng; ig++) {
    ctl->retq_zmin[ig]=scan_ctl(argc, argv, "RETQ_ZMIN", ig, "-999", NULL);
    ctl->retq_zmax[ig]=scan_ctl(argc, argv, "RETQ_ZMAX", ig, "-999", NULL);
  }
  for(iw=0; iw<ctl->nw; iw++) {
    ctl->retk_zmin[iw]=scan_ctl(argc, argv, "RETK_ZMIN", iw, "-999", NULL);
    ctl->retk_zmax[iw]=scan_ctl(argc, argv, "RETK_ZMAX", iw, "-999", NULL);
  }
  ctl->retnn=(int)scan_ctl(argc, argv, "RETNN", -1, "0", NULL);
  ctl->retrr=(int)scan_ctl(argc, argv, "RETRR", -1, "0", NULL);
  ctl->retss=(int)scan_ctl(argc, argv, "RETSS", -1, "0", NULL);

  ctl->retnn_zmin=(int)scan_ctl(argc, argv, "RETNN_ZMIN", -1, "-999", NULL);
  ctl->retnn_zmax=(int)scan_ctl(argc, argv, "RETNN_ZMAX", -1, "-999", NULL);

  /* Output flags... */
  ctl->write_bbt=(int)scan_ctl(argc, argv, "WRITE_BBT", -1, "0", NULL);
  ctl->write_matrix=(int)scan_ctl(argc, argv, "WRITE_MATRIX", -1, "0", NULL);
}

/*****************************************************************************/

double scan_ctl(int argc,
                char *argv[],
                const char *varname,
                int arridx,
                const char *defvalue,
                char *value) {
  
  FILE *in=NULL;
  
  char dummy[LEN], fullname1[LEN], fullname2[LEN], line[LEN],
    msg[2*LEN], rvarname[LEN], rval[LEN];
  
  int contain=0, i;
  
  /* Open file... */
  if(argv[1][0]!='-')
    if(!(in=fopen(argv[1], "r")))
      ERRMSG("Cannot open file!");
  
  /* Set full variable name... */
  if(arridx>=0) {
    sprintf(fullname1, "%s[%d]", varname, arridx);
    sprintf(fullname2, "%s[*]", varname);
  } else{
    sprintf(fullname1, "%s", varname);
    sprintf(fullname2, "%s", varname);
  }
  
  /* Read data... */
  if(in!=NULL)
    while(fgets(line, LEN, in))
      if(sscanf(line, "%s %s %s", rvarname, dummy, rval)==3)
        if(strcasecmp(rvarname, fullname1)==0 || 
           strcasecmp(rvarname, fullname2)==0) {
          contain=1;
          break;
        }
  for(i=1; i<argc-1; i++)
    if(strcasecmp(argv[i], fullname1)==0 || 
       strcasecmp(argv[i], fullname2)==0) {
      sprintf(rval, "%s", argv[i+1]);
      contain=1;
      break;
    }
  
  /* Close file... */
  if(in!=NULL)
    fclose(in);
  
  /* Check for missing variables... */
  if(!contain) {
    if(strlen(defvalue)>0)
      sprintf(rval, "%s", defvalue);
    else {
      sprintf(msg, "Missing variable %s!\n", fullname1);
      ERRMSG(msg);
    }
  }
  
  /* Write info... */
  printf("%s = %s\n", fullname1, rval);
  
  /* Return values... */
  if(value!=NULL)
    sprintf(value, "%s", rval);
  return atof(rval);
}
