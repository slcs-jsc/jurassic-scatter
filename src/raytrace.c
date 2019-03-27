#include "jurassic.h"
#include "control.h"
#include "lineofsight.h"

int main(int argc, char *argv[]) {
  
  static atm_t atm;
  
  static ctl_t ctl;
  
  static los_t los;
  
  static obs_t obs;
  
  static aero_t aero;
  
  FILE *out;
  
  char filename[LEN], aerofile[LEN];
    
  int id, ig, ip, ir, iw;
  
  /* Check arguments... */
  if(argc<4)
    ERRMSG("Give parameters: <ctl> <obs> <atm>");
  
  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  
  /* Get aero... */
  scan_ctl(argc, argv, "AEROFILE", -1, "-", aerofile);
  
  /* Read observation geometry... */
  read_obs(NULL, argv[2], &ctl, &obs);
  
  /* Read atmospheric data... */
  read_atm(NULL, argv[3], &ctl, &atm);
  
  /* Read aerosol and cloud data */
  if(aerofile[0]!='-' && ctl.sca_n>0) {
    read_aero(NULL, aerofile, &ctl, &aero);
    /* Get aerosol/cloud optical properties */
    get_opt_prop(&ctl, &aero);
  } else if (aerofile[0]=='-' && ctl.sca_n>0) {
    ERRMSG("Please give aerosol file name or set SCA_N=0 for clear air simulation!");
  }
  
  /* Loop over rays... */
  for(ir=0; ir<obs.nr; ir++) {
    
    /* Raytracing... */
    raytrace(&ctl, &atm, &obs, &aero, &los, ir);
    
    /* Create file... */
    sprintf(filename, "los.%d", ir);
    if(!(out=fopen(filename, "w")))
      ERRMSG("Cannot create los.tab!");
    
    /* Write header... */
    fprintf(out,
	    "# $1 = time (seconds since 2000-01-01T00:00Z)\n"
	    "# $2 = LOS point altitude [km]\n"
	    "# $3 = LOS point longitude [deg]\n"
	    "# $4 = LOS point latitude [deg]\n"
	    "# $5 = LOS point pressure [hPa]\n"
	    "# $6 = LOS point temperature [K]\n");
    for(ig=0; ig<ctl.ng; ig++)
      fprintf(out, "# $%d = LOS point %s volume mixing ratio \n",
	      7+ig, ctl.emitter[ig]);
    for(iw=0; iw<ctl.nw; iw++)
      fprintf(out, "# $%d = LOS point window %d extinction [1/km]\n",
	      7+ctl.ng+iw, iw);
    for(ig=0; ig<ctl.ng; ig++)
      fprintf(out, "# $%d = LOS point %s column density [molec/cm^2] \n",
	      7+ctl.ng+ctl.nw+ig, ctl.emitter[ig]);
    for(id=0; id<ctl.nd; id++)
      fprintf(out, "# $%d = LOS point beta_e(%.4f cm-1) [km-1] \n", 7+ctl.ng+ctl.nw+ctl.ng+id, ctl.nu[id]);
    for(id=0; id<ctl.nd; id++)
      fprintf(out, "# $%d = LOS point beta_s(%.4f cm-1) [km-1] \n", 7+ctl.ng+ctl.nw+ctl.ng+ctl.nd+id, ctl.nu[id]);
    for(id=0; id<ctl.nd; id++)
      fprintf(out, "# $%d = LOS point beta_a(%.4f cm-1) [km-1] \n", 7+ctl.ng+ctl.nw+ctl.ng+2*ctl.nd+id, ctl.nu[id]);
    fprintf(out, "# $%d = LOS segement length [km] \n", 7+ctl.ng+ctl.nw+ctl.ng+3*ctl.nd);
    
    fprintf(out, "\n");
    
    /* Loop over LOS points... */
    for(ip=0; ip<los.np; ip++) {
      fprintf(out, "%.2f %g %g %g %g %g", obs.time[ir],
	      los.z[ip], los.lon[ip], los.lat[ip],
	      los.p[ip], los.t[ip]);
      for(ig=0; ig<ctl.ng; ig++)
	fprintf(out, " %g", los.q[ip][ig]);
      for(iw=0; iw<ctl.nw; iw++)
	fprintf(out, " %g", los.k[ip][iw]);
      for(ig=0; ig<ctl.ng; ig++)
	fprintf(out, " %g", los.u[ip][ig]);
      for(id=0; id<ctl.nd; id++)
	fprintf(out, " %g", aero.beta_e[los.aeroi[ip]][id]*los.aerofac[ip]);
      for(id=0; id<ctl.nd; id++)
	fprintf(out, " %g", aero.beta_s[los.aeroi[ip]][id]*los.aerofac[ip]);
      for(id=0; id<ctl.nd; id++)
	fprintf(out, " %g", aero.beta_a[los.aeroi[ip]][id]*los.aerofac[ip]);
      fprintf(out, " %g", los.ds[ip]);      

      fprintf(out, "\n");
    }
    
    /* Close file... */
    fclose(out);

    printf("Wrote output to %s \n",filename);
  }

  return EXIT_SUCCESS;
}
