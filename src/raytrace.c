#include "jurassic.h"
#include "control.h"
#include "lineofsight.h"

int main(int argc, char *argv[]) {
  
  static atm_t atm, atm2;
  
  static ctl_t ctl;
  
  static los_t los;
  
  static obs_t obs;

  static aero_i aeroin;

  static aero_t aero;
  
  FILE *out;
  
  char filename[LEN], losbase[LEN], massbase[LEN], aerofile[LEN];
  
  double cgu[NLOS][NGMAX], s;
  
  int ig, ip, ir, iw;
  
  /* Check arguments... */
  if(argc<4)
    ERRMSG("Give parameters: <ctl> <obs> <atm>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  
  /* Get basenames... */
  scan_ctl(argc, argv, "LOSBASE", -1, "los", losbase);
  scan_ctl(argc, argv, "MASSBASE", -1, "-", massbase);

  /* Get aero... */
  scan_ctl(argc, argv, "AEROFILE", -1, "-", aerofile);
  
  /* Read observation geometry... */
  read_obs(NULL, argv[2], &ctl, &obs);
  
  /* Read atmospheric data... */
  read_atm(NULL, argv[3], &ctl, &atm);
  
  /* Read aerosol and cloud data */
  if(aerofile[0]!='-' && ctl.sca_n>0) {
    read_aero(NULL, aerofile, &ctl, &aeroin);
    /* Get aerosol/cloud optical properties */
    get_opt_prop(&ctl, &aeroin, &aero);
  } 
  else if (aerofile[0]=='-' && ctl.sca_n>0) {
    ERRMSG("Please give aerosol file name or set SCA_N=0 for clear air simulation!");
  }

  /* Create output file... */
  if(!(out=fopen("raytrace.tab", "w")))
    ERRMSG("Cannot create raytrace.tab!");
  
  /* Write header... */
  fprintf(out,
	  "# $1 = time (seconds since 2000-01-01T00:00Z)\n"
	  "# $2 = observer altitude [km]\n"
	  "# $3 = observer longitude [deg]\n"
	  "# $4 = observer latitude [deg]\n"
	  "# $5 = view point altitude [km]\n"
	  "# $6 = view point longitude [deg]\n"
	  "# $7 = view point latitude [deg]\n"
	  "# $8 = tangent point altitude [km]\n"
	  "# $9 = tangent point longitude [deg]\n"
	  "# $10 = tangent point latitude [deg]\n"
	  "# $11 = ray path index\n"
	  "# $12 = ray path length [km]\n");
  for(ig=0; ig<ctl.ng; ig++)
    fprintf(out, "# $%d = %s column density [molec/cm^2]\n",
	    13+ig, ctl.emitter[ig]);
  fprintf(out, "\n");
  
  /* Loop over rays... */
  for(ir=0; ir<obs.nr; ir++) {
    
    /* Raytracing... */
    raytrace(&ctl, &atm, &obs, &aero, &los, ir);
    
    /* Get ray path length... */
    s=0;
    for(ip=0; ip<los.np; ip++)
      s+=los.ds[ip];
    
    /* Get column densities... */
    for(ig=0; ig<ctl.ng; ig++) {
      cgu[0][ig]=los.u[0][ig];
      for(ip=1; ip<los.np; ip++)
	cgu[ip][ig]=cgu[ip-1][ig]+los.u[ip][ig];
    }
    
    /* Write to file... */
    fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g %d %g",
	    obs.time[ir], obs.obsz[ir], obs.obslon[ir], obs.obslat[ir],
	    obs.vpz[ir], obs.vplon[ir], obs.vplat[ir], obs.tpz[ir],
	    obs.tplon[ir], obs.tplat[ir], ir, s);
    for(ig=0; ig<ctl.ng; ig++)
      fprintf(out, " %g", cgu[los.np-1][ig]);
    fprintf(out, "\n");
    
    /* Atmospheric data... */
    if(losbase[0]!='-') {
      
      /* Copy data... */
      atm2.np=los.np;
      for(ip=0; ip<los.np; ip++) {
	atm2.time[ip]=obs.time[ir];
	atm2.z[ip]=los.z[ip];
	atm2.lon[ip]=los.lon[ip];
	atm2.lat[ip]=los.lat[ip];
	atm2.p[ip]=los.p[ip];
	atm2.t[ip]=los.t[ip];
	for(ig=0; ig<ctl.ng; ig++)
	  atm2.q[ig][ip]=los.q[ip][ig];
	for(iw=0; iw<ctl.nw; iw++)
	  atm2.k[iw][ip]=los.k[ip][iw];
      }
      
      /* Save data... */
      sprintf(filename, "%s.%d", losbase, ir);
      write_atm(NULL, filename, &ctl, &atm2);
    }

    /* Number and column densities... */
    if(massbase[0]!='-') {
      
      /* Copy data... */
      for(ip=0; ip<los.np; ip++) {
	atm2.p[ip]=atm2.t[ip]=0;
	for(ig=0; ig<ctl.ng; ig++)
	  atm2.q[ig][ip]=los.u[ip][ig]/(los.ds[ip]*1e5);
	for(iw=0; iw<ctl.nw; iw++)
	  atm2.k[iw][ip]=0;
      }
      
      /* Save data... */
      sprintf(filename, "%s.numdens.%d", massbase, ir);
      write_atm(NULL, filename, &ctl, &atm2);
      
      /* Copy data... */
      for(ip=0; ip<los.np; ip++)
	for(ig=0; ig<ctl.ng; ig++)
	  atm2.q[ig][ip]=cgu[ip][ig];
      
      /* Save data... */
      sprintf(filename, "%s.coldens.%d", massbase, ir);
      write_atm(NULL, filename, &ctl, &atm2);
    }
  }
  
  /* Close file... */
  fclose(out);
  
  return EXIT_SUCCESS;
}
