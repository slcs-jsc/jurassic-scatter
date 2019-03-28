#include "jurassic.h"
#include "control.h"
#include "atmosphere.h"

/* Little Helper: Program to define an atmosphere grid and to 
   add atmospheric parameters. One can either interpolate the 
   Remedios climatology for the polar winter, polar summer, 
   midlatitude, or equatorial region onto the defined grid, 
   or provide a single external profile to be interpolated 
   onto the defined grid. The midlatitude climatology and a 
   1 km are the defaults if no values are given. */

/* ------------------------------------------------------------
   Structs...
   ------------------------------------------------------------ */

/* Climate data struct. */
typedef struct {
  
  /* Altitude */
  double z[121];

  /* Pressure */
  double pre[121];

  /* Temperature */
  double tem[121];

  /* And all the long gas list ... */
  double c2h2[121];
  double c2h6[121];
  double ccl4[121];
  double ch4[121];
  double clo[121];
  double clono2[121];
  double co[121];
  double co2[121];
  double cof2[121];
  double f11[121];
  double f12[121];
  double f13[121];
  double f113[121];
  double f114[121];
  double f14[121];
  double f22[121];
  double h2o[121];
  double h2o2[121];
  double hcn[121];
  double hno3[121];
  double hno4[121];
  double hocl[121];
  double n2o[121];
  double n2o5[121];
  double nh3[121];
  double no[121];
  double no2[121];
  double o3[121];
  double ocs[121];
  double pan[121];
  double sf6[121];
  double so2[121];

} clim_t;

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Interpolate climatological data. */
void climatology(ctl_t *ctl,
		 clim_t *clim,
 		 atm_t *atm_mean); 

void get_clim_data(int czone,
		   clim_t *clim,
		   char *climpath);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(int argc, char *argv[]) {
  
  static atm_t atm;

  static ctl_t ctl;
  
  static clim_t clim;

  double dz, z, zmax, zmin;
  
  char zone[LEN], climpath[LEN];

  int czone=0;

  /* Check arguments... */
  if(argc<4)
    ERRMSG("Give parameters: <ctl> <atm_pts> <atm_mean>");
  
  /* Read forward model control parameters... */
  read_ctl(argc, argv, &ctl);
  
  /* Read geolocations... */
  zmin=scan_ctl(argc, argv, "ZMIN", -1, "0", NULL);
  zmax=scan_ctl(argc, argv, "ZMAX", -1, "0", NULL);
  dz=scan_ctl(argc, argv, "DZ", -1, "1", NULL);
  
  /* Set geolocations... */
  if(argv[2][0]!='-'){
    read_atm(NULL, argv[2], &ctl, &atm);
  } else {
    atm.np=0;
    for(z=zmin; z<=zmax; z+=dz) {
      atm.z[atm.np]=z;
      if((++atm.np)>NPMAX)
	ERRMSG("Too many data points!");
    }
    if(atm.np<=1)
      ERRMSG("Could not set atmospheric grid!");
  }
  
  /* Get climate zone */
  scan_ctl(argc, argv, "CLIMZONE", -1, "midl", zone);

  if (strcmp(zone, "midln")==0)
    czone=0;
  if (strcmp(zone, "pwin")==0)
    czone=1;
  if (strcmp(zone, "psum")==0)
    czone=2;
  if (strcmp(zone, "equn")==0)
    czone=3;
  if (strcmp(zone, "profile")==0)
    czone=4;

  /* get path to climatology file */
  scan_ctl(argc, argv, "CLIMPATH", -1, "", climpath);

  /* Get climate zone climatology */
  get_clim_data(czone, &clim, climpath);

  /* Interpolate climatological data... */
  climatology(&ctl, &clim, &atm);
  
  /* Write data to disk... */
  write_atm(NULL, argv[3], &ctl, &atm);
  
  return EXIT_SUCCESS;
}

/*****************************************************************************/

void get_clim_data(int czone,
		   clim_t *clim,
		   char *climpath){

  FILE *in;
  char file[LEN], line[LEN], *tok, name[LEN];
  int lat, month, zangle, lat1, lat2, month1, month2, zang1, zang2;
  double alt, mean, var;
  
  /* Set climate zone identifyer */
  /* midlatitude night - default */
  lat = 20;
  month = 1;
  zangle = 0;
  /* pwin */
  if (czone == 1){
    lat = 65;
    month = 1;
  }
  /* psum */
  if (czone == 2){
    lat = 65;
    month = 4;
  }
  /* equ */
  if (czone == 3)
    lat = -20;
      
  /* Read in climatology file */
  /* Set filename... */
  /* sprintf(file, "%s%s", climpath, "clim_remedios.tab"); */
  sprintf(file, "%s", climpath);
  
  /* Write info... */
  printf("Read climatology: %s\n", file);
  
  /* Open file... */
  if(!(in=fopen(file, "r")))
    ERRMSG("Cannot open file!");

  /* Read data... */
  while(fgets(line, LEN, in)) {

    /* Read data... */
    TOK(line, tok, "%s", *name);
    TOK(NULL, tok, "%d", lat1);
    TOK(NULL, tok, "%d", lat2);
    TOK(NULL, tok, "%d", month1);
    TOK(NULL, tok, "%d", month2);
    TOK(NULL, tok, "%d", zang1);
    TOK(NULL, tok, "%d", zang2);
    TOK(NULL, tok, "%lg", alt);
    TOK(NULL, tok, "%lg", mean);
    TOK(NULL, tok, "%lg", var);
    
    /* Select climate zone */
    if ((lat==lat1 && month==month1 && zangle==zang1) || czone== 4){
      
      if (strcmp(name, "TEMPERATURE")==0){
	clim->tem[(int)alt] = mean;
	clim->z[(int)alt] = alt;
      }
      if (strcmp(name, "PRESSURE")==0)
	clim->pre[(int)alt] = mean;
      if (strcmp(name, "C2H2")==0)
	clim->c2h2[(int)alt] = mean;
      if (strcmp(name, "C2H6")==0)
	clim->c2h6[(int)alt] = mean;
      if (strcmp(name, "CCL4")==0)
	clim->ccl4[(int)alt] = mean;
      if (strcmp(name, "CH4")==0)
	clim->ch4[(int)alt] = mean;
      if (strcmp(name, "CLO")==0)
	clim->clo[(int)alt] = mean;
      if (strcmp(name, "CLONO2")==0)
	clim->clono2[(int)alt] = mean;
      if (strcmp(name, "CO")==0)
	clim->co[(int)alt] = mean;
      if (strcmp(name, "CO2")==0)
      	clim->co2[(int)alt] = mean;
      if (strcmp(name, "COF2")==0)
	clim->cof2[(int)alt] = mean;
      if (strcmp(name, "F11")==0)
	clim->f11[(int)alt] = mean;
      if (strcmp(name, "F12")==0)
	clim->f12[(int)alt] = mean;
      if (strcmp(name, "F14")==0)
	clim->f14[(int)alt] = mean;
      if (strcmp(name, "F22")==0)
	clim->f22[(int)alt] = mean;
      if (strcmp(name, "H2O")==0)
	clim->h2o[(int)alt] = mean;
      if (strcmp(name, "H2O2")==0)
	clim->h2o2[(int)alt] = mean;
      if (strcmp(name, "HCN")==0)
	clim->hcn[(int)alt] = mean;
      if (strcmp(name, "HNO3")==0)
	clim->hno3[(int)alt] = mean;
      if (strcmp(name, "HNO4")==0)
	clim->hno4[(int)alt] = mean;
      if (strcmp(name, "HOCL")==0)
	clim->hocl[(int)alt] = mean;
      if (strcmp(name, "N2O")==0)
	clim->n2o[(int)alt] = mean;
      if (strcmp(name, "N2O5")==0)
	clim->n2o5[(int)alt] = mean;
      if (strcmp(name, "NH3")==0)
	clim->nh3[(int)alt] = mean;
      if (strcmp(name, "NO")==0)
	clim->no[(int)alt] = mean;
      if (strcmp(name, "NO2")==0)
	clim->no2[(int)alt] = mean;
      if (strcmp(name, "O3")==0)
	clim->o3[(int)alt] = mean;
      if (strcmp(name, "OCS")==0)
	clim->ocs[(int)alt] = mean;
      if (strcmp(name, "PAN")==0)
	clim->pan[(int)alt] = mean;
      if (strcmp(name, "SF6")==0)
	clim->sf6[(int)alt] = mean;
      if (strcmp(name, "SO2")==0)
	clim->so2[(int)alt] = mean; 
    }
    if (strcmp(name, "F113")==0)
      clim->f113[(int)alt] = mean;
    if (strcmp(name, "F114")==0)
      clim->f114[(int)alt] = mean;
  }
  
  /* Close file... */
  fclose(in);
}

/*****************************************************************************/

void climatology(ctl_t *ctl,
		 clim_t *clim,
		 atm_t *atm) {

  /* static int ig_co2=-999; */
  
  double *q[NGMAX]={NULL}; /* co2, */
  
  int ig, ip, iw, iz;

  /* /\* Find emitter index of CO2... *\/ */
  /* if(ig_co2==-999) */
  /*   ig_co2=find_emitter(ctl, "CO2"); */
  
  /* Identify variable... */
  for(ig=0; ig<ctl->ng; ig++) {
    q[ig]=NULL;
    if(strcasecmp(ctl->emitter[ig], "C2H2")==0) q[ig]=clim->c2h2;
    if(strcasecmp(ctl->emitter[ig], "C2H6")==0) q[ig]=clim->c2h6;
    if(strcasecmp(ctl->emitter[ig], "CCl4")==0) q[ig]=clim->ccl4;
    if(strcasecmp(ctl->emitter[ig], "CH4")==0) q[ig]=clim->ch4;
    if(strcasecmp(ctl->emitter[ig], "ClO")==0) q[ig]=clim->clo;
    if(strcasecmp(ctl->emitter[ig], "ClONO2")==0) q[ig]=clim->clono2;
    if(strcasecmp(ctl->emitter[ig], "CO")==0) q[ig]=clim->co;
    if(strcasecmp(ctl->emitter[ig], "CO2")==0) q[ig]=clim->co2;
    if(strcasecmp(ctl->emitter[ig], "COF2")==0) q[ig]=clim->cof2;
    if(strcasecmp(ctl->emitter[ig], "F11")==0) q[ig]=clim->f11;
    if(strcasecmp(ctl->emitter[ig], "F113")==0) q[ig]=clim->f113;
    if(strcasecmp(ctl->emitter[ig], "F114")==0) q[ig]=clim->f114;
    if(strcasecmp(ctl->emitter[ig], "F12")==0) q[ig]=clim->f12;
    if(strcasecmp(ctl->emitter[ig], "F14")==0) q[ig]=clim->f14;
    if(strcasecmp(ctl->emitter[ig], "F22")==0) q[ig]=clim->f22;
    if(strcasecmp(ctl->emitter[ig], "H2O")==0) q[ig]=clim->h2o;
    if(strcasecmp(ctl->emitter[ig], "H2O2")==0) q[ig]=clim->h2o2;
    if(strcasecmp(ctl->emitter[ig], "HCN")==0) q[ig]=clim->hcn;
    if(strcasecmp(ctl->emitter[ig], "HNO3")==0) q[ig]=clim->hno3;
    if(strcasecmp(ctl->emitter[ig], "HNO4")==0) q[ig]=clim->hno4;
    if(strcasecmp(ctl->emitter[ig], "HOCl")==0) q[ig]=clim->hocl;
    if(strcasecmp(ctl->emitter[ig], "N2O")==0) q[ig]=clim->n2o;
    if(strcasecmp(ctl->emitter[ig], "N2O5")==0) q[ig]=clim->n2o5;
    if(strcasecmp(ctl->emitter[ig], "NH3")==0) q[ig]=clim->nh3;
    if(strcasecmp(ctl->emitter[ig], "NO")==0) q[ig]=clim->no;
    if(strcasecmp(ctl->emitter[ig], "NO2")==0) q[ig]=clim->no2;
    if(strcasecmp(ctl->emitter[ig], "O3")==0) q[ig]=clim->o3;
    if(strcasecmp(ctl->emitter[ig], "OCS")==0) q[ig]=clim->ocs;
    if(strcasecmp(ctl->emitter[ig], "PAN")==0) q[ig]=clim->pan;
    if(strcasecmp(ctl->emitter[ig], "SF6")==0) q[ig]=clim->sf6;
    if(strcasecmp(ctl->emitter[ig], "SO2")==0) q[ig]=clim->so2;
  }
  
  /* Loop over atmospheric data points... */
  for(ip=0; ip<atm->np; ip++) {
    
    /* Get altitude index... */
    iz=locate(clim->z, 121, atm->z[ip]);
    
    printf("%f %f \n", atm->z[ip], atm->p[ip] );

    /* Interpolate pressure... */
    atm->p[ip]=EXP(clim->z[iz], clim->pre[iz], clim->z[iz+1], clim->pre[iz+1], atm->z[ip]);
    
    /* Interpolate temperature... */
    atm->t[ip]=LIN(clim->z[iz], clim->tem[iz], clim->z[iz+1], clim->tem[iz+1], atm->z[ip]);
    
    /* Interpolate trace gases... */
    for(ig=0; ig<ctl->ng; ig++)
      if(q[ig]!=NULL)
	atm->q[ig][ip]=LIN(clim->z[iz], q[ig][iz], clim->z[iz+1], q[ig][iz+1], atm->z[ip]);
      else
	atm->q[ig][ip]=0;

    /* /\* Set CO2... *\/ */
    /* if(ig_co2>=0) { */
    /*   co2=371.92429e-6+1.840618e-6*(atm->time[ip]-63158400.)/31557600.; */
    /*   /\* co2=co2*1; *\/ */
    /*   atm->q[ig_co2][ip]=co2; */
    /* } */
    
    /* Set extinction to zero... */
    for(iw=0; iw<ctl->nw; iw++)
      atm->k[iw][ip]=0;
  }
}

