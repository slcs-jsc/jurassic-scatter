#include "misc.h"

/*****************************************************************************/

void cart2geo(double *x,
	      double *z,
	      double *lon,
	      double *lat) {
  
  double radius;
  
  radius=NORM(x);
  *lat=asin(x[2]/radius)*180/M_PI;
  *lon=atan2(x[1], x[0])*180/M_PI;
  *z=radius-RE;
}

/*****************************************************************************/

void copy_obs(ctl_t *ctl,
	      obs_t *obs_dest,
	      obs_t *obs_src,
	      int init) {
  
  int id, ir;
  
  size_t s;
  
  /* Data size... */
  s=(size_t)obs_src->nr*sizeof(double);
  
  /* Copy data... */
  obs_dest->nr=obs_src->nr;
  memcpy(obs_dest->time, obs_src->time, s);
  memcpy(obs_dest->obsz, obs_src->obsz, s);
  memcpy(obs_dest->obslon, obs_src->obslon, s);
  memcpy(obs_dest->obslat, obs_src->obslat, s);
  memcpy(obs_dest->vpz, obs_src->vpz, s);
  memcpy(obs_dest->vplon, obs_src->vplon, s);
  memcpy(obs_dest->vplat, obs_src->vplat, s);
  memcpy(obs_dest->tpz, obs_src->tpz, s);
  memcpy(obs_dest->tplon, obs_src->tplon, s);
  memcpy(obs_dest->tplat, obs_src->tplat, s);
  for(id=0; id<ctl->nd; id++) {
    memcpy(obs_dest->rad[id], obs_src->rad[id], s);
    memcpy(obs_dest->tau[id], obs_src->tau[id], s);
  }
  
  /* Initialize... */
  if(init)
    for(ir=0; ir<obs_dest->nr; ir++)
      for(id=0; id<ctl->nd; id++)
	if(gsl_finite(obs_dest->rad[id][ir])) {
	  obs_dest->rad[id][ir]=0;
	  obs_dest->tau[id][ir]=0;
	}
}

/*****************************************************************************/

void geo2cart(double z, 
	      double lon,
	      double lat,
	      double *x) {
  
  double radius;
  
  radius=z+RE;
  x[0]=radius*cos(lat/180*M_PI)*cos(lon/180*M_PI);
  x[1]=radius*cos(lat/180*M_PI)*sin(lon/180*M_PI);
  x[2]=radius*sin(lat/180*M_PI);
}

/*****************************************************************************/

void init_tbl(ctl_t *ctl,
	      tbl_t *tbl) {
  
  FILE *in;
  
  char filename[LEN], line[LEN];
  
  double eps, eps_old, press, press_old, temp, temp_old, u, u_old,
    f[NSHAPE], fsum, nu[NSHAPE], tmin=100, tmax=400;
  
  int i, id, ig, ip, it, n;
  
  /* Loop over trace gases and channels... */
  for(ig=0; ig<ctl->ng; ig++)
    for(id=0; id<ctl->nd; id++) {
      
      /* Set filename... */
      sprintf(filename, "%s_%.4f_%s.bin",
	      ctl->tblbase, ctl->nu[id], ctl->emitter[ig]);
      
      /* Try to open binary file... */
      if((in=fopen(filename, "r"))) {
	
	/* Write info... */
	LOGMSG(2, printf("Read emissivity table: %s\n", filename));
	
	/* Read data... */
	FREAD(&tbl->np[ig][id], int, 1, in);
	if(tbl->np[ig][id]>TBLNPMAX)
	  ERRMSG("Too many pressure levels!");
	FREAD(&tbl->p[ig][id], double, tbl->np[ig][id], in);
	FREAD(tbl->nt[ig][id], int, tbl->np[ig][id], in);
	for(ip=0; ip<tbl->np[ig][id]; ip++) {
	  if(tbl->nt[ig][id][ip]>TBLNTMAX)
	    ERRMSG("Too many temperatures!");
	  FREAD(&tbl->t[ig][id][ip], double, tbl->nt[ig][id][ip], in);
	  FREAD(tbl->nu[ig][id][ip], int, tbl->nt[ig][id][ip], in);
	  for(it=0; it<tbl->nt[ig][id][ip]; it++) {
	    FREAD(&tbl->u[ig][id][ip][it], float,
		  GSL_MIN(tbl->nu[ig][id][ip][it], TBLNUMAX), in);
	    FREAD(&tbl->eps[ig][id][ip][it], float,
		  GSL_MIN(tbl->nu[ig][id][ip][it], TBLNUMAX), in);
	  }
	}
	
	/* Close file... */
	fclose(in);
      }
      
      /* Try to read ASCII file... */
      else {
	
	/* Initialize... */
	tbl->np[ig][id]=-1;
	eps_old=-999;
	press_old=-999;
	temp_old=-999;
	u_old=-999;
	
	/* Try to open file... */
	sprintf(filename, "%s_%.4f_%s.tab",
		ctl->tblbase, ctl->nu[id], ctl->emitter[ig]);

	if(!(in=fopen(filename, "r"))) {
	  LOGMSG(2, printf("Missing emissivity table: %s\n", filename));
	  continue;
	}
	LOGMSG(2, printf("Read emissivity table: %s\n", filename));
	
	/* Read data... */
	while(fgets(line, LEN, in)) {
	  
	  /* Parse line... */
	  if(sscanf(line,"%lg %lg %lg %lg", &press, &temp, &u, &eps)!=4)
	    continue;
	  
	  /* Determine pressure index... */
	  if(press!=press_old) {
	    press_old=press;
	    if((++tbl->np[ig][id])>=TBLNPMAX)
	      ERRMSG("Too many pressure levels!");
	    tbl->nt[ig][id][tbl->np[ig][id]]=-1;
	  }
	  
	  /* Determine temperature index... */
	  if(temp!=temp_old) {
	    temp_old=temp;
	    if((++tbl->nt[ig][id][tbl->np[ig][id]])>=TBLNTMAX)
	      ERRMSG("Too many temperatures!");
	    tbl->nu[ig][id][tbl->np[ig][id]]
	      [tbl->nt[ig][id][tbl->np[ig][id]]]=-1;
	  }
	  
	  /* Determine column density index... */
	  if((eps>eps_old && u>u_old) || tbl->nu[ig][id][tbl->np[ig][id]]
	     [tbl->nt[ig][id][tbl->np[ig][id]]]<0) {
	    eps_old=eps;
	    u_old=u;
	    if((++tbl->nu[ig][id][tbl->np[ig][id]]
		[tbl->nt[ig][id][tbl->np[ig][id]]])>=TBLNUMAX) {
	      tbl->nu[ig][id][tbl->np[ig][id]]
		[tbl->nt[ig][id][tbl->np[ig][id]]]--;
	      continue;
	    }
	  }
	  
	  /* Store data... */
	  tbl->p[ig][id][tbl->np[ig][id]]=press;
	  tbl->t[ig][id][tbl->np[ig][id]][tbl->nt[ig][id][tbl->np[ig][id]]]
	    =temp;
	  tbl->u[ig][id][tbl->np[ig][id]][tbl->nt[ig][id][tbl->np[ig][id]]]
	    [tbl->nu[ig][id][tbl->np[ig][id]]
	     [tbl->nt[ig][id][tbl->np[ig][id]]]]=(float)u;
	  tbl->eps[ig][id][tbl->np[ig][id]][tbl->nt[ig][id][tbl->np[ig][id]]]
	    [tbl->nu[ig][id][tbl->np[ig][id]]
	     [tbl->nt[ig][id][tbl->np[ig][id]]]]=(float)eps;
	}
	
	/* Increment counters... */
	tbl->np[ig][id]++;
	for(ip=0; ip<tbl->np[ig][id]; ip++) {
	  tbl->nt[ig][id][ip]++;
	  for(it=0; it<tbl->nt[ig][id][ip]; it++)
	    tbl->nu[ig][id][ip][it]++;
	}
	
	/* Close file... */
	fclose(in);
      }
    }
  
  /* Write info... */
  LOGMSG(2, printf("Initialize source function table...\n"));
  
  /* Loop over channels... */
  for(id=0; id<ctl->nd; id++) {
    
    /* Read filter function... */
    sprintf(filename, "%s_%.4f.filt", ctl->tblbase, ctl->nu[id]);
    read_shape(filename, nu, f, &n);
    
    /* Compute source function table... */
    for(it=0; it<TBLNSMAX; it++) {
      
      /* Set temperature... */
      tbl->st[it]=LIN(0.0, tmin, TBLNSMAX-1.0, tmax, (double)it);
      
      /* Integrate Planck function... */
      fsum=0;
      tbl->sr[id][it]=0;
      for(i=0; i<n; i++) {
	fsum+=f[i];
	tbl->sr[id][it]+=f[i]*planck(tbl->st[it], nu[i]);
      }
      tbl->sr[id][it]/=fsum;
    }
  }
}

/*****************************************************************************/

int locate(double *xx,
	   int n,
	   double x) {
  
  int i, ilo, ihi;
  
  ilo=0;
  ihi=n-1;
  i=(ihi+ilo)>>1;
  
  if(xx[i]<xx[i+1])
    while(ihi>ilo+1) {
      i=(ihi+ilo)>>1;
      if(xx[i]>x)
	ihi=i;
      else
	ilo=i;
    }
  else
    while(ihi>ilo+1) {
      i=(ihi+ilo)>>1;
      if(xx[i]<=x)
	ihi=i;
      else
	ilo=i;
    }
  
  return ilo;
}

/*****************************************************************************/

void read_obs(const char *dirname,
	      const char *filename,
	      ctl_t *ctl,
	      obs_t *obs) {
  
  FILE *in;
  
  char file[LEN], line[LEN], *tok;
  
  int id;
  
  /* Init... */
  obs->nr=0;
  
  /* Set filename... */
  if(dirname!=NULL)
    sprintf(file, "%s/%s", dirname, filename);
  else
    sprintf(file, "%s", filename);
  
  /* Write info... */
  printf("Read observation data: %s\n", file);
  
  /* Open file... */
  if(!(in=fopen(file, "r")))
    ERRMSG("Cannot open file!");
  
  /* Read data... */
  while(fgets(line, LEN, in)) {
    
    /* Read data... */
    TOK(line, tok, "%lg", obs->time[obs->nr]);
    TOK(NULL, tok, "%lg", obs->obsz[obs->nr]);
    TOK(NULL, tok, "%lg", obs->obslon[obs->nr]);
    TOK(NULL, tok, "%lg", obs->obslat[obs->nr]);
    TOK(NULL, tok, "%lg", obs->vpz[obs->nr]);
    TOK(NULL, tok, "%lg", obs->vplon[obs->nr]);
    TOK(NULL, tok, "%lg", obs->vplat[obs->nr]);
    TOK(NULL, tok, "%lg", obs->tpz[obs->nr]);
    TOK(NULL, tok, "%lg", obs->tplon[obs->nr]);
    TOK(NULL, tok, "%lg", obs->tplat[obs->nr]);
    for(id=0; id<ctl->nd; id++)
      TOK(NULL, tok, "%lg", obs->rad[id][obs->nr]);
    for(id=0; id<ctl->nd; id++)
      TOK(NULL, tok, "%lg", obs->tau[id][obs->nr]);
    
    /* Increment counter... */
    if((++obs->nr)>NRMAX)
      ERRMSG("Too many rays!");
  }
  
  /* Close file... */
  fclose(in);
  
  /* Check number of points... */
  if(obs->nr<1)
    ERRMSG("Could not read any data!");
}

/*****************************************************************************/

void write_obs(const char *dirname,
	       const char *filename,
	       ctl_t *ctl,
	       obs_t *obs) {
  
  FILE *out;

  char file[LEN];
  
  int id, ir, n=10;

  /* Set filename... */
  if(dirname!=NULL)
    sprintf(file, "%s/%s", dirname, filename);
  else
    sprintf(file, "%s", filename);
  
  /* Write info... */
  printf("Write observation data: %s\n", file);
  
  /* Create file... */
  if(!(out=fopen(file, "w")))
    ERRMSG("Cannot create file!");
  
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
	  "# $10 = tangent point latitude [deg]\n");
  for(id=0; id<ctl->nd; id++)
    fprintf(out, "# $%d = channel %g: radiance [W/(m^2 sr cm^-1)]\n",
	    ++n, ctl->nu[id]);
  for(id=0; id<ctl->nd; id++)
    fprintf(out, "# $%d = channel %g: transmittance\n",
	    ++n, ctl->nu[id]);
  
  /* Write data... */
  for(ir=0; ir<obs->nr; ir++) {
    if(ir==0 || obs->time[ir]!=obs->time[ir-1])
      fprintf(out, "\n");
    fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g", obs->time[ir],
	    obs->obsz[ir], obs->obslon[ir], obs->obslat[ir],
	    obs->vpz[ir], obs->vplon[ir], obs->vplat[ir],
	    obs->tpz[ir], obs->tplon[ir], obs->tplat[ir]);
    for(id=0; id<ctl->nd; id++)
      fprintf(out, " %g", obs->rad[id][ir]);
    for(id=0; id<ctl->nd; id++)
      fprintf(out, " %g", obs->tau[id][ir]);
    fprintf(out, "\n");
  }
  
  /* Close file... */
  fclose(out);
}

