#include "lineofsight.h"

/*****************************************************************************/

void raytrace(ctl_t *ctl,
	      atm_t *atm,
	      obs_t *obs,
	      aero_t *aero,
	      los_t *los,
	      int ir) {
  
  double cosa, d, dmax, dmin=0, ds, ex0[3], ex1[3], h=0.02, k[NWMAX],
    lat, lon, n, naux, ng[3], norm, p, q[NGMAX], 
    t, x[3], xh[3], xobs[3], xvp[3], z=1e99, zmax, zmin, zrefrac=60;
  /*zrefrac=25 for CRISTA-NF*/
  int i, ig, ip, iw, stop=0;

  /* Initialize... */
  los->np=0;
  los->tsurf=-999;
  obs->tpz[ir]=obs->vpz[ir];
  obs->tplon[ir]=obs->vplon[ir];
  obs->tplat[ir]=obs->vplat[ir];

  /* Get altitude range of atmospheric data... */
  zmin=gsl_stats_min(atm->z, 1, (size_t)atm->np);
  zmax=gsl_stats_max(atm->z, 1, (size_t)atm->np);
  
  /* Check observer altitude... */
  if(obs->obsz[ir]<zmin)
    ERRMSG("Observer below surface!");
  
  /* Check view point altitude... */
  if(obs->vpz[ir]>zmax-0.001)
    return;
  
  /* Determine Cartesian coordinates for observer and view point... */
  geo2cart(obs->obsz[ir], obs->obslon[ir], obs->obslat[ir], xobs);
  geo2cart(obs->vpz[ir], obs->vplon[ir], obs->vplat[ir], xvp);
  
  /* Determine initial tangent vector... */
  for(i=0; i<3; i++)
    ex0[i]=xvp[i]-xobs[i];
  norm=NORM(ex0);
  for(i=0; i<3; i++)
    ex0[i]/=norm;
  
  /* Observer within atmosphere... */
  for(i=0; i<3; i++)
    x[i]=xobs[i];
  
  /* Observer above atmosphere (search entry point)... */
  if(obs->obsz[ir]>zmax) {
    dmax=norm;
    while(fabs(dmin-dmax)>0.001) {
      d=(dmax+dmin)/2;
      for(i=0; i<3; i++)
	x[i]=xobs[i]+d*ex0[i];
      cart2geo(x, &z, &lon, &lat);
      if(z<=zmax && z>zmax-0.001)
	break;
      if(z<zmax-0.0005)
	dmax=d;
      else
	dmin=d;
    }
  }
  
  /* Ray-tracing... */
  while(1) {
    
    /* Set step length... */
    ds=ctl->rayds;
    if(ctl->raydz>0) {
      norm=NORM(x);
      for(i=0; i<3; i++)
	xh[i]=x[i]/norm;
      cosa=fabs(DOTP(ex0, xh));
      if(cosa!=0)
	ds=GSL_MIN(ctl->rayds, ctl->raydz/cosa);
    }
    
    /* Determine geolocation... */
    cart2geo(x, &z, &lon, &lat);
    
    /* Check if LOS hits the ground or has left atmosphere and save last los point. */
    if(z<zmin+0.001 || z>zmax+0.001) {
      stop=(z<zmin+0.001 ? 2 : 1);
      los->z[los->np] = z;
      los->lon[los->np]=lon;
      los->lat[los->np]=lat;
      intersection_point(ctl, atm, (z<zmin+0.001 ? &zmin : &zmax), los, los->np, los, los->np);
      los->ds[los->np]=0.;
    }
    
    /* Save first and middle los points. */
    if(stop==0) {
      intpol_atm_geo(ctl, atm, z, lon, lat, &p, &t, q, k);
      
      los->lon[los->np]=lon;
      los->lat[los->np]=lat;
      los->z[los->np]=z;
      los->p[los->np]=p;
      los->t[los->np]=t;
      for(ig=0; ig<ctl->ng; ig++)
	los->q[los->np][ig]=q[ig];
      for(iw=0; iw<ctl->nw; iw++)
	los->k[los->np][iw]=k[iw];
      los->ds[los->np]=ds;
    }

    /* Increment and check number of LOS points... */
    if((++los->np)>NLOS)
      ERRMSG("Too many LOS points!");

    /* Check stop flag... */
    if(stop) {
      los->tsurf=(stop==2 ? t : -999);
      break;
    }
    
    /* Determine refractivity... */
    if(ctl->refrac && z<=zrefrac)
      n=1+refractivity(p, t);
    else
      n=1;
    
    /* Construct new tangent vector (first term)... */
    for(i=0; i<3; i++)
      ex1[i]=ex0[i]*n;
    
    /* Compute gradient of refractivity... */
    if(ctl->refrac && z<=zrefrac) {
      for(i=0; i<3; i++)
	xh[i]=x[i]+0.5*ds*ex0[i];
      cart2geo(xh, &z, &lon, &lat);
      intpol_atm_geo(ctl, atm, z, lon, lat, &p, &t, q, k);
      n=refractivity(p, t);
      for(i=0; i<3; i++) {
	xh[i]+=h;
	cart2geo(xh, &z, &lon, &lat);
	intpol_atm_geo(ctl, atm, z, lon, lat, &p, &t, q, k);
	naux=refractivity(p, t);
	ng[i]=(naux-n)/h;
	xh[i]-=h;
      }
    } else
      for(i=0; i<3; i++)
	ng[i]=0;
    
    /* Construct new tangent vector (second term)... */
    for(i=0; i<3; i++)
      ex1[i]+=ds*ng[i];
    
    /* Normalize new tangent vector... */
    norm=NORM(ex1);
    for(i=0; i<3; i++)
      ex1[i]/=norm;
    
    /* Determine next point of LOS... */
    for(i=0; i<3; i++)
      x[i]+=0.5*ds*(ex0[i]+ex1[i]);
    
    /* Copy tangent vector... */
    for(i=0; i<3; i++)
      ex0[i]=ex1[i];
  }
  
  /* Check length of last segment... */
  if(los->ds[los->np-2]<1e-3 && los->np-1>1)
    los->np--;
  
  /* Get tangent point (to be done before changing segment lengths!)... */
  tangent_point(los, &obs->tpz[ir], &obs->tplon[ir], &obs->tplat[ir]);
  
  /* Change segment lengths according to trapezoid rule... */
  los->ds[0]=0.5*los->ds[0];
  for(ip=1; ip<los->np; ip++)
    los->ds[ip]=0.5*(los->ds[ip-1]+los->ds[ip]);

  /* Compute column density... */
  for(ip=0; ip<los->np; ip++)
    for(ig=0; ig<ctl->ng; ig++)
      los->u[ip][ig]=10*los->q[ip][ig]*los->p[ip]
	/(GSL_CONST_MKSA_BOLTZMANN*los->t[ip])*los->ds[ip];

  /* Add additional los points for aerosol layers and add aerosol data */
  if (ctl->sca_n > 0)
    add_aerosol_layers(ctl,atm,los,aero);
}

/*****************************************************************************/
void add_aerosol_layers(ctl_t *ctl,
			atm_t *atm,
			los_t *los,
			aero_t *aero){

  los_t *los_aero;

  double alti[4*NLMAX], altimax, altimin, x1[3], x2[3], x3[3], tt=0., epsilon=0.005; 
  /* deltatop=10., deltabot=10., */

  int il, ig, iw, jl=0, ip, it;

  size_t s;

  /* Allocate extended los... */
  ALLOC(los_aero, los_t, 1);

  /* Create altitudes to sample aerosol edges */
  for (il=0; il<aero->nl;il++){
    alti[jl] = aero->top[il] + epsilon;
    alti[jl+1] = aero->top[il] - epsilon;
    alti[jl+2] = aero->bottom[il] + epsilon;
    alti[jl+3] = aero->bottom[il] - epsilon;
    jl = jl+4;

    /* Create altitudes to sample transition layers */
    if (aero->trans[il] > ctl->transs) {
      tt = aero->trans[il] / ctl->transs;
      
      alti[jl] = aero->top[il] + aero->trans[il] + epsilon;
      alti[jl+1] = aero->top[il] + aero->trans[il] - epsilon;
      alti[jl+2] = aero->bottom[il] - aero->trans[il] + epsilon;
      alti[jl+3] = aero->bottom[il] - aero->trans[il] - epsilon;
      jl = jl+4;
      for (it=1; it<(int)tt; it++){
	alti[jl] = aero->top[il] + aero->trans[il] - epsilon - it*ctl->transs;
	jl++;
	alti[jl] = aero->bottom[il] - aero->trans[il] + epsilon + it*ctl->transs;
	jl++;
      }
    }  
  }

  /* Sort all altitudes from top-down */
  for (il=0; il<jl;il++)
    alti[il]=alti[il]*(-1.);
  gsl_sort(alti,1,(size_t)jl);
  for (il=0; il<jl;il++)
    alti[il]=alti[il]*(-1.);
  
  altimax = gsl_stats_max(alti, 1, (size_t)jl);
  altimin = gsl_stats_min(alti, 1, (size_t)jl);

  /* Copy los to new los and add additional points */
  los_aero->tsurf = los->tsurf;
  los_aero->z[0] = los->z[0];
  los_aero->lat[0] = los->lat[0];
  los_aero->lon[0] = los->lon[0];
  los_aero->ds[0] = los->ds[0];  
  los_aero->t[0] = los->t[0];
  los_aero->p[0] = los->p[0];
  for(ig=0; ig<ctl->ng; ig++)
    los_aero->q[0][ig] = los->q[0][ig];
  for(iw=0; iw<ctl->nw; iw++)
    los_aero->k[0][iw] = los->k[0][iw];
  los_aero->np = 1;

  for (ip=1; ip<los->np; ip++){
    
    /* add new los points around cloud edges */
    if ( (los->z[ip-1] < altimax || los->z[ip] < altimax) &&
	 (los->z[ip-1] > altimin || los->z[ip] > altimin) ) { 
      for (il=0; il<jl;il++){ /* loop over cloud edges */
	/* von oben */
	if(los->z[ip-1] > alti[il] && los->z[ip] < alti[il]){
	  intersection_point(ctl, atm, &alti[il], los, ip, los_aero, los_aero->np);
	  los_aero->np++; 
	}
	/* von unten */
	if(los->z[ip-1] < alti[jl-il-1] && los->z[ip] > alti[jl-il-1]){
	  intersection_point(ctl, atm, &alti[jl-il-1], los, ip, los_aero, los_aero->np);
	  los_aero->np++;
	}
      }
    }

    /* /\* check if current altitude is closer than 2m *\/ */
    /* /\* to any cloud top or bottom *\/ */
    /* deltatop = 10; */
    /* deltabot = 10; */
    /* for (il=0; il<aero->nl;il++){ */
    /*   deltatop = fabs(los->z[ip] - aero->top[il]); */
    /*   deltabot = fabs(los->z[ip] - aero->bottom[il]); */
    /*   if ( deltatop < epsilon*2. || deltabot < epsilon*2. ){ */
    /* 	continue; */
    /*   } */
    /* } */
    
    /* only copy old los points, if they are outside top||bottom +-2m */ 
    /* if ( deltatop >= epsilon*2. && deltabot >= epsilon*2. ) {  */
    /* copy old los points */
    los_aero->z[los_aero->np] = los->z[ip];
    los_aero->lat[los_aero->np] = los->lat[ip];
    los_aero->lon[los_aero->np] = los->lon[ip];
    los_aero->t[los_aero->np] = los->t[ip];
    los_aero->p[los_aero->np] = los->p[ip];
    for(ig=0; ig<ctl->ng; ig++)
      los_aero->q[los_aero->np][ig] = los->q[ip][ig];
    for(iw=0; iw<ctl->nw; iw++)
      los_aero->k[los_aero->np][iw] = los->k[ip][iw];
    
    /* Increment and check number of new LOS points */
    if((los_aero->np++)>NLOS)
      ERRMSG("Too many LOS points!");
    /* } */
  }

  /* Compute segment length following trapezoidal rule */
  geo2cart(los_aero->z[0], los_aero->lon[0],los_aero->lat[0], x1);
  geo2cart(los_aero->z[1], los_aero->lon[1],los_aero->lat[1], x2);
  los_aero->ds[0]= 0.5 * (DIST(x1,x2));
  for(ip=1; ip<los_aero->np-1; ip++){
    geo2cart(los_aero->z[ip-1], los_aero->lon[ip-1],los_aero->lat[ip-1], x1);
    geo2cart(los_aero->z[ip], los_aero->lon[ip],los_aero->lat[ip], x2);
    geo2cart(los_aero->z[ip+1], los_aero->lon[ip+1],los_aero->lat[ip+1], x3);
    los_aero->ds[ip] = 0.5 * (DIST(x1,x2) + DIST(x2,x3));
  }
  geo2cart(los_aero->z[los_aero->np-1], los_aero->lon[los_aero->np-1], 
	   los_aero->lat[los_aero->np-1], x1);
  geo2cart(los_aero->z[los_aero->np-2], los_aero->lon[los_aero->np-2], 
	   los_aero->lat[los_aero->np-2], x2);
  los_aero->ds[los_aero->np-1] = 0.5 * (DIST(x1,x2));

  /* add aerosol/cloud information and column density u to new los  */
  for (ip=0; ip<los_aero->np; ip++){

    /* Compute column density... */
    for(ig=0; ig<ctl->ng; ig++)
      los_aero->u[ip][ig]=10*los_aero->q[ip][ig]*los_aero->p[ip]
	/(GSL_CONST_MKSA_BOLTZMANN*los_aero->t[ip])*los_aero->ds[ip];

    /* Get aerosol/cloud layer id and factor */
    los_aero->aeroi[ip] = -999;
    los_aero->aerofac[ip] = 0.;
    if ( (los_aero->z[ip-1] < altimax || los_aero->z[ip] < altimax) &&
	 (los_aero->z[ip-1] > altimin || los_aero->z[ip] > altimin) ) { 
      for (il=0; il<aero->nl;il++){
        /* Aerosol info within layer centre */
	if (los_aero->z[ip] <= aero->top[il] && 
	    los_aero->z[ip] >= aero->bottom[il]){
	  los_aero->aeroi[ip] = il;
	  los_aero->aerofac[ip] = 1.;	  
	}
	/* Aerosol info in transition region */
	if (aero->trans[il] > ctl->transs &&
	    los_aero->z[ip] <= (aero->top[il] + aero->trans[il]) &&
	    los_aero->z[ip] > aero->top[il]){
	  los_aero->aeroi[ip] = il;
	  los_aero->aerofac[ip] = (aero->top[il] + aero->trans[il] - los_aero->z[ip])/
	    aero->trans[il];
	}
	if (aero->trans[il] > ctl->transs &&
	    los_aero->z[ip] < aero->bottom[il] &&
	    los_aero->z[ip] >= (aero->bottom[il] - aero->trans[il])){
	  los_aero->aeroi[ip] = il;
	  los_aero->aerofac[ip] = fabs(aero->bottom[il]-aero->trans[il]-los_aero->z[ip])/
	    aero->trans[il];
	}
      }
    }
  }

  /* Copy los */
  /*   *los = *los_aero; */
  s=sizeof(los_t);
  memcpy(los, los_aero, s);
   
  /* Free help los... */
  free(los_aero);
}

/*****************************************************************************/

double refractivity(double p,
		    double t) {
  
  /* Refractivity of air at 4 to 15 micron... */
  return 7.753e-05*p/t;
}

/*****************************************************************************/

void intersection_point(ctl_t *ctl,
			atm_t *atm,
			double *znew,
			los_t *los,
			int ip,
			los_t *los_aero,
			int jp){
  
  double frac, x1[3], x2[3];
  int i ;
  
  frac = (los->z[ip-1] - *znew) / (los->z[ip-1] - los->z[ip]);
  geo2cart(los->z[ip-1], los->lon[ip-1], los->lat[ip-1], x1);
  geo2cart(los->z[ip], los->lon[ip], los->lat[ip], x2);

  for(i=0; i<3; i++)
    x2[i]=x1[i]+frac*(x2[i]-x1[i]);
  /* get new coordinates */
  cart2geo(x2, &los_aero->z[jp], &los_aero->lon[jp], &los_aero->lat[jp]);
  /* get atmosphere parameters */
  intpol_atm_geo(ctl, atm, los_aero->z[jp], los_aero->lon[jp], los_aero->lat[jp],
		 &los_aero->p[jp], &los_aero->t[jp], los_aero->q[jp], los_aero->k[jp]);
}

/*****************************************************************************/

void tangent_point(los_t *los,
		   double *tpz,
		   double *tplon,
		   double *tplat) {
  
  double a, b, c, dummy, v[3], v0[3], v2[3], x, x1, x2, yy0, yy1, yy2;
  
  size_t i, ip;
  
  /* Find minimum altitude... */
  ip=gsl_stats_min_index(los->z, 1, (size_t)los->np);
  
  /* Nadir or zenith... */
  if(ip<=0 || ip>=(size_t)los->np-1) {
    *tpz=los->z[los->np-1];
    *tplon=los->lon[los->np-1];
    *tplat=los->lat[los->np-1];
  }
  
  /* Limb... */
  else {
    
    /* Determine interpolating polynomial y=a*x^2+b*x+c... */
    yy0=los->z[ip-1];
    yy1=los->z[ip];
    yy2=los->z[ip+1];
    x1=sqrt(gsl_pow_2(los->ds[ip])-gsl_pow_2(yy1-yy0));
    x2=x1+sqrt(gsl_pow_2(los->ds[ip+1])-gsl_pow_2(yy2-yy1));
    a=1/(x1-x2)*(-(yy0-yy1)/x1+(yy0-yy2)/x2);
    b=-(yy0-yy1)/x1-a*x1;
    c=yy0;
    
    /* Get tangent point location... */
    x=-b/(2*a);
    *tpz=a*x*x+b*x+c;
    geo2cart(los->z[ip-1], los->lon[ip-1], los->lat[ip-1], v0);
    geo2cart(los->z[ip+1], los->lon[ip+1], los->lat[ip+1], v2);
    for(i=0; i<3; i++)
      v[i]=LIN(0.0, v0[i], x2, v2[i], x);
    cart2geo(v, &dummy, tplon, tplat);
  }
}
