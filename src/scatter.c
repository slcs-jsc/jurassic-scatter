#include "scatter.h"

/*****************************************************************************/

void bascoord(double *dz,
	      double *dy,
	      double *ex,
	      double *ey,
	      double *ez) {
  
  double dotp, norm;
  
  int i;
  
  /* Set first orthonormal vector... */
  norm=NORM(dz);
  for(i=0; i<3; i++)
    ez[i]=dz[i]/norm;
  
  /* Normalize second direction... */
  norm=NORM(dy);
  for(i=0; i<3; i++)
    ey[i]=dy[i]/norm;
  
  /* Avoid linear dependences... */
  dotp=DOTP(ey, ez);
  if(dotp==1) {
    ey[0]=1;
    ey[1]=0;
    ey[2]=0;
    dotp=ez[0];
    if(dotp==1) {
      ey[0]=0;
      ey[1]=1;
      dotp=ez[1];
    }
  }
  
  /* Set second orthonormal vector (Gram-Schmidt)... */
  for(i=0; i<3; i++)
    ey[i]-=dotp*ez[i];
  
  norm=NORM(ey);
  for(i=0; i<3; i++)
    ey[i]/=norm;
  
  /* Get third orthonormal vector (cross-product)... */
  ex[0]=ey[1]*ez[2]-ey[2]*ez[1];
  ex[1]=ey[2]*ez[0]-ey[0]*ez[2];
  ex[2]=ey[0]*ez[1]-ey[1]*ez[0];
}

/*****************************************************************************/

void bhmie(double x,
	   double n_real,
	   double n_imag,
	   double *phase,
	   double *qext,
	   double *qsca) {
  
  gsl_complex cxan, cxbn, cxxi, cxy, cxxi1, cxtemp, cxref, 
    cxs1[NTHETA+2], cxs2[NTHETA+2], cxd[10000];
  
  double apsi, apsi1, chi, chi0, chi1, dang, fn, p, rn, t, theta2,
    xstop, ymod, dn, dx, psi, psi0, psi1, amu[NTHETA+2], pi[NTHETA+2],
    pi0[NTHETA+2], pi1[NTHETA+2], tau[NTHETA+2];
  
  int i, j, jj, n, nmx, nn, nstop, ntheta;
  
  /* Set scattering angles, ntheta=(NTHETA+1)/2... */
  if((NTHETA+1)%2!=0)
    ERRMSG("NTHETA needs to be odd!");
  ntheta=(NTHETA+1)/2;
  
  /* Bohren-Huffman Mie code... */
  cxref=gsl_complex_rect(n_real, n_imag);
  dx=x;
  cxy=gsl_complex_mul(gsl_complex_rect(x, 0.0), cxref);
  
  /* Series expansion terminated after NSTOP terms */  
  xstop=x+4.E0*pow(x,0.3333)+2.0;
  nstop=(int)xstop;
  ymod = gsl_complex_abs(cxy);
  nmx=(int)((xstop>ymod ? xstop : ymod)+15);
  if(nmx>10000)
    ERRMSG("Too many Mie terms!");
  dang=.5E0*M_PI/(double)(ntheta-1);
  for(j=1; j<=ntheta; j++) {
    theta2=(double)(j-1)*dang;
    amu[j]=cos(theta2);
  }

  /* Logarithmic derivative D(J) calculated by downward recurrence */
  /*  beginning with initial value (0.,0.) at J=NMX */ 
  cxd[nmx] = gsl_complex_rect(0.E0, 0.E0);
  nn=nmx-1;
  for(n=1; n<= nn; n++) {
    rn=nmx-n+1;
    cxtemp=gsl_complex_add(cxd[nmx-n+1],
			   gsl_complex_div(gsl_complex_rect(rn, 0.0), cxy));
    cxtemp=gsl_complex_div(gsl_complex_rect(1.0, 0.0), cxtemp);
    cxd[nmx-n]
      =gsl_complex_sub(gsl_complex_div(gsl_complex_rect(rn, 0.0), cxy), cxtemp);
  }
  
  for(j=1; j<=ntheta; j++) {
    pi0[j]=0.E0;
    pi1[j]=1.E0;
  }
  
  nn=2*ntheta-1;
  for(j=1; j<=nn; j++) {
    cxs1[j] = gsl_complex_rect(0.E0, 0.E0);
    cxs2[j] = gsl_complex_rect(0.E0, 0.E0);
  }

  /* Riccati-Bessel functions with real argument X */
  /*  calculated by upward recurrence */
  psi0=cos(dx);
  psi1=sin(dx);
  chi0=-sin(x);
  chi1=cos(x);
  apsi1=psi1;
  cxxi1 = gsl_complex_rect(apsi1,-chi1);
  *qsca=0.E0;

  for(n=1; n<=nstop; n++) {  
    dn=n;
    rn=n;
    fn=(2.E0*rn+1.E0)/(rn*(rn+1.E0));
    psi=(2.E0*dn-1.E0)*psi1/dx-psi0;
    apsi=psi;
    chi=(2.E0*rn-1.E0)*chi1/x-chi0;
    cxxi=gsl_complex_rect(apsi,-chi);

    cxan=gsl_complex_div(cxd[n], cxref);
    cxan=gsl_complex_add(cxan, gsl_complex_rect(rn/x, 0.0));
    cxan=gsl_complex_mul(cxan, gsl_complex_rect(apsi, 0.0));
    cxan=gsl_complex_sub(cxan, gsl_complex_rect(apsi1, 0.0));

    cxtemp=gsl_complex_div(cxd[n], cxref);
    cxtemp=gsl_complex_add(cxtemp, gsl_complex_rect(rn/x, 0.0));
    cxtemp=gsl_complex_mul(cxtemp, cxxi);
    cxtemp=gsl_complex_sub(cxtemp, cxxi1);
    cxan=gsl_complex_div(cxan, cxtemp);

    cxbn=gsl_complex_mul(cxref, cxd[n]);
    cxbn=gsl_complex_add(cxbn, gsl_complex_rect(rn/x, 0.0));
    cxbn=gsl_complex_mul(cxbn, gsl_complex_rect(apsi, 0.0));
    cxbn=gsl_complex_sub(cxbn, gsl_complex_rect(apsi1, 0.0));

    cxtemp=gsl_complex_mul(cxref, cxd[n]);
    cxtemp=gsl_complex_add(cxtemp, gsl_complex_rect(rn/x, 0.0));
    cxtemp=gsl_complex_mul(cxtemp, cxxi);
    cxtemp=gsl_complex_sub(cxtemp, cxxi1);
    cxbn=gsl_complex_div(cxbn, cxtemp);
 
    *qsca=*qsca+(2.*rn+1.)
      *(gsl_complex_abs(cxan)*gsl_complex_abs(cxan)
	+gsl_complex_abs(cxbn)*gsl_complex_abs(cxbn));
    for(j=1; j<=ntheta; j++) {
      jj=2*ntheta-j;
      pi[j]=pi1[j];
      tau[j]=rn*amu[j]*pi[j]-(rn+1.E0)*pi0[j];
      p=pow(-1.0,n-1);
      
      cxtemp=gsl_complex_mul(cxan, gsl_complex_rect(pi[j], 0.0));
      cxtemp
	=gsl_complex_add(cxtemp, gsl_complex_mul(cxbn, gsl_complex_rect(tau[j],
									0.0)));
      cxtemp=gsl_complex_mul(gsl_complex_rect(fn, 0.0), cxtemp);
      cxs1[j]=gsl_complex_add(cxs1[j], cxtemp);
      
      t=pow(-1.0,n);
      
      cxtemp=gsl_complex_mul(cxan, gsl_complex_rect(tau[j], 0.0));
      cxtemp
	=gsl_complex_add(cxtemp, gsl_complex_mul(cxbn, gsl_complex_rect(pi[j],
									0.0)));
      cxtemp=gsl_complex_mul(gsl_complex_rect(fn, 0.0), cxtemp);
      cxs2[j]=gsl_complex_add(cxs2[j], cxtemp);
      
      if(j!=jj) {
	cxtemp=gsl_complex_mul(cxan, gsl_complex_rect(pi[j]*p, 0.0));
	cxtemp
	  =gsl_complex_add(cxtemp,
			   gsl_complex_mul(cxbn,
					   gsl_complex_rect(tau[j]*t, 0.0)));
	cxtemp=gsl_complex_mul(gsl_complex_rect(fn, 0.0), cxtemp);
	cxs1[jj]=gsl_complex_add(cxs1[jj], cxtemp);
	cxtemp=gsl_complex_mul(cxan, gsl_complex_rect(tau[j]*t, 0.0));
	cxtemp
	  =gsl_complex_add(cxtemp,
			   gsl_complex_mul(cxbn,
					   gsl_complex_rect(pi[j]*p, 0.0)));
	cxtemp=gsl_complex_mul(gsl_complex_rect(fn, 0.0), cxtemp);
	cxs2[jj]=gsl_complex_add(cxs2[jj], cxtemp);
      }
    }
    
    psi0=psi1;
    psi1=psi;
    apsi1=psi1;
    chi0=chi1;
    chi1=chi;
    cxxi1=gsl_complex_rect(apsi1, -chi1);
    
    /* For each angle J, compute pi_n+1 from PI = pi_n , PI0 = pi_n-1 */
    for(j=1; j<=ntheta; j++) {
      pi1[j]=((2.*rn+1.)*amu[j]*pi[j]-(rn+1.)*pi0[j])/rn;
      pi0[j]=pi[j];
    }
  }
  
  /* Compute efficiencies... */
  *qsca=(2.E0/(x*x))*(*qsca);
  *qext=(4.E0/(x*x))*cxs1[1].dat[0];
  
  /* Compute phase function from scattering amplitudes... */
  /* calculate phase function following Liou: p192 eq 5.2.111a */
  /* P_11 = 4*PI*(i_1+i_2)/(2*k^2*sigma_s)) */
  /* i_1(theta), i_2(theta) = abs(S_1(theta))^2, abs(S_2(theta))^2 */
  /* intensity functions for perpendicular and parallel components */
  /* Q_s = sigma_s/(PI *rad^2) - scattering efficiency (here qsca) */
  /* sigma_s = scattering cross section */
  /* x = k*rad - size parameter (x=2*PI*rad/lambda) */

  for(i=0; i<2*ntheta-1; i++)
    phase[i]=2/(x*x*(*qsca))
      *(cxs1[i+1].dat[0]*cxs1[i+1].dat[0]+cxs1[i+1].dat[1]*cxs1[i+1].dat[1]
	+cxs2[i+1].dat[0]*cxs2[i+1].dat[0]+cxs2[i+1].dat[1]*cxs2[i+1].dat[1]);
}

/*****************************************************************************/

void gauher(double *x,
	    double *w){

  /* Calculate abcissas (x) and weights (w) for Gauss-Hermite quadrature. */
  /* This routine basically follows the Numerical Recipes. */
  /* W. H. Press and S. A. Teukolsky and W. T. Vetterling and B. P. Flannery: */
  /* Numerical Recipes in C, 3rd edition, Cambridge University Press, 2007 */

  const double EPS=1.0e-14, PIM4=0.7511255444649425;

  const int MAXIT=10000;

  int i,its,j,m;

  double p1,p2,p3,pp,z=0,z1;

  int n=NRAD;

  m = (n+1)/2;
  for (i=0; i<m; i++) {
    if (i == 0) {
      z = sqrt((double)(2*n+1))-1.85575*pow((double)(2*n+1),-0.16667);
    } else if (i == 1) {
      z -= 1.14*pow((double)n,0.426)/z;
    } else if (i == 2) {
      z = 1.86*z-0.86*x[0];
    } else if (i == 3) {
      z = 1.91*z-0.91*x[1];
    } else {
      z = 2.0*z-x[i-2];
    }
    for (its=0; its<MAXIT; its++) {
      p1 = PIM4;
      p2 = 0.0;
      for (j=0; j<n; j++) {
	p3 = p2;
	p2 = p1;
	p1 = z*sqrt(2.0/((double)j+1))*p2-sqrt((double)j/(j+1.))*p3;
      }
      pp = sqrt(2.*(double)n)*p2;
      z1 = z;
      z = z1-p1/pp;
      if (fabs(z-z1) <= EPS) break;
    }
    if (its >= MAXIT) {
      ERRMSG("Too many iterations in gauher.");
    }
    x[i] = z;
    x[n-1-i] = -z;
    w[i] = 2.0/(pp*pp);
    w[n-1-i] = w[i];
  }
}

/*****************************************************************************/
void get_opt_prop(ctl_t *ctl,
		  aero_i *aeroin,
		  aero_t *aero){

  int nl=1, nm=1, count=0;
  int ii,jj, id, itheta;  
  double mbeta_e[NDMAX], mbeta_s[NDMAX], mp[NDMAX][NTHETA];
 
  /* check input data */
  /* ToDo: improve check and sort data */
  if(aeroin->top[0]<=aeroin->bottom[0])
    ERRMSG("Aerosol top altitude is smaller than bottom altitude. Please check aero.tab.");

  aero->top[0] = aeroin->top[0];
  aero->bottom[0] = aeroin->bottom[0];
  aero->trans[0] = aeroin->trans[0];
  
  for (ii=1; ii<aeroin->nm; ii++){
    if(aeroin->top[ii]>aeroin->top[ii-1] ||
       aeroin->top[ii]<=aeroin->bottom[ii] ||
       aeroin->bottom[ii]>aeroin->bottom[ii-1])
      ERRMSG("Aerosol/Cloud altitudes and/or transition layers are wrong. Please check aero.tab.");
    /* Identify number of aerosol/cloud layers and set top, */
    /* bottom and transition layer */
    /* Check if input is already sorted list from top to bottom; */
    if(aeroin->top[ii]!=aeroin->top[ii-1]){
      aero->top[nl] = aeroin->top[ii];
      aero->bottom[nl] = aeroin->bottom[ii];
      aero->trans[nl] = aeroin->trans[ii];
      aero->nm[nl-1] = nm;
      nl++;
      nm=0;
    }
    nm++;
  }
  aero->nm[nl-1] = nm;
  aero->nl = nl;
  
  /* Get optical properties for each layer. */
  count=0;
  for (ii=0; ii<aero->nl; ii++){

    /* Initialise each layer */
    for(id=0; id<ctl->nd; id++){
      aero->beta_e[ii][id] = 0.;
      aero->beta_a[ii][id] = 0.;
      aero->beta_s[ii][id] = 0.;
      for(itheta=0; itheta<NTHETA; itheta++){
	aero->p[ii][id][itheta] = 0.;
	mp[id][itheta] = 0.;
      }
      mbeta_e[id] = 0.;
      mbeta_s[id] = 0.;
    }

    /* Get optical properties for each mode. */
    for (jj=0; jj<aero->nm[ii]; jj++){

      if(strcasecmp(aeroin->type[count], "MIE")==0){
    	/* Get optical properties for log-normal mode using Mie theory. */ 
	/* Gauss-Hermite integration */
	opt_prop_mie_log(ctl, aeroin, count, mbeta_e, mbeta_s, mp);
      } 
      else if(strcasecmp(aeroin->type[count], "Ext")==0){ 
	/* Get optical properties from external data base. Selects properties from closest wavenumber in data base file. */
	opt_prop_external(ctl, aeroin, count, mbeta_e, mbeta_s, mp);
      } else if(strcasecmp(aeroin->type[count], "Const")==0){ 
    	printf("Using constant extinction [1/km]: %g\n", aeroin->nn[count]);
    	ERRMSG("Implement me!");
      }
      else {
 	ERRMSG("Please give valid scattering model (MIE, Ext, Const)!");
      }
      
      /* Sum up optical properties for each layer */
      for(id=0; id<ctl->nd; id++){
    	aero->beta_e[ii][id] += mbeta_e[id];
    	aero->beta_s[ii][id] += mbeta_s[id];
     	aero->beta_a[ii][id] += (mbeta_e[id] - mbeta_s[id]);
    	for(itheta=0; itheta<NTHETA; itheta++){
    	  aero->p[ii][id][itheta] += mp[id][itheta];
	  mp[id][itheta] = 0.;
	}
	mbeta_e[id] = 0.;
	mbeta_s[id] = 0.;
      }
      count++;
    }
  }
}

/*****************************************************************************/

void opt_prop_mie_log(ctl_t *ctl,
		     aero_i *aeroin,
		     int count,
		     double *beta_ext,
		     double *beta_sca,
		     double phase[NDMAX][NTHETA]){

  FILE *in;
  
  char line[LEN];

  static int init=0;

  static double nu[REFMAX], nr[REFMAX], ni[REFMAX], n_imag[NDMAX], n_real[NDMAX], 
    rad_min=0.001, rad_max=1000., weights[NRAD], zs[NRAD];

  int npts=0, id, idx, nn, jj;

  double K1, rad, lambda, x, qext, qsca, qphase[NTHETA];

  /* Read and interpolate refractive indices... */
  /* Check if previous mode has the same refractive index */
  if(count==0 || strcmp(aeroin->filepath[count], aeroin->filepath[count-1])!=0) { 
    
    /* Read data... */
    printf("Read refractive indices: %s\n", aeroin->filepath[count]);
    if(!(in=fopen(aeroin->filepath[count], "r")))
      ERRMSG("Cannot open file!");
    while(fgets(line, LEN, in))
      if(sscanf(line, "%lg %lg %lg", &nu[npts], &nr[npts], &ni[npts])==3)
	if((++npts)>REFMAX)
	  ERRMSG("Too many data points!");
    fclose(in);
    
    /* Interpolate... */
    for(id=0; id<ctl->nd; id++) {
      idx=locate(nu, npts, ctl->nu[id]);
      n_real[id]=LIN(nu[idx], nr[idx], nu[idx+1], nr[idx+1], ctl->nu[id]);
      n_imag[id]=LIN(nu[idx], ni[idx], nu[idx+1], ni[idx+1], ctl->nu[id]);
    }
  } 

  /* Check log-normal parameters... */
  if(aeroin->nn[count]<=0 || aeroin->rr[count]<=0 || aeroin->ss[count]<=1)
    ERRMSG("The log-normal parameters are nonsense. ((?_?)) ");
    
  /* Integrate Mie parameters over log-normal mode */
  if(!init) {
    init=1;  
    /* get abcissas and weights for Gauss-Hermite quadrature */
    gauher(zs, weights); 
  }

  /* set coefficient */
  K1 = aeroin->nn[count] * 1e-3 * sqrt(M_PI);

  /* sum up Gaussian nodes */
  for (nn=0; nn<NRAD; ++nn) {
    rad = exp(sqrt(2) * log(aeroin->ss[count]) * zs[nn] + log(aeroin->rr[count]));
    if (rad >= rad_min && rad <= rad_max && K1 > 0.) {
 
     for(id=0; id<ctl->nd; id++){

	/* size parameter */
	lambda = 1./(ctl->nu[id])* pow(10,4.);
	x = 2*M_PI*rad/lambda;

	/* evaluate Mie Code at the nodes */
	bhmie(x, n_real[id], n_imag[id], qphase, &qext, &qsca);
	/* bhmie(rad, wavn, nang); */

	beta_ext[id] += K1 * pow(rad,2) *  qext * weights[nn];
	beta_sca[id] += K1 * pow(rad,2) *  qsca * weights[nn];
	    
	for (jj=0; jj<NTHETA; ++jj)
	  phase[id][jj] += K1 * qsca * pow(rad,2) * qphase[jj] * weights[nn];

      } 
    } 
  }
  
  /* Weight phase function with beta_s */
  for(id=0; id<ctl->nd; id++){
    for (jj=0; jj<NTHETA; ++jj)
      phase[id][jj] /= beta_sca[id];
  }
}

/*****************************************************************************/

void opt_prop_external(ctl_t *ctl,
		       aero_i *aeroin,
		       int count,
		       double *beta_ext,
		       double *beta_sca,
		       double phase[NDMAX][NTHETA]){

  FILE *in;
  
  char line[LEN], *tok; 

  static int init=0;

  static double nu[REFMAX], n_ext[REFMAX], n_sca[REFMAX], n_phase[REFMAX][NTHETA];

  int npts=0, ia, id, im;

  /* evtl. später Interpolation zw. Wellenlängen als Option einbauen */

  /* Read optical properties and find closest match to each wavenumber */
  if(aeroin->nn[count]==0) {
  
    /* Check if previous mode has the same optical properties */
    if(count==0 || strcmp(aeroin->filepath[count], aeroin->filepath[count-1])!=0) {
      
      /* Check for file... */
      printf("Read non-spherical optical properties: %s\n", aeroin->filepath[count]);
      if(!(in=fopen(aeroin->filepath[count], "r")))
	ERRMSG("Cannot open file!");
      
      /* Read data... */
      while(fgets(line, LEN, in)) {
    	
	TOK(line, tok, "%lg", nu[npts]);  
	TOK(NULL, tok, "%lg", n_ext[npts]);
	TOK(NULL, tok, "%lg", n_sca[npts]);
	for(ia=0; ia<NTHETA; ia++)
	  TOK(NULL, tok, "%lg", n_phase[npts][ia]);
	
	if((++npts)>REFMAX)
	  ERRMSG("Too many data points!");
      }
      
      /* Close file... */
      fclose(in);
      
      /* Check number of points... */
      if(npts<1)
	ERRMSG("Could not read any data!");
    }

    /* Find closest match in wavenumber for each spectral point */
    for(id=0; id<ctl->nd; id++){
      im=locate(nu, npts, ctl->nu[id]);
      if(im != npts && abs(nu[im] - ctl->nu[id]) > abs(nu[im+1] - ctl->nu[id]))
	im=im+1;
      
      beta_ext[id] = n_ext[im];
      beta_sca[id] = n_sca[im];
      for (ia=0; ia<NTHETA; ++ia)
      	phase[id][ia] = n_phase[im][ia];
    }
  }

  /* Read optical properties and interpolate to wavenumber */
  if(aeroin->nn[count]==1) {
    ERRMSG("Implement me!");
  }
}

/*****************************************************************************/

void read_aero(const char *dirname,
	       const char *filename,
	       ctl_t *ctl,
	       aero_i *aeroin){

  FILE *in;
  
  char file[LEN], line[LEN], *tok;
  
  /* Init... */
  aeroin->nm=0;
  
  /* Set filename... */
  if(dirname!=NULL)
    sprintf(file, "%s/%s", dirname, filename);
  else
    sprintf(file, "%s", filename);
  
  /* Write info... */
  printf("Read aerosol data: %s\n", file);
  
  /* Open file... */
  if(!(in=fopen(file, "r")))
    ERRMSG("Cannot open file!");
  
  /* Read data... */
  while(fgets(line, LEN, in)) {
    
    /* Read data... */
    TOK(line, tok, "%lg", aeroin->top[aeroin->nm]);
    TOK(NULL, tok, "%lg", aeroin->bottom[aeroin->nm]);
    TOK(NULL, tok, "%lg", aeroin->trans[aeroin->nm]);
    TOK(NULL, tok, "%s",  aeroin->type[aeroin->nm][0]);
    TOK(NULL, tok, "%s",  aeroin->filepath[aeroin->nm][0]); 
    TOK(NULL, tok, "%lg", aeroin->nn[aeroin->nm]); 
    TOK(NULL, tok, "%lg", aeroin->rr[aeroin->nm]);
    TOK(NULL, tok, "%lg", aeroin->ss[aeroin->nm]);
    
    /* Increment counter... */
    if((++aeroin->nm)>SCAMOD)
      ERRMSG("Too many aerosol models!");
  }
  
  /* Close file... */
  fclose(in);
  
  /* Check number of points... */
  if(aeroin->nm<1)
    ERRMSG("Could not read any data!");

  /* Check consistency */
  if(ctl->sca_n != aeroin->nm)
    ERRMSG("Number of scattering models in control file and aerosol file does not match.");

}

/*****************************************************************************/

void srcfunc_sca(ctl_t *ctl,
		 atm_t *atm,
		 aero_t *aero,
		 double sec,
		 double *x,
		 double *dx,
		 int il,
		 double *src_sca,
		 int scattering) {
  
  /* Compute scattering of thermal radiation... */
  if(ctl->ip==1)
    srcfunc_sca_1d(ctl, atm, aero, x, dx, il, src_sca, scattering);
  else
    srcfunc_sca_3d(ctl, atm, aero, x, dx, il, src_sca, scattering);
  
  /* Compute scattering of solar radiation... */
  if(TSUN>0)
    srcfunc_sca_sun(ctl, atm, aero, sec, x, dx, il, src_sca);
}

/*****************************************************************************/

void srcfunc_sca_1d(ctl_t *ctl,
		    atm_t *atm,
		    aero_t *aero,
		    double *x,
		    double *dx,
		    int il,
		    double *src_sca,
		    int scattering) {
  
  obs_t *obs2;
  
  double alpha[NTHETA], alpha2, dnorth[3], ek[3], lx[3], ly[3], lz[3], phi,
    phase2, rad, sx[3], sy[3], sz[3], theta[NTHETA], theta2, w=0, wsum[NDMAX], xv[3];
  
  int i, id, idp, idx, iphi, ir, itheta, nalpha=28, nphi=180, ntheta2=180;

  int n1=1, n2=2;
  double midang=83, up=92, down=81, step=1+2;
  
  /* Allocate... */
  ALLOC(obs2, obs_t, 1);
  
  /* Set scattering phase function angles... */
  for(itheta=0; itheta<NTHETA; itheta++)
    theta[itheta]=(double)itheta*M_PI/180.;
  
  /* Get local coordinate system... */
  dnorth[0]=-x[0];
  dnorth[1]=-x[1];
  dnorth[2]=2*RE-x[2];
  bascoord(x, dnorth, lx, ly, lz);
  
  /* Set angles - tested version with nalpha=28 and Fibonacci Numbers */
  alpha[0] = 0;
  for (ir=nalpha/2-4; ir<nalpha/2+4; ++ir){
    alpha[ir] = midang*M_PI/180.;
    midang++;
  }
  for (ir=0; ir<nalpha/2-5; ++ir){
    alpha[nalpha/2+ir+4] = up*M_PI/180.;
    alpha[nalpha/2-ir-5] = down*M_PI/180.;
    up+=step;
    down-=step;
    if (step < 13){
      step = (double)n1 + (double)n2;
      n1 = n2;
      n2 = (int)step;
    } 
  }
  alpha[nalpha-1] = M_PI;

  /* Get incident radiation... */
  /* nalpha=181; */
  for(ir=0; ir<nalpha; ir++) {
    
    /* Set angle... */
    /* traditional 0-180 deg in 1 deg steps with nalpha=181 */
    /* alpha[ir] = ir*M_PI/180.; */

    /* initial version with nalpha=21 */
    /* alpha[ir]=acos(2*(double)ir/(nalpha-1.0)-1.0); */
  
    /* Set view point... */
    for(i=0; i<3; i++)
      /* xv[i]=x[i]+10.*cos(alpha[ir])*lz[i]+10.*sin(alpha[ir])*ly[i];  */
      xv[i]=x[i]+10.*cos(alpha[ir])*(-1)*lz[i]+10.*sin(alpha[ir])*ly[i];
    
    /* Set observation geometry... */
    obs2->nr=nalpha;
    cart2geo(x, &obs2->obsz[ir], &obs2->obslon[ir], &obs2->obslat[ir]);
    cart2geo(xv, &obs2->vpz[ir], &obs2->vplon[ir], &obs2->vplat[ir]);
    
    /* Get pencil beam radiance... */
    formod_pencil(ctl, atm, obs2, aero, scattering-1, ir);
  }
  
  /* Get orthonormal basis (with respect to LOS)... */
  bascoord(dx, x, sx, sy, sz);  
  
  /* Initialize... */
  for(id=0; id<ctl->nd; id++){
    src_sca[id]=0;
    wsum[id]=0;
  }
  
  /* Loop over phase function angles... */
  for(itheta=0; itheta<ntheta2; itheta++) {
    
    /* Set phase function angle in 1° steps (avoid 0 and 180°)... */
    theta2=(0.5+itheta)*M_PI/180.;
    
    /* Loop over azimuth angles... */
    for(iphi=0; iphi<nphi; iphi++) {
      
      /* Set azimuth angle in 2° steps... */
      phi=2.*(0.5+iphi)*M_PI/180.;
      
      /* Get unit vector on sphere... */
      for(i=0; i<3; i++)
	ek[i]
	  =sin(theta2)*sin(phi)*sx[i]
	  +sin(theta2)*cos(phi)*sy[i]
	  +cos(theta2)*sz[i];

      /* Get phase function index */
      idp=locate(theta, NTHETA, theta2);

      /* Get zenith angle... */
      alpha2=ANGLE(-1.*lz, ek);

      /* Get source ray angle index */
      idx=locate(alpha, nalpha, alpha2);
	
      /* Loop over channels... */
      for(id=0; id<ctl->nd; id++) {
	
	/* Interpolate phase function... */
	phase2=LIN(theta[idp], aero->p[il][id][idp],
		   theta[idp+1], aero->p[il][id][idp+1], theta2);

	/* Get weighting factor (area of surface element * phase function)... */
	w=sin(theta2)*phase2;
	
	/* Interpolate radiance to particular angle... */
	rad=LIN(alpha[idx], obs2->rad[id][idx],
	 	alpha[idx+1], obs2->rad[id][idx+1], alpha2);
	
	/* Integrate... */
	src_sca[id]+=w*rad;
	wsum[id]+=w;
      }
    }
  }
  
  /* Normalize... */
  for(id=0; id<ctl->nd; id++)
    src_sca[id]/=wsum[id]; 
  
  /* Free... */
  free(obs2);
}

/*****************************************************************************/

void srcfunc_sca_3d(ctl_t *ctl,
		    atm_t *atm,
		    aero_t *aero,
		    double *x,
		    double *dx,
		    int il,
		    double *src_sca,
		    int scattering) {
  
  obs_t *obs2;
  
  double phi, phase2, sx[3], sy[3], sz[3], theta[NTHETA], theta2,
    w, wsum=0, xv[3];
  
  int i, id, idx, iphi, itheta, nphi=180, ntheta2=180;
  
  /* Allocate... */
  ALLOC(obs2, obs_t, 1);
  
  /* Set scattering phase function angles... */
  for(itheta=0; itheta<NTHETA; itheta++)
    theta[itheta]=M_PI*(double)itheta/(NTHETA-1);
  
  /* Get orthonormal basis (with respect to LOS)... */
  bascoord(dx, x, sx, sy, sz);  
  
  /* Initialize... */
  for(id=0; id<ctl->nd; id++)
    src_sca[id]=0;
  
  /* Loop over phase function angles... */
  for(itheta=0; itheta<ntheta2; itheta++) {
    
    /* Set phase function angle... */
    theta2=(0.5+itheta)/ntheta2*M_PI;
    
    /* Loop over azimuth angles... */
    for(iphi=0; iphi<nphi; iphi++) {
      
      /* Set azimuth angle... */
      phi=(0.5+iphi)/nphi*2*M_PI;
      
      /* Set view point... */
      for(i=0; i<3; i++)
	xv[i]=x[i]
	  +10*sin(theta2)*sin(phi)*sx[i]
	  +10*sin(theta2)*cos(phi)*sy[i]
	  +10*cos(theta2)*sz[i];
      
      /* Set observation geometry... */
      obs2->nr=1;
      cart2geo(x, &obs2->obsz[0], &obs2->obslon[0], &obs2->obslat[0]);
      cart2geo(xv, &obs2->vpz[0], &obs2->vplon[0], &obs2->vplat[0]);
      
      /* Get incident radiation... */
      formod_pencil(ctl, atm, obs2, aero, scattering-1, 0);
      
      /* Get phase function index */
      idx=locate(theta, NTHETA, theta2);

      /* Loop over channels... */
      for(id=0; id<ctl->nd; id++) {
	
	/* Interpolate phase function... */
	phase2=LIN(theta[idx], aero->p[il][id][idx],
		   theta[idx+1], aero->p[il][id][idx+1], theta2);
	
	/* Get weighting factor (surface element area * phase function)... */
	w=M_PI/ntheta2*2*M_PI*sin(theta2)/nphi*phase2;
	
	/* Integrate... */
	src_sca[id]+=w*obs2->rad[id][0];
	wsum+=w;
      }
    }
  }
  
  /* Normalize... */
  for(id=0; id<ctl->nd; id++)
    src_sca[id]/=wsum;
  
  /* Free... */
  free(obs2);
}

/*****************************************************************************/

void srcfunc_sca_sun(ctl_t *ctl,
		     atm_t *atm,
		     aero_t *aero,
		     double sec,
		     double *x,
		     double *dx,
		     int il,
		     double *src_sca) {

  los_t *los;
  
  obs_t *obs;
  
  double azi, dnorth[3], dout[3], ek[3], lx[3], ly[3], lz[3], phase2, sza,
    sza_beam, sza_cor, theta[NTHETA], theta2, x0[3], x1[3];
  
  int i, i2, id, idx, itheta;
  
  /* Allocate... */
  ALLOC(los, los_t, 1);
  ALLOC(obs, obs_t, 1);

  /* Set scattering phase function angles... */
  for(itheta=0; itheta<NTHETA; itheta++)
    theta[itheta]=M_PI*(double)itheta/(NTHETA-1);
  
  /* Get local coordinate system... */
  dnorth[0]=-x[0];
  dnorth[1]=-x[1];
  dnorth[2]=2*RE-x[2];
  bascoord(x, dnorth, lx, ly, lz);
  
  /* Get geometric coordinates of the Sun... */
  obs->nr=1;
  cart2geo(x, &obs->obsz[0], &obs->obslon[0], &obs->obslat[0]);
  suncoord(sec, obs->obslon[0], obs->obslat[0], &azi, &sza);
  
  /* Find true elevation angle of Sun... */
  sza_cor=sza;
  for(i2=0; i2<10; i2++) {
      
    /* Set observation geometry... */
    for(i=0; i<3; i++) {
      ek[i]
	=sin(sza_cor*M_PI/180)*sin(azi*M_PI/180)*lx[i]
	+sin(sza_cor*M_PI/180)*cos(azi*M_PI/180)*ly[i]
	+cos(sza_cor*M_PI/180)*lz[i];
      x1[i]=x[i]+10*ek[i];
    }
    cart2geo(x1, &obs->vpz[0], &obs->vplon[0], &obs->vplat[0]);
    
    /* Get zenith angle at end of beam... */
    raytrace(ctl, atm, obs, aero, los, 0);
    if(los->np<2)
      break;
    geo2cart(los->z[los->np-2], los->lon[los->np-2], los->lat[los->np-2], x0);
    geo2cart(los->z[los->np-1], los->lon[los->np-1], los->lat[los->np-1], x1);
    for(i=0; i<3; i++)
      dout[i]=x1[i]-x0[i];
    sza_beam=ANGLE(x, dout)*180/M_PI;
    
    /* Test for convergence... */
    if(fabs(sza_beam-sza)<0.01)
      break;
    
    /* Adapt geometric solar zenith angle (0.61803 golden ratio)... */
    sza_cor-=0.61803*(sza_beam-sza);
    sza_cor=GSL_MIN(GSL_MAX(sza_cor, 0), 180);
  }
  
  /* Check that LOS doesn't hit the ground... */
  if(los->tsurf<0) {
    
    /* Compute path transmittance... */
    formod_pencil(ctl, atm, obs, aero, 0, 0);
    
    /* Get phase function position... */
    theta2=ANGLE(ek, dx);
    idx=locate(theta, NTHETA, theta2);

    /* Loop over channels... */
    for(id=0; id<ctl->nd; id++) {
      
      /* Get phase function... */
      phase2=LIN(theta[idx], aero->p[il][id][idx],
		 theta[idx+1], aero->p[il][id][idx+1], theta2);
      
      /* Add solar radiance (6.764e-5 solid angle of the sun)... */
      src_sca[id]+=6.764e-5*phase2
	*planck(TSUN, ctl->nu[id])*obs->tau[id][0];
    }
  }

  /* Free... */
  free(los);
  free(obs);
}

/*****************************************************************************/

void suncoord(double sec,
	      double lon,
	      double lat,
	      double *azi,
	      double *sza) {
  
  double D, dec, e, g, GMST, h, L, LST, q, ra;
  
  /* Number of days and fraction with respect to 2000-01-01T12:00Z... */
  D=sec/86400-0.5;
  
  /* Geocentric apparent ecliptic longitude [rad]... */
  g=(357.529+0.98560028*D)*M_PI/180;
  q=280.459+0.98564736*D;
  L=(q+1.915*sin(g)+0.020*sin(2*g))*M_PI/180;
  
  /* Mean obliquity of the ecliptic [rad]... */
  e=(23.439-0.00000036*D)*M_PI/180;
  
  /* Declination [rad]... */
  dec=asin(sin(e)*sin(L));
  
  /* Right ascension [rad]... */
  ra=atan2(cos(e)*sin(L), cos(L));
  
  /* Greenwich Mean Sidereal Time [h]... */
  GMST=18.697374558+24.06570982441908*D;
  
  /* Local Sidereal Time [h]... */
  LST=GMST+lon/15;
  
  /* Hour angle [rad]... */
  h=LST/12*M_PI-ra;
  
  /* Convert latitude... */
  lat*=M_PI/180;
  
  /* Azimuth [deg]... */
  *azi=atan2(-sin(h), cos(lat)*tan(dec)-sin(lat)*cos(h))*180/M_PI;
  
  /* Solar zenith angle (90 deg - elevation) [deg]... */
  *sza=acos(sin(lat)*sin(dec)+cos(lat)*cos(dec)*cos(h))*180/M_PI;
}
