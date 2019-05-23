#include "retrievalmodel.h"

/*****************************************************************************/

void analyze_avk(ret_t *ret,
		 ctl_t *ctl,
		 atm_t *atm,
		 aero_t *aero,
		 int *iqa,
		 int *ipa,
		 gsl_matrix *avk) {
  
  static atm_t atm_cont, atm_res;
  static aero_t aero_cont, aero_res;
  
  int ig, iq, iw;
  
  size_t i, n, n0[NQMAX], n1[NQMAX];
  
  /* Get sizes... */
  n=avk->size1;
  
  /* Find sub-matrices for different quantities... */
  for(iq=0; iq<NQMAX; iq++) {
    n0[iq]=NMAX;
    for(i=0; i<n; i++) {
      if(iqa[i]==iq && n0[iq]==NMAX)
	n0[iq]=i;
      if(iqa[i]==iq)
	n1[iq]=i-n0[iq]+1;
    }
  }
  
  /* Initialize... */
  copy_atm(ctl, &atm_cont, atm, 1);
  copy_atm(ctl, &atm_res, atm, 1);
  copy_aero(ctl, &aero_cont, aero, 1);
  copy_aero(ctl, &aero_res, aero, 1);
  
  /* Analyze quantities... */
  analyze_avk_quantity(avk, IDXP, ipa, n0, n1, atm_cont.p, atm_res.p);
  analyze_avk_quantity(avk, IDXT, ipa, n0, n1, atm_cont.t, atm_res.t);
  analyze_avk_quantity(avk, IDXNN, ipa, n0, n1, aero_cont.nn, aero_res.nn);
  analyze_avk_quantity(avk, IDXRR, ipa, n0, n1, aero_cont.rr, aero_res.rr);
  analyze_avk_quantity(avk, IDXSS, ipa, n0, n1, aero_cont.ss, aero_res.ss);
  for(ig=0; ig<ctl->ng; ig++)
    analyze_avk_quantity(avk, IDXQ(ig), ipa, n0, n1,
			 atm_cont.q[ig], atm_res.q[ig]);
  for(iw=0; iw<ctl->nw; iw++)
    analyze_avk_quantity(avk, IDXK(iw), ipa, n0, n1,
			 atm_cont.k[iw], atm_res.k[iw]);
  
  /* analyze_avk_quantity(avk, IDXNN, ipa, n0, n1, */
  /* 			 aero_cont.nn, aero_res.nn); */
  
  
  /* Write results to disk... */
  write_atm(ret->dir, "atm_cont.tab", ctl, &atm_cont); /* contribution from averaging kernel */
  write_atm(ret->dir, "atm_res.tab", ctl, &atm_res);
  if(ctl->retnn || ctl->retrr || ctl->retss)
    write_aero(ret->dir, "aero_cont.tab", &aero_cont);
}

/*****************************************************************************/

void analyze_avk_quantity(gsl_matrix *avk,
			  int iq,
			  int *ipa,
			  size_t *n0,
			  size_t *n1,
			  double *cont,
			  double *res) {
  
  size_t i, j;
  
  /* Loop over state vector elements... */
  if(n0[iq]<NMAX)
    for(i=0; i<n1[iq]; i++) {
      
      /* Get area of averagig kernel... */
      for(j=0; j<n1[iq]; j++)
	cont[ipa[n0[iq]+i]]+=gsl_matrix_get(avk, n0[iq]+i, n0[iq]+j);
      
      /* Get information density... */
      res[ipa[n0[iq]+i]]=1/gsl_matrix_get(avk, n0[iq]+i, n0[iq]+i);
    }
}

/*****************************************************************************/

double corr_function(double z0,
		     double lon0,
		     double lat0,
		     double z1,
		     double lon1,
		     double lat1,
		     double cz,
		     double ch) {
  
  double x0[3], x1[3];
  
  /* Get Cartesian coordinates... */
  geo2cart(0, lon0, lat0, x0);
  geo2cart(0, lon1, lat1, x1);
  
  /* Compute correlations... */
  return exp(-DIST(x0, x1)/ch-fabs(z0-z1)/cz);
}

/*****************************************************************************/

double cost_function(FILE *out,
                     int it,
                     gsl_vector *dx,
                     gsl_vector *dy,
                     gsl_matrix *s_a_inv,
                     gsl_vector *sig_eps_inv) {
  
  gsl_vector *x_aux, *y_aux;
  
  double chisq, chisq_a, chisq_m=0;
  
  size_t i, m, n;
  
  /* Get sizes... */
  m=dy->size;
  n=dx->size;
  
  /* Allocate... */
  x_aux=gsl_vector_alloc(n);
  y_aux=gsl_vector_alloc(m);
  
  /* Determine normalized cost function...
     (chi^2 = 1/m * [dy^T * S_eps^{-1} * dy + dx^T * S_a^{-1} * dx]) */
  for(i=0; i<m; i++)
    chisq_m+=gsl_pow_2(gsl_vector_get(dy, i)*gsl_vector_get(sig_eps_inv, i));
  chisq_m/=(double)m;
  gsl_blas_dgemv(CblasNoTrans, 1.0, s_a_inv, dx, 0.0, x_aux);
  gsl_blas_ddot(dx, x_aux, &chisq_a);
  chisq_a/=(double)m;
  chisq=chisq_m+chisq_a;
  
  /* Write info... */
  printf("it= %d / chi^2/m= %g (meas: %g / apr: %g)\n",
         it, chisq, chisq_m, chisq_a);
  
  /* Write header... */
  if(it==0)
    fprintf(out,
            "# $1 = iteration number\n"
            "# $2 = normalized cost function: total\n"
            "# $3 = normalized cost function: measurements\n"
            "# $4 = normalized cost function: a priori\n"
            "# $5 = number of measurements\n"
            "# $6 = number of state vector elements\n\n");
  
  /* Write data... */
  fprintf(out, "%d %g %g %g %d %d\n",
          it, chisq, chisq_m, chisq_a, (int)m, (int)n);
  
  /* Free... */
  gsl_vector_free(x_aux);
  gsl_vector_free(y_aux);
  
  /* Return cost function value... */
  return chisq;
}

/*****************************************************************************/

void idx2name(ctl_t *ctl,
	      int idx,
	      char *quantity) {
  
  int ig, iw;
  
  if(idx==IDXP)
    sprintf(quantity, "PRESSURE");
  
  if(idx==IDXT)
    sprintf(quantity, "TEMPERATURE");
  
  for(ig=0; ig<ctl->ng; ig++)
    if(idx==IDXQ(ig))
      sprintf(quantity, "%s", ctl->emitter[ig]);
  
  for(iw=0; iw<ctl->nw; iw++)
    if(idx==IDXK(iw))
      sprintf(quantity, "EXTINCT_WINDOW%d", iw);

  if(idx==IDXNN)
    sprintf(quantity, "CONCENTRATION");
  if(idx==IDXRR)
    sprintf(quantity, "RADIUS");
  if(idx==IDXSS)
    sprintf(quantity, "PSD_WIDTH");
}

/*****************************************************************************/

void kernel(ctl_t *ctl,
	    atm_t *atm,
	    obs_t *obs,
	    aero_t *aero,
	    gsl_matrix *k) {
  
  static atm_t atm1;
  
  static obs_t obs1;
  
  static aero_t aero1;
  
  static int ipa[NMAX], iqa[NMAX], ira[MMAX];
  
  gsl_vector *x0, *x1, *yy0, *yy1;
  
  double h;
  
  size_t i, j, n, m;
  
  /* Get sizes... */
  m=k->size1;
  n=k->size2;
  
  /* Allocate... */
  x0=gsl_vector_alloc(n);
  x1=gsl_vector_alloc(n);
  yy0=gsl_vector_alloc(m);
  yy1=gsl_vector_alloc(m);
  /* Compute radiance for undisturbed atmospheric data... */
  formod(ctl, atm, obs, aero);
  /* Compose vectors... */
  atm2x(ctl, atm, aero, x0, iqa, ipa);
  obs2y(ctl, obs, yy0, NULL, ira);

  /* Initialize kernel matrix... */
  gsl_matrix_set_zero(k);
  /* Loop over state vector elements... */
  for(j=0; j<n; j++) {
    
    /* Set perturbation size... */
    h=-999;
    
    /* Pressure... */
    if(iqa[j]==IDXP)
      h=GSL_MAX(fabs(0.01*gsl_vector_get(x0, j)), 1e-7);
    
    /* Temperature... */
    if(iqa[j]==IDXT)
      h=1;
    
    /* Volume mixing ratios... */
    if(iqa[j]>=IDXQ(0) && iqa[j]<IDXQ(ctl->ng))
      h=GSL_MAX(fabs(0.01*gsl_vector_get(x0, j)), 1e-15);
    
    /* Extinction... */
    if(iqa[j]>=IDXK(0) && iqa[j]<IDXK(ctl->nw))
      h=1e-4;
    /* Particle concentration... */
    if(iqa[j]==IDXNN)
      h=GSL_MAX(fabs(0.5*gsl_vector_get(x0, j)), 0.1);
    /* Particle size... */
    if(iqa[j]==IDXRR)
      h=GSL_MAX(fabs(0.25*gsl_vector_get(x0, j)), 0.1);
    /* Particle size distribution width... */
    if(iqa[j]==IDXSS)
      h=GSL_MAX(fabs(0.1*gsl_vector_get(x0, j)), 0.01);

    /* Check perturbation size... */
    if(h<=0)
      ERRMSG("Cannot set perturbation size!");
      
    /* Disturb state vector element... */
    gsl_vector_memcpy(x1, x0);
    gsl_vector_set(x1, j, gsl_vector_get(x1, j)+h);
    copy_atm(ctl, &atm1, atm, 0);
    copy_aero(ctl, &aero1, aero, 0);
    copy_obs(ctl, &obs1, obs, 0);
    x2atm(ctl, x1, &atm1, &aero1);

    /* Compute radiance for disturbed atmospheric data... */
    formod(ctl, &atm1, &obs1, &aero1);

    /* Compose measurement vector for disturbed radiance data... */
    obs2y(ctl, &obs1, yy1, NULL, NULL);
    
    /* Compute derivatives... */
    for(i=0; i<m; i++)
      gsl_matrix_set(k, i, j,
		     (gsl_vector_get(yy1, i)-gsl_vector_get(yy0, i))/h);
  }
  
  /* Free... */
  gsl_vector_free(x0);
  gsl_vector_free(x1);
  gsl_vector_free(yy0);
  gsl_vector_free(yy1);
}

/*****************************************************************************/

void matrix_invert(gsl_matrix *a) {
  
  size_t diag=1, i, j, n;
  
  /* Get size... */
  n=a->size1;
  
  /* Check if matrix is diagonal... */
  for(i=0; i<n && diag; i++)
    for(j=i+1; j<n; j++)
      if(gsl_matrix_get(a, i, j)!=0) {
        diag=0;
        break;
      }
  
  /* Quick inversion of diagonal matrix... */
  if(diag)
    for(i=0; i<n; i++)
      gsl_matrix_set(a, i, i, 1/gsl_matrix_get(a, i, i));
  
  /* Matrix inversion by means of Cholesky decomposition... */
  else {
    gsl_linalg_cholesky_decomp(a);
    gsl_linalg_cholesky_invert(a);
  }
}

/*****************************************************************************/

void matrix_product(gsl_matrix *a,
		    gsl_vector *b,
		    int transpose,
		    gsl_matrix *c) {
  
  gsl_matrix *aux;
  
  size_t i, j, m, n;
  
  /* Set sizes... */
  m=a->size1;
  n=a->size2;
  
  /* Allocate... */
  aux=gsl_matrix_alloc(m, n);
  
  /* Compute A^T B A... */
  if(transpose==1) {
    
    /* Compute B^1/2 A... */
    for(i=0; i<m; i++)
      for(j=0; j<n; j++)
	gsl_matrix_set(aux, i, j,
		       gsl_vector_get(b, i)*gsl_matrix_get(a, i, j));
    
    /* Compute A^T B A = (B^1/2 A)^T (B^1/2 A)... */
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, aux, aux, 0.0, c);
  }
  
  /* Compute A B A^T... */
  else if(transpose==2) {

    /* Compute A B^1/2... */
    for(i=0; i<m; i++)
      for(j=0; j<n; j++)
	gsl_matrix_set(aux, i, j,
		       gsl_matrix_get(a, i, j)*gsl_vector_get(b, j));
    
    /* Compute A B A^T = (A B^1/2) (A B^1/2)^T... */
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, aux, aux, 0.0, c);
  }

  /* Free... */
  gsl_matrix_free(aux);
}

/*****************************************************************************/

size_t obs2y(ctl_t *ctl,
	     obs_t *obs,
	     gsl_vector *y,
	     int *ida,
	     int *ira) {
  
  int id, ir;
  
  size_t m=0;
  
  /* Determine measurement vector... */
  for(ir=0; ir<obs->nr; ir++)
    for(id=0; id<ctl->nd; id++)
      if(gsl_finite(obs->rad[id][ir])) {
	if(y!=NULL)
	  gsl_vector_set(y, m, obs->rad[id][ir]);
	if(ida!=NULL)
	  ida[m]=id;
	if(ira!=NULL)
	  ira[m]=ir;
	m++;
      }
  
  return m;
}

/*****************************************************************************/

void optimal_estimation(ret_t *ret,
			ctl_t *ctl,
			obs_t *obs_meas,
			obs_t *obs_i,
			atm_t *atm_apr,
			atm_t *atm_i,
			aero_t *aero_apr,
			aero_t *aero_i) {
  
  static int ipa[NMAX], iqa[NMAX];
  
  gsl_matrix *a, *auxnm, *cov, *gain, *k_i, *s_a_inv;
  
  gsl_vector *b, *dx, *dy, *sig_eps_inv, *sig_formod, *sig_noise,
    *x_a, *x_i, *x_step, *y_aux, *y_i, *y_m;
  
  FILE *out;
  
  char filename[2*LEN];
  
  double chisq, chisq_old, disq=0, lmpar=0.001;
  
  int ig, ip, it, it2, iw, im;
  
  size_t i, j, m, n;
  
  /* ------------------------------------------------------------
     Initialize...
     ------------------------------------------------------------ */
  
  /* Get sizes... */
  m=obs2y(ctl, obs_meas, NULL, NULL, NULL);
  n=atm2x(ctl, atm_apr, aero_apr, NULL, iqa, ipa);
  if(m<=0 || n<=0)
    ERRMSG("Check problem definition!");
  
  /* Write info... */
  printf("Problem size: m= %d / n= %d (alloc= %.4g MB / stat= %.4g MB)\n",
         (int)m, (int)n,
         (double)(3*m*n+3*n*n+8*m+8*n)*sizeof(double)/1024./1024.,
         (double)(5*sizeof(atm_t)+3*sizeof(obs_t)
                  +2*NMAX*sizeof(int)+3*sizeof(aero_t))/1024./1024.);
  
  /* Allocate... */
  a=gsl_matrix_alloc(n, n);
  cov=gsl_matrix_alloc(n, n);
  k_i=gsl_matrix_alloc(m, n);
  s_a_inv=gsl_matrix_alloc(n, n);
  
  b=gsl_vector_alloc(n);
  dx=gsl_vector_alloc(n);
  dy=gsl_vector_alloc(m);
  sig_eps_inv=gsl_vector_alloc(m);
  sig_formod=gsl_vector_alloc(m);
  sig_noise=gsl_vector_alloc(m);
  x_a=gsl_vector_alloc(n);
  x_i=gsl_vector_alloc(n);
  x_step=gsl_vector_alloc(n);
  y_aux=gsl_vector_alloc(m);
  y_i=gsl_vector_alloc(m);
  y_m=gsl_vector_alloc(m);
  
  /* Set initial state... */
  copy_atm(ctl, atm_i, atm_apr, 0);
  copy_obs(ctl, obs_i, obs_meas, 0);
  copy_aero(ctl, aero_i, aero_apr, 0);
  formod(ctl, atm_i, obs_i, aero_i);

  /* Set state vectors and observation vectors... */
  atm2x(ctl, atm_apr, aero_apr, x_a, NULL, NULL);
  atm2x(ctl, atm_i, aero_i, x_i, NULL, NULL);
  obs2y(ctl, obs_meas, y_m, NULL, NULL);
  obs2y(ctl, obs_i, y_i, NULL, NULL);
  
  /* Set inverse a priori covariance S_a^-1... */
  set_cov_apr(ret, ctl, atm_apr, aero_apr, iqa, ipa, s_a_inv);
  write_matrix(ret->dir, "matrix_cov_apr.tab", ctl, s_a_inv,
	       atm_i, aero_i, obs_i, "x", "x", "r");
  matrix_invert(s_a_inv);

  /* Get measurement errors... */
  set_cov_meas(ret, ctl, obs_meas, sig_noise, sig_formod, sig_eps_inv);

  /* Create cost function file... */
  sprintf(filename, "%s/costs.tab", ret->dir);
  if(!(out=fopen(filename, "w")))
    ERRMSG("Cannot create cost function file!");
  
  /* Determine dx = x_i - x_a and dy = y - F(x_i) ... */
  gsl_vector_memcpy(dx, x_i);
  gsl_vector_sub(dx, x_a);
  gsl_vector_memcpy(dy, y_m);
  gsl_vector_sub(dy, y_i);

  /* Compute and check cost function... */
  if(!gsl_finite(chisq=cost_function(out, 0, dx, dy, s_a_inv, sig_eps_inv))) {
    printf("Retrieval failed!\n");
    return;
  }
  
  /* Compute initial kernel... */
  kernel(ctl, atm_i, obs_i, aero_i, k_i);

  /* ------------------------------------------------------------
     Levenberg-Marquardt minimization...
     ------------------------------------------------------------ */
  
  /* Outer loop... */
  for(it=1; it<=ret->conv_itmax; it++) {
    
    /* Store current state... */
    chisq_old=chisq;
    
    /* Compute kernel matrix K_i... */
    if(it%ret->kernel_recomp==0 && !(ret->kernel_recomp==1 && it==1))
      kernel(ctl, atm_i, obs_i, aero_i, k_i);
    
    /* Compute K_i^T * S_eps^{-1} * K_i ... */
    if(it%ret->kernel_recomp==0 || it==1)
      matrix_product(k_i, sig_eps_inv, 1, cov);
    
    /* Determine dx = x_i - x_a ... */
    gsl_vector_memcpy(dx, x_i);
    gsl_vector_sub(dx, x_a);

    /* Determine dy = y - F(x_i) ... */
    gsl_vector_memcpy(dy, y_m);
    gsl_vector_sub(dy, y_i);

    /* Determine b = K_i^T * S_eps^{-1} * dy - S_a^{-1} * dx ...*/
    for(i=0; i<m; i++)
      gsl_vector_set(y_aux, i, gsl_vector_get(dy, i)
		     *gsl_pow_2(gsl_vector_get(sig_eps_inv, i)));
    gsl_blas_dgemv(CblasTrans, 1.0, k_i, y_aux, 0.0, b);
    gsl_blas_dgemv(CblasNoTrans, -1.0, s_a_inv, dx, 1.0, b);
    
    /* Inner loop... */
    for(it2=0; it2<20; it2++) {
      
      /* Compute A = (1 + lmpar) * S_a^{-1} + K_i^T * S_eps^{-1} * K_i ... */
      gsl_matrix_memcpy(a, s_a_inv);
      gsl_matrix_scale(a, 1+lmpar);
      gsl_matrix_add(a, cov);
      
      /* Solve A * x_step = b by means of Cholesky decomposition... */
      gsl_linalg_cholesky_decomp(a);
      gsl_linalg_cholesky_solve(a, b, x_step);
      
      /* Update atmospheric state... */
      gsl_vector_add(x_i, x_step);
      copy_atm(ctl, atm_i, atm_apr, 0);
      copy_obs(ctl, obs_i, obs_meas, 0);
      copy_aero(ctl, aero_i, aero_apr, 0);
      x2atm(ctl, x_i, atm_i, aero_i);

      /* Check atmospheric state... */
      for(ip=0; ip<atm_i->np; ip++) {
	atm_i->p[ip]=GSL_MIN(GSL_MAX(atm_i->p[ip], 5e-7), 5e4);
	atm_i->t[ip]=GSL_MIN(GSL_MAX(atm_i->t[ip], 100), 400);
	for(ig=0; ig<ctl->ng; ig++)
	  atm_i->q[ig][ip]=GSL_MIN(GSL_MAX(atm_i->q[ig][ip], 0), 1);
	for(iw=0; iw<ctl->nw; iw++)
	  atm_i->k[iw][ip]=GSL_MAX(atm_i->k[iw][ip], 0);
      }

      /* check aerosol parameters... */
      for(im=0; im<aero_i->nm; im++) {
	aero_i->nn[im]=GSL_MIN(GSL_MAX(aero_i->nn[im], 1e-6), 1e4);
      }
      for(im=0; im<aero_i->nm; im++) {
	aero_i->rr[im]=GSL_MIN(GSL_MAX(aero_i->rr[im], 1e-2), 100);
      }
      for(im=0; im<aero_i->nm; im++) {
	aero_i->ss[im]=GSL_MIN(GSL_MAX(aero_i->ss[im], 1.1), 3);
      }
      

      /* printf("  n= %g\n", aero_i->nn[0]); */
      /* printf("  r= %g\n", aero_i->rr[0]); */


      /* Forward calculation... */
      formod(ctl, atm_i, obs_i, aero_i);
      obs2y(ctl, obs_i, y_i, NULL, NULL);

      /* Determine dx = x_i - x_a and dy = y - F(x_i) ... */
      gsl_vector_memcpy(dx, x_i);
      gsl_vector_sub(dx, x_a);
      gsl_vector_memcpy(dy, y_m);
      gsl_vector_sub(dy, y_i);
      
      /* Compute cost function... */
      chisq=cost_function(out, it, dx, dy, s_a_inv, sig_eps_inv);
      
      /* Modify Levenberg-Marquardt parameter... */
      if(chisq>chisq_old) {
	lmpar*=10;
	gsl_vector_sub(x_i, x_step);
      } else {
	lmpar/=10;
	break;
      }
    }
    
    /* Get normalized step size in state space... */
    gsl_blas_ddot(x_step, b, &disq);
    disq/=(double)n;
    
    /* Convergence test... */
    if(disq<ret->conv_dmin)
      break;
  }
  
  /* Close cost function file... */
  fclose(out);
  
  /* Store results... */
  write_obs(ret->dir, "obs_final.tab", ctl, obs_i);
  write_atm(ret->dir, "atm_final.tab", ctl, atm_i);
  write_aero(ret->dir, "aero_final.tab", aero_i);
  write_matrix(ret->dir, "matrix_kernel.tab", ctl, k_i,
	       atm_i, aero_i, obs_i, "y", "x", "r");
  
  /* ------------------------------------------------------------
     Analysis of retrieval results...
     ------------------------------------------------------------ */
  
  /* Check if error analysis is requested... */
  if(ret->err_ana) {
    
    /* Allocate... */
    auxnm=gsl_matrix_alloc(n, m);
    gain=gsl_matrix_alloc(n, m);
    
    /* Compute inverse retrieval covariance...
       cov^{-1} = S_a^{-1} + K_i^T * S_eps^{-1} * K_i */
    matrix_product(k_i, sig_eps_inv, 1, cov);
    gsl_matrix_add(cov, s_a_inv);
    
    /* Compute retrieval covariance... */
    matrix_invert(cov);
    write_matrix(ret->dir, "matrix_cov_ret.tab", ctl, cov,
		 atm_i, aero_i, obs_i, "x", "x", "r");

    write_stddev("total", ret, ctl, atm_i, aero_i, cov);
    
    /* Compute gain matrix...
       G = cov * K^T * S_eps^{-1} */
    for(i=0; i<n; i++)
      for(j=0; j<m; j++)
	gsl_matrix_set(auxnm, i, j, gsl_matrix_get(k_i, j, i)
		       *gsl_pow_2(gsl_vector_get(sig_eps_inv, j)));
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, cov, auxnm, 0.0, gain);
    write_matrix(ret->dir, "matrix_gain.tab", ctl, gain,
		 atm_i, aero_i, obs_i, "x", "y", "c");
    
    /* Compute retrieval error due to noise... */
    matrix_product(gain, sig_noise, 2, a);
    write_stddev("noise", ret, ctl, atm_i, aero_i, a);
    
    /* Compute retrieval error  due to forward model errors... */
    matrix_product(gain, sig_formod, 2, a);
    write_stddev("formod", ret, ctl, atm_i, aero_i, a);
    
    /* Compute averaging kernel matrix
       A = G * K ... */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, gain, k_i, 0.0, a);
    write_matrix(ret->dir, "matrix_avk.tab", ctl, a,
		 atm_i, aero_i, obs_i, "x", "x", "r");
    
    /* Analyze averaging kernel matrix... */
    analyze_avk(ret, ctl, atm_i, aero_i, iqa, ipa, a);

    /* Free... */
    gsl_matrix_free(auxnm);
    gsl_matrix_free(gain);
  }
  
  /* ------------------------------------------------------------
     Finalize...
     ------------------------------------------------------------ */
  
  gsl_matrix_free(a);
  gsl_matrix_free(cov);
  gsl_matrix_free(k_i);
  gsl_matrix_free(s_a_inv);
  
  gsl_vector_free(b);
  gsl_vector_free(dx);
  gsl_vector_free(dy);
  gsl_vector_free(sig_eps_inv);
  gsl_vector_free(sig_formod);
  gsl_vector_free(sig_noise);
  gsl_vector_free(x_a);
  gsl_vector_free(x_i);
  gsl_vector_free(x_step);
  gsl_vector_free(y_aux);
  gsl_vector_free(y_i);
  gsl_vector_free(y_m);
}

/*****************************************************************************/

void read_ret(int argc,
	      char *argv[],
	      ctl_t *ctl,
	      ret_t *ret) {
  
  int id, ig, iw;
  
  /* Iteration control... */
  ret->kernel_recomp=(int)scan_ctl(argc, argv, "KERNEL_RECOMP", -1, "1", NULL);
  ret->conv_itmax=(int)scan_ctl(argc, argv, "CONV_ITMAX", -1, "20", NULL);
  ret->conv_dmin=scan_ctl(argc, argv, "CONV_DMIN", -1, "0.1", NULL);
  
  /* Filtering of bad observations... */
  ret->resmax=scan_ctl(argc, argv, "RESMAX", -1, "-999", NULL);
  
  /* Error analysis... */
  ret->err_ana=(int)scan_ctl(argc, argv, "ERR_ANA", -1, "1", NULL);
  
  for(id=0; id<ctl->nd; id++)
    ret->err_formod[id]=scan_ctl(argc, argv, "ERR_FORMOD", id, "0", NULL);

  for(id=0; id<ctl->nd; id++)
    ret->err_noise[id]=scan_ctl(argc, argv, "ERR_NOISE", id, "0", NULL);
  
  ret->err_press=scan_ctl(argc, argv, "ERR_PRESS", -1, "0", NULL);
  ret->err_press_cz=scan_ctl(argc, argv, "ERR_PRESS_CZ", -1, "-999", NULL);
  ret->err_press_ch=scan_ctl(argc, argv, "ERR_PRESS_CH", -1, "-999", NULL);

  ret->err_temp=scan_ctl(argc, argv, "ERR_TEMP", -1, "0", NULL);
  ret->err_temp_cz=scan_ctl(argc, argv, "ERR_TEMP_CZ", -1, "-999", NULL);
  ret->err_temp_ch=scan_ctl(argc, argv, "ERR_TEMP_CH", -1, "-999", NULL);

  ret->err_nn=scan_ctl(argc, argv, "ERR_NN", -1, "0", NULL);
  ret->err_rr=scan_ctl(argc, argv, "ERR_RR", -1, "0", NULL);
  ret->err_ss=scan_ctl(argc, argv, "ERR_SS", -1, "0", NULL);

  for(ig=0; ig<ctl->ng; ig++) {
    ret->err_q[ig]=scan_ctl(argc, argv, "ERR_Q", ig, "0", NULL);
    ret->err_q_cz[ig]=scan_ctl(argc, argv, "ERR_Q_CZ", ig, "-999", NULL);
    ret->err_q_ch[ig]=scan_ctl(argc, argv, "ERR_Q_CH", ig, "-999", NULL);
  }

  for(iw=0; iw<ctl->nw; iw++) {
    ret->err_k[iw]=scan_ctl(argc, argv, "ERR_K", iw, "0", NULL);
    ret->err_k_cz[iw]=scan_ctl(argc, argv, "ERR_K_CZ", iw, "-999", NULL);
    ret->err_k_ch[iw]=scan_ctl(argc, argv, "ERR_K_CH", iw, "-999", NULL);
  }
}

/*****************************************************************************/

void set_cov_apr(ret_t *ret,
		 ctl_t *ctl,
		 atm_t *atm,
		 aero_t *aero,
		 int *iqa,
		 int *ipa,
		 gsl_matrix *s_a) {
  
  gsl_vector *x_a;
  
  double ch, cz;
  
  int ig, iw;
  
  size_t i, j, n;
  
  /* Get sizes... */
  n=s_a->size1;
  
  /* Allocate... */
  x_a=gsl_vector_alloc(n);
  
  /* Get sigma vector... */
  atm2x(ctl, atm, aero, x_a, NULL, NULL);
  for(i=0; i<n; i++) {
    if(iqa[i]==IDXP)
      gsl_vector_set(x_a, i, ret->err_press/100*gsl_vector_get(x_a, i));
    if(iqa[i]==IDXT)
      gsl_vector_set(x_a, i, ret->err_temp);
    for(ig=0; ig<ctl->ng; ig++)
      if(iqa[i]==IDXQ(ig))
	gsl_vector_set(x_a, i, ret->err_q[ig]/100*gsl_vector_get(x_a, i));
    for(iw=0; iw<ctl->nw; iw++)
      if(iqa[i]==IDXK(iw))
	gsl_vector_set(x_a, i, ret->err_k[iw]);
    if(iqa[i]==IDXNN)
      gsl_vector_set(x_a, i, ret->err_nn);
    if(iqa[i]==IDXRR)
      gsl_vector_set(x_a, i, ret->err_rr);
    if(iqa[i]==IDXSS)
      gsl_vector_set(x_a, i, ret->err_ss);
  }
  
  /* Initialize diagonal covariance... */
  gsl_matrix_set_zero(s_a);
  for(i=0; i<n; i++)
    gsl_matrix_set(s_a, i, i, gsl_pow_2(gsl_vector_get(x_a, i)));
  
  /* Loop over matrix elements... */
  for(i=0; i<n; i++)
    for(j=0; j<n; j++)
      if(i!=j && iqa[i]==iqa[j]) {
	
	/* Initialize... */
	cz=ch=0;
	
	/* Set correlation lengths for pressure... */
	if(iqa[i]==IDXP) {
	  cz=ret->err_press_cz;
	  ch=ret->err_press_ch;
	}
	
	/* Set correlation lengths for temperature... */
	if(iqa[i]==IDXT) {
	  cz=ret->err_temp_cz;
	  ch=ret->err_temp_ch;
	}
	
	/* Set correlation lengths for volume mixing ratios... */
	for(ig=0; ig<ctl->ng; ig++)
	  if(iqa[i]==IDXQ(ig)) {
	    cz=ret->err_q_cz[ig];
	    ch=ret->err_q_ch[ig];
	  }
	
	/* Set correlation lengths for extinction... */
	for(iw=0; iw<ctl->nw; iw++)
	  if(iqa[i]==IDXK(iw)) {
	    cz=ret->err_k_cz[iw];
	    ch=ret->err_k_ch[iw];
	  }
	
	/* Compute correlations... */
	if(cz>0 && ch>0)
	  gsl_matrix_set(s_a, i, j,
			 gsl_vector_get(x_a, i)*gsl_vector_get(x_a, j)
			 *corr_function(atm->z[ipa[i]],
					atm->lon[ipa[i]],
					atm->lat[ipa[i]],
					atm->z[ipa[j]],
					atm->lon[ipa[j]],
					atm->lat[ipa[j]], cz, ch));
      }
  
  /* Free... */
  gsl_vector_free(x_a);
}

/*****************************************************************************/

void set_cov_meas(ret_t *ret,
		  ctl_t *ctl,
		  obs_t *obs,
		  gsl_vector *sig_noise,
		  gsl_vector *sig_formod,
		  gsl_vector *sig_eps_inv) {
  
  static obs_t obs_err;
  
  int id, ir;
  
  size_t i, m;
  
  /* Get size... */
  m=sig_eps_inv->size;
  
  /* Noise error (always considered in retrieval fit)... */
  copy_obs(ctl, &obs_err, obs, 1);
  for(ir=0; ir<obs_err.nr; ir++)
    for(id=0; id<ctl->nd; id++)
      obs_err.rad[id][ir]
	=(gsl_finite(obs->rad[id][ir]) ? ret->err_noise[id] : GSL_NAN);
  obs2y(ctl, &obs_err, sig_noise, NULL, NULL);
  
  /* Forward model error (always considered in retrieval fit)... */
  copy_obs(ctl, &obs_err, obs, 1);
  for(ir=0; ir<obs_err.nr; ir++)
    for(id=0; id<ctl->nd; id++)
      obs_err.rad[id][ir]
	=fabs(ret->err_formod[id]/100*obs->rad[id][ir]);
  obs2y(ctl, &obs_err, sig_formod, NULL, NULL);
  
  /* Total error... */
  for(i=0; i<m; i++)
    gsl_vector_set(sig_eps_inv, i,
		   1/sqrt(gsl_pow_2(gsl_vector_get(sig_noise, i))
			  +gsl_pow_2(gsl_vector_get(sig_formod, i))));
}

/*****************************************************************************/

void write_stddev(const char *quantity,
		  ret_t *ret,
		  ctl_t *ctl,
		  atm_t *atm,
		  aero_t *aero,
		  gsl_matrix *s) {
  
  static atm_t atm_aux;

  static aero_t aero_aux;
  
  gsl_vector *x_aux;
  
  char filename[LEN];
  
  size_t i, n;
  
  /* Get sizes... */
  n=s->size1;
  
  /* Allocate... */
  x_aux=gsl_vector_alloc(n);
  
  /* Compute standard deviation... */
  for(i=0; i<n; i++)
    gsl_vector_set(x_aux, i, sqrt(gsl_matrix_get(s, i, i)));
  
  /* Write to disk... */
  copy_atm(ctl, &atm_aux, atm, 1);
  copy_aero(ctl, &aero_aux, aero, 1);
  x2atm(ctl, x_aux, &atm_aux, &aero_aux);
  sprintf(filename, "atm_err_%s.tab", quantity);
  write_atm(ret->dir, filename, ctl, &atm_aux);
  if (ctl->retnn || ctl->retrr || ctl->retss){
    sprintf(filename, "aero_err_%s.tab", quantity);
    write_aero(ret->dir, filename, &aero_aux);
  }

  /* Free... */
  gsl_vector_free(x_aux);
}

/*****************************************************************************/

void write_matrix(const char *dirname,
		  const char *filename,
		  ctl_t *ctl,
		  gsl_matrix *matrix,
		  atm_t *atm,
		  aero_t *aero,
		  obs_t *obs,
		  const char *rowspace,
		  const char *colspace,
		  const char *sort) {
  
  static int cida[MMAX], ciqa[NMAX], cipa[NMAX], cira[MMAX],
    rida[MMAX], riqa[NMAX], ripa[NMAX], rira[MMAX];
  
  FILE *out;
  
  char file[LEN], quantity[LEN];
  
  size_t i, j, nc, nr;
  
  /* Check output flag... */
  if(!ctl->write_matrix)
    return;
  
  /* Set filename... */
  if(dirname!=NULL)
    sprintf(file, "%s/%s", dirname, filename);
  else
    sprintf(file, "%s", filename);
  
  /* Write info... */
  printf("Write matrix: %s\n", file);
  
  /* Create file... */
  if(!(out=fopen(file, "w")))
    ERRMSG("Cannot create file!");
  
  /* Write header (row space)... */
  if(rowspace[0]=='y') {
    
    fprintf(out,
	    "# $1 = Row: index (measurement space)\n"
	    "# $2 = Row: channel wavenumber [cm^-1]\n"
	    "# $3 = Row: time (seconds since 2000-01-01T00:00Z)\n"
	    "# $4 = Row: view point altitude [km]\n"
	    "# $5 = Row: view point longitude [deg]\n"
	    "# $6 = Row: view point latitude [deg]\n");
    
    /* Get number of rows... */
    nr=obs2y(ctl, obs, NULL, rida, rira);
      
  } else {
    
    fprintf(out,
	    "# $1 = Row: index (state space)\n"
	    "# $2 = Row: name of quantity\n"
	    "# $3 = Row: time (seconds since 2000-01-01T00:00Z)\n"
	    "# $4 = Row: altitude [km]\n"
	    "# $5 = Row: longitude [deg]\n"
	    "# $6 = Row: latitude [deg]\n");
    
    /* Get number of rows... */
    nr=atm2x(ctl, atm, aero, NULL, riqa, ripa);
  }
  
  /* Write header (column space)... */  
  if(colspace[0]=='y') {
    
    fprintf(out,
	    "# $7 = Col: index (measurement space)\n"
	    "# $8 = Col: channel wavenumber [cm^-1]\n"
	    "# $9 = Col: time (seconds since 2000-01-01T00:00Z)\n"
	    "# $10 = Col: view point altitude [km]\n"
	    "# $11 = Col: view point longitude [deg]\n"
	    "# $12 = Col: view point latitude [deg]\n");
    
    /* Get number of columns... */
    nc=obs2y(ctl, obs, NULL, cida, cira);
    
    } else {
    
    fprintf(out,
	    "# $7 = Col: index (state space)\n"
	    "# $8 = Col: name of quantity\n"
	    "# $9 = Col: time (seconds since 2000-01-01T00:00Z)\n"
	    "# $10 = Col: altitude [km]\n"
	    "# $11 = Col: longitude [deg]\n"
	    "# $12 = Col: latitude [deg]\n");
    
    /* Get number of columns... */
    nc=atm2x(ctl, atm, aero, NULL, ciqa, cipa);
  }
  
  /* Write header entry... */
  fprintf(out, "# $13 = Matrix element\n\n");
  
  /* Write matrix data... */
  i=j=0;
  while(i<nr && j<nc) {
    
    /* Check matrix value... */
    if(gsl_matrix_get(matrix, i, j)!=0) {
      
      /* Write info about the row... */
      if(rowspace[0]=='y')
	fprintf(out, "%d %g %.2f %g %g %g",
		(int)i, ctl->nu[rida[i]],
		obs->time[rira[i]], obs->vpz[rira[i]],
		obs->vplon[rira[i]], obs->vplat[rira[i]]);
      else {
	idx2name(ctl, riqa[i], quantity);
	fprintf(out, "%d %s %.2f %g %g %g", (int)i, quantity,
		atm->time[ripa[i]], atm->z[ripa[i]],
		atm->lon[ripa[i]], atm->lat[ripa[i]]);
	}
      
      /* Write info about the column... */
      if(colspace[0]=='y')
	fprintf(out, " %d %g %.2f %g %g %g",
		(int)j, ctl->nu[cida[j]],
		obs->time[cira[j]], obs->vpz[cira[j]],
		obs->vplon[cira[j]], obs->vplat[cira[j]]);
      else {
	idx2name(ctl, ciqa[j], quantity);
	fprintf(out, " %d %s %.2f %g %g %g", (int)j, quantity,
		atm->time[cipa[j]], atm->z[cipa[j]],
		atm->lon[cipa[j]], atm->lat[cipa[j]]);
      }
      
      /* Write matrix entry... */
      fprintf(out, " %g\n", gsl_matrix_get(matrix, i, j));
    }
    
    /* Set matrix indices... */
    if(sort[0]=='r') {
      j++;
      if(j>=nc) {
	j=0;
	i++;
	fprintf(out, "\n");
      }
    } else {
      i++;
      if(i>=nr) {
	i=0;
	j++;
	fprintf(out, "\n");
      }
    }
  }
  
  /* Close file... */
  fclose(out);
}

/*****************************************************************************/

void x2atm(ctl_t *ctl,
	   gsl_vector *x,
	   atm_t *atm,
	   aero_t *aero) {
  
  int ig, iw; 
  
  size_t n=0;
  
  /* Set pressure... */
  x2atm_help(atm, ctl->retp_zmin, ctl->retp_zmax, atm->p, x, &n);
  
  /* Set temperature... */
  x2atm_help(atm, ctl->rett_zmin, ctl->rett_zmax, atm->t, x, &n);
  
  /* Set volume mixing ratio... */
  for(ig=0; ig<ctl->ng; ig++)
    x2atm_help(atm, ctl->retq_zmin[ig], ctl->retq_zmax[ig], atm->q[ig], x, &n);
  
  /* Set extinction... */
  for(iw=0; iw<ctl->nw; iw++)
    x2atm_help(atm, ctl->retk_zmin[iw], ctl->retk_zmax[iw], atm->k[iw], x, &n);
 
  /* Set particle concentration... */
  if(ctl->retnn) {
    aero->nn[0]=gsl_vector_get(x, n);
    n++;
  }
 
  /* Set particle mode radius... */
  if(ctl->retrr) {
    aero->rr[0]=gsl_vector_get(x, n);
    n++;
  }
 
  /* Set particle size distribution width... */
  if(ctl->retss) {
    aero->ss[0]=gsl_vector_get(x, n);
    n++;
  }
}

/*****************************************************************************/

void x2atm_help(atm_t *atm,
		double zmin,
		double zmax,
		double *value,
		gsl_vector *x,
		size_t *n) {
  
  int ip;

  /* Extract state vector elements... */
  for(ip=0; ip<atm->np; ip++){
    if(atm->z[ip]>=zmin && atm->z[ip]<=zmax) {
      value[ip]=gsl_vector_get(x, *n);
      (*n)++;
    }
  }
}

/*****************************************************************************/

void y2obs(ctl_t *ctl,
	   gsl_vector *y,
	   obs_t *obs) {
  
  int id, ir;
  
  size_t m=0;
  
  /* Decompose measurement vector... */
  for(ir=0; ir<obs->nr; ir++)
    for(id=0; id<ctl->nd; id++)
      if(gsl_finite(obs->rad[id][ir])) {
	obs->rad[id][ir]=gsl_vector_get(y, m);
	m++;
      }
}
