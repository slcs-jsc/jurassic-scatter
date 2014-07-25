#include "retrievalmodel.h"

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
  atm2x(ctl, atm, x0, iqa, ipa);
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
    
    /* Check perturbation size... */
    if(h<=0)
      ERRMSG("Cannot set perturbation size!");
      
    /* Disturb state vector element... */
    gsl_vector_memcpy(x1, x0);
    gsl_vector_set(x1, j, gsl_vector_get(x1, j)+h);
    copy_atm(ctl, &atm1, atm, 0);
    copy_obs(ctl, &obs1, obs, 0);
    x2atm(ctl, x1, &atm1);
    
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

void write_matrix(const char *dirname,
		  const char *filename,
		  ctl_t *ctl,
		  gsl_matrix *matrix,
		  atm_t *atm,
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
    nr=atm2x(ctl, atm, NULL, riqa, ripa);
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
    nc=atm2x(ctl, atm, NULL, ciqa, cipa);
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
	   atm_t *atm) {
  
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
  for(ip=0; ip<atm->np; ip++)
    if(atm->z[ip]>=zmin && atm->z[ip]<=zmax) {
      value[ip]=gsl_vector_get(x, *n);
      (*n)++;
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
