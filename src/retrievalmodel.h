#ifndef RETRIEVALMODEL_H
#define RETRIEVALMODEL_H

#include "jurassic.h"
#include "forwardmodel.h"

/* Determine name of state vector quantity for given index. */
void idx2name(ctl_t *ctl,
	      int idx,
	      char *quantity);

/* Compute Jacobians. */
void kernel(ctl_t *ctl,
	    atm_t *atm,
	    obs_t *obs,
	    aero_t *aero,
	    gsl_matrix *k);

/* Invert symmetric matrix. */
void matrix_invert(gsl_matrix *a);

/* Compute matrix product A^TBA or ABA^T for diagonal matrix B. */
void matrix_product(gsl_matrix *a,
		    gsl_vector *b,
		    int transpose,
		    gsl_matrix *c);

/* Compose measurement vector. */
size_t obs2y(ctl_t *ctl,
	     obs_t *obs,
	     gsl_vector *y,
	     int *ida,
	     int *ira);

/* Write matrix. */
void write_matrix(const char *dirname,
		  const char *filename,
		  ctl_t *ctl,
		  gsl_matrix *matrix,
		  atm_t *atm,
		  obs_t *obs,
		  const char *rowspace,
		  const char *colspace,
		  const char *sort);

/* Decompose parameter vector or state vector. */
void x2atm(ctl_t *ctl,
	   gsl_vector *x,
	   atm_t *atm);

/* Extract elements from state vector. */
void x2atm_help(atm_t *atm,
		double zmin,
		double zmax,
		double *value,
		gsl_vector *x,
		size_t *n);

/* Decompose measurement vector. */
void y2obs(ctl_t *ctl,
	   gsl_vector *y,
	   obs_t *obs);


#endif
