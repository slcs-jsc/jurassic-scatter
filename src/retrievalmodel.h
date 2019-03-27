#ifndef RETRIEVALMODEL_H
#define RETRIEVALMODEL_H

/* #include "jurassic.h" */
#include "control.h"
#include "forwardmodel.h"

/* Compute information content and resolution. */
void analyze_avk(ret_t *ret,
		 ctl_t *ctl,
		 atm_t *atm,
		 aero_t *aero,
		 int *iqa,
		 int *ipa,
		 gsl_matrix *avk);

/* Analyze averaging kernels for individual retrieval target. */
void analyze_avk_quantity(gsl_matrix *avk,
			  int iq,
			  int *ipa,
			  size_t *n0,
			  size_t *n1,
			  double *cont,
			  double *res);

/* Compute correlations based on spatial distance. */
double corr_function(double z0,
		     double lon0,
		     double lat0,
		     double z1,
		     double lon1,
		     double lat1,
		     double cz,
		     double ch);

/* Compute cost function. */
double cost_function(FILE *out,
                     int it,
                     gsl_vector *dx,
                     gsl_vector *dy,
                     gsl_matrix *s_a_inv,
                     gsl_vector *sig_eps_inv);

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

/* Carry out optimal estimation retrieval. */
void optimal_estimation(ret_t *ret,
			ctl_t *ctl,
			obs_t *obs_meas,
			obs_t *obs_i,
			atm_t *atm_apr,
			atm_t *atm_i,
			aero_t *aero_apr,
			aero_t *aero_i);

/* Read retrieval control parameters. */
void read_ret(int argc,
	      char *argv[],
	      ctl_t *ctl,
	      ret_t *ret);

/* Set a priori covariance. */
void set_cov_apr(ret_t *ret,
		 ctl_t *ctl,
		 atm_t *atm,
		 aero_t *aero,
		 int *iqa,
		 int *ipa,
		 gsl_matrix *s_a);

/* Set measurement errors. */
void set_cov_meas(ret_t *ret,
		  ctl_t *ctl,
		  obs_t *obs,
		  gsl_vector *sig_noise,
		  gsl_vector *sig_formod,
		  gsl_vector *sig_eps_inv);

/* Write retrieval error to file. */
void write_stddev(const char *quantity,
		  ret_t *ret,
		  ctl_t *ctl,
		  atm_t *atm,
		  aero_t *aero,
		  gsl_matrix *s);

/* Write matrix. */
void write_matrix(const char *dirname,
		  const char *filename,
		  ctl_t *ctl,
		  gsl_matrix *matrix,
		  atm_t *atm,
		  aero_t *aero,
		  obs_t *obs,
		  const char *rowspace,
		  const char *colspace,
		  const char *sort);

/* Decompose parameter vector or state vector. */
void x2atm(ctl_t *ctl,
	   gsl_vector *x,
	   atm_t *atm,
	   aero_t *aero);

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
