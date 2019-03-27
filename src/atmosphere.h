#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

/* Module for atmosphere related functions. */

#include "jurassic.h"
#include "misc.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Compose state vector or parameter vector. */
size_t atm2x(ctl_t *ctl,
	     atm_t *atm,
	     aero_t *aero,
	     gsl_vector *x,
	     int *iqa,
	     int *ipa);

/* Add elements to state vector. */
void atm2x_help(atm_t *atm,
		double zmin,
		double zmax,
		double *value,
		int val_iqa,
		gsl_vector *x,
		int *iqa,
		int *ipa,
		size_t *n);

/* Copy and initialize atmospheric data. */
void copy_atm(ctl_t *ctl,
	      atm_t *atm_dest,
	      atm_t *atm_src,
	      int init);

/* Find index of an emitter. */
int find_emitter(ctl_t *ctl,
		 const char *emitter);

/* Determine gravity of Earth. */
double gravity(double z, 
	       double lat);

/* Set hydrostatic equilibrium. */
void hydrostatic(ctl_t *ctl,
		 atm_t *atm);

/* Set hydrostatic equilibrium for individual profile. */
void hydrostatic_1d(ctl_t *ctl,
		    atm_t *atm,
		    int ip0,
		    int ip1);
/* Interpolate complete atmospheric data set. */
void intpol_atm(ctl_t *ctl,
		atm_t *atm_dest,
		atm_t *atm_src);

/* Interpolate atmospheric data for given geolocation. */
void intpol_atm_geo(ctl_t *ctl,
		    atm_t *atm,
		    double z0,
		    double lon0,
		    double lat0,
		    double *p,
		    double *t,
		    double *q,
		    double *k);

/* Interpolate 1D atmospheric data (vertical profile). */
void intpol_atm_1d(ctl_t *ctl,
		   atm_t *atm,
		   int idx0,
		   int n,
		   double z0,
		   double *p,
		   double *t,
		   double *q,
		   double *k);

/* Interpolate 2D atmospheric data (satellite track). */
void intpol_atm_2d(ctl_t *ctl,
		   atm_t *atm,
		   double z0,
		   double lon0,
		   double lat0,
		   double *p,
		   double *t,
		   double *q,
		   double *k);

/* Interpolate 3D atmospheric data (Lagrangian grid). */
void intpol_atm_3d(ctl_t *ctl,
		   atm_t *atm,
		   double z0,
		   double lon0,
		   double lat0,
		   double *p,
		   double *t,
		   double *q,
		   double *k);

/* Read atmospheric data. */
void read_atm(const char *dirname,
	      const char *filename,
	      ctl_t *ctl,
	      atm_t *atm);

/* Write atmospheric data. */
void write_atm(const char *dirname,
	       const char *filename,
	       ctl_t *ctl,
	       atm_t *atm);

#endif
