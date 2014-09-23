#ifndef SONSTIGE_H
#define SONSTIGE_H

#include "jurassic.h"
#include "scatter.h"
#include "forwardmodel.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Convert Cartesian coordinates to geolocation. */
void cart2geo(double *x,
	      double *z,
	      double *lon,
	      double *lat);

/* Copy and initialize observation data. */
void copy_obs(ctl_t *ctl,
	      obs_t *obs_dest,
	      obs_t *obs_src,
	      int init);

/* Convert geolocation to Cartesian coordinates. */
void geo2cart(double z,
	      double lon,
	      double lat,
	      double *x);

/* Initialize look-up tables. */
void init_tbl(ctl_t *ctl,
	      tbl_t *tbl);

/* Find array index. */
int locate(double *xx,  /* array */
	   int n,       /* array size */ 
	   double x);   /* value */

/* Read observation data. */
/* Reads observations e.g for retrieval */
void read_obs(const char *dirname,
	      const char *filename,
	      ctl_t *ctl,
	      obs_t *obs);

/* Write observation data. */
void write_obs(const char *dirname,
	       const char *filename,
	       ctl_t *ctl,
	       obs_t *obs);

#endif
