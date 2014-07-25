#ifndef LINEOFSIGHT_H
#define LINEOFSIGHT_H

#include "jurassic.h"
#include "scatter.h"

/* Do ray-tracing to determine LOS. */
void raytrace(ctl_t *ctl,
	      atm_t *atm,
	      obs_t *obs,
	      aero_t *aero,
	      los_t *los,
	      int ir);

/* Compute refractivity (return value is n - 1). */
double refractivity(double p,
		    double t);

/* Find tangent point of a given LOS. */
void tangent_point(los_t *los,
		   double *tpz,
		   double *tplon,
		   double *tplat);

/* Find ground or TOA intersection point of a LOS. */
void intersection_point(ctl_t *ctl,
			atm_t *atm,
			double *znew,
			los_t *los,
			int ip,
			los_t *los_aero,
			int jp);

/* Add points to LOS for fine sampling of aerosol/cloud layers. */
void add_aerosol_layers(ctl_t *ctl,
			atm_t *atm,
			los_t *los,
			aero_t *aero);

#endif
