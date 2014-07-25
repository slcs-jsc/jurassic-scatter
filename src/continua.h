#ifndef CONTINUA_H
#define CONTINUA_H

#include "jurassic.h"
#include "atmosphere.h"

/* Compute carbon dioxide continuum (optical depth). */
double ctmco2(ctl_t *ctl,
	      double nu,
	      double p,
	      double t,
	      double *u);

/* Compute water vapor continuum (optical depth). */
double ctmh2o(ctl_t *ctl,
	      double nu,
	      double p,
	      double t,
	      double *q,
	      double *u);

/* Compute nitrogen continuum (absorption coefficient). */
double ctmn2(double nu,
	     double p,
	     double t);

/* Compute oxygen continuum (absorption coefficient). */
double ctmo2(double nu,
	     double p,
	     double t);

#endif
