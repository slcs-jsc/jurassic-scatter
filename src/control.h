#ifndef CONTROL_H
#define CONTROL_H

#include "jurassic.h"

/* Read in any kind of input data from control file. */

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Read forward model control parameters. */
void read_ctl(int argc,
	      char *argv[],
	      ctl_t *ctl);

/* Search control parameter file for variable entry. */
double scan_ctl(int argc,
		char *argv[],
		const char *varname,
		int arridx,
		const char *defvalue,
		char *value);


#endif
