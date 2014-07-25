#ifndef JURASSIC_H
#define JURASSIC_H

#include <gsl/gsl_blas.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_complex_math.h>
#include <math.h>
#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------
   Macros...
   ------------------------------------------------------------ */

/* Allocate memory. */
#define ALLOC(ptr, type, n)				\
  if((ptr=malloc((size_t)(n)*sizeof(type)))==NULL)	\
    ERRMSG("Out of memory!");

/* Compute angle between two vectors. */
#define ANGLE(a, b)						\
  acos(GSL_MIN(GSL_MAX(DOTP(a, b)/NORM(a)/NORM(b), -1), 1))

/* Compute Cartesian distance between two vectors. */
#define DIST(a, b) sqrt(DIST2(a, b))

/* Compute squared distance between two vectors. */
#define DIST2(a, b)							\
  ((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))

/* Compute dot product of two vectors. */
#define DOTP(a, b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

/* Print error message and quit program. */
#define ERRMSG(msg) {                                                   \
    printf("\nError (%s, %s, l%d): %s\n\n",                             \
           __FILE__, __FUNCTION__, __LINE__, msg);                      \
    exit(EXIT_FAILURE);                                                 \
  }

/* Compute exponential interpolation. */
#define EXP(x0, y0, x1, y1, x)					\
  (((y0)>0 && (y1)>0)						\
   ? ((y0)*exp(log((y1)/(y0))/((x1)-(x0))*((x)-(x0))))          \
   : LIN(x0, y0, x1, y1, x))

/* Read binary data. */
#define FREAD(ptr, type, nmemb, stream) {                               \
    if(fread(ptr, sizeof(type), (size_t)nmemb, stream)!=(size_t)nmemb)  \
      ERRMSG("Error while reading!");                                   \
  }

/* Write binary data. */
#define FWRITE(ptr, type, nmemb, stream) {				\
    if(fwrite(ptr, sizeof(type), (size_t)nmemb, stream)!=(size_t)nmemb)	\
      ERRMSG("Error while writing!");					\
  }

/* Compute linear interpolation. */
#define LIN(x0, y0, x1, y1, x)			\
  ((y0)+((y1)-(y0))/((x1)-(x0))*((x)-(x0)))

/* Write log message. */
#define LOGMSG(lev, cmd) {if(lev<=VERBLEV) cmd;}

/* Execute netCDF library command and check result. */
#define NC(cmd) {				     \
    if((cmd)!=NC_NOERR)				     \
      ERRMSG(nc_strerror(cmd));			     \
  }

/* Compute norm of a vector. */
#define NORM(a) sqrt(DOTP(a, a))

/* Print macro for debugging. */
#define PRINT(format, var)                                              \
  printf("Print (%s, %s, l%d): %s= "format"\n",                         \
         __FILE__, __FUNCTION__, __LINE__, #var, var);

/* Read string tokens. */
#define TOK(line, tok, format, var) {			\
    if(((tok)=strtok((line), " \t"))) {			\
      if(sscanf(tok, format, &(var))!=1) continue;	\
    } else ERRMSG("Error while reading!");		\
  }

/* ------------------------------------------------------------
   Redefinable constants...
   ------------------------------------------------------------ */

/* Verbosity level. */
#ifndef VERBLEV
#define VERBLEV 2
#endif

/* ------------------------------------------------------------
   Constants...
   ------------------------------------------------------------ */

/* First spectroscopic constant (c_1 = 2 h c^2) [W/(m^2 sr cm^-4)]. */
#define C1 1.19104259e-8

/* Second spectroscopic constant (c_2 = h c / k) [K/cm^-1]. */
#define C2 1.43877506

/* Standard pressure [hPa]. */
#define P0 1013.25

/* Mean radius of Earth [km]. */
#define RE 6367.421

/* Mass of Earth [kg]. */
#define ME 5.976e24

/* Temperature of the Sun [K]. */
#define TSUN 5780.

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/* Maximum length of ASCII data lines. */
#define LEN 5000

/* Maximum size of measurement vector. */
#define MMAX (NRMAX*NDMAX)

/* Maximum size of state vector. */
#define NMAX (NQMAX*NPMAX)

/* Maximum number of quantities. */
#define NQMAX (2+NGMAX+NWMAX)

/* Maximum number of spectral windows. */
#define NWMAX 5

/* Maximum number of radiance channels. */
#define NDMAX 50

/* Maximum number of emitters. */
#define NGMAX 15

/* Maximum number of LOS points. */
#define NLOS 1000

/* Maximum number of atmospheric data points. */
#define NPMAX 1000

/* Maximum number of ray paths. */
#define NRMAX 1000

/* Maximum number of shape function grid points. */
#define NSHAPE 10000

/* Number of ray paths used for FOV calculations. */
#define NFOV 5

/* Maximum number of pressure levels in emissivity tables. */
#define TBLNPMAX 40

/* Maximum number of source function temperature levels. */
#define TBLNSMAX 1201

/* Maximum number of temperatures in emissivity tables. */
#define TBLNTMAX 30

/* Maximum number of column densities in emissivity tables. */
#define TBLNUMAX 320

/* Maximum number of scattering models. */
#define SCAMOD 15

/* Maximum number of aerosol/cloud layers */
#define NLMAX 10

/* Number of scattering angles (from 0 to 180 deg). */
#define NTHETA 181

/* Number of points for Gauss-Hermite integration. */
#define NRAD 170

/* Maximum number of refractive indices. */
#define REFMAX 5000

/* ------------------------------------------------------------
   Quantity indices...
   ------------------------------------------------------------ */

/* Pressure. */
#define IDXP 0

/* Temperature. */
#define IDXT 1

/* Volume mixing ratios. */
#define IDXQ(ig) (2+ig)

/* Extinction. */
#define IDXK(iw) (2+ctl->ng+iw)

/* ------------------------------------------------------------
   Global Structs...
   ------------------------------------------------------------ */
/* Forward model control parameters. */
typedef struct {
  
  /* Number of emitters. */
  int ng;
  
  /* Name of each emitter. */
  char emitter[NGMAX][LEN];
  
  /* Number of radiance channels. */
  int nd;

  /* Number of spectral windows. */
  int nw;
  
  /* Centroid wavenumber of each channel [cm^-1]. */
  double nu[NDMAX];

  /* Window index of each channel. */
  int window[NDMAX];
  
  /* Basename for table files and filter function files. */
  char tblbase[LEN];
  
  /* Reference height for hydrostatic pressure profile (-999 to skip) [km]. */
  double hydz;
  
  /* Compute CO2 continuum (0=no, 1=yes). */
  int ctm_co2;
  
  /* Compute H2O continuum (0=no, 1=yes). */
  int ctm_h2o;
  
  /* Compute N2 continuum (0=no, 1=yes). */
  int ctm_n2;
  
  /* Compute O2 continuum (0=no, 1=yes). */
  int ctm_o2;
  
  /* Number of scattering models. */
  int sca_n;
  
  /* Number of recursions for multiple scattering. */
  /* (0=no scattering, 1=single scattering, 2<=multiple scattering) */
  int sca_mult;

  /* Extinction coefficient type if sca_mult=0 */
  char sca_ext[LEN];

  /* Interpolation method (1=profile, 2=satellite track, 3=Lagrangian grid). */
  int ip;

  /* Influence length for vertical interpolation [km]. */
  double cz;
  
  /* Influence length for horizontal interpolation [km]. */
  double cx;
  
  /* Take into account refractivity (0=no, 1=yes). */
  int refrac;
  
  /* Maximum step length for raytracing [km]. */
  double rayds;
  
  /* Vertical step length for raytracing [km]. */
  double raydz;

  /* Sampling step for transition layers [km]. */
  double transs;
  
  /* Field-of-view data file. */
  char fov[LEN];
  
  /* Minimum altitude for pressure retrieval [km]. */
  double retp_zmin;
  
  /* Maximum altitude for pressure retrieval [km]. */
  double retp_zmax;
  
  /* Minimum altitude for temperature retrieval [km]. */
  double rett_zmin;
  
  /* Maximum altitude for temperature retrieval [km]. */
  double rett_zmax;
  
  /* Minimum altitude for volume mixing ratio retrieval [km]. */
  double retq_zmin[NGMAX];
  
  /* Maximum altitude for volume mixing ratio retrieval [km]. */
  double retq_zmax[NGMAX];
  
  /* Minimum altitude for extinction retrieval [km]. */
  double retk_zmin[NWMAX];
  
  /* Maximum altitude for extinction retrieval [km]. */
  double retk_zmax[NWMAX];
  
  /* Use brightness temperature instead of radiance (0=no, 1=yes). */
  int write_bbt;
  
  /* Write matrix data (0=no, 1=yes). */
  int write_matrix;
  
} ctl_t;
/* ------------------------------------------------------------*/

/* Line-of-sight data. */
typedef struct {
  
  /* Number of LOS points. */
  int np;
  
  /* Altitude [km]. */
  double z[NLOS];
  
  /* Longitude [deg]. */
  double lon[NLOS];
  
  /* Latitude [deg]. */
  double lat[NLOS];
  
  /* Pressure [hPa]. */
  double p[NLOS];
  
  /* Temperature [K]. */
  double t[NLOS];
  
  /* Volume mixing ratio. */
  double q[NLOS][NGMAX];

  /* Extinction [1/km]. */
  double k[NLOS][NWMAX];

  /* Aerosol/cloud layer index */
  int aeroi [NLOS];

  /* Aerosol/cloud layer scaling factor for transition layer */
  double aerofac[NLOS];
 
  /* Surface temperature [K]. */
  double tsurf;
  
  /* Segment length [km]. */
  double ds[NLOS];
  
  /* Column density [molecules/cm^2]. */
  double u[NLOS][NGMAX];
  
} los_t;
/* ------------------------------------------------------------*/

/* Observation geometry and radiance data. */
typedef struct {
  
  /* Number of ray paths. */
  int nr;
  
  /* Time (seconds since 2000-01-01T00:00Z). */
  double time[NRMAX];
  
  /* Observer altitude [km]. */
  double obsz[NRMAX];
  
  /* Observer longitude [deg]. */
  double obslon[NRMAX];
  
  /* Observer latitude [deg]. */
  double obslat[NRMAX];
  
  /* View point altitude [km]. */
  double vpz[NRMAX];
  
  /* View point longitude [deg]. */
  double vplon[NRMAX];
  
  /* View point latitude [deg]. */
  double vplat[NRMAX];
  
  /* Tangent point altitude [km]. */
  double tpz[NRMAX];
  
  /* Tangent point longitude [deg]. */
  double tplon[NRMAX];
  
  /* Tangent point latitude [deg]. */
  double tplat[NRMAX];
  
  /* Transmittance of ray path. */
  double tau[NDMAX][NRMAX];
  
  /* Radiance [W/(m^2 sr cm^-1)]. */
  double rad[NDMAX][NRMAX];
  
} obs_t;
/* ------------------------------------------------------------*/

/* Emissivity look-up tables. */
typedef struct {
  
  /* Number of pressure levels. */
  int np[NGMAX][NDMAX];
  
  /* Number of temperatures. */
  int nt[NGMAX][NDMAX][TBLNPMAX];
  
  /* Number of column densities. */
  int nu[NGMAX][NDMAX][TBLNPMAX][TBLNTMAX];
  
  /* Pressure [hPa]. */
  double p[NGMAX][NDMAX][TBLNPMAX];
  
  /* Temperature [K]. */
  double t[NGMAX][NDMAX][TBLNPMAX][TBLNTMAX];
  
  /* Column density [molecules/cm^2]. */
  float u[NGMAX][NDMAX][TBLNPMAX][TBLNTMAX][TBLNUMAX];
  
  /* Emissivity. */
  float eps[NGMAX][NDMAX][TBLNPMAX][TBLNTMAX][TBLNUMAX];
  
  /* Source function temperature [K]. */
  double st[TBLNSMAX];
  
  /* Source function radiance [W/(m^2 sr cm^-1)]. */
  double sr[NDMAX][TBLNSMAX];

} tbl_t;
/* ------------------------------------------------------------*/

/* Aerosol and Cloud optical properties. */
typedef struct { 

  /* Number of aerosol/cloud layers */
  int nl;

  /* Number of modes per layer */
  int nm[NLMAX];

  /* Layer top altitude [km] */
  double top[NLMAX];

  /* Layer bottom altitude [km] */
  double bottom[NLMAX];

  /* Transition layer thickness [km] */
  double trans[NLMAX];

  /* Extinction coefficient [1/km]*/
  double beta_e[NLMAX][NDMAX];

  /* Scattering coefficient [1/km]*/
  double beta_s[NLMAX][NDMAX];

  /* Absorption coefficient [1/km]*/
  double beta_a[NLMAX][NDMAX];
  
  /* Phase function for each layer, angle and wave number */
  double p[NLMAX][NDMAX][NTHETA];

 } aero_t;
/* ------------------------------------------------------------*/

 /* Aerosol and cloud file input data. */
  typedef struct{
    
    /* Number of aerosol/cloud models */
    int nm;

    /* Layer top altitude [km] */
    double top[SCAMOD];
    
    /* Layer bottom altitude [km] */
    double bottom[SCAMOD];
    
    /* Transition layer thickness [km] */
    double trans[SCAMOD];
    
    /* Optical properties source */
    char type[SCAMOD][LEN];

    /* Refractive index file or optical properties file */
    char filepath[SCAMOD][LEN];

    /* Number concentration [cm-3] or extinction coefficient [km-1] */
    double nn[SCAMOD];

    /* Median radius of log-normal size distribution [mum] */
    double rr[SCAMOD];
    
    /* Width of log-normal size distribution */
    double ss[SCAMOD];
    
  } aero_i;
/* ------------------------------------------------------------*/

/* Atmospheric data. */
typedef struct { 

   /* Number of data points. */
   int np; 

   /* Time (seconds since 2000-01-01T00:00Z). */
   double time[NPMAX]; 

   /* Altitude [km]. */
   double z[NPMAX];

   /* Longitude [deg]. */
   double lon[NPMAX]; 

   /* Latitude [deg]. */
   double lat[NPMAX]; 

   /* Pressure [hPa]. */
   double p[NPMAX];

   /* Temperature [K]. */
   double t[NPMAX];

   /* Volume mixing ratio. */
   double q[NGMAX][NPMAX];

   /* Extinction [1/km]. */
   double k[NWMAX][NPMAX];

   /* Init flag for interpolation (0=no, 1=yes). */
   int init;

 } atm_t;

#endif
