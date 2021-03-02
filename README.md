# jurassic-scatter

The JUelich RApid Spectral SImulation Code with scattering is based on JURASSIC 
(https://github.com/slcs-jsc/jurassic), but was a bit reorganized to include modules 
that allow simulations with scattering on aerosol and cloud particles.

![GitHub tag (latest SemVer)](https://img.shields.io/github/tag/slcs-jsc/jurassic-scatter.svg)
![GitHub top language](https://img.shields.io/github/languages/top/slcs-jsc/jurassic-scatter.svg)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/slcs-jsc/jurassic-scatter.svg)
![GitHub last commit](https://img.shields.io/github/last-commit/slcs-jsc/jurassic-scatter.svg)
![GitHub](https://img.shields.io/github/license/slcs-jsc/jurassic-scatter.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4572986.svg)](https://doi.org/10.5281/zenodo.4572986)

## Short Description

The Juelich Rapid Spectral Simulation Code (JURASSIC) is a fast radiative transfer
model for the mid-infrared spectral region (Hoffmann, 2006). It was used in several
studies for the infrared limb sounder Michelson Interferometer for Passive Atmospheric
Sounding (MIPAS) (Hoffmann et al., 2005, 2008), Cryogenic Infrared Spectrometers
and Telescopes for the Atmosphere - New Frontiers (CRISTA-NF) (Hoffmann et al.,
2009; Weigel et al., 2010), and Gimballed Limb Observer for Radiance Imaging of
the Atmosphere (GLORIA) (Ungermann et al., 2010) and the nadir instrument Atmospheric
Infrared Sounder (AIRS) (Hoffmann and Alexander, 2009; Grimsdell et al.,
2010; Hoffmann et al., 2013).

For fast simulations, it applies pre-calculated look-up tables of spectral emissivities
and approximations to radiative transfer calculations, such as the emissivity growth
approximation (EGA) (Weinreb and Neuendorffer, 1973; Gordley and Russell, 1981;
Marshall et al., 1994).The look-up-tables were calculated with the Reference Forward
Model (RFM) (Dudhia et al., 2002; Dudhia, 2014), which is an exact line-by-line model
specifically developed for MIPAS. For selected spectral windows, JURASSIC-scatter has been
compared to the line-by-line models RFM and Karlsruhe Optimized and Precise Radiative
transfer Algorithm (KOPRA) (Stiller, 2000; Stiller et al., 2002; Höpfner and
Emde, 2005) and shows good agreement (Griessbach et al., 2013).

JURASSIC-scatter contains a scattering module that allows for radiative transfer simulations
including single scattering on aerosol and cloud particles (Grießbach, 2012; Griessbach et al., 
2013). Forward simulations with scattering on volcanic ash, ice and sulfate aerosol have been 
used to develop and characterise a volcanic ash detection method for MIPAS (Griessbach et al., 
2012, 2014). 

Retrieval of large satellite data sets require plenty of computing time that can be provided 
by supercomputers. JURASSIC-scatter shows best performance using a hybrid MPI/openMP 
parallelization on JURECA and JUWELS at the Jülich Supercomputing Centre (JSC), Forschungszentrum 
Jülich GmbH. It can also be used on workstations or laptops using pure MPI or openMP parallelization.

## Documentation

In the documentation you will find a quick start, examples of how to run simulations, a description of 
the input and output file formats, a list of control flags and their options, as well as examples 
of how to use the modules. The documentation is neither 
ready nor perfect, but it is a start to introduce you to JURASSIC-scatter and to enable you 
to work and to do science with this code package.

## Contributing

We are interested in sharing JURASSIC-scatter for operational or research applications.

Please do not hesitate to contact us, if you have any further questions or need support.

## License

JURASSIC-scatter is distributed under the GNU GPL v3.

## Contact

Dr. Sabine Grießbach

Dr. Lars Hoffmann  

Jülich Supercomputing Centre, Forschungszentrum Jülich
