FSPS: Flexible Stellar Population Synthesis
=====
![Version Badge](https://img.shields.io/badge/version-v3.2-blue)

References
---------
When using this code please cite the following papers:
 * [Conroy, Gunn, & White 2009, ApJ, 699, 486](https://ui.adsabs.harvard.edu/abs/2009ApJ...699..486C)
 * [Conroy & Gunn 2010, ApJ, 712, 833](https://ui.adsabs.harvard.edu/abs/2010ApJ...712..833C)

Installation
----------
If you have Git installed, FSPS can be obtained with the following commands:

```sh
cd /path/to/desired/location/
git clone https://github.com/cconroy20/fsps
```
Otherwise download a gzipped tarball from [here](https://github.com/cconroy20/fsps/releases). Then follow the instructions at [`doc/INSTALL`](doc/INSTALL).

You should not need to update the Git repository until an update is announced (which is why you need to be on the mailing list - see [`doc/INSTALL`](doc/INSTALL)).  If you've obtained FSPS using Git then when an update is announced you will need to simply type `cd $SPS_HOME; git pull` and then recompile (just type `make`).  If you have made your own edits to the FSPS files, Git will attempt to gracefully merge your local version with the repository version.

Documentation
------
See the [Manual](doc/MANUAL.pdf)

## C Driver API

FSPS ships a C driver intended for use by language bindings (e.g., Python-FSPS refactor or a Julia wrapper).

- Header: include/fsps.h
- Shared library: build/libfsps.so (build via `make shared`)
- C driver test: `make test_c`

For pkg-config users, a sample file is provided at fsps.pc (update `prefix` as needed).

See [doc/FSPS_C_API.md](doc/FSPS_C_API.md) for the API reference, array layout, and error handling.

Contents
---------
Below is a brief description of the contents of the directories in the
fsps root directory:

 * [`OUTPUTS`](OUTPUTS): Contains the outputs of a few example calls of the routines
autosps and simple.  You may wish to use this directory for all
outputs of the fsps routines.

 * `build`: Contains the compiled object files (`.o`) and Fortran modules (`.mod`).
This directory is automatically created when you run `make`.

 * [`data`](data): Contains model inputs and lookup tables. Subdirectories include:
	 - `data/isochrones`: isochrone libraries (BaSTI, Padova, MIST, etc.)
	 - `data/spectra`: spectral libraries and reference spectra (A0V, Sun, Hot_spectra)
	 - `data/nebular`: nebular emission tables (continuum + lines)
	 - `data/dust`: dust attenuation/emission and AGN torus models
	 - `data/`: filter lists, indices, IMFs, SFHs, and other small tables

 * [`doc`](doc): Contains the manual, revision history, and installation
instructions.

* [`pro`](pro): Contains IDL files for reading in the .mag, .indx, and .spec
output files

* [`src`](src): Contains the Fortran sources organized by subsystem (core, physics,
spectra, SFH, math, cosmology, ABI, and program entry points).

* [`tests`](tests): Contains the regression test suite and scripts for generating
reference comparison data.
