[![DOI](https://zenodo.org/badge/324291710.svg)](http://doi.org/10.5281/zenodo.4528985)
[![pipeline status](https://gitlab.cern.ch/rest-for-physics/axionlib/badges/master/pipeline.svg)](https://gitlab.cern.ch/rest-for-physics/axionlib/-/commits/master)
[![website](https://img.shields.io/badge/user-guide-E8B6FF.svg)](https://rest-for-physics.github.io)
[![api](https://img.shields.io/badge/user-API-FFCA78.svg)](https://sultan.unizar.es/rest/)
[![forum](https://img.shields.io/badge/user-forum-AAFF90.svg)](https://rest-forum.unizar.es/)
[![github pipeline](https://github.com/rest-for-physics/axionlib/actions/workflows/validation.yml/badge.svg)](https://github.com/rest-for-physics/axionlib/commits/master)

This is a REST-for-Physics library used to generate solar axions and obtain the detection probability function. The library allows to define a generic helioscope setup - buffer gas, magnetic field, optics response, photon transmission - through dedicated REST-for-Physics metadata definitions found in this library.

This library is distributed together with the REST-for-Physics framework as a submodule. In order to use it within the framework it is necessary to pull the library submodule. The main user website guide, at the [downloads page](https://rest-for-physics.github.io/downloading.html), will provide accurate instructions on how to download the axionlib and other modules. In short, the following recipe should pull the latest official version:

```
cd rest-framework
python3 pull-submodules.py --onlylibs
```

### Prerequisites

As any REST library, this library requires a running installion of REST-for-Physics Framework. But on top of that, some calculations require higher precision arithmetics and we need to use a external library named `mpfr`.  We use a c++ wrapper that is available at the [following site](http://www.holoborodko.com/pavel/mpfr/#intro). This wrapper (mpreal.h) is already available/uploaded to this repository.

All you need to be able to compile RestAxionLib is to install the [mpfr](https://www.mpfr.org) and [mpir](http://mpir.org) libraries. Download the source, compile, and install. 
Usually as simple as running `./configure`, `make` and `make install` at the downloaded source directory.

#### mpfr and mpir installation:

Download the file `mpfr-4.0.2.tar.gz` from the following [site](https://www.mpfr.org/mpfr-current/#download) compile and install executing the following commands.

```
wget https://www.mpfr.org/mpfr-current/mpfr-4.1.0.tar.gz
tar -zxvf mpfr-4.1.0.tar.gz
cd mpfr-4.1.0/
./configure
make
sudo make install (to do a system installation)
```

In order to install mpir, download the file `mpir-3.0.0.tar.bz2` from the following [link](http://mpir.org/downloads.html) into the `~/mpir-3.0.0` directory

2) Executing the following commands should do the job:

```
wget http://mpir.org/mpir-3.0.0.tar.bz2
tar xvjf mpir-3.0.0.tar.bz2
cd mpir-3.0.0
./configure
make
sudo make install (to do a system installation)
```

### RestAxionLib installation

Once you have all the prerequisites installed you need to add the library at the REST-for-Physics framework compilation stage. Go to the main framework build directory and add the axion library as a compilation option.

```
cd framework/build
cmake -DRESTLIB_AXION=ON -DREST_MPFR=ON ../
make -j4 install
```

It is necessary to enable the `REST_MPFR` cmake option in order to benefit from the MPFR routines integration. If not, classes that use MPFR will not be compiled. `REST_MPFR` might be optionally disabled in case MPFR libraries are not available. However, this will lead to a non-complete installation of axionlib.

### Publications

This repository makes use of the following published codes:
- K. Altenmuller et al, REST-for-Physics, a ROOT-based framework for event oriented data analysis and combined Monte Carlo response, [Computer Physics Communications 273, April 2022, 108281](https://doi.org/10.1016/j.cpc.2021.108281).
- S.Hoof, J.Jaeckel, T.J.Lennert, Quantifying uncertainties in the solar axion flux and their impact on determining axion model parameters, [JCAP09(2021)006](https://doi.org/10.1088/1475-7516/2021/09/006).
- T.Kittelmann, E.Klinkby, E.B.Knudsen, P.Willendrup, X.X.Cai, K.Kanaki, Monte Carlo Particle Lists: MCPL, Computer Physics Communications 218 (2017) 17–42](http://dx.doi.org/10.17632/cby92vsv5g.1).
