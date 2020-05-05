# RestAxion

A REST library used to generate solar axions and obtain the detection probability function. The library allows to define a generic helioscope setup - buffer gas, magnetic field, optics response, photon transmission - through REST metadata structures.

### Retrieving RestAxionLib with sub-modules

The data directory inside RestAxionLib is an idependent projec repository installed as a submodule. In order to download the RestAxionLib project, and all the sub-modules we need to clone using the following.

```
git clone --recursive git@lfna.unizar.es:iaxo/RestAxionLib.git
```

If `--recursive` is not used the directory `/data` will not be populated. The `data` repository will not be updated as often as the code repository, but to make sure we get the latest data we should pull the code using the following command.

```
git pull --recurse-submodules
```

If the repository was cloned without `--recursive` option, then it is possible to download the modules later using

```
git submodule update --init --recursive
```


### Prerequisites

As any REST library, RestAxionLib requires a running installion of REST Framework. But on top of that, some calculations require higher precision arithmetics and we need to use a external library named `mpfr`. 
We use a c++ wrapper that is available at the [following site](http://www.holoborodko.com/pavel/mpfr/#intro). This wrapper (mpreal.h) is already available/uploaded to this repository.

All you need to be able to compile RestAxionLib is to install the [mpfr](https://www.mpfr.org) and [mpir](http://mpir.org) libraries. Download the source, compile, and install. 
Usually as simple as running `./configure`, `make` and `make install` at the downloaded source directory.

#### mpfr and mpir installation:

Download the file `mpfr-4.0.2.tar.gz` from the following [site](https://www.mpfr.org/mpfr-current/#download) compile and install executing the following commands.

```
tar -zxvf mpfr-4.0.2.tar.gz
cd mpfr-4.0.2/
./configure
make
sudo make install (to do a system installation)
```

In order to install mpir, download the file `mpir-3.0.0.tar.bz2` from the following [link](http://mpir.org/downloads.html) into the `~/mpir-3.0.0` directory

2) Then inside `~/mpir` directory execute the following commands:
```
cd mpir-3.0.0
tar xvjf mpir-3.0.0.tar.bz2
cd mpir-3.0.0/
./configure
make
sudo make install
```

### Installation

REST libraries must be installed in the system.

You may check REST libraries are available in the system by doing.
```
rest-config --welcome
```

Then, just proceed to create the build directory, compile and install.

We assume here that you are now at the RestAxionLib directory.
```
mkdir build
cd build
cmake ../
make -j4 install
```

This will install the library, headers and data at your REST_PATH 
install directory. If you have no write access permissions to that directory
you may change the install directory at the cmake command.

By executing the following commands at your build directory

```
cmake -DCMAKE_INSTALL_PREFIX=YOUR_INSTALL_DIR ../
make -j4 install
```



Adding new metadata or processes to this library
================================================
TOBE written

