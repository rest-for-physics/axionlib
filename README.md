# RestAxion

A REST library used to generate solar axions and obtain the detection probability function. The library allows to define a generic helioscope setup - buffer gas, magnetic field, optics response, photon transmission - through REST metadata structures.

### Prerequisites

As any REST library, RestAxionLib requires a running installion of REST Framework. But on top of that, some calculations require higher precision arithmetics and we need to use a external library named `mpfr`. We use a c++ wrapper that is available at the [following site](http://www.holoborodko.com/pavel/mpfr/#intro).

All you need to be able to compile RestAxionLib is to install the [mpfr](https://www.mpfr.org) and [mpir](http://mpir.org) libraries. Download the source, compile, and install. Usually as simple as running `./configure`, `make` and `make install`.

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

