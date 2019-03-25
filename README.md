# RestAxion

A REST library used to generate solar axions and obtain the detection probability function. The library allows to define a generic helioscope setup - buffer gas, magnetic field, optics response, photon transmission - through REST metadata structures.

Installation
============

REST libraries must be installed and loaded in the system.

You may check REST libraries are available in the system by doing.
```
rest-config --welcome
```

Then, just proceed to create the build directory, compile and install.

We assume here that you are now at the RestAxionLib directory.
```
mkdir build
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

