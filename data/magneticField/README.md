
File : README.md

Creation Date : 13rd June 2019

---

Description
===========

This data directory contains different definitions of magnetic fields to be read by a future `TRestAxionMagneticField` metadata in order to fill regularly a parallelepiped domain defining a magnetic bore.


Conventions
===========

- `File contents` : Each .dat file contains 6 columns. The three first columns are the position X,Y,Z in m and the three last columns the values of the magnetic field in these positions in T. These files are for regular mesh.

- `Filename` : magneticField_sizeMesh_AUTHOR_YYYYMM.dat. For example, if points are defined every 0.04 meters the file name will be magneticField_004_AUTHOR_YYYMM.dat.

- `X, Y, Z definitions` : The origin is set in the center of the domain. The Z direction is along the magnetic bore and the axis increases from the entry to the end of the magnetic bore. The X and Y axis are defined such that (X,Y,Z) is a direct system.


List of magnetic fields available 
=================================

- `MagneticField_Bykovskiy_201906.dat` 
