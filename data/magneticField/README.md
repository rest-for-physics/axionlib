
File : README.md

Creation Date : 13rd June 2019

---

Description
===========

This data directory contains different definitions of magnetic fields to be read by a future `TRestAxionMagneticField` metadata in order to fill regularly a parallelepiped domain defining a magnetic bore.


Conventions
===========

- `File contents` : Each .dat file contains a first line with 4 numbers (x_max, y_max, z_max and meshSize in mm) and then 6 columns with many lines. The three first columns are the position X,Y,Z in mm and the three last columns the values of the magnetic field in these positions in T. These files are for regular mesh and for rectanglar magnetic field.
All positions are in mm wich are the units by default of REST !

- `Filename` : AUTHOR_YYYYMM.dat.

- `X, Y, Z definitions` : The origin is set in the center of the domain, so x_max=-x_min (the same for y and z). The Z direction is along the magnetic bore and the axis increases from the entry to the end of the magnetic bore. The X and Y axis are defined such that (X,Y,Z) is a direct system.


List of magnetic fields available 
=================================

- `Bykovskiy_201906.dat` 
