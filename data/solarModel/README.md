
File : README.md

Creation Date : 4th April 2019

Database manager : Maurizio Gianotti

---

Description
===========

This data directory contains the description of different solar axion models that can be accessed using the `TRestAxionSolarModel` object. To learn on how to use `TRestAxionSolarModel` visit the file containing the source code, TRestAxionSolarModel.cxx.

Each axion solar model (or production mechanism) provides the pre-integrated(+) solar disc axion flux as a function of the solar radius and energy.

(+) The effective solar axion flux is integrated to the corresponding solar disk (ring) area.

Basic conventions
=================

- `File contents` : Each value inside the file corresponds to the calculated axion flux on Earth for a given Earth-Sun distance. Sun-Earth distance = XXX units. 

- `Filename` : MODEL_AUTHOR_YYYYMM.dat, where MODEL stands for the solar axion production mechanism.

- `Units` : The solar flux given in the tables are measured in cm-2 s-1 keV-1

- `File format` : The format of the tables is fixed to 100 lines. Each line inside the file corresponds to one integrated solar ring, the first line corresponding to the inner ring (Rtop = 0.01 x Rsun), and the last line corresponding to the outer ring (Rbottom=0.99 x Rsun). Each column value corresponds to the integrated energy spectrum in steps of 100 eV. The first value in each row corresponds to the contribution to the solar axion flux from 0 to 100 eV. The energy range goes to 20 keV.


List of solar models available
==============================

- `Primakoff_Gianotti_201904.dat` : TOBE WRITTEN. References and notes.

