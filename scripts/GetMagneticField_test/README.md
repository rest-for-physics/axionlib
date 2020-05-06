The scripts in this directory are used to test the implementation of the trilinear interpolation in the `TRestAxionMagneticField::GetMagneticField` method
which is used to return the values of the magnetic field components at the specific point in the magnetic field volume given as the input parameter.

### Description of files in this directory

1) **create_magnetic_field.cxx** - a restRoot script that produces the magnetic field map and stores it in the file "Magnetic_field.dat". The map contains the values of the magnetic field compoments calculated at a number of grid points in a volume defined by:
```
-350 <= x <= 350 mm
-350 <= y <= 350 mm
-5000 <= z <= 5000 mm
```
The distance between the neighbouring grid points (i.e. the mesh size) is 50 mm in x, y and z directions.
At each grid point described by coordinates (x, y, z) the values of the magnetic field are given by:

```
Bx = 5x -2y +2z
By = 8x +5y -3z
Bz = -4x -4y +z
```

2) **Magnetic_field.dat**: this file contains the magnetic field map produced by the "create_magnetic_field.cxx" script. 
Each row in the file contains data for one grid point. First three columns contain `x`, `y` and `z` coordinate of the grid point, respectively,
while other three columns are Bx, By and Bz values for the corresponding grid point.

3) **GetMagneticField_test.cxx**: is a restRoot script which tests the `TRestAxionMagneticField::GetMagneticField` method. 
The test is performed by comparing the values of the magnetic field components obtained by the `GetMagneticField` method 
with the expected values, i.e., values calculated by the equations for `Bx`, `By` and `Bz` given above. 
This comparison is performed on a set of 31 points. 
The `TVector3` "offset" variable represents the offset of the magnetic field volume with the respect to the laboratory frame. 
The output of the script is shown on screen and written in the output file `GetMagneticField_test_output.txt`. 

4) **my_metadata.rml** is a config file used to load the magnetic field map by the `GetMagneticField_test.cxx` script. 
The parameters in this file are described at the beginning of the `TRestAxionMagneticField.cxx` file in `RestAxionLib/src` directory. 
Note that the `position` parameter represents the offset of the magnetic field volume with the respect to the laboratory frame 
and it should be equal to the value of the `offset` variable in the `GetMagneticField_test.cxx` script.
