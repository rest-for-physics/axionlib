
### Description

The scripts in this directory are used to test the implementation of the `TRestAxionMagneticField::GetVolumeBoundaries` and `TRestAxionMagneticField::GetFieldBoundaries` methods. These methods are used to find the coordinates of the points where the trajectory of the particle intersects the boundary planes of magnetic field regions through which the particle passes. In general, the magnetic field volume can consist of several regions where each region is defined in a separate `<addMagnetVolume...>` line in the configuration file. Also, in some (or all) regions the magnetic field can be zero in the outer parts of the region, i.e., near the borders. In this case, it is possible that the particle, right after entering the region, passes through the part where `B=0` before traversing the section where `B>0`. Also, just before exiting the region, the particle can again pass through the region where `B=0`. Therefore, two methods are used to find the boundary points of the particle trajectory in each region. The `TRestAxionMagneticField::GetVolumeBoundaries` method searches for the points where the trajectory intersects the boundary planes of the region. If two such points (entry point and exit point) are found, their coordinates are returned. The `TRestAxionMagneticField::GetFieldBoundaries` method checks if the particle, right after entering the region, first passes through the section of the region where `B=0` and determines the coordinates of the first point where `B>0` on the particle trajectory after entering the region. In a similar manner, the method checks if just before exiting the region the particle passes through the section where `B=0` and determines the last point on the trajectory where `B>0` before exiting the region. Thus these two points are actually the boundary points of one segment of the particle trajectory along which `B>0`.

To perform the test just execute the command `restRoot -b -q Boundaries_test_interactive.C`.

There is also a non-interactive version of this test `Boundaries_test.C` with four predefined inputs suitable to launch the test as a part of the gitlab pipeline chain.

### Description of the magnetic field volume used in this test

The magntic field volume consists of two identical regions placed one behind the other along the z-direction which is aligned with the axis of the magnet bore. Each region has dimensions:

```
-70 <= x <= +70 mm
-70 <= y <= +70 mm
-1000 <= z <= +1000 mm
```
(The coordinates are given in the region's local reference frame with origin (0, 0, 0) placed in the center of the region.)

The first region (`Volume #0`) is centered in the laboratory frame at (0, 0, 0), while the second region (`Volume #1`) is centered at (0, 0, 2000.01).

Both regions have identical magnetic field maps given in the file `B_Field_boundary_test.dat` which is created with the script `create_magnetic_field.C`.
(The coordinates in the field map are also given in the region's local reference frame with origin (0, 0, 0) in the center of the region.)

The field map contains values of `Bx`, `By` and `Bz` on a 3D set of grid points. The mesh size of the grid is 10 mm in `x`, `y` and `z` direction.
The `Bz` component is equal to 0 in all grid points, while the `Bx` and `By` components are equal to 0 in all grid points except in the grid points in the subvolume defined by:

```
-40 <= x <= +40 mm
-40 <= y <= +40 mm
-500 <= z <= +600 mm
```
where `Bx` = 4 and `By` = 3.

However, since the magnetic field is given on a discrete set of grid points separated by 10 mm in x, y and z directions, the `Bx` and `By` components are not equal to 0 in the following ranges of points:

```
-50 < x < +50 mm
-50 < y < +50 mm
-510 < z < +610 mm
```
as can be seen in the `Bx_vs_x_y.jpg` and `By_vs_x_y_SURF3.jpg` plots.

The reason is that the `TRestAxionMagneticField::GetMagneticField` method (which returns `Bx`, `By` and `Bz` values at the point specified as the input parameter) uses trilinear interpolation to determine the magnetic field values at the given point. Therefore, for the points with `x` coordinate, e.g.,
-50 <= x <= -40, the value of `Bx` will linearly increase from `Bx` = 0 (which is the value at `x` = -50 grid point) to `Bx` = 4 at `x` = -40 grid point.

### Description of files in this directory

1) **create_magnetic_field.C** - a restRoot script that produces the magnetic field map described above and stores it in the file `B_Field_boundary_test.dat`.

2) **B_Field_boundary_test.dat**: this file contains the magnetic field map produced by the `create_magnetic_field.C` script.
Each row in the file contains data for one grid point. First three columns contain `x`, `y` and `z` coordinate of the grid point, respectively,
while other three columns are `Bx`, `By` and `Bz` values for the corresponding grid point.

3) **Boundaries_test_interactive.C** - a restRoot script which tests the `TRestAxionMagneticField::GetVolumeBoundaries` and `TRestAxionMagneticField::GetFieldBoundaries` methods. The test is interactive in a way that the user should, after starting the script, enter the coordinates of the particle's initial position and the components of the particle's direction vector. Then the script calls `GetVolumeBoundaries` and `GetFieldBoundaries` methods to determine the boudary points for each of the two magnetic field regions (`Volume #0` and `Volume #1`) and prints them on screen as well as in the file `Boundaries_test_interactive_output.txt`. The example of one such test is shown in `example.jpeg`.

4) **Boundaries_test.C** - a non-interactive version of the script that uses four predefined inputs to perform the test and is suitable to launch in the gitlab pipeline chain.

5) **Bx_vs_x_y.jpeg**: this file shows the profile of the `Bx` component in the x-y plane in the middle of the region, i.e., at `z=0`.

6) **By_vs_x_y_SURF3.jpeg**: this file shows the profile of the `Bx` component in the x-y plane in the middle of the region, i.e., at `z=0`.

7) **example.jpeg**: this file shows the example of the screen-output produced by the script `Boundaries_test_interactive.C`.

8) **Boundaries_test_interactive_output.txt** : this text file contains the output produced by the script `Boundaries_test_interactive.C` that can be useful for detailed studies.
