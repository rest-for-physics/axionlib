void create_magnetic_field(){

ofstream Bfield_file;
Double_t x, y, z;
Double_t Bx, By, Bz;

Bfield_file.open("B_Field_boundary_test.dat");
Bfield_file << std::setprecision(1) << std::fixed;

for (int i = 0; i < 15; i++) {
     x = -70.0 + i * 10; 
     for (int j = 0; j < 15; j++) {
          y = -70.0 + j * 10;
          for (int k = 0; k < 201; k++) {
               z = -1000.0 + k * 10;
               Bx = 0.0; By = 0.0;
               if ((x >= -40.0) && (x <= 40.0) && (y >= -40.0) && (y <= 40.0) && (z >= -500.0) && (z <= 600.0)) {
                   Bx = 4.0;
                   By = 3.0;
                   Bz = 0.0;
               }
               else {
                   Bx = 0.0;
                   By = 0.0;
                   Bz = 0.0;
               }
               Bfield_file << x << "\t" << y << "\t" << z << "\t" << Bx << "\t" << By << "\t" << Bz << "\n";
          }
    }
}
Bfield_file.close();
}
