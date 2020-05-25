#include <iomanip>
void create_magnetic_field() {
    ofstream Bfield_file;
    Double_t x, y, z;
    Double_t Bx, By, Bz;

    Bfield_file.open("Magnetic_field.dat");
    Bfield_file << std::setprecision(1) << std::fixed;
    for (int i = 0; i < 15; i++) {
        x = -350.0 + i * 50.0;
        for (int j = 0; j < 15; j++) {
            y = -350.0 + j * 50.0;
            for (int k = 0; k < 201; k++) {
                z = -5000.0 + k * 50.0;
                Bx = 5.0 * x - 2.0 * y + 2.0 * z;
                By = 8.0 * x + 5.0 * y - 3.0 * z;
                Bz = -4.0 * x - 4.0 * y + z;
                Bfield_file << x << "\t" << y << "\t" << z << "\t" << Bx << "\t" << By << "\t" << Bz << "\n";
            }
        }
    }
    Bfield_file.close();
}
