#include <fstream>
#include <iomanip>
#include <iostream>
using namespace std;

TRestAxionMagneticField* myField;
ofstream myfile;
TVector3 offset = TVector3(800.0, 800.0, 8000.0);

void PrintPoint(TVector3 point);

void GetMagneticField_test() {
    TVector3 coordinates;
    myfile.open("GetMagneticField_test_output.txt");
    myfile << std::setprecision(1) << std::fixed;
    cout << std::setprecision(1) << std::fixed;
    myField = new TRestAxionMagneticField("../fields.rml", "bField");

    // changing x, while y and z are constant
    coordinates = TVector3(-349.8, -325.0, -4975.0);
    PrintPoint(coordinates);
    coordinates = TVector3(-325.0, -325.0, -4975.0);
    PrintPoint(coordinates);
    coordinates = TVector3(-300.2, -325.0, -4975.0);
    PrintPoint(coordinates);
    coordinates = TVector3(-300.0, -325.0, -4975.0);
    PrintPoint(coordinates);
    coordinates = TVector3(-299.8, -325.0, -4975.0);
    PrintPoint(coordinates);
    coordinates = TVector3(-275.0, -325.0, -4975.0);
    PrintPoint(coordinates);
    coordinates = TVector3(-250.2, -325.0, -4975.0);
    PrintPoint(coordinates);
    coordinates = TVector3(-250.0, -325.0, -4975.0);
    PrintPoint(coordinates);
    coordinates = TVector3(-249.8, -325.0, -4975.0);
    PrintPoint(coordinates);
    coordinates = TVector3(-200.2, -325.0, -4975.0);
    PrintPoint(coordinates);

    // changing y, while x and z are constant
    coordinates = TVector3(-200.2, -300.2, -4975.0);
    PrintPoint(coordinates);
    coordinates = TVector3(-200.2, -300.0, -4975.0);
    PrintPoint(coordinates);
    coordinates = TVector3(-200.2, -299.8, -4975.0);
    PrintPoint(coordinates);
    coordinates = TVector3(-200.2, -275.0, -4975.0);
    PrintPoint(coordinates);
    coordinates = TVector3(-200.2, -250.2, -4975.0);
    PrintPoint(coordinates);
    coordinates = TVector3(-200.2, -250.0, -4975.0);
    PrintPoint(coordinates);
    coordinates = TVector3(-200.2, -249.8, -4975.0);
    PrintPoint(coordinates);
    coordinates = TVector3(-200.2, -200.1, -4975.0);
    PrintPoint(coordinates);

    // changing z, while x and y are constant
    coordinates = TVector3(-200.2, -200.1, -4950.2);
    PrintPoint(coordinates);
    coordinates = TVector3(-200.2, -200.1, -4950.0);
    PrintPoint(coordinates);
    coordinates = TVector3(-200.2, -200.1, -4949.8);
    PrintPoint(coordinates);
    coordinates = TVector3(-200.2, -200.1, -4900.2);
    PrintPoint(coordinates);

    // various x, y, z combinations
    coordinates = TVector3(50.0, 50.0, 50.0);
    PrintPoint(coordinates);
    coordinates = TVector3(350.0, 200.1, 2000.0);
    PrintPoint(coordinates);
    coordinates = TVector3(225.0, 350.0, 4950.2);
    PrintPoint(coordinates);
    coordinates = TVector3(140.6, 71.0, 5000.0);
    PrintPoint(coordinates);
    coordinates = TVector3(350.0, 350.0, 5000.0);
    PrintPoint(coordinates);
    // points out of boundaries
    coordinates = TVector3(350.1, 200.0, 1001.0);
    PrintPoint(coordinates);
    coordinates = TVector3(-200.2, -351.0, -1950.2);
    PrintPoint(coordinates);
    coordinates = TVector3(155.0, -275.1, 5001.0);
    PrintPoint(coordinates);
    coordinates = TVector3(-350.1, 350.1, 5001.0);
    PrintPoint(coordinates);

    myfile.close();
}

void PrintPoint(TVector3 point) {
    TVector3 B = myField->GetMagneticField(point + offset);
    TVector3 true_value = TVector3(0.0, 0.0, 0.0);

    true_value[0] = 5.0 * point.X() - 2.0 * point.Y() + 2.0 * point.Z();
    true_value[1] = 8.0 * point.X() + 5.0 * point.Y() - 3.0 * point.Z();
    true_value[2] = -4.0 * point.X() - 4.0 * point.Y() + point.Z();
    cout << "x = " << point.X() + offset.X() << "  y = " << point.Y() + offset.Y()
         << "  z = " << point.Z() + offset.Z() << endl;
    cout << "Bx = " << B.X() << "  By = " << B.Y() << "  Bz = " << B.Z() << "   True values: "
         << "TBx = " << true_value.X() << "  TBy = " << true_value.Y() << "  TBz = " << true_value.Z()
         << endl;
    cout << "Difference = (" << true_value.X() - B.X() << ", " << true_value.Y() - B.Y() << ", "
         << true_value.Z() - B.Z() << ")" << endl
         << endl;

    myfile << "x = " << point.X() + offset.X() << "  y = " << point.Y() + offset.Y()
           << "  z = " << point.Z() + offset.Z() << "\n";
    myfile << "Bx = " << B.X() << "  By = " << B.Y() << "  Bz = " << B.Z() << "   True values: "
           << "TBx = " << true_value.X() << "  TBy = " << true_value.Y() << "  TBz = " << true_value.Z()
           << "\n";
    myfile << "Difference = (" << true_value.X() - B.X() << ", " << true_value.Y() - B.Y() << ", "
           << true_value.Z() - B.Z() << ")"
           << "\n\n";
}
