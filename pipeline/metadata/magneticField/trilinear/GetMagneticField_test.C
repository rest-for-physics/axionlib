#include <fstream>
#include <iomanip>
#include <iostream>
using namespace std;

TRestAxionMagneticField* myField;
TVector3 offset = TVector3(800.0, 800.0, 8000.0);
Int_t internal_points_check = 0;
Int_t external_points_check = 0;

Int_t CheckPoint(TVector3 point);

Int_t GetMagneticField_test() {
    TVector3 coordinates;
    myField = new TRestAxionMagneticField("fields.rml", "bField");

    // changing x, while y and z are constant
    coordinates = TVector3(-349.8, -325.0, -4975.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(-325.0, -325.0, -4975.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(-300.2, -325.0, -4975.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(-300.0, -325.0, -4975.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(-299.8, -325.0, -4975.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(-275.0, -325.0, -4975.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(-250.2, -325.0, -4975.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(-250.0, -325.0, -4975.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(-249.8, -325.0, -4975.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(-200.2, -325.0, -4975.0);
    if (CheckPoint(coordinates)) return 5;

    // changing y, while x and z are constant
    coordinates = TVector3(-200.2, -300.2, -4975.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(-200.2, -300.0, -4975.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(-200.2, -299.8, -4975.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(-200.2, -275.0, -4975.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(-200.2, -250.2, -4975.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(-200.2, -250.0, -4975.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(-200.2, -249.8, -4975.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(-200.2, -200.1, -4975.0);
    if (CheckPoint(coordinates)) return 5;

    // changing z, while x and y are constant
    coordinates = TVector3(-200.2, -200.1, -4950.2);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(-200.2, -200.1, -4950.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(-200.2, -200.1, -4949.8);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(-200.2, -200.1, -4900.2);
    if (CheckPoint(coordinates)) return 5;

    // various x, y, z combinations
    coordinates = TVector3(50.0, 50.0, 50.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(350.0, 200.1, 2000.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(225.0, 350.0, 4950.2);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(140.6, 71.0, 5000.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(350.0, 350.0, 5000.0);
    if (CheckPoint(coordinates)) return 5;

    cout << "Checking points evaluated outside magnetic regions." << endl;
    cout << " -	4 warning messages should appear from TRestAxionMagneticField" << endl;
    // points out of boundaries
    coordinates = TVector3(350.1, 200.0, 1001.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(-200.2, -351.0, -1950.2);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(155.0, -275.1, 5001.0);
    if (CheckPoint(coordinates)) return 5;
    coordinates = TVector3(-350.1, 350.1, 5001.0);
    if (CheckPoint(coordinates)) return 5;

    if ((internal_points_check == 27) && (external_points_check == 4)) {
        cout << "GetMagneticField test successful." << endl;
        return 0;
    } else {
        cout << "Errors in GetMagneticField test." << endl;
        return 5;
    }
}

Int_t CheckPoint(TVector3 point) {
    TVector3 B = myField->GetMagneticField(point + offset);
    TVector3 true_value = TVector3(0.0, 0.0, 0.0);

    true_value[0] = 5.0 * point.X() - 2.0 * point.Y() + 2.0 * point.Z();
    true_value[1] = 8.0 * point.X() + 5.0 * point.Y() - 3.0 * point.Z();
    true_value[2] = -4.0 * point.X() - 4.0 * point.Y() + point.Z();
    if ((fabs(true_value.X() - B.X()) <= 0.0001 * std::max(fabs(true_value.X()), fabs(B.X()))) &&
        (fabs(true_value.Y() - B.Y()) <= 0.0001 * std::max(fabs(true_value.Y()), fabs(B.Y()))) &&
        (fabs(true_value.Z() - B.Z()) <= 0.0001 * std::max(fabs(true_value.Z()), fabs(B.Z())))) {
        internal_points_check++;
        return 0;
    } else if (B.X() == 0.0 && B.Y() == 0.0 && B.Z() == 0.0) {
        external_points_check++;
        return 0;
    } else {
        cout << "Error occured when evaluating point (" << point.X() + offset.X() << ", "
             << point.Y() + offset.Y() << ", " << point.Z() + offset.Z() << ")" << endl;
        return 5;
    }
}
