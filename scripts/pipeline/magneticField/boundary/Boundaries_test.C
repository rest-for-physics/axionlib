#include <fstream>
#include <iomanip>
#include <iostream>
using namespace std;

TRestAxionMagneticField* myField;

Int_t CheckValues(TVector3 real_value, TVector3 true_value);

void Boundaries_test() {
    TVector3 position;
    TVector3 direction;
    std::vector<TVector3> boundaries;
    std::vector<TVector3> field_boundaries;
    std::vector<TVector3> true_boundaries;
    std::vector<TVector3> true_field_boundaries;
    myField = new TRestAxionMagneticField("../fields.rml", "bField_2_regions");

    // Example 1
    position.SetXYZ(10.0, 10.0, -1010.0);
    direction.SetXYZ(0.0, 0.0, 1.0);
    // testing boundaries in Volume #0
    boundaries.clear();
    boundaries = myField->GetVolumeBoundaries(0, position, direction);
    true_boundaries.resize(2);
    true_boundaries[0].SetXYZ(10.0, 10.0, -1000.0);
    true_boundaries[1].SetXYZ(10.0, 10.0, 1000.0);
    if ((CheckValues(boundaries[0], true_boundaries[0])) ||
        (CheckValues(boundaries[1], true_boundaries[1]))) {
        cout << "Error occured in GetVolumeBoundaries method when evaluating in Volume #0 for initial "
                "position  ("
             << position.X() << ", " << position.Y() << ", " << position.Z() << ") and direction ("
             << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << endl;
        return 5;
    }
    field_boundaries.clear();
    field_boundaries = myField->GetFieldBoundaries(0, position, direction);
    true_field_boundaries.resize(2);
    true_field_boundaries[0].SetXYZ(10.0, 10.0, -505.0);
    true_field_boundaries[1].SetXYZ(10.0, 10.0, 605.0);
    if ((CheckValues(field_boundaries[0], true_field_boundaries[0])) ||
        (CheckValues(field_boundaries[1], true_field_boundaries[1]))) {
        cout << "Error occured in GetFieldBoundaries method when evaluating in Volume #0 for initial "
                "position  ("
             << position.X() << ", " << position.Y() << ", " << position.Z() << ") and direction ("
             << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << endl;
        return 5;
    }
    // testing boundaries in Volume #1
    boundaries.clear();
    boundaries = myField->GetVolumeBoundaries(1, position, direction);
    true_boundaries.clear();
    true_boundaries[0].SetXYZ(10.0, 10.0, 1000.01);
    true_boundaries[1].SetXYZ(10.0, 10.0, 3000.01);
    if ((CheckValues(boundaries[0], true_boundaries[0])) ||
        (CheckValues(boundaries[1], true_boundaries[1]))) {
        cout << "Error occured in GetVolumeBoundaries method when evaluating in Volume #1 for initial "
                "position  ("
             << position.X() << ", " << position.Y() << ", " << position.Z() << ") and direction ("
             << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << endl;
        return 5;
    }
    field_boundaries.clear();
    field_boundaries = myField->GetFieldBoundaries(1, position, direction);
    true_field_boundaries.clear();
    true_field_boundaries[0].SetXYZ(10.0, 10.0, 1495.01);
    true_field_boundaries[1].SetXYZ(10.0, 10.0, 2605.01);
    if ((CheckValues(field_boundaries[0], true_field_boundaries[0])) ||
        (CheckValues(field_boundaries[1], true_field_boundaries[1]))) {
        cout << "Error occured in GetFieldBoundaries method when evaluating in Volume #1 for initial "
                "position  ("
             << position.X() << ", " << position.Y() << ", " << position.Z() << ") and direction ("
             << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << endl;
        return 5;
    }

    // Example 2
    position.SetXYZ(20.0, 20.0, 3010.0);
    direction.SetXYZ(0.0, 0.0, -1.0);
    // testing boundaries in Volume #0
    boundaries.clear();
    boundaries = myField->GetVolumeBoundaries(0, position, direction);
    true_boundaries.resize(2);
    true_boundaries[0].SetXYZ(20.0, 20.0, 1000.0);
    true_boundaries[1].SetXYZ(20.0, 20.0, -1000.0);
    if ((CheckValues(boundaries[0], true_boundaries[0])) ||
        (CheckValues(boundaries[1], true_boundaries[1]))) {
        cout << "Error occured in GetVolumeBoundaries method when evaluating in Volume #0 for initial "
                "position  ("
             << position.X() << ", " << position.Y() << ", " << position.Z() << ") and direction ("
             << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << endl;
        return 5;
    }
    field_boundaries.clear();
    field_boundaries = myField->GetFieldBoundaries(0, position, direction);
    true_field_boundaries.resize(2);
    true_field_boundaries[0].SetXYZ(20.0, 20.0, 605.0);
    true_field_boundaries[1].SetXYZ(20.0, 20.0, -505.0);
    if ((CheckValues(field_boundaries[0], true_field_boundaries[0])) ||
        (CheckValues(field_boundaries[1], true_field_boundaries[1]))) {
        cout << "Error occured in GetFieldBoundaries method when evaluating in Volume #0 for initial "
                "position  ("
             << position.X() << ", " << position.Y() << ", " << position.Z() << ") and direction ("
             << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << endl;
        return 5;
    }
    // testing boundaries in Volume #1
    boundaries.clear();
    boundaries = myField->GetVolumeBoundaries(1, position, direction);
    true_boundaries.clear();
    true_boundaries[0].SetXYZ(20.0, 20.0, 3000.01);
    true_boundaries[1].SetXYZ(20.0, 20.0, 1000.01);
    if ((CheckValues(boundaries[0], true_boundaries[0])) ||
        (CheckValues(boundaries[1], true_boundaries[1]))) {
        cout << "Error occured in GetVolumeBoundaries method when evaluating in Volume #1 for initial "
                "position  ("
             << position.X() << ", " << position.Y() << ", " << position.Z() << ") and direction ("
             << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << endl;
        return 5;
    }
    field_boundaries.clear();
    field_boundaries = myField->GetFieldBoundaries(1, position, direction);
    true_field_boundaries.clear();
    true_field_boundaries[0].SetXYZ(20.0, 20.0, 2605.01);
    true_field_boundaries[1].SetXYZ(20.0, 20.0, 1495.01);
    if ((CheckValues(field_boundaries[0], true_field_boundaries[0])) ||
        (CheckValues(field_boundaries[1], true_field_boundaries[1]))) {
        cout << "Error occured in GetFieldBoundaries method when evaluating in Volume #1 for initial "
                "position  ("
             << position.X() << ", " << position.Y() << ", " << position.Z() << ") and direction ("
             << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << endl;
        return 5;
    }

    // Example 3
    position.SetXYZ(-60.0, 0.0, -1010.0);
    direction.SetXYZ(1.0, 0.0, 35.0);
    // testing boundaries in Volume #0
    boundaries.clear();
    boundaries = myField->GetVolumeBoundaries(0, position, direction);
    true_boundaries.resize(2);
    true_boundaries[0].SetXYZ(-59.7143, 0.0, -1000.0);
    true_boundaries[1].SetXYZ(-2.57143, 0.0, 1000.0);
    if ((CheckValues(boundaries[0], true_boundaries[0])) ||
        (CheckValues(boundaries[1], true_boundaries[1]))) {
        cout << "Error occured in GetVolumeBoundaries method when evaluating in Volume #0 for initial "
                "position  ("
             << position.X() << ", " << position.Y() << ", " << position.Z() << ") and direction ("
             << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << endl;
        return 5;
    }
    field_boundaries.clear();
    field_boundaries = myField->GetFieldBoundaries(0, position, direction);
    true_field_boundaries.resize(2);
    true_field_boundaries[0].SetXYZ(-45.5772, 0.0, -505.202);
    true_field_boundaries[1].SetXYZ(-13.8525, 0.0, 605.161);
    if ((CheckValues(field_boundaries[0], true_field_boundaries[0])) ||
        (CheckValues(field_boundaries[1], true_field_boundaries[1]))) {
        cout << "Error occured in GetFieldBoundaries method when evaluating in Volume #0 for initial "
                "position  ("
             << position.X() << ", " << position.Y() << ", " << position.Z() << ") and direction ("
             << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << endl;
        return 5;
    }
    // testing boundaries in Volume #1
    boundaries.clear();
    boundaries = myField->GetVolumeBoundaries(1, position, direction);
    true_boundaries.clear();
    true_boundaries[0].SetXYZ(-2.57114, 0.0, 1000.01);
    true_boundaries[1].SetXYZ(54.5717, 0.0, 3000.01);
    if ((CheckValues(boundaries[0], true_boundaries[0])) ||
        (CheckValues(boundaries[1], true_boundaries[1]))) {
        cout << "Error occured in GetVolumeBoundaries method when evaluating in Volume #1 for initial "
                "position  ("
             << position.X() << ", " << position.Y() << ", " << position.Z() << ") and direction ("
             << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << endl;
        return 5;
    }
    field_boundaries.clear();
    field_boundaries = myField->GetFieldBoundaries(1, position, direction);
    true_field_boundaries.clear();
    true_field_boundaries[0].SetXYZ(11.5659, 0.0, 1494.81);
    true_field_boundaries[1].SetXYZ(43.2906, 0.0, 2605.17);
    if ((CheckValues(field_boundaries[0], true_field_boundaries[0])) ||
        (CheckValues(field_boundaries[1], true_field_boundaries[1]))) {
        cout << "Error occured in GetFieldBoundaries method when evaluating in Volume #1 for initial "
                "position  ("
             << position.X() << ", " << position.Y() << ", " << position.Z() << ") and direction ("
             << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << endl;
        return 5;
    }

    // Example 4
    position.SetXYZ(-70.0, 0.0, -1020.0);
    direction.SetXYZ(0.7, 0.0, 40.0);
    // testing boundaries in Volume #0
    boundaries.clear();
    boundaries = myField->GetVolumeBoundaries(0, position, direction);
    true_boundaries.resize(2);
    true_boundaries[0].SetXYZ(-69.65, 0.0, -1000.0);
    true_boundaries[1].SetXYZ(-34.65, 0.0, 1000.0);
    if ((CheckValues(boundaries[0], true_boundaries[0])) ||
        (CheckValues(boundaries[1], true_boundaries[1]))) {
        cout << "Error occured in GetVolumeBoundaries method when evaluating in Volume #0 for initial "
                "position  ("
             << position.X() << ", " << position.Y() << ", " << position.Z() << ") and direction ("
             << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << endl;
        return 5;
    }
    field_boundaries.clear();
    field_boundaries = myField->GetFieldBoundaries(0, position, direction);
    true_field_boundaries.resize(2);
    true_field_boundaries[0].SetXYZ(-49.9655, 0.0, 124.828);
    true_field_boundaries[1].SetXYZ(-41.5614, 0.0, 605.06);
    if ((CheckValues(field_boundaries[0], true_field_boundaries[0])) ||
        (CheckValues(field_boundaries[1], true_field_boundaries[1]))) {
        cout << "Error occured in GetFieldBoundaries method when evaluating in Volume #0 for initial "
                "position  ("
             << position.X() << ", " << position.Y() << ", " << position.Z() << ") and direction ("
             << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << endl;
        return 5;
    }
    // testing boundaries in Volume #1
    boundaries.clear();
    boundaries = myField->GetVolumeBoundaries(1, position, direction);
    true_boundaries.clear();
    true_boundaries[0].SetXYZ(-34.6498, 0.0, 1000.01);
    true_boundaries[1].SetXYZ(0.350175, 0.0, 3000.01);
    if ((CheckValues(boundaries[0], true_boundaries[0])) ||
        (CheckValues(boundaries[1], true_boundaries[1]))) {
        cout << "Error occured in GetVolumeBoundaries method when evaluating in Volume #1 for initial "
                "position  ("
             << position.X() << ", " << position.Y() << ", " << position.Z() << ") and direction ("
             << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << endl;
        return 5;
    }
    field_boundaries.clear();
    field_boundaries = myField->GetFieldBoundaries(1, position, direction);
    true_field_boundaries.clear();
    true_field_boundaries[0].SetXYZ(-25.9887, 0.0, 1494.93);
    true_field_boundaries[1].SetXYZ(-6.56127, 0.0, 2605.07);
    if ((CheckValues(field_boundaries[0], true_field_boundaries[0])) ||
        (CheckValues(field_boundaries[1], true_field_boundaries[1]))) {
        cout << "Error occured in GetFieldBoundaries method when evaluating in Volume #1 for initial "
                "position  ("
             << position.X() << ", " << position.Y() << ", " << position.Z() << ") and direction ("
             << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << endl;
        return 5;
    }

    cout << "GetVolumeBoundaries and GetFieldBoundaries test successful." << endl;
    return 0;
}

Int_t CheckValues(TVector3 real_value, TVector3 true_value) {
    if ((fabs(true_value.X() - real_value.X()) <=
         0.0001 * std::max(fabs(true_value.X()), fabs(real_value.X()))) &&
        (fabs(true_value.Y() - real_value.Y()) <=
         0.0001 * std::max(fabs(true_value.Y()), fabs(real_value.Y()))) &&
        (fabs(true_value.Z() - real_value.Z()) <=
         0.0001 * std::max(fabs(true_value.Z()), fabs(real_value.Z()))))
        return 0;
    else
        return 5;
}
