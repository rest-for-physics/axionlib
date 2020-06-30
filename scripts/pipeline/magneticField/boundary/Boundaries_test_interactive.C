#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

TRestAxionMagneticField *myField;
ofstream myfile;

void PrintVolumeBoundaries (Int_t, TVector3, TVector3, std::vector<TVector3>);
void PrintFieldBoundaries (std::vector<TVector3>);

void Boundaries_test_interactive() {
    
    Int_t volume_ID;
    TVector3 position;
    TVector3 direction;
    std::vector<TVector3> boundaries;
    std::vector<TVector3> field_boundaries;
    myfile.open ("Boundaries_test_interactive_output.txt");
    myField = new TRestAxionMagneticField("../fields.rml","bField_2_regions");
    char choice = 'y';
    Double_t x, y, z;
    Double_t dir_x, dir_y, dir_z;
    while (choice == 'y') {
        cout << "Enter x, y and z coordinate of the particle initial position : ";
        cin >> x >> y >> z;
        cout << "Enter x, y and z coordinate of the particle direction : ";
        cin >> dir_x >> dir_y >> dir_z;
        position.SetXYZ(x, y, z);
        direction.SetXYZ(dir_x, dir_y, dir_z);
        for (volume_ID = 0; volume_ID < myField->GetNumberOfVolumes(); volume_ID++) {
            boundaries.clear();
            boundaries = myField->GetVolumeBoundaries(volume_ID, position, direction);
            PrintVolumeBoundaries(volume_ID, position, direction, boundaries);
            field_boundaries.clear();
            field_boundaries = myField->GetFieldBoundaries(volume_ID, position, direction);
            PrintFieldBoundaries(field_boundaries);
        }
        cout << "continue? (y/n) : ";
        cin >> choice;
    }
    cout << "Test finished" << endl;
    myfile.close();
}

void PrintVolumeBoundaries (Int_t volume_ID, TVector3 position, TVector3 direction, std::vector<TVector3> boundaries) {
    if (volume_ID == 0) {
        cout << "Initial position : (" << position.X() << ", " << position.Y() << ", " << position.Z() << ")" << endl;
        cout << "Direction : (" << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << endl << endl;        
    }
    cout << "Volume #" << volume_ID << endl;
    cout << "Number of volume boundary points found : " << boundaries.size() << endl;
    for (int i = 0; i < boundaries.size(); i++) 
        cout << "boundaries[" << i <<"] = (" << boundaries[i].X() << ", " << boundaries[i].Y() << ", " << boundaries[i].Z() << ")" << endl;
    cout << endl << endl;
    
    if (volume_ID == 0) {
        myfile << "Initial position : (" << position.X() << ", " << position.Y() << ", " << position.Z() << ")" << "\n";
        myfile << "Direction : (" << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << "\n\n";
    }
    myfile << "Volume #" << volume_ID << "\n";
    myfile << "Number of volume boundary points found : " << boundaries.size() << "\n";
    for (int i = 0; i < boundaries.size(); i++) 
        myfile << "boundaries[" << i <<"] = (" << boundaries[i].X() << ", " << boundaries[i].Y() << ", " << boundaries[i].Z() << ")" << "\n";
    myfile << "\n\n";
}

void PrintFieldBoundaries (std::vector<TVector3> fieldboundaries) {
    cout << "Number of field boundary points found : " << fieldboundaries.size() << endl;
    for (int i = 0; i < fieldboundaries.size(); i++) 
        cout << "fieldboundaries[" << i <<"] = (" << fieldboundaries[i].X() << ", " << fieldboundaries[i].Y() << ", " << fieldboundaries[i].Z() << ")" << endl;
    cout << "---------------------------------" << endl;
    
    myfile << "Number of field boundary points found : " << fieldboundaries.size() << "\n";
    for (int i = 0; i < fieldboundaries.size(); i++) 
        myfile << "fieldboundaries[" << i <<"] = (" << fieldboundaries[i].X() << ", " << fieldboundaries[i].Y() << ", " << fieldboundaries[i].Z() << ")" << "\n";
    myfile << "---------------------------------" << "\n";
}
