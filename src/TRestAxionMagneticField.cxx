/******************** REST disclaimer ***********************************
 * This file is part of the REST software framework.                     *
 *                                                                       *
 * Copyright (C) 2016 GIFNA/TREX (University of Zaragoza)                *
 * For more information see http://gifna.unizar.es/trex                  *
 *                                                                       *
 * REST is free software: you can redistribute it and/or modify          *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation, either version 3 of the License, or     *
 * (at your option) any later version.                                   *
 *                                                                       *
 * REST is distributed in the hope that it will be useful,               *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have a copy of the GNU General Public License along with   *
 * REST in $REST_PATH/LICENSE.                                           *
 * If not, see http://www.gnu.org/licenses/.                             *
 * For the list of contributors see $REST_PATH/CREDITS.                  *
 *************************************************************************/

/***************** DOXYGEN DOCUMENTATION ********************************
/// TRestAxionMagneticField is a class that allows to load externally
/// defined magnetic fields, and create magnetic volume regions using those
/// pre-generated definitions. Once the field maps have been loaded this
/// class will evaluate if a coordinate (x,y,z) is in a given magnetic
/// region, and it will return a non-zero value of the magnetic field in
/// case (x,y,z) is inside a magnetic region.
///
/// ----
/// THIS SHOULD BE GENERALIZED SOON. For now we consider that the bounds of
/// the volume are X,Y: -0.24 to 0.24 and Z:-5 to 5
/// ----
///
/// We can add any number of magnetic volumes inside the RML definition
/// as shown in the following piece of code.
///
/// <TRestAxionMagneticField>
///         <addMagneticVolume file="magnetic.file" position="(30,0,0)mm" />
///         <addMagneticVolume file="magnetic.file" position="(-30,0,0)mm" />
/// <TRestAxionMagneticField/>
///
/// where we produce 2 magnetic regions, using the same magnetic map provided
/// in file `magnetic.file` and shifted by x=-30mm and x=30mm. The parameters
/// available in the `addMagneticVolume` definition are described in this list.
///
/// - file : This allows to specify the filename that contains the values of the
/// magnetic field. Few files will be found under `data/magneticField`. They all
/// contain 3 columns for the position in the volume and 3 columns to define the
/// magnetic field vector. It is not the full path but only the name of the file.
///
/// - position : By convention, the volume is build using the coordinates provided
/// in the magnetic field file given. However, it is possible to translate the
/// volume using the `position` field.
///
/// \TODO maybe add a section to visualize the magnetic volume (Can be done with
/// ViewField()).
///
/// \warning It seems to be difficult to apply rotations to the field coordinates,
/// and assign it to Garfield methods. It seems also not possible to rotate the field
/// in Garfield routines. Still, we may keep always our magnet in horizontal position
/// to perform our studies.
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2019-June: First concept and implementation of TRestAxionMagneticField class.
///            Eve Pachoud
///
/// \class      TRestAxionMagneticField
/// \author     Eve Pachoud
///
/// <hr>
///
 *************************************************************************/

#include "TRestAxionMagneticField.h"

using namespace std;
#ifdef USE_Garfield
using namespace Garfield;
#endif

ClassImp(TRestAxionMagneticField);

TRestAxionMagneticField::TRestAxionMagneticField() : TRestMetadata() { Initialize(); }

TRestAxionMagneticField::TRestAxionMagneticField(const char* cfgFileName, string name)
    : TRestMetadata(cfgFileName) {
    cout << "Entering TRestAxionMagneticField constructor( cfgFileName, name )" << endl;

    Initialize();

    LoadConfigFromFile(fConfigFileName, name);

    PrintMetadata();
}

TRestAxionMagneticField::~TRestAxionMagneticField() {
    // TRestAxionMagneticField destructor

    debug << "Entering ... TRestAxionMagneticField() destructor." << endl;

#if defined USE_Garfield
    delete fSetOfField;
#endif
}

void TRestAxionMagneticField::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    fNofVolumes = 0;

#if defined USE_Garfield
    fSetOfField = new Garfield::Sensor();
#endif
}

void TRestAxionMagneticField::LoadMagneticVolumes() {
#ifdef USE_Garfield

    
    for( unsigned int n = 0; n < fPositions.size(); n++ ) {

    /*** Read information from the fileName ***/

    fstream file_r;
    file_r.open((string)getenv("REST_PATH") + "/data/magneticField/" + (string)fFileNames[n]);  // original file
    if (file_r.is_open()) {
        ofstream file_w;
        file_w.open("/tmp/tmp_bField_"+(string)getenv("USER")+".txt"); 
        TVector3 coordinates;                   
        double read;
        int i = 0;
        int j = 0;
	double xmax,xmin;
	double ymax,ymin;
	double zmax,zmin;
	double sizeMesh;
        while (file_r >> read) {
            if (i < 4) {
                if (i == 0) {
                    xmax = read;
                    xmin = -xmax;
                    i = i + 1;
                } else if (i == 1) {
                    ymax = read;
                    ymin = -ymax;
                    i = i + 1;
                } else if (i == 2) {
                    zmax = read;
                    zmin = -zmax;
                    i = i + 1;
                } else if (i == 3) {
                    sizeMesh = read;
                    i = i + 1;
                }

            }

            else {
                if (j < 3) {
                    coordinates[j] = read;
                    if (j == 2) {
                        for (int k = 0; k < 3; k++) {
                            coordinates[k] = coordinates[k] + fPositions[n][k];
                            file_w << coordinates[k];
                            file_w << "\t";
                        }
                    }
                    j = j + 1;
                }

                else {
                    if (j == 5) {
                        file_w << read;
                        file_w << "\n";
                        j = 0;
                    } else {
                        file_w << read;
                        file_w << "\t";
                        j = j + 1;
                    }
                }
            }
        }
        file_w.close();
        file_r.close();

        /*** Create the mesh if file is open ***/
        int nx = (int)(2 * xmax / sizeMesh) + 1;
        int ny = (int)(2 * ymax / sizeMesh) + 1;
        int nz = (int)(2 * zmax / sizeMesh) + 1;
      
	Garfield::ComponentVoxel * mesh = new Garfield::ComponentVoxel();
        mesh->SetMesh(nx, ny, nz, xmin + fPositions[n][0], xmax + fPositions[n][0], ymin + fPositions[n][1], ymax + fPositions[n][1],
                       zmin + fPositions[n][2], zmax + fPositions[n][2]);

        /*** Fill the mesh with the magnetic field if file is open ***/
        mesh->LoadMagneticField("/tmp/tmp_bField_"+(string)getenv("USER")+".txt", "XYZ", 1,1); 
        mesh->EnableInterpolation(true);
 
        fSetOfField->AddComponent(mesh);
	fNofVolumes++;
	fXmax.push_back(xmax);	
	fXmin.push_back(xmin);	
	fYmax.push_back(ymax);	
	fYmin.push_back(ymin);	
	fZmax.push_back(zmax);	
	fXmin.push_back(zmin);	
	fSizeMesh.push_back(sizeMesh);
    } else
        cout << " Cannot find the file " << endl;
}
#else
    cout << "This REST is not complied with garfield, it cannot load any magnetic field Volume!" << endl;
#endif
}


void TRestAxionMagneticField::InitFromConfigFile() {
    this->Initialize();
    string bVolume;
    size_t pos = 0;
    while ((bVolume = GetKEYDefinition("addMagneticVolume", pos)) != "")  
    {
        TString filename = GetParameter("file");
	fFileNames.push_back(filename);
        TVector3 position = Get3DVectorParameterWithUnits("position");
        fPositions.push_back(position);
    }
    LoadMagneticVolumes();
}

void TRestAxionMagneticField::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    metadata << " - Number of magnetic volumes : " << fNofVolumes << endl;
    metadata << " ------------------------------------------------ " << endl;
    double x, y, z;
    for (int p = 0; p < fNofVolumes; p++) {
        x = fPositions[p][0];
        y = fPositions[p][1];
        z = fPositions[p][2];
        metadata << "* Volume " << p + 1 << " : "
                 << "  - Set in (" << x << "," << y << "," << z << ")" <<  endl;
	metadata << "  - Bounds : " << endl;
	metadata << "    xmin : " << fXmin[p] << " , xmax : " << fXmax[p]<<endl;
	metadata << "    ymin : " << fYmin[p] << " , ymax : " << fYmax[p]<<endl;      
	metadata << "    zmin : " << fZmin[p] << " , zmax : " << fZmax[p]<<endl;
	metadata << "  - Size of the mesh : "<< fSizeMesh[p] << endl;            
	metadata << "  - File loaded : " << fFileNames[p] << endl;
    }  
    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}
