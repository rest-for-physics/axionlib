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

TVector3 TRestAxionMagneticField::GetMagneticField(Double_t x, Double_t y, Double_t z) {
    Double_t bX = 0, bY = 0, bZ = 0;
    Int_t st;
#if defined USE_Garfield
    fSetOfField->MagneticField(x, y, z, bX, bY, bZ, st);

    debug << "Bx: " << bX << " T" << endl;
    debug << "By: " << bY << " T" << endl;
    debug << "Bz: " << bZ << " T" << endl;
#else
    cout << "This REST is not compiled with garfield, it cannot get field values using ComponentVoxel!"
         << endl;
#endif

    return TVector3(bX, bY, bZ);
}

void TRestAxionMagneticField::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    fNofVolumes = 0;
    fCanvas = NULL;
    fHisto = NULL;

#if defined USE_Garfield
    fSetOfField = new Sensor();
#endif
}

TCanvas* TRestAxionMagneticField::DrawHistogram(TString projection, TString Bcomp, Int_t VIndex,
                                                Double_t step) {
    if (fCanvas != NULL) {
        delete fCanvas;
        fCanvas = NULL;
    }

    if (fHisto != NULL) {
        delete fHisto;
        fHisto = NULL;
    }

    if (VIndex == -1) VIndex = 0;

    if (step == -1) step = fSizeMesh[0];

    if (VIndex >= fNofVolumes) ferr << VIndex << " corresponds to none volume index " << endl;

    Double_t xmax, xmin, ymax, ymin, zmax, zmin;
    xmax = fXmax[VIndex];
    xmin = fXmin[VIndex];
    ymax = fYmax[VIndex];
    ymin = fYmin[VIndex];
    zmax = fZmax[VIndex];
    zmin = fZmin[VIndex];
    Int_t nBinsX = (xmax - xmin) / step;
    Int_t nBinsY = (ymax - ymin) / step;
    Int_t nBinsZ = (zmax - zmin) / step;

    Double_t x, y, z;
    Double_t B;
    TVector3 Bvec;

#if defined USE_Garfield

    if (projection == "XY") {
        fCanvas = new TCanvas("fCanvas", "");
        fHisto = new TH2D("", "", nBinsX, xmin, xmax, nBinsY, ymin, ymax);

        z = (zmin + zmax) / 2.0;
        x = xmin;

        for (Int_t i = 0; i < nBinsX; i++) {
            y = ymin;
            for (Int_t j = 0; j < nBinsY; j++) {
                Bvec = GetMagneticField(x, y, z);
                if (Bcomp == "X")
                    B = Bvec[0];
                else {
                    if (Bcomp == "Y")
                        B = Bvec[1];
                    else {
                        if (Bcomp == "Z")
                            B = Bvec[2];
                        else
                            ferr << "You entered : " << Bcomp
                                 << " as a B component but you have to choose X, Y or Z" << endl;
                    }
                }
                fHisto->Fill(x, y, B);
                y = y + step;
            }
            x = x + step;
        }

        fCanvas->cd();
        fHisto->SetBit(TH1::kNoStats);
        fHisto->GetXaxis()->SetTitle("x (mm)");
        fHisto->GetYaxis()->SetTitle("y (mm)");

        if (Bcomp == "X") {
            fHisto->SetTitle("B_{x} against x and y");
            fCanvas->SetTitle("B_{x} against x and y");
        }
        if (Bcomp == "Y") {
            fHisto->SetTitle("B_{y} against x and y");
            fCanvas->SetTitle("B_{y} against x and y");
        }
        if (Bcomp == "Z") {
            fHisto->SetTitle("B_{z} against x and y");
            fCanvas->SetTitle("B_{z} against x and y");
        }

        fHisto->Draw("COLZ0");
        return fCanvas;
    }

    else {
        if (projection == "XZ") {
            TCanvas* fCanvas = new TCanvas("fCanvas", "");
            fHisto = new TH2D("", "", nBinsX, xmin, xmax, nBinsZ, zmin, zmax);

            y = (ymin + ymax) / 2.0;
            x = xmin;

            for (Int_t i = 0; i < nBinsX; i++) {
                z = zmin;
                for (Int_t j = 0; j < nBinsZ; j++) {
                    Bvec = GetMagneticField(x, y, z);
                    if (Bcomp == "X")
                        B = Bvec[0];
                    else {
                        if (Bcomp == "Y")
                            B = Bvec[1];
                        else {
                            if (Bcomp == "Z")
                                B = Bvec[2];
                            else
                                ferr << "You entered : " << Bcomp
                                     << " as a B component but you have to choose X, Y or Z" << endl;
                        }
                    }
                    fHisto->Fill(x, z, B);
                    z = z + step;
                }
                x = x + step;
            }

            fCanvas->cd();
            fHisto->SetBit(TH1::kNoStats);
            fHisto->GetXaxis()->SetTitle("x (mm)");
            fHisto->GetYaxis()->SetTitle("z (mm)");

            if (Bcomp == "X") {
                fHisto->SetTitle("B_{x} against x and z");
                fCanvas->SetTitle("B_{x} against x and z");
            }
            if (Bcomp == "Y") {
                fHisto->SetTitle("B_{y} against x and z");
                fCanvas->SetTitle("B_{y} against x and z");
            }
            if (Bcomp == "Z") {
                fHisto->SetTitle("B_{z} against x and z");
                fCanvas->SetTitle("B_{z} against x and z");
            }

            fHisto->Draw("COLZ0");
            return fCanvas;
        }

        else {
            if (projection == "YZ") {
                TCanvas* fCanvas = new TCanvas("fCanvas", "");
                fHisto = new TH2D("", "", nBinsY, ymin, ymax, nBinsZ, zmin, zmax);

                x = (xmin + xmax) / 2.0;
                y = ymin;

                for (Int_t i = 0; i < nBinsY; i++) {
                    z = zmin;
                    for (Int_t j = 0; j < nBinsZ; j++) {
                        Bvec = GetMagneticField(x, y, z);
                        if (Bcomp == "X")
                            B = Bvec[0];
                        else {
                            if (Bcomp == "Y")
                                B = Bvec[1];
                            else {
                                if (Bcomp == "Z")
                                    B = Bvec[2];
                                else
                                    ferr << "You entered : " << Bcomp
                                         << " as a B component but you have to choose X, Y or Z" << endl;
                            }
                        }
                        fHisto->Fill(y, z, B);
                        z = z + step;
                    }
                    y = y + step;
                }

                fCanvas->cd();
                fHisto->SetBit(TH1::kNoStats);
                fHisto->GetXaxis()->SetTitle("y (mm)");
                fHisto->GetYaxis()->SetTitle("z (mm)");

                if (Bcomp == "X") {
                    fHisto->SetTitle("B_{x} against y and z");
                    fCanvas->SetTitle("B_{x} against y and z");
                }
                if (Bcomp == "Y") {
                    fHisto->SetTitle("B_{y} against y and z");
                    fCanvas->SetTitle("B_{y} against y and z");
                }
                if (Bcomp == "Z") {
                    fHisto->SetTitle("B_{z} against y and z");
                    fCanvas->SetTitle("B_{z} against y and z");
                }

                fHisto->Draw("COLZ0");
                return fCanvas;
            }

            else
                ferr << "You entered : " << projection
                     << " as a projection but you have to choose XY, XY or XZ" << endl;
        }
    }

#else

    cout << "This REST is not compiled with garfield, it cannot get field values using Sensor !" << endl;
#endif
    return fCanvas;
}

void TRestAxionMagneticField::LoadMagneticVolumes() {
#ifdef USE_Garfield
    for (unsigned int n = 0; n < fPositions.size(); n++) {
        /*** Read information from the fileName ***/

        fstream file_r;
        file_r.open((string)getenv("REST_PATH") + "/data/magneticField/" +
                    (string)fFileNames[n]);  // original file
        if (file_r.is_open()) {
            ofstream file_w;
            file_w.open("/tmp/tmp_bField_" + (string)getenv("USER") + ".txt");
            TVector3 coordinates;
            double read;
            int i = 0;
            int j = 0;
            double xmax, xmin;
            double ymax, ymin;
            double zmax, zmin;
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

            ComponentVoxel* mesh = new ComponentVoxel();
            cout << "Setting Garfield mesh : Positions size : " << fPositions.size() << endl;
            mesh->SetMesh(nx, ny, nz, xmin + fPositions[n][0], xmax + fPositions[n][0],
                          ymin + fPositions[n][1], ymax + fPositions[n][1], zmin + fPositions[n][2],
                          zmax + fPositions[n][2]);

            /*** Fill the mesh with the magnetic field if file is open ***/
            mesh->LoadMagneticField("/tmp/tmp_bField_" + (string)getenv("USER") + ".txt", "XYZ", 1, 1);
            mesh->EnableInterpolation(true);

            fSetOfField->AddComponent(mesh);
            fNofVolumes++;

            fXmin.push_back(xmin + fPositions[n][0]);
            fXmax.push_back(xmax + fPositions[n][0]);

            fYmin.push_back(ymin + fPositions[n][1]);
            fYmax.push_back(ymax + fPositions[n][1]);

            fZmin.push_back(zmin + fPositions[n][2]);
            fZmax.push_back(zmax + fPositions[n][2]);
            fSizeMesh.push_back(sizeMesh);
        } else
            cout << " Cannot find the file " << endl;
    }
#else
    cout << "This REST is not compiled with garfield, it cannot load any magnetic field Volume!" << endl;
#endif
}

void TRestAxionMagneticField::InitFromConfigFile() {
    this->Initialize();

    string bVolume;
    size_t pos = 0;
    while ((bVolume = GetKEYDefinition("addMagneticVolume", pos)) != "") {
        TString filename = GetFieldValue("fileName", bVolume);
        fFileNames.push_back(filename);

        TVector3 position = Get3DVectorFieldValueWithUnits("position", bVolume);
        fPositions.push_back(position);
    }

    LoadMagneticVolumes();

    PrintMetadata();
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
        metadata << "* Volume " << p << " : "
                 << "  - Set in (" << x << "," << y << "," << z << ")"
                 << " mm" << endl;
        metadata << "  - Bounds : " << endl;
        metadata << "    xmin : " << fXmin[p] << " mm , xmax : " << fXmax[p] << " mm" << endl;
        metadata << "    ymin : " << fYmin[p] << " mm, ymax : " << fYmax[p] << " mm" << endl;
        metadata << "    zmin : " << fZmin[p] << " mm, zmax : " << fZmax[p] << " mm" << endl;
        metadata << "  - Size of the mesh : " << fSizeMesh[p] << " mm" << endl;
        metadata << "  - File loaded : " << fFileNames[p] << endl;
    }
    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}
