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

//////////////////////////////////////////////////////////////////////////
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
/// \code
/// <TRestAxionMagneticField>
///         <addMagneticVolume file="magnetic.file" position="(30,0,0)mm" />
///         <addMagneticVolume file="magnetic.file" position="(-30,0,0)mm" />
/// <TRestAxionMagneticField/>
///
/// \endcode
///
/// where we produce 2 magnetic regions, using the same magnetic map provided
/// in file `magnetic.file` and shifted by x=-30mm and x=30mm. The parameters
/// available in the `addMagneticVolume` definition are described in this list.
///
/// - *file* : This allows to specify the filename that contains the values of the
/// magnetic field. Few files will be found under `data/magneticField`. They all
/// contain 3 columns for the position in the volume and 3 columns to define the
/// magnetic field vector. It is not the full path but only the name of the file.
///
/// - *position* : By convention, the volume is build using the coordinates provided
/// in the magnetic field file given. However, it is possible to translate the
/// volume using the `position` field.
///
/// \todo Magnetic field volume rotations are not implemented. Are they necessary?
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
/// 2020-February: Reviewing TRestAxionMagneticField class.
///            Javier Galan
///
/// \class      TRestAxionMagneticField
/// \author     Eve Pachoud
/// \author     Javier Galan
///
/// <hr>
///

#include "TRestAxionMagneticField.h"

using namespace std;

ClassImp(TRestAxionMagneticField);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionMagneticField::TRestAxionMagneticField() : TRestMetadata() { Initialize(); }

///////////////////////////////////////////////
/// \brief Constructor loading data from a config file
///
/// If no configuration path is defined using TRestMetadata::SetConfigFilePath
/// the path to the config file must be specified using full path, absolute or
/// relative.
///
/// The default behaviour is that the config file must be specified with
/// full path, absolute or relative.
///
/// \param cfgFileName A const char* giving the path to an RML file.
/// \param name The name of the specific metadata. It will be used to find the
/// corresponding TRestG4Metadata section inside the RML.
///
TRestAxionMagneticField::TRestAxionMagneticField(const char* cfgFileName, string name)
    : TRestMetadata(cfgFileName) {
    debug << "Entering TRestAxionMagneticField constructor( cfgFileName, name )" << endl;

    Initialize();

    LoadConfigFromFile(fConfigFileName, name);

    PrintMetadata();
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionMagneticField::~TRestAxionMagneticField() {
    debug << "Entering ... TRestAxionMagneticField() destructor." << endl;
}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionMagneticField members
///
void TRestAxionMagneticField::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    fCanvas = NULL;
    fHisto = NULL;
}

///////////////////////////////////////////////
/// \brief A method that creates a canvas where magnetic field map is drawn
///
TCanvas* TRestAxionMagneticField::DrawHistogram(TString projection, TString Bcomp, Int_t VIndex,
                                                Double_t step) {
    /* Quarantined
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

if (VIndex >= GetNumberOfVolumes()) ferr << VIndex << " corresponds to none volume index " << endl;

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
    */

    return fCanvas;
}

///////////////////////////////////////////////
/// \brief A method to help loading magnetic field data, as x,y,z,Bx,By,Bz into a magnetic volume definition
/// using its corresponding mesh.
///
/// This method will be made private since it will only be used internally.
///
void TRestAxionMagneticField::LoadMagneticFieldData(MagneticFieldVolume& mVol,
                                                    std::vector<std::vector<Double_t>> data) {
    mVol.field.resize(mVol.mesh.GetNodesX());
    for (int n = 0; n < mVol.field.size(); n++) {
        mVol.field[n].resize(mVol.mesh.GetNodesY());
        for (int m = 0; m < mVol.field[n].size(); m++) mVol.field[n][m].resize(mVol.mesh.GetNodesZ());
    }

    debug << "TRestAxionMagneticField::LoadMagneticFieldData. Printing first 5 data rows" << endl;
    for (Int_t n = 0; n < data.size(); n++) {
        Int_t nX = mVol.mesh.GetNodeX(data[n][0]);
        Int_t nY = mVol.mesh.GetNodeY(data[n][1]);
        Int_t nZ = mVol.mesh.GetNodeZ(data[n][2]);

        if (n < 5) {
            debug << "X: " << data[n][0] << " Y: " << data[n][1] << " Z: " << data[n][2] << endl;
            debug << "nX: " << nX << " nY: " << nY << " nZ: " << nZ << endl;
            debug << "Bx: " << data[n][3] << " By: " << data[n][4] << " Bz: " << data[n][5] << endl;
        }

        mVol.field[nX][nY][nZ] = TVector3(data[n][3], data[n][4], data[n][5]);
    }
}

///////////////////////////////////////////////
/// \brief It will load the magnetic field data from the data filenames specified at the RML definition.
///
/// This method will be made private since it will only be used internally.
///
void TRestAxionMagneticField::LoadMagneticVolumes() {
    for (unsigned int n = 0; n < fPositions.size(); n++) {
        string fullPathName = SearchFile((string)fFileNames[n]);

        std::vector<std::vector<Double_t>> fieldData;
        TRestTools::ReadASCIITable(fullPathName, fieldData);

        // First line are the limits of field and size of mesh
        Double_t xmax = fieldData[0][0];
        Double_t ymax = fieldData[0][1];
        Double_t zmax = fieldData[0][2];

        Double_t sizeMesh = fieldData[0][3];

        // We keep in the vector only the field data. We remove first row
        fieldData.erase(fieldData.begin());

        if (GetVerboseLevel() >= REST_Debug) {
            cout << "Reading magnetic field map" << endl;
            cout << "--------------------------" << endl;

            cout << "Full path : " << fullPathName << endl;

            cout << "xMax: " << xmax << " yMax: " << ymax << " zMax: " << zmax << endl;
            cout << "Mesh size : " << sizeMesh << endl;

            cout << "Printing beginning of magnetic file table : " << fieldData.size() << endl;
            TRestTools::PrintTable(fieldData, 0, 5);
        }

        // Number of nodes
        Int_t nx = (Int_t)(2 * xmax / sizeMesh) + 1;
        Int_t ny = (Int_t)(2 * ymax / sizeMesh) + 1;
        Int_t nz = (Int_t)(2 * zmax / sizeMesh) + 1;

        // We create an auxiliar mesh helping to initialize the fieldMap
        // The mesh is centered at zero. Absolute position is defined in the Magnetic volume
        TRestMesh restMesh;
        restMesh.SetSize(2 * xmax, 2 * ymax, 2 * zmax);
        restMesh.SetOrigin(TVector3(-xmax, -ymax, -zmax));
        restMesh.SetNodes(nx, ny, nz);

        MagneticFieldVolume mVolume;
        mVolume.position = fPositions[n];
        mVolume.mesh = restMesh;

        LoadMagneticFieldData(mVolume, fieldData);

        fMagneticFieldVolumes.push_back(mVolume);
    }
}

///////////////////////////////////////////////
/// \brief It returns the magnetic field vector at x,y,z
///
TVector3 TRestAxionMagneticField::GetMagneticField(Double_t x, Double_t y, Double_t z) {
    return GetMagneticField(TVector3(x, y, z));
}

///////////////////////////////////////////////
/// \brief It returns the magnetic field vector at TVector3(pos)
///
TVector3 TRestAxionMagneticField::GetMagneticField(TVector3 pos) {
    for (int n = 0; n < fMagneticFieldVolumes.size(); n++) {
        TVector3 relativePosition = pos - fPositions[n];
        if (fMagneticFieldVolumes[n].mesh.IsInside(relativePosition)) {
            TVector3 node = GetMagneticVolumeNode(fMagneticFieldVolumes[n], pos);
            Int_t nX = node.X();
            Int_t nY = node.Y();
            Int_t nZ = node.Z();

            // We just return the field at nX, nY, nZ.
            // Which is the bottom,left,down node.
            // We could interpolate here using the other nodes field nX+1,nY+1,nZ+1
            return fMagneticFieldVolumes[n].field[nX][nY][nZ];
        }
    }
}

///////////////////////////////////////////////
/// \brief It returns the intensity of the transversal magnetic field component for the defined propagation
/// direction at the position given by argument.
///
Double_t TRestAxionMagneticField::GetTransversalComponent(TVector3 position, TVector3 direction) {
    return abs(GetMagneticField(position).Perp(direction));
}

///////////////////////////////////////////////
/// \brief It returns a vector describing the transversal magnetic field component between `from` and `to`
/// positions given by argument.
///
/// The differential element `dl` is by default 1mm, but it can be modified through the third argument of this
/// function.
///
/// The maximum number of divisions (unlimited by default) of the output vector can be fixed by the forth
/// argument. In that case, the differential element `dl` length might be increased to fullfil such condition.
///
std::vector<Double_t> TRestAxionMagneticField::GetTransversalComponentAlongPath(TVector3 from, TVector3 to,
                                                                                Double_t dl, Int_t Nmax) {
    Double_t length = (to - from).Mag();

    Double_t diff = dl;
    if (Nmax > 0) {
        if (length / dl > Nmax) {
            diff = length / Nmax;
            warning << "TRestAxionMagneticField::GetTransversalComponentAlongPath. Nmax reached!" << endl;
            warning << "Nmax = " << Nmax << endl;
            warning << "Adjusting differential step to : " << diff << " mm" << endl;
        }
    }

    TVector3 direction = (to - from).Unit();

    std::vector<Double_t> Bt;
    for (Double_t d = 0; d < length; d += diff) {
        Bt.push_back(GetTransversalComponent(from + d * direction, direction));
    }

    return Bt;
}

///////////////////////////////////////////////
/// \brief It returns the average of the transversal magnetic field intensity between the 3-d coordinates
/// `from` and `to`.
///
/// The differential element `dl` defines the integration step, and it is by default 1mm, but it can be
/// modified through the third argument of this function.
///
/// The maximum number of divisions of the output vector can be fixed by the forth argument. In that case, the
/// differential element `dl` length might be increased to fullfil such condition.
///
Double_t TRestAxionMagneticField::GetTransversalFieldAverage(TVector3 from, TVector3 to, Double_t dl,
                                                             Int_t Nmax) {
    Double_t length = (to - from).Mag();

    Double_t Bavg = 0.;
    std::vector<Double_t> Bt = GetTransversalComponentAlongPath(from, to, dl, Nmax);
    for (auto& b : Bt) Bavg += b;

    if (length > 0) return Bavg / length;

    error << "TRestAxionMagneticField::GetTransversalFieldAverage. Lenght is zero!" << endl;
    return 0.;
}

///////////////////////////////////////////////
/// \brief It returns the corresponding mesh node in the magnetic volume
///
/// The corresponging node to x,y,z is the bottom, down, left node in the cell volume
/// defined by 8-nodes.
///
/// This method will be made private, no reason to use it outside this class.
///
TVector3 TRestAxionMagneticField::GetMagneticVolumeNode(MagneticFieldVolume mVol, TVector3 pos) {
    TVector3 relPosition = pos - mVol.position;

    Int_t nx = mVol.mesh.GetNodeX(relPosition.X());
    Int_t ny = mVol.mesh.GetNodeY(relPosition.Y());
    Int_t nz = mVol.mesh.GetNodeZ(relPosition.Z());
    return TVector3(nx, ny, nz);
}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionMagnetic field members through a RML file
///
void TRestAxionMagneticField::InitFromConfigFile() {
    this->Initialize();

    string bVolume;
    size_t pos = 0;
    while ((bVolume = GetKEYDefinition("addMagneticVolume", pos)) != "") {
        string filename = GetFieldValue("fileName", bVolume);
        fFileNames.push_back(filename);

        TVector3 position = Get3DVectorFieldValueWithUnits("position", bVolume);
        fPositions.push_back(position);
    }

    LoadMagneticVolumes();
}

///////////////////////////////////////////////
/// \brief Prints on screen the details about the Geant4 simulation
/// conditions, stored in TRestG4Metadata.
///
void TRestAxionMagneticField::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    metadata << " - Number of magnetic volumes : " << GetNumberOfVolumes() << endl;
    metadata << " ------------------------------------------------ " << endl;
    for (int p = 0; p < GetNumberOfVolumes(); p++) {
        MagneticFieldVolume vol = fMagneticFieldVolumes[p];

        Double_t centerX = fPositions[p][0];
        Double_t centerY = fPositions[p][1];
        Double_t centerZ = fPositions[p][2];

        Double_t halfSizeX = vol.mesh.GetNetSizeX() / 2.;
        Double_t halfSizeY = vol.mesh.GetNetSizeY() / 2.;
        Double_t halfSizeZ = vol.mesh.GetNetSizeZ() / 2.;

        Double_t xMin = centerX - halfSizeX;
        Double_t yMin = centerY - halfSizeY;
        Double_t zMin = centerZ - halfSizeZ;

        Double_t xMax = centerX + halfSizeX;
        Double_t yMax = centerY + halfSizeY;
        Double_t zMax = centerZ + halfSizeZ;

        metadata << "* Volume " << p << " : "
                 << "  - Set in (" << centerX << "," << centerY << "," << centerZ << ")"
                 << " mm" << endl;
        metadata << "  - Bounds : " << endl;
        metadata << "    xmin : " << xMin << " mm , xmax : " << xMax << " mm" << endl;
        metadata << "    ymin : " << yMin << " mm, ymax : " << yMax << " mm" << endl;
        metadata << "    zmin : " << zMin << " mm, zmax : " << zMax << " mm" << endl;
        metadata << "  - File loaded : " << fFileNames[p] << endl;
    }
    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

