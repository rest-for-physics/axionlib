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
/// defined magnetic fields, and create magnetic volume regions associated to
/// those pre-generated definitions. Once the field maps have been loaded this
/// class will be able to evaluate the field vector at any coordinate (x,y,z).
/// If the coordinate (x,y,z) is ourside any defined region, the returned
/// field will be (0,0,0).
///
/// TODO Description of magnetic field interpolation
///
/// ### RML definition
///
/// We can add any number of magnetic volumes inside the RML definition
/// as shown in the following piece of code,
///
/// \code
/// <TRestAxionMagneticField>
///         <addMagneticVolume file="magnetic.file" position="(30,0,0)cm" />
///         <addMagneticVolume file="magnetic.file" position="(-30,0,0)cm" />
/// <TRestAxionMagneticField/>
/// \endcode
///
/// where we produce 2 magnetic regions, using the same magnetic map provided
/// in file `magnetic.file` and shifted by x=-30cm and x=30cm. The following
/// parameters might be defined at each `addMagneticVolume` entry.
///
/// - *file* : This allows to specify the filename that contains the values of the
/// magnetic field. Few files will be found under `data/magneticField`. They all
/// contain 3 columns for the position in the volume and 3 columns to define the
/// magnetic field vector. It is not the full path but only the name of the file.
/// The default value is "noField.dat".
///
/// - *position* : By convention, the volume is build using the coordinates provided
/// in the magnetic field file given. However, it is possible to translate the
/// volume using the `position` field. The default value is (0,0,0).
///
/// - *field* : A 3D vector given in magnetic field units. If no units are given
/// the vector is assumed to be expressed in Teslas. This field will be added as an
/// offset to the magnetic field map given. The default value is (0,0,0).
///
/// All parameteres are optional, and if not provided they will take their default
/// values.
///
/// \note If no magnetic field map is provided we still need to provide a file
/// defining the boundary box size and mesh size, such as the default file provided
/// `data/magneticFile/noField.dat`. We can still define a constant magnetic field
/// vector for that volume using the `field` parameter. In that case
/// the method TRestMagneticField::IsFieldConstant will return true.
///
/// ### Adding gas properties to each of the magnetic volumes.
///
/// On top of the magnetic field map that we associate to a magnetic region (or
/// volume) we can - optionally - assign gas properties to the volume. This is done
/// using an interface to the TRestAxionBufferGas. And we may initialize as follows
/// at the `addMagneticVolume` entry in the corresponding RML section.
///
/// \code
/// <TRestAxionMagneticField>
///         <addMagneticVolume file="magnetic.file" position="(30,0,0)cm"
///							gasMixture="He" gasDensity="0.1mg/dm3"/>
///         <addMagneticVolume file="magnetic.file" position="(-30,0,0)cm"
///							gasMixture="He+Xe" gasDensity="3g/m3+0.2mg/dm3"/>
/// <TRestAxionMagneticField/>
/// \endcode
///
/// \note If no magnetic field map, neither a constant field is provided, this
/// class might still serve to define different regions distributed in the space
/// that describe different gas properties, even though the magnetic field will
/// be equal to zero.
///
/// ### The magnetic field file format.
///
/// For the moment, the file used for magnetic field description is written in
/// plain-text format. The data file should describe a box delimited by the
/// vertexes `(-xMax, -yMax, -zMax)` and `(xMax,yMax,zMax)`. The field `(Bx,By,Bz)`
/// is given for all the points in a regular grid inside that box with a grid
/// element given by `meshSize`, default units are `T` for the field and `mm` for
/// distances.
///
/// The first row in the file corresponds with the values of `xMax`, `yMax`,
/// `zMax` and `meshSize`. The next rows provide the magnetic field vector
/// at each of the elements of the grid,  where we will find the followoing
/// data components `x`, `y`, `z`  and `Bx`, `By`,`Bz`.
///
/// ### Visualizing the magnetic field
///
/// TODO Review and validate DrawHistogram drawing method and describe its
/// use here.
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
/// 2020-April: Reviewing and validating TRestAxionMagneticField class.
///            Javier Galan and Krešimir Jakovčić
///
/// \class      TRestAxionMagneticField
/// \author     Eve Pachoud
/// \author     Javier Galan <javier.galan@unizar.es>
/// \author     Krešimir Jakovčić <kjakov@irb.hr>
///
/// <hr>
///

#include "TRestAxionMagneticField.h"

using namespace std;

#include "TRestPhysics.h"
using namespace REST_Physics;

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
/// corresponding TRestAxionMagneticField section inside the RML.
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
    for (int n = 0; n < fMagneticFieldVolumes.size(); n++) delete fMagneticFieldVolumes[n].bGas;
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
/// TODO Add detailed documentation here
///
TCanvas* TRestAxionMagneticField::DrawHistogram(TString projection, TString Bcomp, Int_t volIndex,
                                                Double_t step) {
    if (!FieldLoaded()) LoadMagneticVolumes();

    if (fCanvas != NULL) {
        delete fCanvas;
        fCanvas = NULL;
    }

    if (fHisto != NULL) {
        delete fHisto;
        fHisto = NULL;
    }

    if (volIndex < 0) volIndex = 0;

    if (volIndex >= GetNumberOfVolumes()) {
        ferr << volIndex << " corresponds to none volume index " << endl;
        ferr << "Total number of volumes : " << GetNumberOfVolumes() << endl;
        ferr << "Setting volIndex to the first volume" << endl;
        volIndex = 0;
    }

    if (step < 0) step = fMeshSize[volIndex];

    MagneticFieldVolume* vol = GetMagneticVolume(volIndex);
    if (!vol) return fCanvas;

    Double_t centerX = fPositions[volIndex][0];
    Double_t centerY = fPositions[volIndex][1];
    Double_t centerZ = fPositions[volIndex][2];

    Double_t halfSizeX = vol->mesh.GetNetSizeX() / 2.;
    Double_t halfSizeY = vol->mesh.GetNetSizeY() / 2.;
    Double_t halfSizeZ = vol->mesh.GetNetSizeZ() / 2.;

    Double_t xMin = centerX - halfSizeX;
    Double_t yMin = centerY - halfSizeY;
    Double_t zMin = centerZ - halfSizeZ;

    Double_t xMax = centerX + halfSizeX;
    Double_t yMax = centerY + halfSizeY;
    Double_t zMax = centerZ + halfSizeZ;

    Int_t nBinsX = (xMax - xMin) / step;
    Int_t nBinsY = (yMax - yMin) / step;
    Int_t nBinsZ = (zMax - zMin) / step;

    Double_t x, y, z;
    Double_t B;
    TVector3 Bvec;

    if (projection == "XY") {
        fCanvas = new TCanvas("fCanvas", "");
        fHisto = new TH2D("", "", nBinsX, xMin, xMax, nBinsY, yMin, yMax);

        z = (zMin + zMax) / 2.0;
        x = xMin;

        for (Int_t i = 0; i < nBinsX; i++) {
            y = yMin;
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
            fHisto = new TH2D("", "", nBinsX, xMin, xMax, nBinsZ, zMin, zMax);

            y = (yMin + yMax) / 2.0;
            x = xMin;

            for (Int_t i = 0; i < nBinsX; i++) {
                z = zMin;
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
                fHisto = new TH2D("", "", nBinsY, yMin, yMax, nBinsZ, zMin, zMax);

                x = (xMin + xMax) / 2.0;
                y = yMin;

                for (Int_t i = 0; i < nBinsY; i++) {
                    z = zMin;
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
        Int_t nX = mVol.mesh.GetNodeX((Int_t)data[n][0]);
        Int_t nY = mVol.mesh.GetNodeY((Int_t)data[n][1]);
        Int_t nZ = mVol.mesh.GetNodeZ((Int_t)data[n][2]);

        if (n < 5) {
            debug << "X: " << data[n][0] << " Y: " << data[n][1] << " Z: " << data[n][2] << endl;
            debug << "nX: " << nX << " nY: " << nY << " nZ: " << nZ << endl;
            debug << "Bx: " << data[n][3] << " By: " << data[n][4] << " Bz: " << data[n][5] << endl;
        }

        if (mVol.field[nX][nY][nZ] != TVector3(0.0, 0.0, 0.0)) {
            warning << "X: " << data[n][0] << " Y: " << data[n][1] << " Z: " << data[n][2] << endl;
            warning << "nX: " << nX << " nY: " << nY << " nZ: " << nZ << endl;
            warning << "WARNING: field[nX][nY][nZ] element not equal to initial value (0, 0, 0) !!" << endl;
            warning << "It has value: "
                    << "mVol.field[" << nX << "][" << nY << "][" << nZ << "] = ("
                    << mVol.field[nX][nY][nZ].X() << " , " << mVol.field[nX][nY][nZ].Y() << " , "
                    << mVol.field[nX][nY][nZ].Z() << ")" << endl;
            warning << "Values to write: "
                    << "Bx: " << data[n][3] << " By: " << data[n][4] << " Bz: " << data[n][5] << endl
                    << endl;

            this->SetError("There was a problem assigning the field matrix!");
            GetChar();
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
        if (!TRestTools::ReadASCIITable(fullPathName, fieldData)) {
            ferr << "Problem reading file : " << fullPathName << endl;
            exit(1);
        }

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

            if (fieldData.size() > 4) {
                cout << "Printing beginning of magnetic file table : " << fieldData.size() << endl;
                TRestTools::PrintTable(fieldData, 0, 5);
            } else {
                cout << "The data file constains no field map" << endl;
            }
        }

        // Number of nodes
        Int_t nx = (Int_t)(2 * xmax / sizeMesh) + 1;
        Int_t ny = (Int_t)(2 * ymax / sizeMesh) + 1;
        Int_t nz = (Int_t)(2 * zmax / sizeMesh) + 1;

        fMeshSize.push_back(sizeMesh);

        // We create an auxiliar mesh helping to initialize the fieldMap
        // The mesh is centered at zero. Absolute position is defined in the Magnetic volume
        // TODO It would be interesting that TRestMesh could be used in cylindrical coordinates.
        TRestMesh restMesh;
        restMesh.SetSize(2 * xmax, 2 * ymax, 2 * zmax);
        restMesh.SetOrigin(TVector3(-xmax, -ymax, -zmax));
        restMesh.SetNodes(nx, ny, nz);

        MagneticFieldVolume mVolume;
        mVolume.position = fPositions[n];
        mVolume.mesh = restMesh;

        if (mVolume.bGas == NULL) mVolume.bGas = new TRestAxionBufferGas();
        if (fGasMixtures[n] != "vacuum") mVolume.bGas->SetGasMixture(fGasMixtures[n], fGasDensities[n]);

        if (fieldData.size() > 0) LoadMagneticFieldData(mVolume, fieldData);

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
    Int_t id = GetVolumeIndex(pos);

    if (id >= 0) {
        if (IsFieldConstant(id)) return fConstantField[id];
        TVector3 node = GetMagneticVolumeNode(fMagneticFieldVolumes[id], pos);
        Int_t nX = node.X();
        Int_t nY = node.Y();
        Int_t nZ = node.Z();

        // We just return the field at nX, nY, nZ.
        // Which is the bottom,left,down node.
        // We could interpolate here using the other nodes field nX+1,nY+1,nZ+1
        return fMagneticFieldVolumes[id].field[nX][nY][nZ] + fConstantField[id];

    } else {
        warning << "TRestAxionMagneticField::GetMagneticField position is outside any volume" << endl;
    }
    return TVector3(0, 0, 0);
}

///////////////////////////////////////////////
/// \brief It returns the effective photon mass in eV for the given position `(x,y,z)` and energy `en` using
/// the gas properties defined at the corresponding magnetic volume.
///
Double_t TRestAxionMagneticField::GetPhotonMass(Double_t x, Double_t y, Double_t z, Double_t en) {
    return GetPhotonMass(TVector3(x, y, z), en);
}

///////////////////////////////////////////////
/// \brief It returns the effective photon mass in eV for the given position `pos` and energy `en` using the
/// gas properties defined at the corresponding magnetic volume.
///
Double_t TRestAxionMagneticField::GetPhotonMass(TVector3 pos, Double_t en) {
    Int_t id = GetVolumeIndex(pos);

    return GetPhotonMass(id, en);
}

///////////////////////////////////////////////
/// \brief It returns the effective photon mass in eV at the corresponding magnetic volume id.
///
Double_t TRestAxionMagneticField::GetPhotonMass(Int_t id, Double_t en) {
    if (id >= 0) {
        return fMagneticFieldVolumes[id].bGas->GetPhotonMass(en);
    } else {
        warning << "TRestAxionMagneticField::GetPhotonMass position is outside any volume" << endl;
    }
    return -1;
}

///////////////////////////////////////////////
/// \brief It returns the photon absorption length in cm-1 for the given position `(x,y,z)` and energy `en`
/// using the gas properties defined at the corresponding magnetic volume.
///
Double_t TRestAxionMagneticField::GetPhotonAbsorptionLength(Double_t x, Double_t y, Double_t z, Double_t en) {
    return GetPhotonMass(TVector3(x, y, z), en);
}

///////////////////////////////////////////////
/// \brief It returns the photon absorption length in cm-1 for the given position `pos` and energy `en` using
/// the gas properties defined at the corresponding magnetic volume.
///
Double_t TRestAxionMagneticField::GetPhotonAbsorptionLength(TVector3 pos, Double_t en) {
    Int_t id = GetVolumeIndex(pos);

    return GetPhotonAbsorptionLength(id, en);
}

///////////////////////////////////////////////
/// \brief It returns the photon absorption length in cm-1 for the given position `pos` and energy `en` at the
/// given volume id.
///
Double_t TRestAxionMagneticField::GetPhotonAbsorptionLength(Int_t id, Double_t en) {
    if (id >= 0) {
        return fMagneticFieldVolumes[id].bGas->GetPhotonAbsorptionLength(en);
    } else {
        warning << "TRestAxionMagneticField::GetPhotonAbsorptionLength position is outside any volume"
                << endl;
    }
    return -1;
}

///////////////////////////////////////////////
/// \brief it returns the corresponding volume index at the given position. If not found it will return
/// -1.
///
Int_t TRestAxionMagneticField::GetVolumeIndex(TVector3 pos) {
    if (!FieldLoaded()) LoadMagneticVolumes();

    for (int n = 0; n < fMagneticFieldVolumes.size(); n++) {
        TVector3 relativePosition = pos - fPositions[n];
        if (fMagneticFieldVolumes[n].mesh.IsInside(relativePosition)) return n;
    }
    return -1;
}

///////////////////////////////////////////////
/// \brief it returns the volume position (or center) for the given volume `id`.
///
TVector3 TRestAxionMagneticField::GetVolumeCenter(Int_t id) { return GetVolumePosition(id); }

///////////////////////////////////////////////
/// \brief it returns the volume position (or center) for the given volume `id`.
///
TVector3 TRestAxionMagneticField::GetVolumePosition(Int_t id) {
    if (GetNumberOfVolumes() > id)
        return fPositions[id];
    else {
        warning << "TRestAxionMagneticField::GetVolumePosition. Id : " << id << " out of range!" << endl;
        warning << "Number of volumes defined : " << GetNumberOfVolumes() << endl;
        return TVector3(0, 0, 0);
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

    ferr << "TRestAxionMagneticField::GetTransversalFieldAverage. Lenght is zero!" << endl;
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
/// \brief It will return true if the magnetic the regions overlap
///
Bool_t TRestAxionMagneticField::CheckOverlaps() {
    for (int n = 0; n < GetNumberOfVolumes(); n++) {
        for (int m = n + 1; m < GetNumberOfVolumes(); m++) {
            TVector3 b = GetMagneticVolume(m)->mesh.GetVertex(0) - fPositions[m];
            TVector3 t = GetMagneticVolume(m)->mesh.GetVertex(1) - fPositions[m];
            if (GetMagneticVolume(n)->mesh.IsInside(b)) return true;
            if (GetMagneticVolume(n)->mesh.IsInside(t)) return true;
        }
    }
    return false;
}

///////////////////////////////////////////////
/// \brief Finds the in/out particle trajectory boundaries for a particular magnetic region bounding box.
///
/// This method checks if the trajectory defined by the position `pos` and direction `dir` passes through
/// the magnetic ﬁeld region/volume `id` given. If two such points (entry point and exit point) are found,
/// their coordinates are returned. In the example shown in Fig. 1 from TRestAxionFieldPropagationProcess
/// these points are: IN 1 and OUT 1 for the region #1 and IN2 and OUT 2 for the region #2.
///
/// If no intersection is found, or the particle is not moving towards the volume, the returned
/// std::vector
/// will be empty.
///
std::vector<TVector3> TRestAxionMagneticField::GetVolumeBoundaries(Int_t id, TVector3 pos, TVector3 dir) {
    MagneticFieldVolume* vol = GetMagneticVolume(id);

    std::vector<TVector3> boundaries;

    if (vol) boundaries = vol->mesh.GetTrackBoundaries(pos, dir);

    return boundaries;
}

///////////////////////////////////////////////
/// \brief Finds the in/out particle trajectory boundaries for a particular magnetic region, similar to
/// the method TRestAxionMagneticField::GetVolumeBoudaries, but requiring that the in/out points are the
/// first/last points where the **transversal** field intensity is not zero.
///
/// If no precision is given, the mesh size of the corresponding volume will be used as reference. The
/// precision will be meshSize/2.
///
/// If no intersection is found the returned std::vector will be empty.
///
std::vector<TVector3> TRestAxionMagneticField::GetFieldBoundaries(Int_t id, TVector3 pos, TVector3 dir,
                                                                  Double_t precision) {
    std::vector<TVector3> volumeBoundaries = GetVolumeBoundaries(id, pos, dir);
    if (volumeBoundaries.size() != 2) return volumeBoundaries;

    if (IsFieldConstant(id)) return volumeBoundaries;

    MagneticFieldVolume* vol = GetMagneticVolume(id);
    if (!vol) return volumeBoundaries;

    if (precision == 0) precision = fMeshSize[id] / 2.;

    TVector3 unit = dir.Unit();

    TVector3 in = volumeBoundaries[0];
    while (GetTransversalComponent(in, dir) == 0) in = MoveByDistanceFast(in, unit, precision);

    TVector3 out = volumeBoundaries[1];
    while (GetTransversalComponent(out, -dir) == 0) out = MoveByDistanceFast(out, -unit, precision);

    std::vector<TVector3> fieldBoundaries;
    fieldBoundaries.push_back(in);
    fieldBoundaries.push_back(out);

    return fieldBoundaries;
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
        if (filename == "Not defined")
            fFileNames.push_back("noField.dat");
        else
            fFileNames.push_back(filename);

        TVector3 position = Get3DVectorFieldValueWithUnits("position", bVolume);
        if (position == TVector3(-1, -1, -1))
            fPositions.push_back(TVector3(0, 0, 0));
        else
            fPositions.push_back(position);

        TVector3 field = Get3DVectorFieldValueWithUnits("field", bVolume);
        if (field == TVector3(-1, -1, -1))
            fConstantField.push_back(TVector3(0, 0, 0));
        else
            fConstantField.push_back(field);

        TString gasMixture = GetFieldValue("gasMixture", bVolume);
        if (gasMixture == "Not defined") gasMixture = "vacuum";
        fGasMixtures.push_back(gasMixture);

        TString gasDensity = GetFieldValue("gasDensity", bVolume);
        if (gasDensity == "Not defined") gasDensity = "0";
        fGasDensities.push_back(gasDensity);

        debug << "Reading new magnetic volume" << endl;
        debug << "-----" << endl;
        debug << "Filename : " << filename << endl;
        debug << "Position: ( " << position.X() << ", " << position.Y() << ", " << position.Z() << ") mm"
              << endl;
        debug << "Field: ( " << field.X() << ", " << field.Y() << ", " << field.Z() << ") T" << endl;
        debug << "Gas mixture : " << gasMixture << endl;
        debug << "Gas density : " << gasDensity << endl;
        debug << "----" << endl;
    }

    LoadMagneticVolumes();

    // TODO we should check that the volumes do not overlap
}

///////////////////////////////////////////////
/// \brief Prints on screen the information about the metadata members of TRestAxionMagneticField
///
void TRestAxionMagneticField::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    metadata << " - Number of magnetic volumes : " << GetNumberOfVolumes() << endl;
    metadata << " ------------------------------------------------ " << endl;
    for (int p = 0; p < GetNumberOfVolumes(); p++) {
        if (p > 0) metadata << " ------------------------------------------------ " << endl;
        MagneticFieldVolume* vol = GetMagneticVolume(p);

        Double_t centerX = fPositions[p][0];
        Double_t centerY = fPositions[p][1];
        Double_t centerZ = fPositions[p][2];

        Double_t halfSizeX = vol->mesh.GetNetSizeX() / 2.;
        Double_t halfSizeY = vol->mesh.GetNetSizeY() / 2.;
        Double_t halfSizeZ = vol->mesh.GetNetSizeZ() / 2.;

        Double_t xMin = centerX - halfSizeX;
        Double_t yMin = centerY - halfSizeY;
        Double_t zMin = centerZ - halfSizeZ;

        Double_t xMax = centerX + halfSizeX;
        Double_t yMax = centerY + halfSizeY;
        Double_t zMax = centerZ + halfSizeZ;

        metadata << "* Volume " << p << " centered at  (" << centerX << "," << centerY << "," << centerZ
                 << ") mm" << endl;
        metadata << "  - Bounds : " << endl;
        metadata << "    xmin : " << xMin << " mm , xmax : " << xMax << " mm" << endl;
        metadata << "    ymin : " << yMin << " mm, ymax : " << yMax << " mm" << endl;
        metadata << "    zmin : " << zMin << " mm, zmax : " << zMax << " mm" << endl;
        metadata << "  - File loaded : " << fFileNames[p] << endl;
        metadata << "  - Buffer gas mixture : " << fGasMixtures[p] << endl;
        metadata << "  - Buffer gas densities : " << fGasDensities[p] << endl;
    }
    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

