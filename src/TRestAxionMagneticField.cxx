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
/// The default value is "none", that identifies when no input field map is given.
/// If no magnetic field map is provided the field will be constant in the volume
/// definition, and fixed by the `field` parameter.
///
/// - *position* : By convention, the volume is build using the coordinates provided
/// in the magnetic field file given. However, it is possible to translate the
/// volume using the `position` field. The default value is (0,0,0).
///
/// - *field* : A 3D vector given in magnetic field units. If no units are given
/// the vector is assumed to be expressed in Teslas. This field will be added as an
/// offset to the magnetic field map given. The default value is (0,0,0).
///
/// - *boundMax* : A 3D vector, `(xMax,yMax,zMax)` that defines the bounding box
/// size. The box size will be bounded by the vertex `(xMin,yMin,zMin)` and
/// `(xMax,yMax,zMax)". This parameter is required if no field map file is given.
/// If a field map file is provided, the bounding box will be extracted from the
/// field map. In that case, this parameter will be only used for validation.
///
/// - *meshSize* : A 3D vector that defines the size of the cell in a regular
/// grid where the field is defined. This parameter is required if no field map
/// file is given. If a field map file is provided, the 3-dimensional cell element
/// dimensions will be extracted from the field map. In that case, this parameter
/// will be only used for validation.
///
/// - *meshType* : It defines the type of mesh boundary. The default value will
/// be cylinder. It defines a cylinder with its axis sitting on z. The radius of
/// the cylinder will be the first component of the *meshSize* 3D vector. The
/// second component will be ignored, as it is defined inside TRestMesh. This
/// will affect to the identification of field boundaries, and to the evaluation
/// of the field, which will only happen when the evaluated coordinates are
/// inside the bounding volume.
///
/// All parameteres are optional, and if not provided they will take their default
/// values.
///
/// \note If no magnetic field map is provided, i.e. we just want to define a
/// constant magnetic field vector, we still need to provide a boundary box size
/// and mesh size through the `boundMax` and `meshSize` parameters. We still have
/// the possibility to define a constant magnetic field vector for that volume
/// using the `field` parameter. In that case the method
/// TRestMagneticField::IsFieldConstant will return true.
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
/// There are two possible formats to provide a magnetic field map.
///
/// * **Plain-text format (.dat)** : The data file should describe the points of
/// a box delimited by the vertexes `(-xMax, -yMax, -zMax)` and `(xMax,yMax,zMax)`.
/// The field `(Bx,By,Bz)` must be given for all the points in a regular grid.
/// The default units are `T` for the field and `mm` for distances. The size of
/// the grid cells will be deduced from the data. Each row in the text file
/// should provide 6-tabulated values `x`, `y`, `z`, `Bx`, `By`, `Bz`.
///
/// * **Binary format (.bin)** : Similar to the plain-text format. It should
/// contain all the elements of the grid, bounded by (-xMax,-yMax,-zMax) and
/// (xMax,yMax,zMax). Each element is built using the 3-coordinates `x`, `y`, `z`
/// and the 3-field `Bx`, `By`, `Bz` values expressed in 4-bytes size, Float_t.
///
/// ### A more detailed example
///
/// The following example shows different allowed volume definition entries.
///
/// \code
///    <TRestAxionMagneticField name="bFieldBabyIAXO" title="First magnetic field definition"
///    verboseLevel="info" >
///
///		<!-- A volume from a text file centered at (0,0,0) in a Helium+Xenon gas mixture -->
///		<addMagneticVolume fileName="Bykovskiy_201906.dat" position="(0,0,0)mm"
///					gasMixture="He+Xe" gasDensity="1mg/dm3+1e-1mg/dm3" />
///
///		<!-- A volume from binary file, include Bx=1T offset, and boundMax validation -->
///		<addMagneticVolume fileName="Bykovskiy_202004.bin" position="(800,800,800)mm"
///					field="(1,0,0)"	boundMax="(350,350,4000)" />
///
///		<!-- A magnetic volume with a constant field definition (no field map) -->
///		<addMagneticVolume field="(0,0,3)T" position="(-800,-800,-800)mm"
///							boundMax="(100,200,300)" meshSize="(10,20,30)" />
///    </TRestAxionMagneticField>
/// \endcode
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

    if (GetVerboseLevel() >= REST_Info) PrintMetadata();
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

    // CAUTION. This needs revision just fixed adding .X() in order compile!
    if (step < 0) step = fMeshSize[volIndex].X();

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
                                                    std::vector<std::vector<Float_t>> data) {
    mVol.field.resize(mVol.mesh.GetNodesX());
    for (int n = 0; n < mVol.field.size(); n++) {
        mVol.field[n].resize(mVol.mesh.GetNodesY());
        for (int m = 0; m < mVol.field[n].size(); m++) mVol.field[n][m].resize(mVol.mesh.GetNodesZ());
    }

    debug << "TRestAxionMagneticField::LoadMagneticFieldData. Printing first 5 data rows" << endl;
    for (Int_t n = 0; n < data.size(); n++) {
        // The magnetic field map is centered at zero.
        // But the mesh definition contains the offset position
        // We shift the data to match the mesh node network.
        Int_t nX = mVol.mesh.GetNodeX((Int_t)(data[n][0] + mVol.mesh.GetNetSizeX() / 2.), true);
        Int_t nY = mVol.mesh.GetNodeY((Int_t)(data[n][1] + mVol.mesh.GetNetSizeY() / 2.), true);
        Int_t nZ = mVol.mesh.GetNodeZ((Int_t)(data[n][2] + mVol.mesh.GetNetSizeZ() / 2.), true);

        if (n < 5) {
            debug << "X: " << data[n][0] << " Y: " << data[n][1] << " Z: " << data[n][2] << endl;
            debug << "absX: " << data[n][0] + mVol.position.X() << " absY: " << data[n][1] + mVol.position.Y()
                  << " absZ: " << data[n][2] + mVol.position.Z() << endl;
            debug << "nX: " << nX << " nY: " << nY << " nZ: " << nZ << endl;
            debug << "Bx: " << data[n][3] << " By: " << data[n][4] << " Bz: " << data[n][5] << endl;
            if (GetVerboseLevel() >= REST_Extreme) GetChar();
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
            if (GetVerboseLevel() >= REST_Extreme) GetChar();
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
        debug << "Reading file : " << fFileNames[n] << endl;
        debug << "Full path : " << fullPathName << endl;

        std::vector<std::vector<Float_t>> fieldData;
        if (fullPathName.find(".dat") != string::npos) {
            debug << "Reading ASCII format" << endl;
            if (!TRestTools::ReadASCIITable(fullPathName, fieldData)) {
                ferr << "Problem reading file : " << fullPathName << endl;
                exit(1);
            }
        } else if (fullPathName.find(".bin") != string::npos) {
            debug << "Reading binary format" << endl;
            if (!TRestTools::ReadBinaryTable(fullPathName, fieldData, 6)) {
                ferr << "Problem reading file : " << fullPathName << endl;
                exit(2);
            }
        } else if (fFileNames[n] != "none") {
            ferr << "Filename : " << fullPathName << endl;
            ferr << "File format not recognized!" << endl;
            exit(3);
        }

        if (fFileNames[n] != "none" && fieldData.size() < 2) {
            ferr << "Field data size is no more than 2 grid points!" << endl;
            ferr << "Filename : " << fullPathName << endl;
            ferr << "Probably something went wrong loading the file" << endl;
            exit(4);
        }

        Float_t xMin = -fBoundMax[n].X(), yMin = -fBoundMax[n].Y(), zMin = -fBoundMax[n].Z();
        Float_t xMax = fBoundMax[n].X(), yMax = fBoundMax[n].Y(), zMax = fBoundMax[n].Z();
        Float_t meshSizeX = fMeshSize[n].X(), meshSizeY = fMeshSize[n].Y(), meshSizeZ = fMeshSize[n].Z();

        // If a field map is defined we get the boundaries, and mesh size from the volume
        if (fieldData.size() > 0) {
            debug << "Reading max boundary values" << endl;
            xMax = TRestTools::GetMaxValueFromTable(fieldData, 0);
            yMax = TRestTools::GetMaxValueFromTable(fieldData, 1);
            zMax = TRestTools::GetMaxValueFromTable(fieldData, 2);

            if (fBoundMax[n] != TVector3(0, 0, 0)) {
                if (fBoundMax[n] != TVector3(xMax, yMax, zMax)) {
                    warning << "Volume : " << n << endl;
                    warning << "boundMax was defined in RML but does not match the field map boundaries!"
                            << endl;
                    warning << "Max. Field map boundaries : (" << xMax << ", " << yMax << ", " << zMax << ")"
                            << endl;
                }
            }

            debug << "Reading min boundary values" << endl;
            xMin = TRestTools::GetMinValueFromTable(fieldData, 0);
            yMin = TRestTools::GetMinValueFromTable(fieldData, 1);
            zMin = TRestTools::GetMinValueFromTable(fieldData, 2);

            if (fBoundMax[n] != TVector3(0, 0, 0)) {
                if (-fBoundMax[n] != TVector3(xMin, yMin, zMin)) {
                    warning << "Volume : " << n << endl;
                    warning << "boundMax was defined in RML but does not match the field map boundaries"
                            << endl;
                    warning << "Min. Field map boundaries : (" << xMin << ", " << yMin << ", " << zMin << ")"
                            << endl;
                }
            }
            fBoundMax[n] = TVector3(xMax, yMax, zMax);

            debug << "Reading mesh size" << endl;
            meshSizeX = TRestTools::GetLowestIncreaseFromTable(fieldData, 0);
            meshSizeY = TRestTools::GetLowestIncreaseFromTable(fieldData, 1);
            meshSizeZ = TRestTools::GetLowestIncreaseFromTable(fieldData, 2);

            if (fMeshSize[n] != TVector3(0, 0, 0)) {
                if (fMeshSize[n] != TVector3(meshSizeX, meshSizeY, meshSizeZ)) {
                    warning << "Volume : " << n << endl;
                    warning << "MeshSize was defined in RML but does not match the mesh size deduced from "
                               "field map"
                            << endl;
                    warning << "Mesh size : (" << meshSizeX << ", " << meshSizeY << ", " << meshSizeZ << ")"
                            << endl;
                }
            }
            fMeshSize[n] = TVector3(meshSizeX, meshSizeY, meshSizeZ);
        }

        if (GetVerboseLevel() >= REST_Debug) {
            debug << "Reading magnetic field map" << endl;
            debug << "--------------------------" << endl;

            debug << "Full path : " << fullPathName << endl;

            debug << "Boundaries" << endl;
            debug << "xMin: " << xMin << " yMin: " << yMin << " zMin: " << zMin << endl;
            debug << "xMax: " << xMax << " yMax: " << yMax << " zMax: " << zMax << endl;
            debug << "Mesh size" << endl;

            debug << "sX: " << meshSizeX << " sY: " << meshSizeY << " sZ: " << meshSizeZ << endl;

            if (fieldData.size() > 4) {
                debug << "Printing beginning of magnetic file table : " << fieldData.size() << endl;
                TRestTools::PrintTable(fieldData, 0, 5);
            } else {
                debug << "The data file contains no field map" << endl;
            }
        }
        if (GetVerboseLevel() >= REST_Extreme) GetChar();

        // Number of nodes
        Int_t nx = (Int_t)(2 * xMax / meshSizeX) + 1;
        Int_t ny = (Int_t)(2 * yMax / meshSizeY) + 1;
        Int_t nz = (Int_t)(2 * zMax / meshSizeZ) + 1;

        // We create an auxiliar mesh helping to initialize the fieldMap
        // The mesh is centered at zero. Absolute position is defined in the Magnetic volume
        // TODO It would be interesting that TRestMesh could be used in cylindrical coordinates.
        TRestMesh restMesh;
        restMesh.SetSize(2 * xMax, 2 * yMax, 2 * zMax);
        restMesh.SetOrigin(fPositions[n] - TVector3(xMax, yMax, zMax));
        restMesh.SetNodes(nx, ny, nz);
        if (fMeshType[n] == "cylinder")
            restMesh.SetCylindrical(true);
        else
            restMesh.SetCylindrical(false);

        MagneticFieldVolume mVolume;
        mVolume.position = fPositions[n];
        mVolume.mesh = restMesh;

        if (mVolume.bGas == NULL) mVolume.bGas = new TRestAxionBufferGas();
        if (fGasMixtures[n] != "vacuum") mVolume.bGas->SetGasMixture(fGasMixtures[n], fGasDensities[n]);

        if (fieldData.size() > 0) LoadMagneticFieldData(mVolume, fieldData);

        if (fBoundMax[n] == TVector3(0, 0, 0)) {
            ferr << "The bounding box was not defined for volume " << n << "!" << endl;
            ferr << "Please review RML configuration for TRestAxionMagneticField" << endl;
            exit(22);
        } else if (fMeshSize[n] == TVector3(0, 0, 0)) {
            ferr << "The mesh grid size was not defined for volume " << n << "!" << endl;
            ferr << "Please review RML configuration for TRestAxionMagneticField" << endl;
            exit(22);
        }
        fMagneticFieldVolumes.push_back(mVolume);
    }

    if (CheckOverlaps()) {
        ferr << "TRestAxionMagneticField::LoadMagneticVolumes. Volumes overlap!" << endl;
        exit(1);
    }
    debug << "Finished loading magnetic volumes" << endl;
}

///////////////////////////////////////////////
/// \brief It returns the magnetic field vector at x,y,z
///
TVector3 TRestAxionMagneticField::GetMagneticField(Double_t x, Double_t y, Double_t z) {
    return GetMagneticField(TVector3(x, y, z));
}

///////////////////////////////////////////////
/// \brief It returns the magnetic field vector at TVector3(pos) using trilinear interpolation
/// that is implemented following instructions given at https://en.wikipedia.org/wiki/Trilinear_interpolation
///
TVector3 TRestAxionMagneticField::GetMagneticField(TVector3 pos) {
    Int_t id = GetVolumeIndex(pos);

    if (id < 0) {
        warning << "TRestAxionMagneticField::GetMagneticField position is outside any volume" << endl;
        return TVector3(0, 0, 0);
    } else {
        if (IsFieldConstant(id)) return fConstantField[id];
        TVector3 node = GetMagneticVolumeNode(fMagneticFieldVolumes[id], pos);
        Int_t nX = node.X();
        Int_t nY = node.Y();
        Int_t nZ = node.Z();

        Int_t nX_1, nY_1, nZ_1;

        if ((nX + 1) < fMagneticFieldVolumes[id].mesh.GetNodesX())
            nX_1 = nX + 1;
        else
            nX_1 = nX;
        if ((nY + 1) < fMagneticFieldVolumes[id].mesh.GetNodesY())
            nY_1 = nY + 1;
        else
            nY_1 = nY;
        if ((nZ + 1) < fMagneticFieldVolumes[id].mesh.GetNodesZ())
            nZ_1 = nZ + 1;
        else
            nZ_1 = nZ;

        TVector3 C000 = fMagneticFieldVolumes[id].field[nX][nY][nZ] + fConstantField[id];
        TVector3 C100 = fMagneticFieldVolumes[id].field[nX_1][nY][nZ] + fConstantField[id];
        TVector3 C010 = fMagneticFieldVolumes[id].field[nX][nY_1][nZ] + fConstantField[id];
        TVector3 C110 = fMagneticFieldVolumes[id].field[nX_1][nY_1][nZ] + fConstantField[id];
        TVector3 C001 = fMagneticFieldVolumes[id].field[nX][nY][nZ_1] + fConstantField[id];
        TVector3 C101 = fMagneticFieldVolumes[id].field[nX_1][nY][nZ_1] + fConstantField[id];
        TVector3 C011 = fMagneticFieldVolumes[id].field[nX][nY_1][nZ_1] + fConstantField[id];
        TVector3 C111 = fMagneticFieldVolumes[id].field[nX_1][nY_1][nZ_1] + fConstantField[id];

        Double_t x0 = fMagneticFieldVolumes[id].mesh.GetX(nX);
        Double_t x1 = fMagneticFieldVolumes[id].mesh.GetX(nX_1);
        Double_t xd;
        if (x0 == x1)
            xd = 0;
        else
            xd = (pos.X() - x0) / (x1 - x0);
        if ((xd < 0) || (xd > 1))
            warning << "TRestAxionMagneticField::GetMagneticField  Error: xd NOT between 0 and 1" << endl;

        Double_t y0 = fMagneticFieldVolumes[id].mesh.GetY(nY);
        Double_t y1 = fMagneticFieldVolumes[id].mesh.GetY(nY_1);
        Double_t yd;
        if (y0 == y1)
            yd = 0;
        else
            yd = (pos.Y() - y0) / (y1 - y0);
        if ((yd < 0) || (yd > 1))
            warning << "TRestAxionMagneticField::GetMagneticField  Error: yd NOT between 0 and 1" << endl;

        Double_t z0 = fMagneticFieldVolumes[id].mesh.GetZ(nZ);
        Double_t z1 = fMagneticFieldVolumes[id].mesh.GetZ(nZ_1);
        Double_t zd;
        if (z0 == z1)
            zd = 0;
        else
            zd = (pos.Z() - z0) / (z1 - z0);
        if ((zd < 0) || (zd > 1))
            warning << "TRestAxionMagneticField::GetMagneticField  Error: zd NOT between 0 and 1" << endl;

        // first we interpolate along x-axis
        TVector3 C00 = C000 * (1.0 - xd) + C100 * xd;
        TVector3 C01 = C001 * (1.0 - xd) + C101 * xd;
        TVector3 C10 = C010 * (1.0 - xd) + C110 * xd;
        TVector3 C11 = C011 * (1.0 - xd) + C111 * xd;

        // then we interpolate along y-axis
        TVector3 C0 = C00 * (1.0 - yd) + C10 * yd;
        TVector3 C1 = C01 * (1.0 - yd) + C11 * yd;

        // finally we interpolate along z-axis
        TVector3 C = C0 * (1.0 - zd) + C1 * zd;

        debug << "position = (" << pos.X() << ", " << pos.Y() << ", " << pos.Z() << ")       ";
        debug << "nX = " << nX << " nY = " << nY << " nZ = " << nZ << "     nX_1 = " << nX_1
              << "   nY_1 = " << nY_1 << "   nZ_1 = " << nZ_1 << endl
              << endl;
        debug << "C000 = (" << C000.X() << ", " << C000.Y() << ", " << C000.Z() << ")" << endl << endl;
        debug << "C100 = (" << C100.X() << ", " << C100.Y() << ", " << C100.Z() << ")" << endl << endl;
        debug << "C010 = (" << C010.X() << ", " << C010.Y() << ", " << C010.Z() << ")" << endl << endl;
        debug << "C110 = (" << C110.X() << ", " << C110.Y() << ", " << C110.Z() << ")" << endl << endl;
        debug << "C001 = (" << C001.X() << ", " << C001.Y() << ", " << C001.Z() << ")" << endl << endl;
        debug << "C101 = (" << C101.X() << ", " << C101.Y() << ", " << C101.Z() << ")" << endl << endl;
        debug << "C011 = (" << C011.X() << ", " << C011.Y() << ", " << C011.Z() << ")" << endl << endl;
        debug << "C111 = (" << C111.X() << ", " << C111.Y() << ", " << C111.Z() << ")" << endl << endl;
        debug << " -------------------------------------------------------" << endl;
        debug << "C00 = (" << C00.X() << ", " << C00.Y() << ", " << C00.Z() << ")" << endl << endl;
        debug << "C01 = (" << C01.X() << ", " << C01.Y() << ", " << C01.Z() << ")" << endl << endl;
        debug << "C10 = (" << C10.X() << ", " << C10.Y() << ", " << C10.Z() << ")" << endl << endl;
        debug << "C11 = (" << C11.X() << ", " << C11.Y() << ", " << C11.Z() << ")" << endl << endl;
        debug << " -------------------------------------------------------" << endl;
        debug << "C0 = (" << C0.X() << ", " << C0.Y() << ", " << C0.Z() << ")" << endl << endl;
        debug << "C1 = (" << C1.X() << ", " << C1.Y() << ", " << C1.Z() << ")" << endl << endl;
        debug << " -------------------------------------------------------" << endl;
        debug << "C = (" << C.X() << ", " << C.Y() << ", " << C.Z() << ")" << endl << endl;

        return C;
    }
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
    return GetPhotonAbsorptionLength(TVector3(x, y, z), en);
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
        if (fMagneticFieldVolumes[n].mesh.IsInside(pos)) return n;
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
/// `direction` and `position` given by argument.
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
    Int_t nx = mVol.mesh.GetNodeX(pos.X());
    Int_t ny = mVol.mesh.GetNodeY(pos.Y());
    Int_t nz = mVol.mesh.GetNodeZ(pos.Z());
    return TVector3(nx, ny, nz);
}

///////////////////////////////////////////////
/// \brief It will return true if the magnetic the regions overlap
///
Bool_t TRestAxionMagneticField::CheckOverlaps() {
    debug << "Checking overlaps" << endl;
    for (int n = 0; n < GetNumberOfVolumes(); n++) {
        for (int m = 0; m < GetNumberOfVolumes(); m++) {
            if (m == n) continue;
            debug << "Volume : " << m << endl;

            TVector3 b = GetMagneticVolume(m)->mesh.GetVertex(0);
            debug << "Relative bottom vertex : (" << b.X() << ", " << b.Y() << ", " << b.Z() << ")" << endl;
            if (GetMagneticVolume(n)->mesh.IsInside(b)) return true;

            TVector3 t = GetMagneticVolume(m)->mesh.GetVertex(1);
            debug << "Relative top vertex : (" << t.X() << ", " << t.Y() << ", " << t.Z() << ")" << endl;
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

    if (precision == 0) precision = min(fMeshSize[id].X(), min(fMeshSize[id].Y(), fMeshSize[id].Z())) / 2.;

    TVector3 unit = dir.Unit();
    std::vector<TVector3> fieldBoundaries;

    TVector3 in = volumeBoundaries[0];
    while ((GetTransversalComponent(in, dir) == 0) && (((volumeBoundaries[1] - in) * dir) > 0)) in = MoveByDistanceFast(in, unit, precision);
    if (((volumeBoundaries[1] - in) * dir) > 0)
        fieldBoundaries.push_back(in);
    else return fieldBoundaries;

    TVector3 out = volumeBoundaries[1];
    while ((GetTransversalComponent(out, -dir) == 0) && (((volumeBoundaries[0] - out) * dir) < 0) && (((out - in) * dir) > 0)) out = MoveByDistanceFast(out, -unit, precision);
    if ((((volumeBoundaries[0] - out) * dir) < 0) && (((out - in) * dir) > 0))
        fieldBoundaries.push_back(out);
    else return fieldBoundaries;

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
            fFileNames.push_back("none");
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

        TVector3 boundMax = Get3DVectorFieldValueWithUnits("boundMax", bVolume);
        if (boundMax == TVector3(-1, -1, -1))
            fBoundMax.push_back(TVector3(0, 0, 0));
        else
            fBoundMax.push_back(boundMax);

        TVector3 meshSize = Get3DVectorFieldValueWithUnits("meshSize", bVolume);
        if (meshSize == TVector3(-1, -1, -1))
            fMeshSize.push_back(TVector3(0, 0, 0));
        else
            fMeshSize.push_back(meshSize);

        string type = GetFieldValue("meshType", bVolume);
        if (type == "Not defined")
            fMeshType.push_back("cylinder");
        else
            fMeshType.push_back(type);

        // TRestMesh will only consider the first bounding component anyway
        if (fMeshType.back() == "cylinder" && fBoundMax.back().X() != fBoundMax.back().Y()) {
            warning << "Mesh type is cylinder. But X and Y inside boundMax are not the same!" << endl;
            warning << "Making second bound component Y equal to the X bound component!" << endl;
            fBoundMax.back().SetY(fBoundMax.back().X());
        }

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
        debug << "Max bounding box ( " << boundMax.X() << ", " << boundMax.Y() << ", " << boundMax.Z() << ")"
              << endl;
        debug << "Mesh size ( " << meshSize.X() << ", " << meshSize.Y() << ", " << meshSize.Z() << ")"
              << endl;
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
        metadata << "  - Offset field [T] : (" << fConstantField[p].X() << ", " << fConstantField[p].Y()
                 << ", " << fConstantField[p].Z() << ")" << endl;
        metadata << "  - File loaded : " << fFileNames[p] << endl;
        metadata << "  - Buffer gas mixture : " << fGasMixtures[p] << endl;
        metadata << "  - Buffer gas densities : " << fGasDensities[p] << endl;
        metadata << "  - Bounds : " << endl;
        metadata << "    xmin : " << xMin << " mm , xmax : " << xMax << " mm" << endl;
        metadata << "    ymin : " << yMin << " mm, ymax : " << yMax << " mm" << endl;
        metadata << "    zmin : " << zMin << " mm, zmax : " << zMax << " mm" << endl;
    }
    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

