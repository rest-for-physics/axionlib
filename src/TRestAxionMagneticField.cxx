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
/// `(xMax,yMax,zMax)`. This parameter is required if no field map file is given.
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
/// All parameters are optional, and if not provided they will take their default
/// values.
///
/// \note If no magnetic field map is provided, i.e. we just want to define a
/// constant magnetic field vector, we still need to provide a boundary box size
/// and mesh size through the `boundMax` and `meshSize` parameters. We still have
/// the possibility to define a constant magnetic field vector for that volume
/// using the `field` parameter. In that case the method
/// TRestMagneticField::IsFieldConstant will return true.
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
///		<!-- A volume from a text file centered at (0,0,0) -->
///		<addMagneticVolume fileName="Bykovskiy_201906.dat" position="(0,0,0)mm" />
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
/// ### Using this class
///
/// Once we have created an instance of this class we will be able to access the
/// magnetic field directly, at any point by using trilinear interpolation.
///
/// For example, the following code will retrieve the field at a position found
/// at a distance z=-2m at the magnet axis,
///
/// \code
/// TRestAxionMagneticField *mag = new TRestAxionMagneticField("fields.rml", "babyIAXO" );
/// mag->GetMagneticField( 0, 0, -2000 );
/// \endcode
///
/// or in a similar way,
///
/// \code
/// TRestAxionMagneticField *mag = new TRestAxionMagneticField("fields.rml", "babyIAXO" );
/// mag->GetMagneticField( TVector3(0, 0, -2000), true );
/// \endcode
///
/// where the second argument, `true` will enable the warning message system.
///
/// The following code evaluates the transversal field component at
/// 1m distance along the Z-axis, for a vector pointing in the direction (0,1,1).
///
/// \code
/// TRestAxionMagneticField *mag = new TRestAxionMagneticField("fields.rml", "babyIAXO" );
/// mag->GetTransversalComponent( TVector3(0,0,1000), TVector3(0,1,1) );
/// \endcode
///
/// There are also geometric functions that allow to identify the boundaries of
/// the magnetic volume. A test particle will penetrate in the bounding box and identify
/// the moment where the field changes to a value different from (0,0,0) in order
/// to identify the entrance and exit point.
///
/// The following code will return two TVector3 with the magnetic volume entrance
/// and exit coordinates, for a test particle being placed at x=10cm, y=10cm and
/// z=-4m, pointing towards the magnetic field with a slight direction deviation from
/// the magnet axis.
///
/// \code
/// TRestAxionMagneticField *mag = new TRestAxionMagneticField("fields.rml", "babyIAXO" );
/// mag->GetFieldBoundaries(TVector3(100,100, -4000) , TVector3(0.02, 0.03, 1))
/// \endcode
///
/// In the other hand, TRestAxionMagneticField::GetVolumeBoundaries will return the
/// bounding box containing that magnetic field.
///
/// ### Visualizing the magnetic field
///
/// You may visualize the magnetic field profile along tracks towards a vanishing point
/// by using the TRestAxionMagneticField::DrawTracks method, using:
///
/// \code
///   TRestAxionMagneticField *field = new TRestAxionMagneticField( "fields.rml", "babyIAXO" );
///   field->DrawTracks( TVector3(0,0,8000), 100 );
/// \endcode
///
/// That will produce the following plot:
///
/// \htmlonly <style>div.image img[src="trackBprofile.png"]{width:800px;}</style> \endhtmlonly
///
/// ![Tracks through the magnetic field volume and its corresponding T-field component](trackBprofile.png)
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
///             Javier Galan and Krešimir Jakovčić
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

#include "TGraph.h"
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
    RESTDebug << "Entering TRestAxionMagneticField constructor( cfgFileName, name )" << RESTendl;

    Initialize();

    LoadConfigFromFile(fConfigFileName, name);

    if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Info) PrintMetadata();
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionMagneticField::~TRestAxionMagneticField() {
    RESTDebug << "Entering ... TRestAxionMagneticField() destructor." << RESTendl;
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
/// This method is used to create various 2D visualisations of the `Bx`, `By` and `Bz`
/// magnetic field components for a given magnetic field region.
///
/// The following input parameters can be specified :
///
/// - *projection* : it specifies the plane in the magnetic field region from which
/// the plotted values of the magnetic field component were taken. The allowed values
/// are: "XY", "XZ" and "YZ" for `X-Y`, `X-Z` and `Y-Z` plane, respectively.
///
/// - *Bcomp* : it specifies which component of the magnetic field is plotted.
/// The allowed values are: "X", "Y" and "Z" for `Bx`, `By` and `Bz` components,
/// respectively.
///
/// - *volIndex* : it specifies the index of the magnetic field region/volume for
/// which the plot is shown.
///
/// - *step* : it specifies the step/bin size for the 2D histogram that is shown.
/// If this parameter is not specified, or if it has value < = 0, the values are taken
/// from the corresponding mesh size values stored in data member `fMeshSize`.
///
/// - *style* : it specifies the plotting style for the histogram. It can correspond
/// to any draw option for 2D histograms in ROOT. The list of options can be found at:
/// https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html#drawing-histograms
/// The default value is "COLZ0" which draws a box for each cell in the histogram with a
/// color scale varying with values. The other useful option is "SURF3Z" which draws a
/// surface plot with a coloured contour view on the top.
///
/// - *depth* : it specifies the position of the plane in the magnetic field region from which
/// the plotted values of the magnetic field component were taken. If this parameter is not
/// specified, the histogram is plotted for the plane that goes through the middle of the region.
///
/// The following example will plot the values of the `Bx` component of the magnetic field
/// taken in the X-Y plane positioned at z=5000mm in the magnetic field region with index 0.
/// The step/bin size is 50mm, and the plotting style is "SURF3Z".
///
/// \code
///    field->DrawHistogram("XY","X",0,50.0,"SURF3Z",5000.0)
/// \endcode
/// where `field` is a pointer to the TRestAxionMagneticField object that describes the
/// magnetic field.
///
TCanvas* TRestAxionMagneticField::DrawHistogram(TString projection, TString Bcomp, Int_t volIndex,
                                                Double_t step, TString style, Double_t depth) {
    Double_t step_x, step_y, step_z;
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

    if ((unsigned int)volIndex >= GetNumberOfVolumes()) {
        RESTError << volIndex << " corresponds to none volume index " << RESTendl;
        RESTError << "Total number of volumes : " << GetNumberOfVolumes() << RESTendl;
        RESTError << "Setting volIndex to the first volume" << RESTendl;
        volIndex = 0;
    }

    if (step <= 0) {
        step_x = fMeshSize[volIndex].X();
        step_y = fMeshSize[volIndex].Y();
        step_z = fMeshSize[volIndex].Z();
    } else
        step_x = step_y = step_z = step;

    MagneticFieldVolume* vol = GetMagneticVolume(volIndex);
    if (!vol) return fCanvas;

    if (!(projection == "XY" || projection == "XZ" || projection == "YZ")) {
        RESTError << "You entered : " << projection << " as a projection but you have to choose XY, XZ or YZ"
                  << RESTendl;
        return fCanvas;
    }

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

    Int_t nBinsX = (xMax - xMin) / step_x;
    Int_t nBinsY = (yMax - yMin) / step_y;
    Int_t nBinsZ = (zMax - zMin) / step_z;

    Double_t x = -1, y = -1, z = -1;
    Double_t B = 0;
    TVector3 Bvec;

    if (projection == "XY") {
        fCanvas = new TCanvas("fCanvas", "");
        fHisto = new TH2D("", "", nBinsX, xMin, xMax, nBinsY, yMin, yMax);

        if (depth < -100000.0)
            z = (zMin + zMax) / 2.0;
        else if ((depth >= zMin) && (depth <= zMax))
            z = depth;
        else
            RESTError << "You entered depth = " << depth << ", but you have to choose depth between " << zMin
                      << " and " << zMax << RESTendl;
        x = xMin;

        for (Int_t i = 0; i < nBinsX; i++) {
            y = yMin;
            for (Int_t j = 0; j < nBinsY; j++) {
                Bvec = GetMagneticField(TVector3(x, y, z), false);
                if (Bcomp == "X")
                    B = Bvec[0];
                else {
                    if (Bcomp == "Y")
                        B = Bvec[1];
                    else {
                        if (Bcomp == "Z")
                            B = Bvec[2];
                        else
                            RESTError << "You entered : " << Bcomp
                                      << " as a B component but you have to choose X, Y or Z" << RESTendl;
                    }
                }
                fHisto->Fill(x, y, B);
                y = y + step_y;
            }
            x = x + step_x;
        }

        fCanvas->cd();
        fHisto->SetBit(TH1::kNoStats);
        fHisto->GetXaxis()->SetTitle("x (mm)");
        fHisto->GetYaxis()->SetTitle("y (mm)");

        if (Bcomp == "X") {
            TString title = Form("B_{x} against x and y @ z = %.1f", z);
            fHisto->SetTitle(title);
            fCanvas->SetTitle(title);
        }
        if (Bcomp == "Y") {
            TString title = Form("B_{y} against x and y @ z = %.1f", z);
            fHisto->SetTitle(title);
            fCanvas->SetTitle(title);
        }
        if (Bcomp == "Z") {
            TString title = Form("B_{z} against x and y @ z = %.1f", z);
            fHisto->SetTitle(title);
            fCanvas->SetTitle(title);
        }

        fHisto->Draw(style);
        return fCanvas;
    }

    if (projection == "XZ") {
        TCanvas* fCanvas = new TCanvas("fCanvas", "");
        fHisto = new TH2D("", "", nBinsX, xMin, xMax, nBinsZ, zMin, zMax);

        if (depth < -100000.0)
            y = (yMin + yMax) / 2.0;
        else if ((depth >= yMin) && (depth <= yMax))
            y = depth;
        else
            RESTError << "You entered depth = " << depth << ", but you have to choose depth between " << yMin
                      << " and " << yMax << RESTendl;
        x = xMin;

        for (Int_t i = 0; i < nBinsX; i++) {
            z = zMin;
            for (Int_t j = 0; j < nBinsZ; j++) {
                Bvec = GetMagneticField(TVector3(x, y, z), false);
                if (Bcomp == "X")
                    B = Bvec[0];
                else {
                    if (Bcomp == "Y")
                        B = Bvec[1];
                    else {
                        if (Bcomp == "Z")
                            B = Bvec[2];
                        else
                            RESTError << "You entered : " << Bcomp
                                      << " as a B component but you have to choose X, Y or Z" << RESTendl;
                    }
                }
                fHisto->Fill(x, z, B);
                z = z + step_z;
            }
            x = x + step_x;
        }

        fCanvas->cd();
        fHisto->SetBit(TH1::kNoStats);
        fHisto->GetXaxis()->SetTitle("x (mm)");
        fHisto->GetYaxis()->SetTitle("z (mm)");

        if (Bcomp == "X") {
            TString title = Form("B_{x} against x and z @ y = %.1f", y);
            fHisto->SetTitle(title);
            fCanvas->SetTitle(title);
        }
        if (Bcomp == "Y") {
            TString title = Form("B_{y} against x and z @ y = %.1f", y);
            fHisto->SetTitle(title);
            fCanvas->SetTitle(title);
        }
        if (Bcomp == "Z") {
            TString title = Form("B_{z} against x and z @ y = %.1f", y);
            fHisto->SetTitle(title);
            fCanvas->SetTitle(title);
        }

        fHisto->Draw(style);
        return fCanvas;
    }

    if (projection == "YZ") {
        TCanvas* fCanvas = new TCanvas("fCanvas", "");
        fHisto = new TH2D("", "", nBinsY, yMin, yMax, nBinsZ, zMin, zMax);

        if (depth < -100000.0)
            x = (xMin + xMax) / 2.0;
        else if ((depth >= xMin) && (depth <= xMax))
            x = depth;
        else
            RESTError << "You entered depth = " << depth << ", but you have to choose depth between " << xMin
                      << " and " << xMax << RESTendl;
        y = yMin;

        for (Int_t i = 0; i < nBinsY; i++) {
            z = zMin;
            for (Int_t j = 0; j < nBinsZ; j++) {
                Bvec = GetMagneticField(TVector3(x, y, z), false);
                if (Bcomp == "X")
                    B = Bvec[0];
                else {
                    if (Bcomp == "Y")
                        B = Bvec[1];
                    else {
                        if (Bcomp == "Z")
                            B = Bvec[2];
                        else
                            RESTError << "You entered : " << Bcomp
                                      << " as a B component but you have to choose X, Y or Z" << RESTendl;
                    }
                }
                fHisto->Fill(y, z, B);
                z = z + step_z;
            }
            y = y + step_y;
        }

        fCanvas->cd();
        fHisto->SetBit(TH1::kNoStats);
        fHisto->GetXaxis()->SetTitle("y (mm)");
        fHisto->GetYaxis()->SetTitle("z (mm)");

        if (Bcomp == "X") {
            TString title = Form("B_{x} against y and z @ x = %.1f", x);
            fHisto->SetTitle(title);
            fCanvas->SetTitle(title);
        }
        if (Bcomp == "Y") {
            TString title = Form("B_{y} against y and z @ x = %.1f", x);
            fHisto->SetTitle(title);
            fCanvas->SetTitle(title);
        }
        if (Bcomp == "Z") {
            TString title = Form("B_{z} against y and z @ x = %.1f", x);
            fHisto->SetTitle(title);
            fCanvas->SetTitle(title);
        }

        fHisto->Draw(style);
        return fCanvas;
    }

    return fCanvas;
}

///////////////////////////////////////////////
/// \brief A method that creates a canvas where tracks traversing the magnetic volume
/// are drawm together with their corresponding field intensity profile along the Z-axis.
///
TCanvas* TRestAxionMagneticField::DrawTracks(TVector3 vanishingPoint, Int_t divisions, Int_t volId) {
    if (fCanvas != NULL) {
        delete fCanvas;
        fCanvas = NULL;
    }
    fCanvas = new TCanvas("fCanvas", "", 1600, 600);

    TPad* pad1 = new TPad("pad1", "This is pad1", 0.01, 0.02, 0.99, 0.97);
    pad1->Divide(2, 1);
    pad1->Draw();

    pad1->cd(1);

    Double_t genSizeY = fBoundMax[volId].Y() * 3.;
    Double_t genPositionZ = fPositions[volId][2] - fBoundMax[volId].Z() - 2000;
    Double_t finalPositionZ = fPositions[volId][2] + fBoundMax[volId].Z() + 2000;

    // We generate the particle at different Y-highs
    TGraph* bBox = new TGraph();
    bBox->SetPoint(0, fPositions[volId][2] - fBoundMax[volId].Z(),
                   fPositions[volId][1] - fBoundMax[volId].Y());
    bBox->SetPoint(1, fPositions[volId][2] - fBoundMax[volId].Z(),
                   fPositions[volId][1] + fBoundMax[volId].Y());
    bBox->SetPoint(2, fPositions[volId][2] + fBoundMax[volId].Z(),
                   fPositions[volId][1] + fBoundMax[volId].Y());
    bBox->SetPoint(3, fPositions[volId][2] + fBoundMax[volId].Z(),
                   fPositions[volId][1] - fBoundMax[volId].Y());
    bBox->SetPoint(4, fPositions[volId][2] - fBoundMax[volId].Z(),
                   fPositions[volId][1] - fBoundMax[volId].Y());

    RESTDebug << "Gen position : " << genPositionZ << RESTendl;

    bBox->GetXaxis()->SetLimits(genPositionZ - 500, finalPositionZ + 500);
    bBox->GetHistogram()->SetMaximum(genSizeY + 100);
    bBox->GetHistogram()->SetMinimum(-genSizeY - 100);

    bBox->GetXaxis()->SetTitle("Z [mm]");
    bBox->GetXaxis()->SetTitleSize(0.05);
    bBox->GetXaxis()->SetLabelSize(0.05);
    bBox->GetYaxis()->SetTitle("Y [mm]");
    bBox->GetYaxis()->SetTitleOffset(1.3);
    bBox->GetYaxis()->SetTitleSize(0.05);
    bBox->GetYaxis()->SetLabelSize(0.05);
    bBox->SetLineWidth(2);
    bBox->Draw("AL");

    Int_t n = 0;
    for (Double_t y = genSizeY; y >= -genSizeY; y -= fBoundMax[volId].Y()) {
        TVector3 position(0, y, genPositionZ);
        TVector3 direction = (vanishingPoint - position).Unit();

        std::vector<TVector3> trackBounds = this->GetFieldBoundaries(position, direction);
        TGraph* grB = new TGraph();
        grB->SetPoint(0, trackBounds[0].Z(), trackBounds[0].Y());
        grB->SetPoint(1, trackBounds[1].Z(), trackBounds[1].Y());

        RESTDebug << "Initial" << RESTendl;
        RESTDebug << "-------" << RESTendl;
        if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) position.Print();

        RESTDebug << RESTendl;
        RESTDebug << "Start moving along" << RESTendl;
        RESTDebug << "++++++++++++++++++" << RESTendl;

        TGraph* fieldGr = new TGraph();
        Double_t posZ = fPositions[volId][2] - fBoundMax[volId].Z() - 10;
        Double_t delta = fBoundMax[volId][2] * 2. / divisions;

        Double_t field = this->GetTransversalComponent(position, direction);
        fieldGr->SetPoint(fieldGr->GetN(), posZ, field);

        while (posZ <= fPositions[volId][2] + fBoundMax[volId].Z()) {
            TVector3 posAlongAxis = TVector3(fPositions[volId][0], fPositions[volId][1], posZ);

            position = MoveToPlane(position, direction, TVector3(0, 0, 1), posAlongAxis);
            Double_t field = this->GetTransversalComponent(position, direction);

            fieldGr->SetPoint(fieldGr->GetN(), posZ, field);

            posZ += delta;
        }

        TVector3 posAlongAxis = TVector3(fPositions[volId][0], fPositions[volId][1], posZ + 10);
        position = MoveToPlane(position, direction, TVector3(0, 0, 1), posAlongAxis);

        Double_t field2 = this->GetTransversalComponent(position, direction);
        fieldGr->SetPoint(fieldGr->GetN(), posZ, field2);

        pad1->cd(2);
        fieldGr->SetLineWidth(3);
        fieldGr->SetLineColor(38 + n);
        fieldGr->GetXaxis()->SetLimits(genPositionZ - 500, finalPositionZ + 500);
        fieldGr->GetHistogram()->SetMaximum(2.5);
        fieldGr->GetHistogram()->SetMinimum(0);
        fieldGr->GetXaxis()->SetTitle("Z [mm]");
        fieldGr->GetXaxis()->SetTitleSize(0.05);
        fieldGr->GetXaxis()->SetLabelSize(0.05);
        fieldGr->GetYaxis()->SetTitle("B [T]");
        fieldGr->GetYaxis()->SetTitleOffset(1.3);
        fieldGr->GetYaxis()->SetTitleSize(0.05);
        fieldGr->GetYaxis()->SetLabelSize(0.05);
        if (y == genSizeY)
            fieldGr->Draw("AL");
        else
            fieldGr->Draw("L");
        pad1->cd(1);

        TGraph* gr = new TGraph();
        gr->SetPoint(0, genPositionZ, y);
        position = MoveToPlane(position, direction, TVector3(0, 0, 1), TVector3(0, 0, finalPositionZ));
        gr->SetPoint(1, finalPositionZ, position.Y());

        gr->SetLineWidth(1.5);
        gr->Draw("L");
        grB->SetLineColor(38 + n);
        n++;
        grB->SetLineWidth(3);
        grB->Draw("L");
    }
    bBox->Draw("L");

    return fCanvas;
}

///////////////////////////////////////////////
/// \brief A method to help loading magnetic field data, as x,y,z,Bx,By,Bz into a magnetic volume
/// definition using its corresponding mesh.
///
/// This method will be made private since it will only be used internally.
///
void TRestAxionMagneticField::LoadMagneticFieldData(MagneticFieldVolume& mVol,
                                                    std::vector<std::vector<Float_t>> data) {
    mVol.field.resize(mVol.mesh.GetNodesX());
    for (unsigned int n = 0; n < mVol.field.size(); n++) {
        mVol.field[n].resize(mVol.mesh.GetNodesY());
        for (unsigned int m = 0; m < mVol.field[n].size(); m++)
            mVol.field[n][m].resize(mVol.mesh.GetNodesZ());
    }

    RESTDebug << "TRestAxionMagneticField::LoadMagneticFieldData. Printing first 5 data rows" << RESTendl;
    for (unsigned int n = 0; n < data.size(); n++) {
        // The magnetic field map is centered at zero.
        // But the mesh definition contains the offset position
        // We shift the data to match the mesh node network.
        Int_t nX = mVol.mesh.GetNodeX((Int_t)(data[n][0] + mVol.mesh.GetNetSizeX() / 2.), true);
        Int_t nY = mVol.mesh.GetNodeY((Int_t)(data[n][1] + mVol.mesh.GetNetSizeY() / 2.), true);
        Int_t nZ = mVol.mesh.GetNodeZ((Int_t)(data[n][2] + mVol.mesh.GetNetSizeZ() / 2.), true);

        if (n < 5) {
            RESTDebug << "X: " << data[n][0] << " Y: " << data[n][1] << " Z: " << data[n][2] << RESTendl;
            RESTDebug << "absX: " << data[n][0] + mVol.position.X()
                      << " absY: " << data[n][1] + mVol.position.Y()
                      << " absZ: " << data[n][2] + mVol.position.Z() << RESTendl;
            RESTDebug << "nX: " << nX << " nY: " << nY << " nZ: " << nZ << RESTendl;
            RESTDebug << "Bx: " << data[n][3] << " By: " << data[n][4] << " Bz: " << data[n][5] << RESTendl;
            if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Extreme) GetChar();
        }

        if (mVol.field[nX][nY][nZ] != TVector3(0.0, 0.0, 0.0)) {
            RESTWarning << "X: " << data[n][0] << " Y: " << data[n][1] << " Z: " << data[n][2] << RESTendl;
            RESTWarning << "nX: " << nX << " nY: " << nY << " nZ: " << nZ << RESTendl;
            RESTWarning << "WARNING: field[nX][nY][nZ] element not equal to initial value (0, 0, 0) !!"
                        << RESTendl;
            RESTWarning << "It has value: "
                        << "mVol.field[" << nX << "][" << nY << "][" << nZ << "] = ("
                        << mVol.field[nX][nY][nZ].X() << " , " << mVol.field[nX][nY][nZ].Y() << " , "
                        << mVol.field[nX][nY][nZ].Z() << ")" << RESTendl;
            RESTWarning << "Values to write: "
                        << "Bx: " << data[n][3] << " By: " << data[n][4] << " Bz: " << data[n][5] << RESTendl
                        << RESTendl;

            this->SetError("There was a problem assigning the field matrix!");
            if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Extreme) GetChar();
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
        RESTDebug << "Reading file : " << fFileNames[n] << RESTendl;
        RESTDebug << "Full path : " << fullPathName << RESTendl;

        if (fFileNames[n] != "none" && fullPathName == "") {
            RESTError << "TRestAxionMagneticField::LoadMagneticVolumes. File " << fFileNames[n]
                      << " not found!" << RESTendl;
            RESTError
                << "REST will look for this file at any path given by <searchPath at globals definitions"
                << RESTendl;
            exit(5);
        }

        std::vector<std::vector<Float_t>> fieldData;
        if (fFileNames[n] != "none")
            if (fullPathName.find(".dat") != string::npos) {
                RESTDebug << "Reading ASCII format" << RESTendl;
                if (!TRestTools::ReadASCIITable(fullPathName, fieldData)) {
                    RESTError << "Problem reading file : " << fullPathName << RESTendl;
                    exit(1);
                }
            } else {
                if (fullPathName.find(".bin") != string::npos) {
                    RESTDebug << "Reading binary format" << RESTendl;
                    if (!TRestTools::ReadBinaryTable(fullPathName, fieldData, 6)) {
                        RESTError << "Problem reading file : " << fullPathName << RESTendl;
                        exit(2);
                    }
                }
            }
        else if (fFileNames[n] != "none") {
            RESTError << "Filename : " << fullPathName << RESTendl;
            RESTError << "File format not recognized!" << RESTendl;
            exit(3);
        }

        if (fFileNames[n] != "none" && fieldData.size() < 2) {
            RESTError << "Field data size is no more than 2 grid points!" << RESTendl;
            RESTError << "Filename : " << fullPathName << RESTendl;
            RESTError << "Probably something went wrong loading the file" << RESTendl;
            exit(4);
        }

        Float_t xMin = -fBoundMax[n].X(), yMin = -fBoundMax[n].Y(), zMin = -fBoundMax[n].Z();
        Float_t xMax = fBoundMax[n].X(), yMax = fBoundMax[n].Y(), zMax = fBoundMax[n].Z();
        Float_t meshSizeX = fMeshSize[n].X(), meshSizeY = fMeshSize[n].Y(), meshSizeZ = fMeshSize[n].Z();

        // If a field map is defined we get the boundaries, and mesh size from the volume
        if (fieldData.size() > 0) {
            RESTDebug << "Reading max boundary values" << RESTendl;
            xMax = TRestTools::GetMaxValueFromTable(fieldData, 0);
            yMax = TRestTools::GetMaxValueFromTable(fieldData, 1);
            zMax = TRestTools::GetMaxValueFromTable(fieldData, 2);

            if (fBoundMax[n] != TVector3(0, 0, 0)) {
                if (fBoundMax[n] != TVector3(xMax, yMax, zMax)) {
                    RESTWarning << "Volume : " << n << RESTendl;
                    RESTWarning << "boundMax was defined in RML but does not match the field map boundaries!"
                                << RESTendl;
                    RESTWarning << "Max. Field map boundaries : (" << xMax << ", " << yMax << ", " << zMax
                                << ")" << RESTendl;
                }
            }

            RESTDebug << "Reading min boundary values" << RESTendl;
            xMin = TRestTools::GetMinValueFromTable(fieldData, 0);
            yMin = TRestTools::GetMinValueFromTable(fieldData, 1);
            zMin = TRestTools::GetMinValueFromTable(fieldData, 2);

            if (fBoundMax[n] != TVector3(0, 0, 0)) {
                if (-fBoundMax[n] != TVector3(xMin, yMin, zMin)) {
                    RESTWarning << "Volume : " << n << RESTendl;
                    RESTWarning << "boundMax was defined in RML but does not match the field map boundaries"
                                << RESTendl;
                    RESTWarning << "Min. Field map boundaries : (" << xMin << ", " << yMin << ", " << zMin
                                << ")" << RESTendl;
                }
            }
            fBoundMax[n] = TVector3(xMax, yMax, zMax);

            RESTDebug << "Reading mesh size" << RESTendl;
            meshSizeX = TRestTools::GetLowestIncreaseFromTable(fieldData, 0);
            meshSizeY = TRestTools::GetLowestIncreaseFromTable(fieldData, 1);
            meshSizeZ = TRestTools::GetLowestIncreaseFromTable(fieldData, 2);

            if (fMeshSize[n] != TVector3(0, 0, 0)) {
                if (fMeshSize[n] != TVector3(meshSizeX, meshSizeY, meshSizeZ)) {
                    RESTWarning << "Volume : " << n << RESTendl;
                    RESTWarning
                        << "MeshSize was defined in RML but does not match the mesh size deduced from "
                           "field map"
                        << RESTendl;
                    RESTWarning << "Mesh size : (" << meshSizeX << ", " << meshSizeY << ", " << meshSizeZ
                                << ")" << RESTendl;
                }
            }
            fMeshSize[n] = TVector3(meshSizeX, meshSizeY, meshSizeZ);
        }

        if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) {
            RESTDebug << "Reading magnetic field map" << RESTendl;
            RESTDebug << "--------------------------" << RESTendl;

            RESTDebug << "Full path : " << fullPathName << RESTendl;

            RESTDebug << "Boundaries" << RESTendl;
            RESTDebug << "xMin: " << xMin << " yMin: " << yMin << " zMin: " << zMin << RESTendl;
            RESTDebug << "xMax: " << xMax << " yMax: " << yMax << " zMax: " << zMax << RESTendl;
            RESTDebug << "Mesh size" << RESTendl;

            RESTDebug << "sX: " << meshSizeX << " sY: " << meshSizeY << " sZ: " << meshSizeZ << RESTendl;

            if (fieldData.size() > 4) {
                RESTDebug << "Printing beginning of magnetic file table : " << fieldData.size() << RESTendl;
                TRestTools::PrintTable(fieldData, 0, 5);
            } else {
                RESTDebug << "The data file contains no field map" << RESTendl;
            }
        }
        if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Extreme) GetChar();

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

        if (fieldData.size() > 0) LoadMagneticFieldData(mVolume, fieldData);

        if (fBoundMax[n] == TVector3(0, 0, 0)) {
            RESTError << "The bounding box was not defined for volume " << n << "!" << RESTendl;
            RESTError << "Please review RML configuration for TRestAxionMagneticField" << RESTendl;
            exit(22);
        } else if (fMeshSize[n] == TVector3(0, 0, 0)) {
            RESTError << "The mesh grid size was not defined for volume " << n << "!" << RESTendl;
            RESTError << "Please review RML configuration for TRestAxionMagneticField" << RESTendl;
            exit(22);
        }
        fMagneticFieldVolumes.push_back(mVolume);
    }

    if (CheckOverlaps()) {
        RESTError << "TRestAxionMagneticField::LoadMagneticVolumes. Volumes overlap!" << RESTendl;
        exit(1);
    }

    for (size_t n = 0; n < fMagneticFieldVolumes.size(); n++)
        if (fReMap != TVector3(0, 0, 0)) ReMap(n, fReMap);

    RESTDebug << "Finished loading magnetic volumes" << RESTendl;
}

///////////////////////////////////////////////
/// \brief It returns the magnetic field vector at x,y,z
///
TVector3 TRestAxionMagneticField::GetMagneticField(Double_t x, Double_t y, Double_t z) {
    return GetMagneticField(TVector3(x, y, z));
}

///////////////////////////////////////////////
/// \brief It returns the magnetic field vector at TVector3(pos) using trilinear interpolation
/// that is implemented following instructions given at
/// https://en.wikipedia.org/wiki/Trilinear_interpolation
///
/// The warning in case the evaluated point is found outside any volume might be disabled using
/// the `showWarning` argument.
///
TVector3 TRestAxionMagneticField::GetMagneticField(TVector3 pos, Bool_t showWarning) {
    Int_t id = GetVolumeIndex(pos);

    if (id < 0) {
        if (showWarning)
            RESTWarning << "TRestAxionMagneticField::GetMagneticField position (" << pos.X() << ", "
                        << pos.Y() << ", " << pos.Z() << ") is outside any volume" << RESTendl;
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
        if ((xd < -0.00001) || (xd > 1.00001))
            RESTWarning << "TRestAxionMagneticField::GetMagneticField  Error: xd NOT between 0 and 1"
                        << RESTendl;

        Double_t y0 = fMagneticFieldVolumes[id].mesh.GetY(nY);
        Double_t y1 = fMagneticFieldVolumes[id].mesh.GetY(nY_1);
        Double_t yd;
        if (y0 == y1)
            yd = 0;
        else
            yd = (pos.Y() - y0) / (y1 - y0);
        if ((yd < -0.00001) || (yd > 1.00001))
            RESTWarning << "TRestAxionMagneticField::GetMagneticField  Error: yd NOT between 0 and 1"
                        << RESTendl;

        Double_t z0 = fMagneticFieldVolumes[id].mesh.GetZ(nZ);
        Double_t z1 = fMagneticFieldVolumes[id].mesh.GetZ(nZ_1);
        Double_t zd;
        if (z0 == z1)
            zd = 0;
        else
            zd = (pos.Z() - z0) / (z1 - z0);
        if ((zd < -0.00001) || (zd > 1.00001))
            RESTWarning << "TRestAxionMagneticField::GetMagneticField  Error: zd NOT between 0 and 1"
                        << RESTendl;

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

        RESTDebug << "position = (" << pos.X() << ", " << pos.Y() << ", " << pos.Z() << ")       ";
        RESTDebug << "nX = " << nX << " nY = " << nY << " nZ = " << nZ << "     nX_1 = " << nX_1
                  << "   nY_1 = " << nY_1 << "   nZ_1 = " << nZ_1 << RESTendl << RESTendl;
        RESTDebug << "C000 = (" << C000.X() << ", " << C000.Y() << ", " << C000.Z() << ")" << RESTendl
                  << RESTendl;
        RESTDebug << "C100 = (" << C100.X() << ", " << C100.Y() << ", " << C100.Z() << ")" << RESTendl
                  << RESTendl;
        RESTDebug << "C010 = (" << C010.X() << ", " << C010.Y() << ", " << C010.Z() << ")" << RESTendl
                  << RESTendl;
        RESTDebug << "C110 = (" << C110.X() << ", " << C110.Y() << ", " << C110.Z() << ")" << RESTendl
                  << RESTendl;
        RESTDebug << "C001 = (" << C001.X() << ", " << C001.Y() << ", " << C001.Z() << ")" << RESTendl
                  << RESTendl;
        RESTDebug << "C101 = (" << C101.X() << ", " << C101.Y() << ", " << C101.Z() << ")" << RESTendl
                  << RESTendl;
        RESTDebug << "C011 = (" << C011.X() << ", " << C011.Y() << ", " << C011.Z() << ")" << RESTendl
                  << RESTendl;
        RESTDebug << "C111 = (" << C111.X() << ", " << C111.Y() << ", " << C111.Z() << ")" << RESTendl
                  << RESTendl;
        RESTDebug << " -------------------------------------------------------" << RESTendl;
        RESTDebug << "C00 = (" << C00.X() << ", " << C00.Y() << ", " << C00.Z() << ")" << RESTendl
                  << RESTendl;
        RESTDebug << "C01 = (" << C01.X() << ", " << C01.Y() << ", " << C01.Z() << ")" << RESTendl
                  << RESTendl;
        RESTDebug << "C10 = (" << C10.X() << ", " << C10.Y() << ", " << C10.Z() << ")" << RESTendl
                  << RESTendl;
        RESTDebug << "C11 = (" << C11.X() << ", " << C11.Y() << ", " << C11.Z() << ")" << RESTendl
                  << RESTendl;
        RESTDebug << " -------------------------------------------------------" << RESTendl;
        RESTDebug << "C0 = (" << C0.X() << ", " << C0.Y() << ", " << C0.Z() << ")" << RESTendl << RESTendl;
        RESTDebug << "C1 = (" << C1.X() << ", " << C1.Y() << ", " << C1.Z() << ")" << RESTendl << RESTendl;
        RESTDebug << " -------------------------------------------------------" << RESTendl;
        RESTDebug << "C = (" << C.X() << ", " << C.Y() << ", " << C.Z() << ")" << RESTendl << RESTendl;

        return C;
    }
}

///////////////////////////////////////////////
/// \brief It allows to remap the magnetic field to a larger mesh size. The new mesh
/// size granularity must be provided by argument, and each dimension  must be a factor of the
/// present mesh size.
///
void TRestAxionMagneticField::ReMap(const size_t& n, const TVector3& newMapSize) {
    if (newMapSize.X() == 0 || newMapSize.Y() == 0 || newMapSize.Z() == 0) {
        RESTError << "TRestAxionMagneticField::ReMap. The mesh granularity cannot be 0" << RESTendl;
        RESTError << "Remapping will not take effect" << RESTendl;
        return;
    }

    Double_t remainder = std::fmod(newMapSize.X(), fMeshSize[n].X()) +
                         std::fmod(newMapSize.Y(), fMeshSize[n].Y()) +
                         std::fmod(newMapSize.Z(), fMeshSize[n].Z());
    if (remainder != 0) {
        RESTError << "TRestAxionMagneticField::ReMap. The field cannot be remapped." << RESTendl;
        RESTError << "The new mesh granularity must be a multiple of the existing granularity." << RESTendl;
        RESTError << "Present mesh size : (" << fMeshSize[n].X() << ", " << fMeshSize[n].Y() << ", "
                  << fMeshSize[n].Z() << ")" << RESTendl;
        RESTError << "Requested mesh size : (" << newMapSize.X() << ", " << newMapSize.Y() << ", "
                  << newMapSize.Z() << ")" << RESTendl;
        RESTError << "Remapping will not take effect" << RESTendl;
        return;
    }

    Int_t scaleX = (Int_t)(newMapSize.X() / fMeshSize[n].X());
    Int_t scaleY = (Int_t)(newMapSize.Y() / fMeshSize[n].Y());
    Int_t scaleZ = (Int_t)(newMapSize.Z() / fMeshSize[n].Z());

    Int_t newNodesX = (fMagneticFieldVolumes[n].mesh.GetNodesX() - 1) / scaleX + 1;
    Int_t newNodesY = (fMagneticFieldVolumes[n].mesh.GetNodesY() - 1) / scaleY + 1;
    Int_t newNodesZ = (fMagneticFieldVolumes[n].mesh.GetNodesZ() - 1) / scaleZ + 1;

    for (Int_t nx = 0; nx < newNodesX; nx++)
        for (Int_t ny = 0; ny < newNodesY; ny++)
            for (Int_t nz = 0; nz < newNodesZ; nz++)
                fMagneticFieldVolumes[n].field[nx][ny][nz] =
                    fMagneticFieldVolumes[n].field[nx * scaleX][ny * scaleY][nz * scaleZ];

    fMagneticFieldVolumes[n].mesh.SetNodes(newNodesX, newNodesY, newNodesZ);
    fMagneticFieldVolumes[n].field.resize(fMagneticFieldVolumes[n].mesh.GetNodesX());
    for (unsigned int i = 0; i < fMagneticFieldVolumes[n].field.size(); i++) {
        fMagneticFieldVolumes[n].field[n].resize(fMagneticFieldVolumes[n].mesh.GetNodesY());
        for (unsigned int j = 0; j < fMagneticFieldVolumes[n].field[i].size(); j++)
            fMagneticFieldVolumes[n].field[i][j].resize(fMagneticFieldVolumes[n].mesh.GetNodesZ());
    }

    fMeshSize[n] = TVector3(fMeshSize[n].X() * scaleX, fMeshSize[n].Y() * scaleY, fMeshSize[n].Z() * scaleZ);
}

///////////////////////////////////////////////
/// \brief It returns the corresponding volume index at the given position. If not found it will return
/// -1.
///
Int_t TRestAxionMagneticField::GetVolumeIndex(TVector3 pos) {
    if (!FieldLoaded()) LoadMagneticVolumes();

    for (unsigned int n = 0; n < fMagneticFieldVolumes.size(); n++) {
        if (fMagneticFieldVolumes[n].mesh.IsInside(pos)) return n;
    }
    return -1;
}

///////////////////////////////////////////////
/// \brief It returns true if the given position is found inside a magnetic volume. False otherwise.
///
Bool_t TRestAxionMagneticField::IsInside(TVector3 pos) {
    if (GetVolumeIndex(pos) >= 0) return true;
    return false;
}

///////////////////////////////////////////////
/// \brief it returns the volume position (or center) for the given volume `id`.
///
TVector3 TRestAxionMagneticField::GetVolumeCenter(Int_t id) { return GetVolumePosition(id); }

///////////////////////////////////////////////
/// \brief it returns the volume position (or center) for the given volume `id`.
///
TVector3 TRestAxionMagneticField::GetVolumePosition(Int_t id) {
    if (id >= 0 && GetNumberOfVolumes() > (unsigned int)id)
        return fPositions[id];
    else {
        RESTWarning << "TRestAxionMagneticField::GetVolumePosition. Id : " << id << " out of range!"
                    << RESTendl;
        RESTWarning << "Number of volumes defined : " << GetNumberOfVolumes() << RESTendl;
        return TVector3(0, 0, 0);
    }
}

///////////////////////////////////////////////
/// \brief It returns the intensity of the transversal magnetic field component for the defined
/// propagation `direction` and `position` given by argument.
///
Double_t TRestAxionMagneticField::GetTransversalComponent(TVector3 position, TVector3 direction) {
    return abs(GetMagneticField(position, false).Perp(direction));
}

///////////////////////////////////////////////
/// \brief It returns a vector describing the transversal magnetic field component between `from` and `to`
/// positions given by argument.
///
/// The differential element `dl` is by default 1mm, but it can be modified through the third argument of
/// this function.
///
/// The maximum number of divisions (unlimited by default) of the output vector can be fixed by the forth
/// argument. In that case, the differential element `dl` length might be increased to fullfil such
/// condition.
///
std::vector<Double_t> TRestAxionMagneticField::GetTransversalComponentAlongPath(TVector3 from, TVector3 to,
                                                                                Double_t dl, Int_t Nmax) {
    Double_t length = (to - from).Mag();

    Double_t diff = dl;
    if (Nmax > 0) {
        if (length / dl > Nmax) {
            diff = length / Nmax;
            RESTWarning << "TRestAxionMagneticField::GetTransversalComponentAlongPath. Nmax reached!"
                        << RESTendl;
            RESTWarning << "Nmax = " << Nmax << RESTendl;
            RESTWarning << "Adjusting differential step to : " << diff << " mm" << RESTendl;
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
/// \brief It initializes the field track boundaries data members of this class using a
/// track position and direction so that these values can be used later on in parameterization.
///
void TRestAxionMagneticField::SetTrack(const TVector3& position, const TVector3& direction) {
    std::vector<TVector3> trackBounds = GetFieldBoundaries(position, direction);

    if (trackBounds.size() != 2) {
        fTrackStart = TVector3(0, 0, 0);
        fTrackDirection = TVector3(0, 0, 0);
        fTrackLength = 0;
    }

    fTrackStart = trackBounds[0];
    fTrackLength = (trackBounds[1] - trackBounds[0]).Mag() - 1;
    fTrackDirection = (trackBounds[1] - trackBounds[0]).Unit();
}

///////////////////////////////////////////////
/// \brief It will return the transversal magnetic field component evaluated at a parametric
/// distance `x` (given by argument) for the track defined inside the class. The track will
/// be defined by the data members fStartTrack and fEndTrack which should be initialized by
/// external means by invoking the method SetTrack( position, direction );
///
/// This method will be used by the integration method
///
Double_t TRestAxionMagneticField::GetTransversalComponentInParametricTrack(Double_t x) {
    if (x < 0 || x > fTrackLength) return 0;

    return GetTransversalComponent(fTrackStart + x * fTrackDirection, fTrackDirection);
}

///////////////////////////////////////////////
/// \brief It returns the average of the transversal magnetic field intensity between the 3-d coordinates
/// `from` and `to`.
///
/// The differential element `dl` defines the integration step, and it is by default 1mm, but it can be
/// modified through the third argument of this function.
///
/// The maximum number of divisions of the output vector can be fixed by the forth argument. In that case,
/// the differential element `dl` length might be increased to fullfil such condition.
///
Double_t TRestAxionMagneticField::GetTransversalFieldAverage(TVector3 from, TVector3 to, Double_t dl,
                                                             Int_t Nmax) {
    Double_t Bavg = 0.;
    std::vector<Double_t> Bt = GetTransversalComponentAlongPath(from, to, dl, Nmax);
    for (auto& b : Bt) Bavg += b;

    if (Bt.size() > 0) return Bavg / Bt.size();

    RESTError << "TRestAxionMagneticField::GetTransversalFieldAverage. Lenght is zero!" << RESTendl;
    return 0.;
}

///////////////////////////////////////////////
/// \brief It returns the transverse component of the average magnetic field vector calculated
/// along the line that connects the 3-d coordinates `from` and `to` with respect to that line
///
/// The differential element `dl` defines the step, and it is by default 10mm, but it can be
/// modified through the third argument of this function.
///
/// The maximum number of divisions (unlimited by default)  can be fixed by the forth
/// argument. In that case, the differential element `dl` length might be increased to fullfil such
/// condition.
///
TVector3 TRestAxionMagneticField::GetFieldAverageTransverseVector(TVector3 from, TVector3 to, Double_t dl,
                                                                  Int_t Nmax) {
    Double_t length = (to - from).Mag();

    Double_t diff = dl;
    if (Nmax > 0) {
        if (length / dl > Nmax) {
            diff = length / Nmax;
            RESTWarning << "TRestAxionMagneticField::GetFieldAverageTransverseVector Nmax reached!"
                        << RESTendl;
            RESTWarning << "Nmax = " << Nmax << RESTendl;
            RESTWarning << "Adjusting differential step to : " << diff << " mm" << RESTendl;
        }
    }

    TVector3 direction = (to - from).Unit();
    TVector3 Bavg = TVector3(0.0, 0.0, 0.0);
    TVector3 BTavg = TVector3(0.0, 0.0, 0.0);
    Int_t numberofpoints = 0;

    for (Double_t d = 0; d <= length; d += diff) {
        Bavg = Bavg + GetMagneticField(from + d * direction);
        numberofpoints = numberofpoints + 1;
    }

    if ((length > 0) && (numberofpoints > 0)) {
        Bavg = Bavg * (1.0 / numberofpoints);  // calculates the average magnetic field vector
        BTavg =
            Bavg - (Bavg * direction) *
                       direction;  // calculates the transverse component of the average magnetic field vector
        RESTDebug << "B average vector = (" << Bavg.x() << ", " << Bavg.y() << ", " << Bavg.z() << ")"
                  << RESTendl;
        RESTDebug << "Transverse B average vector = (" << BTavg.x() << ", " << BTavg.y() << ", " << BTavg.z()
                  << ")" << RESTendl;
        return BTavg;
    }
    RESTError << "TRestAxionMagneticField::GetTransversalFieldAverage. Lenght is zero!" << RESTendl;
    return TVector3(0.0, 0.0, 0.0);
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
    RESTDebug << "Checking overlaps" << RESTendl;
    for (unsigned int n = 0; n < GetNumberOfVolumes(); n++) {
        for (unsigned int m = 0; m < GetNumberOfVolumes(); m++) {
            if (m == n) continue;
            RESTDebug << "Volume : " << m << RESTendl;

            TVector3 b = GetMagneticVolume(m)->mesh.GetVertex(0);
            RESTDebug << "Relative bottom vertex : (" << b.X() << ", " << b.Y() << ", " << b.Z() << ")"
                      << RESTendl;
            if (GetMagneticVolume(n)->mesh.IsInsideBoundingBox(b)) return true;

            TVector3 t = GetMagneticVolume(m)->mesh.GetVertex(1);
            RESTDebug << "Relative top vertex : (" << t.X() << ", " << t.Y() << ", " << t.Z() << ")"
                      << RESTendl;
            if (GetMagneticVolume(n)->mesh.IsInsideBoundingBox(t)) return true;
        }
    }
    return false;
}

///////////////////////////////////////////////
/// \brief Finds the in/out particle trajectory boundaries for a particular magnetic region bounding box.
///
/// This method checks if the trajectory defined by the position `pos` and direction `dir` passes through
/// the magnetic field region/volume `id` given. If two such points (entry point and exit point) are
/// found, their coordinates are returned. In the example shown in Fig. 1 from
/// TRestAxionFieldPropagationProcess these points are: IN 1 and OUT 1 for the region #1 and IN2 and OUT 2
/// for the region #2.
///
/// If no intersection is found, or the particle is not moving towards the volume, the returned
/// std::vector
/// will be empty.
///
std::vector<TVector3> TRestAxionMagneticField::GetVolumeBoundaries(TVector3 pos, TVector3 dir, Int_t id) {
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
std::vector<TVector3> TRestAxionMagneticField::GetFieldBoundaries(TVector3 pos, TVector3 dir,
                                                                  Double_t precision, Int_t id) {
    std::vector<TVector3> volumeBoundaries = GetVolumeBoundaries(pos, dir, id);
    if (volumeBoundaries.size() != 2) return volumeBoundaries;

    if (IsFieldConstant(id)) return volumeBoundaries;

    MagneticFieldVolume* vol = GetMagneticVolume(id);
    if (!vol) return volumeBoundaries;

    if (precision == 0) precision = min(fMeshSize[id].X(), min(fMeshSize[id].Y(), fMeshSize[id].Z())) / 2.;

    TVector3 unit = dir.Unit();
    std::vector<TVector3> fieldBoundaries;

    TVector3 in = volumeBoundaries[0];
    while ((((volumeBoundaries[1] - in) * dir) > 0) && (GetTransversalComponent(in, dir) == 0))
        in = MoveByDistanceFast(in, unit, precision);
    if (((volumeBoundaries[1] - in) * dir) > 0)
        fieldBoundaries.push_back(in);
    else
        return fieldBoundaries;

    TVector3 out = volumeBoundaries[1];
    while ((((volumeBoundaries[0] - out) * dir) < 0) && (GetTransversalComponent(out, -dir) == 0) &&
           (((out - in) * dir) > 0))
        out = MoveByDistanceFast(out, -unit, precision);
    if ((((volumeBoundaries[0] - out) * dir) < 0) && (((out - in) * dir) > 0))
        fieldBoundaries.push_back(out);
    else
        return fieldBoundaries;

    return fieldBoundaries;
}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionMagnetic field members through a RML file
///
void TRestAxionMagneticField::InitFromConfigFile() {
    TRestMetadata::InitFromConfigFile();

    auto magVolumeDef = GetElement("addMagneticVolume");
    while (magVolumeDef) {
        string filename = GetFieldValue("fileName", magVolumeDef);
        if (filename == "Not defined")
            fFileNames.push_back("none");
        else
            fFileNames.push_back(filename);

        TVector3 position = Get3DVectorParameterWithUnits("position", magVolumeDef);
        if (position == TVector3(-1, -1, -1))
            fPositions.push_back(TVector3(0, 0, 0));
        else
            fPositions.push_back(position);

        TVector3 field = Get3DVectorParameterWithUnits("field", magVolumeDef);
        if (field == TVector3(-1, -1, -1))
            fConstantField.push_back(TVector3(0, 0, 0));
        else
            fConstantField.push_back(field);

        TVector3 boundMax = Get3DVectorParameterWithUnits("boundMax", magVolumeDef);
        if (boundMax == TVector3(-1, -1, -1))
            fBoundMax.push_back(TVector3(0, 0, 0));
        else
            fBoundMax.push_back(boundMax);

        TVector3 meshSize = Get3DVectorParameterWithUnits("meshSize", magVolumeDef);
        if (meshSize == TVector3(-1, -1, -1))
            fMeshSize.push_back(TVector3(0, 0, 0));
        else
            fMeshSize.push_back(meshSize);

        string type = GetParameter("meshType", magVolumeDef);
        if (type == "NO_SUCH_PARA")
            fMeshType.push_back("cylinder");
        else
            fMeshType.push_back(type);

        // TRestMesh will only consider the first bounding component anyway
        if (fMeshType.back() == "cylinder" && fBoundMax.back().X() != fBoundMax.back().Y()) {
            RESTWarning << "Mesh type is cylinder. But X and Y inside boundMax are not the same!" << RESTendl;
            RESTWarning << "Making second bound component Y equal to the X bound component!" << RESTendl;
            fBoundMax.back().SetY(fBoundMax.back().X());
        }

        RESTDebug << "Reading new magnetic volume" << RESTendl;
        RESTDebug << "-----" << RESTendl;
        RESTDebug << "Filename : " << filename << RESTendl;
        RESTDebug << "Position: ( " << position.X() << ", " << position.Y() << ", " << position.Z() << ") mm"
                  << RESTendl;
        RESTDebug << "Field: ( " << field.X() << ", " << field.Y() << ", " << field.Z() << ") T" << RESTendl;
        RESTDebug << "Max bounding box ( " << boundMax.X() << ", " << boundMax.Y() << ", " << boundMax.Z()
                  << ")" << RESTendl;
        RESTDebug << "Mesh size ( " << meshSize.X() << ", " << meshSize.Y() << ", " << meshSize.Z() << ")"
                  << RESTendl;
        RESTDebug << "----" << RESTendl;

        magVolumeDef = GetNextElement(magVolumeDef);
    }

    LoadMagneticVolumes();

    // TODO we should check that the volumes do not overlap
}

///////////////////////////////////////////////
/// \brief Prints on screen the information about the metadata members of TRestAxionMagneticField
///
void TRestAxionMagneticField::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    if (fReMap != TVector3(0, 0, 0)) {
        RESTMetadata << "Fields original mesh size has been remapped" << RESTendl;
        RESTMetadata << " " << RESTendl;
    }

    RESTMetadata << " - Number of magnetic volumes : " << GetNumberOfVolumes() << RESTendl;
    RESTMetadata << " ------------------------------------------------ " << RESTendl;
    for (unsigned int p = 0; p < GetNumberOfVolumes(); p++) {
        if (p > 0) RESTMetadata << " ------------------------------------------------ " << RESTendl;

        Double_t centerX = fPositions[p][0];
        Double_t centerY = fPositions[p][1];
        Double_t centerZ = fPositions[p][2];

        Double_t halfSizeX = fBoundMax[p].X();
        Double_t halfSizeY = fBoundMax[p].Y();
        Double_t halfSizeZ = fBoundMax[p].Z();

        Double_t xMin = centerX - halfSizeX;
        Double_t yMin = centerY - halfSizeY;
        Double_t zMin = centerZ - halfSizeZ;

        Double_t xMax = centerX + halfSizeX;
        Double_t yMax = centerY + halfSizeY;
        Double_t zMax = centerZ + halfSizeZ;

        RESTMetadata << "* Volume " << p << " centered at  (" << centerX << "," << centerY << "," << centerZ
                     << ") mm" << RESTendl;
        RESTMetadata << "  - Grid mesh element size.  X: " << fMeshSize[p].X() << "mm "
                     << " Y: " << fMeshSize[p].Y() << "mm "
                     << " Z: " << fMeshSize[p].Z() << "mm " << RESTendl;
        RESTMetadata << "  - Offset field [T] : (" << fConstantField[p].X() << ", " << fConstantField[p].Y()
                     << ", " << fConstantField[p].Z() << ")" << RESTendl;
        RESTMetadata << "  - File loaded : " << fFileNames[p] << RESTendl;
        RESTMetadata << " " << RESTendl;
        RESTMetadata << "  - Bounds : " << RESTendl;
        RESTMetadata << "    xmin : " << xMin << " mm , xmax : " << xMax << " mm" << RESTendl;
        RESTMetadata << "    ymin : " << yMin << " mm, ymax : " << yMax << " mm" << RESTendl;
        RESTMetadata << "    zmin : " << zMin << " mm, zmax : " << zMax << " mm" << RESTendl;
        RESTMetadata << " " << RESTendl;
    }
    RESTMetadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << RESTendl;
}
