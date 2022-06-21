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
/// TRestAxionOptics is a class that allows to describe the geometrical
/// and mirror properties of an optics focusing device.
/// This metadata class is a pure abstract class used to define common
/// data members and methods to any specific optics class.
///
/// ## Conceptual optics implementation
///
/// This class implements a method, TRestAxionOptics::PropagatePhoton, that
/// is used to generate the tracking of a photon given an initial position
/// and direction. The photon propagation is done in different phases, where
/// the photon position and direction is calculated at different optical
/// interfaces (entrance/middle/exit) as it is shown in the following
/// figure.
///
/// \htmlonly <style>div.image img[src="opticsSchema.png"]{width:950px;}</style> \endhtmlonly
///
/// ![An schematic showing the principle of use of TRestAxionOptics](opticsSchema.png)
///
/// The transport of the photon in between two interfaces is controlled
/// by two pure virtual methods TRestAxionOptics::FirstMirrorReflection and
/// TRestAxionOptics::SecondMirrorReflection, that must be implemented at
/// the inherited optics class, such as is is done at TRestAxionWolterOptics.
///
/// Other pure virtual methods need to be implemented at the inherited class
/// in order to define the limits of the optical device, such as
/// TRestAxionOptics::GetEntrancePositionZ, TRestAxionOptics::GetExitPositionZ
/// or TRestAxionOptics::GetRadialLimits.
///
/// Furthermore, this class will define common methods that can be exploited
/// by any optical device, such as TRestAxionOptics::FindFocal or
/// TRestAxionOptics::CalculateSpotSize, which are MonteCarlo based, and can be
/// optionally re-implemented at the inherited class.
///
/// ## Optics metadata description
///
/// This class defines generic optics metadata members, that can be further
/// extended at the inherited classes as needed.
///
/// We distinguish the following metadata members (which will be stored on
/// disk and that can be initalized through RML).
///
/// - **mirrorLength**: It defines a common mirror length for the optics.
///
/// - **opticsFile**: It is a text file describing a table with values
/// giving geometrical information about the mirror positioning, thickness,
/// or any other information required to describe the mirrors geometry.
/// Each column on the data file will be related to a property of the
/// optics, such as mirror radius, thickness, angle. While each row
/// will correspond to a particular mirror. The meaning for each column
/// must be implemented at the inherited class. TRestAxionOptics
/// implementation just takes care of reading the file and storing it
/// inside the `fOpticsData` member that is accessible to the inherited
/// classes.
///
/// - **mirrorProperties**: It is an instance of TRestAxionOpticsMirror that contains
/// the mirror properties, such as reflectivity as a function of angle
/// and energy.
///
/// ## Optical masks
///
/// Each of the optical interfaces (entrance/middle/optics) is associated
/// with an optical mask where we can identify regions. These regions allow
/// to determine if a photon entered by and exited from the same optical
/// cavity or allowed region. The use of 3 masks at the entrance, middle
/// and exit levels permits to create a 3-dimensional constrain of the
/// movement of the photon.
///
/// Thus, masks must be constructed by the inherited class using a
/// TRestCombinedMask definition that is instantiated at the
/// TRestAxionOptics level as fEntranceMask, fMiddleMask and fExitMask.
///
/// ## Common drawing methods
///
/// This class also defines common drawing methods such as
/// TRestAxionOptics::DrawScatterMaps, TRestAxionOptics::DrawDensityMaps
/// and TRestAxionOptics::DrawParticleTracks. That may help to visualize
/// and debug the tracking of photons by the new implemented optics device.
///
/// A pure virtual method TRestAxionOptics::DrawMirrors must be implemented
/// at each class in order to draw the mirrors position, to be visualized on
/// the TRestAxionOptics::DrawParticleTracks method.
///
/// The following image is generated using TRestAxionOptics::DrawDensityMaps
/// where we can visualize the XY projection of the photons at different
/// Z-positions.
///
/// \htmlonly <style>div.image img[src="XMM_DMaps.png"]{width:750px;}</style> \endhtmlonly
///
/// ![Hitmaps for the XMM definition, using TRestAxionWolterOptics](XMM_DMaps.png)
///
/// This image was generated using the macro `REST_Axion_XMMPlots.C` which
/// is available when entering the ROOT interface using the `restRootMacros`
/// command. It corresponds with the flux of a perfectly aligned photon flux,
/// deviation=0, for a TRestAxionWolterOptics definition of the XMM optics. The
/// RML description used, `xmm.rml,` can be found at the [axion-lib data
/// repository](https://github.com/rest-for-physics/axionlib-data/tree/master/optics).
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2022-February: First concept and implementation of TRestAxionOptics class.
///            	  Javier Galan
///
/// \class      TRestAxionOptics
/// \author     Javier Galan <javier.galan@unizar.es>
///
/// <hr>
///

#include "TRestAxionOptics.h"

#include <TAxis.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>

using namespace std;

#include "TRestPhysics.h"
using namespace REST_Physics;

ClassImp(TRestAxionOptics);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionOptics::TRestAxionOptics() : TRestMetadata() { Initialize(); }

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
/// corresponding TRestAxionOptics section inside the RML.
///
TRestAxionOptics::TRestAxionOptics(const char* cfgFileName, string name) : TRestMetadata(cfgFileName) {
    RESTDebug << "Entering TRestAxionOptics constructor( cfgFileName, name )" << RESTendl;
    RESTDebug << "File: " << cfgFileName << " Name: " << name << RESTendl;
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionOptics::~TRestAxionOptics() {
    if (fEntranceMask) delete fEntranceMask;
    if (fExitMask) delete fExitMask;
    if (fMiddleMask) delete fMiddleMask;
}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionOptics members
///
void TRestAxionOptics::Initialize() {
    SetLibraryVersion(LIBRARY_VERSION);

    if (fOpticsFile != "") {
        std::string fullPathFileName = SearchFile(fOpticsFile);

        TRestTools::ReadASCIITable(fullPathFileName, fOpticsData, 3);

        // std::cout << "Reading table" << std::endl;
        // TRestTools::PrintTable(fOpticsData, 6, 9);
    }

    if (fEntranceMask != nullptr) {
        delete fEntranceMask;
        fEntranceMask = nullptr;
    }
    fEntranceMask = new TRestCombinedMask();
    fEntranceMask->SetName("Entrance");
    fEntranceMask->SetTitle("Optics entrance mask");
    fEntranceMask->SetVerboseLevel(this->GetVerboseLevel());

    if (fExitMask != nullptr) {
        delete fExitMask;
        fExitMask = nullptr;
    }
    fExitMask = new TRestCombinedMask();
    fExitMask->SetName("Exit");
    fExitMask->SetTitle("Optics exit mask");
    fExitMask->SetVerboseLevel(this->GetVerboseLevel());

    if (fMiddleMask != nullptr) {
        delete fMiddleMask;
        fMiddleMask = nullptr;
    }
    fMiddleMask = new TRestCombinedMask();
    fMiddleMask->SetName("Middle");
    fMiddleMask->SetTitle("Optics middle mask");
    fMiddleMask->SetVerboseLevel(this->GetVerboseLevel());
}

///////////////////////////////////////////////
/// \brief It returns the total number of reflections. Considering maximum
/// 1-reflection per mirror.
///
Int_t TRestAxionOptics::GetNumberOfReflections() {
    if (fFirstInteraction && fSecondInteraction) return 2;
    if (fFirstInteraction) return 1;
    if (fSecondInteraction) return 1;
    return 0;
}

///////////////////////////////////////////////
/// \brief It moves the incoming particle at the entrance of the optics plane and
/// returns the region number where the particle is entering.
///
/// It updates the value of fEntrancePosition and fEntranceDirection.
///
/// If the region is 0, the particle has hit an opaque region and the photon will
/// be lost.
///
Int_t TRestAxionOptics::TransportToEntrance(const TVector3& pos, const TVector3& dir) {
    RESTDebug << "TRestAxionOptics::TransportToEntrance" << RESTendl;
    if (pos.Z() > GetEntrancePositionZ()) {
        RESTWarning << "TRestAxionOptics::TransportToEntrance" << RESTendl;
        RESTWarning << "The particle should be placed before the entrance!" << RESTendl;
        return 0;
    }
    if (dir.Z() <= 0) {
        RESTWarning << "TRestAxionOptics::TransportToEntrance" << RESTendl;
        RESTWarning << "Photon is not moving on the positive Z-direction!" << RESTendl;
        return 0;
    }
    fEntrancePosition =
        REST_Physics::MoveToPlane(pos, dir, TVector3(0, 0, 1), TVector3(0, 0, GetEntrancePositionZ()));

    return fEntranceMask->GetRegion(fEntrancePosition.X(), fEntrancePosition.Y());
}

///////////////////////////////////////////////
/// \brief It moves the incoming particle to the middle of the optics plane and
/// returns the region number where the particle is entering.
///
/// It updates the value of fMiddlePosition and fMiddleDirection.
///
/// If the region is 0, the particle has hit an opaque region and the photon will
/// be lost.
///
Int_t TRestAxionOptics::TransportToMiddle(const TVector3& pos, const TVector3& dir) {
    RESTDebug << "TRestAxionOptics::TransportToMiddle" << RESTendl;
    if (pos.Z() > 0 || pos.Z() < GetEntrancePositionZ()) {
        RESTWarning << "TRestAxionOptics::TransportToMiddle" << RESTendl;
        RESTWarning << "The particle should be placed between entrance and middle!" << RESTendl;
        return 0;
    }
    if (dir.Z() <= 0) {
        RESTWarning << "TRestAxionOptics::TransportToMiddle" << RESTendl;
        RESTWarning << "Photon is not moving on the positive Z-direction!" << RESTendl;
        RESTWarning << "Direction : ( " << dir.X() << ", " << dir.Y() << ", " << dir.Z() << ")" << RESTendl;
        return 0;
    }

    fMiddlePosition = REST_Physics::MoveToPlane(pos, dir, TVector3(0, 0, 1), TVector3(0, 0, 0));
    fMiddleDirection = dir;

    return fMiddleMask->GetRegion(fMiddlePosition.X(), fMiddlePosition.Y());
}

///////////////////////////////////////////////
/// \brief It moves the incoming particle to the exit of the optics plane and
/// returns the region number where the particle is entering.
///
/// It updates the value of fExitPosition and fExitDirection.
///
/// If the region is 0, the particle has hit an opaque region and the photon will
/// be lost.
///
Int_t TRestAxionOptics::TransportToExit(const TVector3& pos, const TVector3& dir) {
    RESTDebug << "TRestAxionOptics::TransportToExit" << RESTendl;
    if (pos.Z() < 0 || pos.Z() > GetExitPositionZ()) {
        RESTWarning << "TRestAxionOptics::TransportToExit" << RESTendl;
        RESTWarning << "The particle should be placed between middle and exit!" << RESTendl;
        return 0;
    }
    if (dir.Z() <= 0) {
        RESTWarning << "TRestAxionOptics::TransportToExit" << RESTendl;
        RESTWarning << "Photon is not moving on the positive Z-direction!" << RESTendl;
        return 0;
    }

    fExitPosition =
        REST_Physics::MoveToPlane(pos, dir, TVector3(0, 0, 1), TVector3(0, 0, GetExitPositionZ()));
    fExitDirection = dir;

    return fExitMask->GetRegion(fExitPosition.X(), fExitPosition.Y());
}

///////////////////////////////////////////////
/// \brief It reinitializes particle positions and directions at the different optical regions
///
void TRestAxionOptics::ResetPositions() {
    fOriginPosition = TVector3(0, 0, 0);

    fEntrancePosition = TVector3(0, 0, 0);
    fEntranceDirection = TVector3(0, 0, 0);

    fFirstInteractionPosition = TVector3(0, 0, 0);

    fMiddlePosition = TVector3(0, 0, 0);
    fMiddleDirection = TVector3(0, 0, 0);

    fSecondInteractionPosition = TVector3(0, 0, 0);

    fExitPosition = TVector3(0, 0, 0);
    fExitDirection = TVector3(0, 0, 0);

    fCurrentMirror = -1;
}

///////////////////////////////////////////////
/// \brief Propagating photon
///
Double_t TRestAxionOptics::PropagatePhoton(const TVector3& pos, const TVector3& dir, Double_t energy) {
    RESTDebug << " --> Entering TRestAxionOptics::PropagatePhoton" << RESTendl;
    Double_t reflectivity = 1;

    RESTDebug << "Reseting positions" << RESTendl;
    ResetPositions();

    fOriginPosition = pos;
    fEntranceDirection = dir;

    /// We move the particle to the entrance optics plane. We update fEntrancePosition
    Int_t entranceRegion = TransportToEntrance(fOriginPosition, fEntranceDirection);
    RESTDebug << "Entrance region : " << entranceRegion << RESTendl;

    if (entranceRegion == 0) return 0.0;

    /// Now that we are placed at the optics entrance plane.
    /// We define the current active mirror (same array index for front and back)
    SetMirror();

    /// We update the position and direction at the first mirror. We update
    /// fFirstInteractionPosition and fMiddleDirection
    FirstMirrorReflection(fEntrancePosition, fEntranceDirection);

    /// We move the particle to the entrance optics plane. We update fEntrancePosition
    Int_t middleRegion = TransportToMiddle(fFirstInteractionPosition, fMiddleDirection);
    RESTDebug << "Middle region : " << middleRegion << RESTendl;
    if (middleRegion != entranceRegion) return 0.0;

    /// We update the position and direction at the second mirror. We update
    /// fSecondInteractionPosition and fExitDirection
    SecondMirrorReflection(fMiddlePosition, fMiddleDirection);

    //// TODO obtain incidence angle and reflectivity

    Int_t exitRegion = TransportToExit(fSecondInteractionPosition, fExitDirection);
    RESTDebug << "Exit region : " << exitRegion << RESTendl;
    if (exitRegion != middleRegion) return 0.0;

    return GetPhotonReflectivity(energy);
}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionOptics field members through a RML file
///
void TRestAxionOptics::InitFromConfigFile() {
    if (fMirrorProperties) {
        delete fMirrorProperties;
        fMirrorProperties = nullptr;
    }

    fMirrorProperties = (TRestAxionOpticsMirror*)this->InstantiateChildMetadata("TRestAxionOpticsMirror");

    TRestMetadata::InitFromConfigFile();

    // If we recover the metadata class from a ROOT file we will need to call Initialize ourselves
    this->Initialize();
}

///////////////////////////////////////////////
/// \brief Prints on screen the information about the metadata members of TRestAxionOptics
///
void TRestAxionOptics::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    RESTMetadata << " - Optics file : " << fOpticsFile << RESTendl;
    RESTMetadata << "---------" << RESTendl;
    RESTMetadata << "Entrance position in Z : " << GetEntrancePositionZ() << " mm" << RESTendl;
    RESTMetadata << "Exit position in Z : " << GetExitPositionZ() << " mm" << RESTendl;
    RESTMetadata << "---------" << RESTendl;
    RESTMetadata << " " << RESTendl;
    RESTMetadata << " Use \"this->PrintMasks()\" to get masks info" << RESTendl;
    RESTMetadata << " Use \"this->PrintMirror()\" to get mirror info" << RESTendl;
    RESTMetadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << RESTendl;
}

///////////////////////////////////////////////
/// \brief Prints on screen the 3-optical masks used on the optics planes
///
void TRestAxionOptics::PrintMasks() {
    if (fEntranceMask) fEntranceMask->PrintMetadata();
    if (fMiddleMask) fMiddleMask->PrintMetadata();
    if (fExitMask) fExitMask->PrintMetadata();
}

///////////////////////////////////////////////
/// \brief Prints on screen the 3-optical masks used on the optics planes
///
void TRestAxionOptics::PrintMirror() {
    if (fMirrorProperties) fMirrorProperties->PrintMetadata();
}

///////////////////////////////////////////////
/// \brief Prints on screen the mask used on the entrance optics plane
///
void TRestAxionOptics::PrintEntranceMask() {
    if (fEntranceMask)
        fEntranceMask->PrintMetadata();
    else
        RESTWarning << "TRestAxionOptics::PrintEntranceMask. Not available" << RESTendl;
}

///////////////////////////////////////////////
/// \brief Prints on screen the mask used on the middle optics plane
///
void TRestAxionOptics::PrintMiddleMask() {
    if (fMiddleMask)
        fMiddleMask->PrintMetadata();
    else
        RESTWarning << "TRestAxionOptics::PrintMiddleMask. Not available" << RESTendl;
}

///////////////////////////////////////////////
/// \brief Prints on screen the mask used on the exit optics plane
///
void TRestAxionOptics::PrintExitMask() {
    if (fExitMask)
        fExitMask->PrintMetadata();
    else
        RESTWarning << "TRestAxionOptics::PrintExitMask. Not available" << RESTendl;
}

///////////////////////////////////////////////
/// \brief Prints the positions taken by the photon after ray-tracing
/// that should have been updated using the method PropagatePhoton
///
void TRestAxionOptics::PrintPhotonTrackingSummary() {
    std::cout << "Photon tracking summary" << std::endl;
    std::cout << "-----------------------" << std::endl;
    std::cout << "Origin : ( " << fOriginPosition.X() << ", " << fOriginPosition.Y() << ", "
              << fOriginPosition.Z() << ")" << std::endl;

    std::cout << "Entrance Direction : ( " << fEntranceDirection.X() << ", " << fEntranceDirection.Y() << ", "
              << fEntranceDirection.Z() << ")" << std::endl;

    std::cout << "Entrance : ( " << fEntrancePosition.X() << ", " << fEntrancePosition.Y() << ", "
              << fEntrancePosition.Z() << ")" << std::endl;

    std::cout << "Entrance radius : "
              << TMath::Sqrt(fEntrancePosition.X() * fEntrancePosition.X() +
                             fEntrancePosition.Y() * fEntrancePosition.Y())
              << std::endl;

    std::cout << "FirstInteraction : ( " << fFirstInteractionPosition.X() << ", "
              << fFirstInteractionPosition.Y() << ", " << fFirstInteractionPosition.Z() << ")" << std::endl;

    std::cout << "Middle Direction : ( " << fMiddleDirection.X() << ", " << fMiddleDirection.Y() << ", "
              << fMiddleDirection.Z() << ")" << std::endl;

    std::cout << "Middle : ( " << fMiddlePosition.X() << ", " << fMiddlePosition.Y() << ", "
              << fMiddlePosition.Z() << ")" << std::endl;

    std::cout << "Middle radius : "
              << TMath::Sqrt(fMiddlePosition.X() * fMiddlePosition.X() +
                             fMiddlePosition.Y() * fMiddlePosition.Y())
              << std::endl;

    std::cout << "SecondInteraction : ( " << fSecondInteractionPosition.X() << ", "
              << fSecondInteractionPosition.Y() << ", " << fSecondInteractionPosition.Z() << ")" << std::endl;

    std::cout << "Exit Direction : ( " << fExitDirection.X() << ", " << fExitDirection.Y() << ", "
              << fExitDirection.Z() << ")" << std::endl;

    std::cout << "Exit : ( " << fExitPosition.X() << ", " << fExitPosition.Y() << ", " << fExitPosition.Z()
              << ")" << std::endl;

    std::cout << "Exit radius : "
              << TMath::Sqrt(fExitPosition.X() * fExitPosition.X() + fExitPosition.Y() * fExitPosition.Y())
              << std::endl;
    std::cout << "-----------------------" << std::endl;
}

///////////////////////////////////////////////
/// \brief A prototype method to be implemented by specific optics to draw an schematic
/// including the mirrors geometry.
///
TPad* TRestAxionOptics::CreatePad(Int_t nx, Int_t ny) {
    if (fPad != nullptr) {
        delete fPad;
        fPad = nullptr;
    }

    fPad = new TPad("optics_pad", "This is the optics drawing pad", 0.01, 0.02, 0.99, 0.97);
    if (nx > 1 || ny > 1) fPad->Divide(nx, ny);
    fPad->Draw();
    fPad->SetRightMargin(0.09);
    fPad->SetLeftMargin(0.2);
    fPad->SetBottomMargin(0.15);

    return fPad;
}

///////////////////////////////////////////////
/// \brief A prototype method to be implemented by specific optics to draw an schematic
/// including the mirrors geometry.
///
Double_t TRestAxionOptics::GetPhotonReflectivity(Double_t energy) {
    if (!fMirrorProperties) return 0;

    Double_t reflectivity = 1.;
    if (IsFirstMirrorReflection()) {
        Double_t angle = REST_Physics::GetVectorsAngle(fEntranceDirection, fMiddleDirection);

        angle = angle * units("deg") / units("rad") / 2.;  // Incidence angle in degrees. We divide by 2

        reflectivity *= fMirrorProperties->GetReflectivity(angle, energy);
    }
    if (IsSecondMirrorReflection()) {
        Double_t angle = REST_Physics::GetVectorsAngle(fMiddleDirection, fExitDirection);

        angle = angle * units("deg") / units("rad") / 2.;  // Incidence angle in degrees. We divide by 2

        reflectivity *= fMirrorProperties->GetReflectivity(angle, energy);
    }
    return reflectivity;
}

///////////////////////////////////////////////
/// \brief A method to draw an optics schematic including the mirrors geometry, and few photon
/// tracks. This method is intended for debugging the photon tracking implementation.
///
/// Two optional parameters are allowed:
/// - **deviation**: It controls the maximum random X and Y direction. 0 by default.
/// - **particles**: The number of particles to be launched. 10 by default.
///
TPad* TRestAxionOptics::DrawParticleTracks(Double_t deviation, Int_t particles) {
    DrawMirrors();

    for (unsigned int n = 0; n < particles; n++) {
        Double_t r = fRandom->Uniform(GetRadialLimits().first, GetRadialLimits().second);
        Double_t angle = fRandom->Uniform(0, 2 * TMath::Pi());
        TVector3 origin(r * TMath::Cos(angle), r * TMath::Sin(angle), -3 * fMirrorLength);
        TVector3 direction(fRandom->Uniform(-deviation, deviation), fRandom->Uniform(-deviation, deviation),
                           1);
        direction = direction.Unit();

        Double_t reflectivity = PropagatePhoton(origin, direction, 1);

        PrintPhotonTrackingSummary();

        TGraph* gr = new TGraph();

        Double_t rO = TMath::Sqrt(origin.X() * origin.X() + origin.Y() * origin.Y());
        gr->SetPoint(gr->GetN(), origin.Z(), rO);

        ////

        Double_t rEntrance = TMath::Sqrt(fEntrancePosition.X() * fEntrancePosition.X() +
                                         fEntrancePosition.Y() * fEntrancePosition.Y());
        if (rEntrance > 0) gr->SetPoint(gr->GetN(), fEntrancePosition.Z(), rEntrance);

        ////

        Double_t rFirst = TMath::Sqrt(fFirstInteractionPosition.X() * fFirstInteractionPosition.X() +
                                      fFirstInteractionPosition.Y() * fFirstInteractionPosition.Y());

        if (rFirst > 0) gr->SetPoint(gr->GetN(), fFirstInteractionPosition.Z(), rFirst);

        ////

        Double_t rMiddle = TMath::Sqrt(fMiddlePosition.X() * fMiddlePosition.X() +
                                       fMiddlePosition.Y() * fMiddlePosition.Y());
        if (rMiddle > 0) gr->SetPoint(gr->GetN(), fMiddlePosition.Z(), rMiddle);

        ////

        Double_t rSecond = TMath::Sqrt(fSecondInteractionPosition.X() * fSecondInteractionPosition.X() +
                                       fSecondInteractionPosition.Y() * fSecondInteractionPosition.Y());

        if (rSecond > 0) gr->SetPoint(gr->GetN(), fSecondInteractionPosition.Z(), rSecond);

        ////

        Double_t rExit =
            TMath::Sqrt(fExitPosition.X() * fExitPosition.X() + fExitPosition.Y() * fExitPosition.Y());
        if (rExit > 0) gr->SetPoint(gr->GetN(), fExitPosition.Z(), rExit);

        ////

        if (reflectivity > 0) {
            TVector3 end = REST_Physics::MoveToPlane(fExitPosition, fExitDirection, TVector3(0, 0, 1),
                                                     TVector3(0, 0, 3 * fMirrorLength));
            Double_t rEnd = TMath::Sqrt(end.X() * end.X() + end.Y() * end.Y());
            if (rEnd > 0) gr->SetPoint(gr->GetN(), end.Z(), rEnd);
        }

        ////

        gr->SetLineWidth(1);
        if (GetMirror() < 0)
            gr->SetLineColor(kBlack);
        else
            gr->SetLineColor(20 + GetMirror() % 20);

        gr->Draw("L");
    }
    return fPad;
}

///////////////////////////////////////////////
/// \brief It implements a generic method to identify the optimum focal point. It can
/// be reimplemented at each specific optics class.
///
/// It is a statistical method, thus the result contains an implicit statistical error.
///
/// It receives 4 arguments.
/// - **from** and **to**: They define the range where the focal will be searched for.
/// - **precision**: It defines the accuracy required
/// - **recalculate**: If `false` it will reuse a previous focal point calculation. If `true`
/// it will force to recalculate the focal point.
/// - **particles**: Number of particles to be launched in order to calculate the spot size.
///
Double_t TRestAxionOptics::FindFocal(Double_t from, Double_t to, Double_t energy, Double_t precision,
                                     Bool_t recalculate, Int_t particles) {
    RESTDebug << "Entering TRestAxionOptics::FindFocal" << RESTendl;

    if (fFocal > 0 && recalculate == false) return fFocal;

    Double_t focal = (from + to) / 2.;
    Double_t spotSize = CalculateSpotSize(energy, focal);

    Double_t step = (to - from) / 10.;

    if (step < 0) {
        step = from;
        from = to;
        to = step;
        step = (to - from) / 10.;
    }

    RESTInfo << "Calculating focal : " << focal << RESTendl;
    for (Double_t f = from; f < to; f += step) {
        Double_t size = CalculateSpotSize(energy, f);
        if (size < spotSize) {
            spotSize = size;
            focal = f;
        }
    }

    if (step > precision) return FindFocal(focal - step, focal + step, energy, precision, recalculate);

    fFocal = focal;
    return fFocal;
}

///////////////////////////////////////////////
/// \brief It measures the spot size through Monte Carlo at a given plane given by z.
/// If z=0 this method will check for the spot size at the focal point, which is the
/// default behaviour.
///
/// The spot is assumed to be centered at (0,0). We are describing a perfectly aligned
/// optics device.
///
/// The Monte Carlo generated photons are also considered to be perfectly aligned
/// with the z-axis.
///
Double_t TRestAxionOptics::CalculateSpotSize(Double_t energy, Double_t z, Int_t particles) {
    Double_t sum = 0;
    Int_t nSum = 0;
    Double_t deviation = 0;
    for (unsigned int n = 0; n < particles; n++) {
        Double_t reflectivity = PropagateMonteCarloPhoton(energy, deviation);

        if (reflectivity == 0) continue;

        TVector3 posZ =
            REST_Physics::MoveToPlane(fExitPosition, fExitDirection, TVector3(0, 0, 1), TVector3(0, 0, z));
        Double_t x = posZ.X();
        Double_t y = posZ.Y();

        sum += reflectivity * x * x + y * y;
        nSum++;
    }

    if (nSum > 0) sum = TMath::Sqrt(sum / nSum);

    return sum;
}

///////////////////////////////////////////////
/// \brief It will produce a MonteCarlo photon spatially distributed in XY as defined by
/// the GetRadialLimits method (extended by 50%), and with direction along the Z-axis
/// with a maximum deviation angle fixed by the `deviation` input parameter. If
/// `deviation=0` the photons will always be parallel to the z-axis. The photons will
/// be launched from z=-3*fMirrorLength.
///
Double_t TRestAxionOptics::PropagateMonteCarloPhoton(Double_t energy, Double_t deviation) {
    Double_t x = fRandom->Uniform(0, 1.5 * GetRadialLimits().second);
    Double_t y = fRandom->Uniform(0, 1.5 * GetRadialLimits().second);
    Double_t r2 = x * x + y * y;
    while (r2 > 1.5 * GetRadialLimits().second || r2 < 0.5 * GetRadialLimits().second) {
        x = fRandom->Uniform(0, 1.5 * GetRadialLimits().second);
        y = fRandom->Uniform(0, 1.5 * GetRadialLimits().second);
        r2 = x * x + y * y;
    }

    Double_t r = fRandom->Uniform(0.5 * GetRadialLimits().first, 1.5 * GetRadialLimits().second);
    Double_t angle = fRandom->Uniform(0, 2 * TMath::Pi());
    TVector3 origin(r * TMath::Cos(angle), r * TMath::Sin(angle), -3 * fMirrorLength);

    Double_t theta = fRandom->Uniform(0, deviation);
    Double_t phi = fRandom->Uniform(0, 2 * TMath::Pi());

    TVector3 direction(0, 0, 1);
    direction.SetTheta(theta);
    direction.SetPhi(phi);
    direction.SetMag(1.);

    return PropagatePhoton(origin, direction, energy);
}

///////////////////////////////////////////////
/// \brief It implements a generic method to identify the optimum focal point. It can
/// be reimplemented at each specific optics class.
///
/// Focal position will be searched around focalHint within half a meter range.
///
TPad* TRestAxionOptics::DrawScatterMaps(Double_t z, Double_t energy, Double_t deviation, Int_t particles,
                                        Double_t focalHint) {
    TRestAxionOptics::CreatePad(2, 2);

    fPad->cd();

    TGraph* grEntrance = new TGraph();
    TGraph* grExit = new TGraph();
    TGraph* grZ = new TGraph();
    TGraph* grFocal = new TGraph();

    Double_t focal = FindFocal(focalHint - 500, focalHint + 500, energy, 1);

    for (unsigned int n = 0; n < particles; n++) {
        Double_t reflectivity = PropagateMonteCarloPhoton(energy, deviation);

        if (fFirstInteractionPosition.Z() == 0) continue;  // The photon hits the entrance mask
        grEntrance->SetPoint(grEntrance->GetN(), fEntrancePosition.X(), fEntrancePosition.Y());

        if (reflectivity == 0) continue;  // The photon hits any mask
        grExit->SetPoint(grExit->GetN(), fExitPosition.X(), fExitPosition.Y());

        TVector3 posZ =
            REST_Physics::MoveToPlane(fExitPosition, fExitDirection, TVector3(0, 0, 1), TVector3(0, 0, z));
        grZ->SetPoint(grZ->GetN(), posZ.X(), posZ.Y());

        TVector3 posFocal = REST_Physics::MoveToPlane(fExitPosition, fExitDirection, TVector3(0, 0, 1),
                                                      TVector3(0, 0, focal));
        grFocal->SetPoint(grFocal->GetN(), posFocal.X(), posFocal.Y());
    }

    fPad->cd(1);

    grEntrance->GetXaxis()->SetLimits(-GetRadialLimits().second * 1.15, GetRadialLimits().second * 1.15);
    grEntrance->GetHistogram()->SetMaximum(GetRadialLimits().second * 1.15);
    grEntrance->GetHistogram()->SetMinimum(-GetRadialLimits().second * 1.15);

    grEntrance->SetTitle("Entrance plane");
    grEntrance->GetXaxis()->SetTitle("X [mm]");
    grEntrance->GetXaxis()->SetTitleSize(0.04);
    grEntrance->GetXaxis()->SetLabelSize(0.04);
    grEntrance->GetXaxis()->SetNdivisions(5);
    grEntrance->GetYaxis()->SetTitle("Y [mm]");
    grEntrance->GetYaxis()->SetTitleOffset(1.4);
    grEntrance->GetYaxis()->SetTitleSize(0.04);
    grEntrance->GetYaxis()->SetLabelSize(0.04);

    grEntrance->SetMarkerStyle(1);
    grEntrance->Draw("AP");

    fPad->cd(2);

    grExit->SetTitle("Exit plane");
    grExit->GetXaxis()->SetTitle("X [mm]");
    grExit->GetXaxis()->SetTitleSize(0.04);
    grExit->GetXaxis()->SetLabelSize(0.04);
    grExit->GetXaxis()->SetNdivisions(5);
    grExit->GetYaxis()->SetTitle("Y [mm]");
    grExit->GetYaxis()->SetTitleOffset(1.4);
    grExit->GetYaxis()->SetTitleSize(0.04);
    grExit->GetYaxis()->SetLabelSize(0.04);

    grExit->SetMarkerStyle(1);
    grExit->Draw("AP");

    fPad->cd(3);

    Double_t xMin = TMath::MinElement(grZ->GetN(), grZ->GetX());
    Double_t xMax = TMath::MaxElement(grZ->GetN(), grZ->GetX());
    Double_t yMin = TMath::MinElement(grZ->GetN(), grZ->GetY());
    Double_t yMax = TMath::MaxElement(grZ->GetN(), grZ->GetY());

    std::string zTitle = "User plane. Z = " + DoubleToString(z) + "mm";
    grZ->SetTitle(zTitle.c_str());
    grZ->GetXaxis()->SetLimits(1.5 * xMin, 1.5 * xMax);
    grZ->GetHistogram()->SetMaximum(1.5 * yMax);
    grZ->GetHistogram()->SetMinimum(1.5 * yMin);

    grZ->GetXaxis()->SetTitle("X [mm]");
    grZ->GetXaxis()->SetTitleSize(0.04);
    grZ->GetXaxis()->SetLabelSize(0.04);
    grZ->GetXaxis()->SetNdivisions(5);
    grZ->GetYaxis()->SetTitle("Y [mm]");
    grZ->GetYaxis()->SetTitleOffset(1.4);
    grZ->GetYaxis()->SetTitleSize(0.04);
    grZ->GetYaxis()->SetLabelSize(0.04);

    grZ->SetMarkerStyle(1);
    grZ->Draw("AP");

    fPad->cd(4);

    std::string focalTitle = "Focal plane. Z = " + DoubleToString(focal) + "mm";
    grFocal->SetTitle(focalTitle.c_str());
    grFocal->GetXaxis()->SetLimits(-10, 10);
    grFocal->GetHistogram()->SetMaximum(10);
    grFocal->GetHistogram()->SetMinimum(-10);

    grFocal->GetXaxis()->SetTitle("X [mm]");
    grFocal->GetXaxis()->SetTitleSize(0.04);
    grFocal->GetXaxis()->SetLabelSize(0.04);
    grFocal->GetXaxis()->SetNdivisions(5);
    grFocal->GetYaxis()->SetTitle("Y [mm]");
    grFocal->GetYaxis()->SetTitleOffset(1.4);
    grFocal->GetYaxis()->SetTitleSize(0.04);
    grFocal->GetYaxis()->SetLabelSize(0.04);

    grFocal->SetMarkerStyle(1);
    grFocal->Draw("AP");

    return fPad;
}

///////////////////////////////////////////////
/// \brief It implements a generic method to identify the optimum focal point. It can
/// be reimplemented at each specific optics class.
///
/// Focal position will be searched around focalHint within half a meter range.
///
TPad* TRestAxionOptics::DrawDensityMaps(Double_t z, Double_t energy, Double_t deviation, Int_t particles,
                                        Double_t focalHint) {
    Double_t focal = FindFocal(focalHint - 500, focalHint + 500, energy, 1);

    TRestAxionOptics::CreatePad(2, 2);

    fPad->cd();

    Double_t lowL = -1.15 * GetRadialLimits().second;
    Double_t highL = 1.15 * GetRadialLimits().second;
    TH2F* hEntrance = new TH2F("entranceH", "Entrance plane", 500, lowL, highL, 500, lowL, highL);
    TH2F* hExit = new TH2F("exitH", "Exit plane", 500, lowL, highL, 500, lowL, highL);

    std::string zTitle = "User plane. Z = " + DoubleToString(z) + "mm";
    TH2F* hZ = new TH2F("zH", zTitle.c_str(), 500, lowL, highL, 500, lowL, highL);

    std::string focalTitle = "Focal plane. Z = " + DoubleToString(focal) + "mm";
    TH2F* hFocal = new TH2F("focalH", focalTitle.c_str(), 500, -10, 10, 500, -10, 10);

    for (unsigned int n = 0; n < particles; n++) {
        Double_t reflectivity = PropagateMonteCarloPhoton(energy, deviation);

        if (fFirstInteractionPosition.Z() == 0) continue;  // The photon hits the entrance mask
        hEntrance->Fill(fEntrancePosition.X(), fEntrancePosition.Y());

        if (reflectivity == 0) continue;  // The photon hits any mask
        hExit->Fill(fExitPosition.X(), fExitPosition.Y());

        TVector3 posZ =
            REST_Physics::MoveToPlane(fExitPosition, fExitDirection, TVector3(0, 0, 1), TVector3(0, 0, z));
        hZ->Fill(posZ.X(), posZ.Y());

        TVector3 posFocal = REST_Physics::MoveToPlane(fExitPosition, fExitDirection, TVector3(0, 0, 1),
                                                      TVector3(0, 0, focal));
        hFocal->Fill(posFocal.X(), posFocal.Y());
    }

    fPad->cd(1);

    hEntrance->GetXaxis()->SetTitle("X [mm]");
    hEntrance->GetXaxis()->SetTitleSize(0.04);
    hEntrance->GetXaxis()->SetLabelSize(0.04);
    hEntrance->GetXaxis()->SetNdivisions(5);
    hEntrance->GetYaxis()->SetTitle("Y [mm]");
    hEntrance->GetYaxis()->SetTitleOffset(1.4);
    hEntrance->GetYaxis()->SetTitleSize(0.04);
    hEntrance->GetYaxis()->SetLabelSize(0.04);

    hEntrance->Draw("colz");

    fPad->cd(2);

    hExit->GetXaxis()->SetTitle("X [mm]");
    hExit->GetXaxis()->SetTitleSize(0.04);
    hExit->GetXaxis()->SetLabelSize(0.04);
    hExit->GetXaxis()->SetNdivisions(5);
    hExit->GetYaxis()->SetTitle("Y [mm]");
    hExit->GetYaxis()->SetTitleOffset(1.4);
    hExit->GetYaxis()->SetTitleSize(0.04);
    hExit->GetYaxis()->SetLabelSize(0.04);

    hExit->Draw("colz");

    fPad->cd(3);

    hZ->GetXaxis()->SetTitle("X [mm]");
    hZ->GetXaxis()->SetTitleSize(0.04);
    hZ->GetXaxis()->SetLabelSize(0.04);
    hZ->GetXaxis()->SetNdivisions(5);
    hZ->GetYaxis()->SetTitle("Y [mm]");
    hZ->GetYaxis()->SetTitleOffset(1.4);
    hZ->GetYaxis()->SetTitleSize(0.04);
    hZ->GetYaxis()->SetLabelSize(0.04);

    hZ->Draw("colz");

    fPad->cd(4);

    hFocal->GetXaxis()->SetTitle("X [mm]");
    hFocal->GetXaxis()->SetTitleSize(0.04);
    hFocal->GetXaxis()->SetLabelSize(0.04);
    hFocal->GetXaxis()->SetNdivisions(5);
    hFocal->GetYaxis()->SetTitle("Y [mm]");
    hFocal->GetYaxis()->SetTitleOffset(1.4);
    hFocal->GetYaxis()->SetTitleSize(0.04);
    hFocal->GetYaxis()->SetLabelSize(0.04);

    hFocal->Draw("colz");

    return fPad;
}
