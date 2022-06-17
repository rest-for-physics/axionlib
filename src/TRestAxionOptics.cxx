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
/// TRestAxionOptics is a class that allows to load externally
/// defined optics response files. This metadata class will be a generic,
/// abstract, class that will be inherited by other more specific metadata
/// classes. This class will define few common metadata members helping to
/// describe the optics alignment, position, and basic geometry specifications,
/// such as number of mirror rings, or additional entrance masks, such as
/// spider mask.
///
/// The derived metadata classes, such as TRestAxionGenericOptics or
/// TRestAxionMCPLOptics must implement the following virtual methods
/// TRestAxionOptics::GetPositionAtExit, TRestAxionOptics::GetDirectionAtExit
/// and TRestAxionOptics::GetEfficiency.
///
/// The following metadata parameters define the optics position, size and
/// alignment:
/// * **center**: It defines the center of the optics, entrance and exit
/// optics planes will be defined using the half lenght and the `center`
/// position.
/// * **axis**: It defines the optical axis direction.
/// * **length**: It defines the size of the optics, used to calculate
/// the optics plane entrance, and the optics plane exit.
///
/// A relevant parameter are the radii that define the mirror rings of
/// the optics. In practice we define the iner and outer radius of a ring
/// or corona (which is the space between two rings or mirrors). Thus, the
/// inner and outer radius of each ring defines the region where photons are
/// capable to go through. Photons hitting other regions should be ignored.
/// The method TRestAxionOptics::GetEntranceRing will return the ring
/// number where the photon entered. If the photon did not enter any
/// ring, then it will return -1.
///
/// The following metadata parameters define the ring entrances:
/// * **ringMinRadii**: It contains a list of lower radius values for
/// each ring.
/// * **ringMaxRadii**: It contains a list of higher radius values for
/// each ring.
///
/// On top of that we may define a spider mask which is usually present
/// at integrated x-ray optics as an structure to keep the rings in
/// position. The spider mask will prevent photons from entering inside
/// the optics mirroring system.
///
/// The following parameters are used to define the spider mask geometry:
/// * **spiderArmsSeparationAngle**: It defines the angular distance, measured in
/// radians, between two consecutive spider arms. If this parameter is equal
/// to 0, the spider net mask will be disabled. It is disabled by default.
/// * **spiderOffsetAngle**: It defines the angle at which the first arm
/// is located. The default is 0, being the first arm located at the
/// positive y-axis. It cannot take a negative value.
/// * **spiderWidth**: The width of each specific spider arm. Measured in
/// radians. Default is 2.5 degrees.
///
/// The number of arms will be determined by those parameters.
///
/// The following image is generated as a validation or a way to visualize the
/// TRestAxionOptics::GetEntranceRing method. Each color represents a particle
/// hitting in a different ring. The position is drawn at both, the generator
/// plane and the optics entrance plane. This creates an effect of diffusion at
/// the generator plane since the generator random direction is slightly tilted
/// respect to the optical axis.
///
/// \htmlonly <style>div.image img[src="opticsBasic.png"]{width:750px;}</style> \endhtmlonly
///
/// ![Basic optics validation for method TRestAxionOptics::GetEntranceRing](opticsBasic.png)
///
/// This image was generated using the pipeline/metadata/optics/basic.py script.
/// And the rings description can be found at the corresponding basic.rml file.
///
/// This class is an abstract class, to see an example of its implementation
/// inside a RML configuration file, check for TRestAxionGenericOptics or
/// TRestAxionMCPLOptics.
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
    if (pos.Z() > GetEntranceZPosition()) {
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
        REST_Physics::MoveToPlane(pos, dir, TVector3(0, 0, 1), TVector3(0, 0, GetEntranceZPosition()));

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
    if (pos.Z() > 0 || pos.Z() < GetEntranceZPosition()) {
        RESTWarning << "TRestAxionOptics::TransportToMiddle" << RESTendl;
        RESTWarning << "The particle should be placed between entrance and middle!" << RESTendl;
        return 0;
    }
    if (dir.Z() <= 0) {
        RESTWarning << "TRestAxionOptics::TransportToMiddle" << RESTendl;
        RESTWarning << "Photon is not moving on the positive Z-direction!" << RESTendl;
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
    if (pos.Z() < 0 || pos.Z() > GetExitZPosition()) {
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
        REST_Physics::MoveToPlane(pos, dir, TVector3(0, 0, 1), TVector3(0, 0, GetExitZPosition()));
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

    //// TODO obtain incidence angle and reflectivity

    /// We move the particle to the entrance optics plane. We update fEntrancePosition
    Int_t middleRegion = TransportToMiddle(fFirstInteractionPosition, fMiddleDirection);
    RESTDebug << "Middle region : " << middleRegion << RESTendl;
    if (middleRegion != entranceRegion) return 0.0;

    /// We update the position and direction at the second mirror
    /// We update the position and direction at the second mirror. We update
    /// fSecondInteractionPosition and fExitDirection
    SecondMirrorReflection(fMiddlePosition, fMiddleDirection);

    //// TODO obtain incidence angle and reflectivity

    Int_t exitRegion = TransportToExit(fSecondInteractionPosition, fExitDirection);
    RESTDebug << "Exit region : " << exitRegion << RESTendl;
    if (exitRegion != middleRegion) return 0.0;

    return reflectivity;
}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionOptics field members through a RML file
///
void TRestAxionOptics::InitFromConfigFile() {
    TRestMetadata::InitFromConfigFile();

    /*
std::vector<Double_t> rMax = StringToElements(GetParameter("ringMaxRadii", "-1"), ",");
std::vector<Double_t> rMin = StringToElements(GetParameter("ringMinRadii", "-1"), ",");

if (rMax.size() != rMin.size())
    SetError(
        "TRestAxionOptics. The number of ring radii definitions do not match! Rings will not be "
        "initialized!");
else {
    fRingsRadii.clear();
    for (unsigned int n = 0; n < rMax.size(); n++) {
        std::pair<Double_t, Double_t> p(rMin[n], rMax[n]);
        fRingsRadii.push_back(p);
    }
}
    */

    // If we recover the metadata class from ROOT file we will need to call Initialize ourselves
    this->Initialize();
}

///////////////////////////////////////////////
/// \brief Prints on screen the information about the metadata members of TRestAxionOptics
///
void TRestAxionOptics::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    RESTMetadata << "Entrance position in Z : " << GetEntranceZPosition() << " mm" << RESTendl;
    RESTMetadata << "Exit position in Z : " << GetExitZPosition() << " mm" << RESTendl;
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
        RESTWarning << "TRestAxionOptics::PrintEntranceMask. Not available" << RESTendl;
}

///////////////////////////////////////////////
/// \brief Prints on screen the mask used on the exit optics plane
///
void TRestAxionOptics::PrintExitMask() {
    if (fExitMask)
        fExitMask->PrintMetadata();
    else
        RESTWarning << "TRestAxionOptics::PrintEntranceMask. Not available" << RESTendl;
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
TPad* TRestAxionOptics::DrawMirrors() {
    if (fPad != nullptr) {
        delete fPad;
        fPad = nullptr;
    }

    fPad = new TPad("optics_pad", "This is the optics drawing pad", 0.01, 0.02, 0.99, 0.97);
    fPad->Draw();
    fPad->SetRightMargin(0.09);
    fPad->SetLeftMargin(0.2);
    fPad->SetBottomMargin(0.15);

    return fPad;
}
