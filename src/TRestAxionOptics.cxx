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
/// hitting in a different ring. The position is drawn at the generator plane,
/// and not at the optics plane. This creates an effect of diffusion since the
/// generator random direction is slightly tilted respect to the optical axis.
///
/// /// \htmlonly <style>div.image img[src="opticsBasic.png"]{width:750px;}</style> \endhtmlonly
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
/// corresponding TRestAxionMagneticField section inside the RML.
///
TRestAxionOptics::TRestAxionOptics(const char* cfgFileName, string name) : TRestMetadata(cfgFileName) {
    debug << "Entering TRestAxionOptics constructor( cfgFileName, name )" << endl;

    Initialize();

    LoadConfigFromFile(fConfigFileName, name);

    if (GetVerboseLevel() >= REST_Info) PrintMetadata();
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionOptics::~TRestAxionOptics() {}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionOptics members
///
void TRestAxionOptics::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    fEntrance = fCenter - 0.5 * fLength * fAxis;
    fExit = fCenter + 0.5 * fLength * fAxis;

    SetMaxAndMinRingRadius();

    if (fSpiderOffsetAngle < 0) fSpiderOffsetAngle = 0;
    if (fSpiderArmsSeparationAngle > 0) InitializeSpiderAngles();

    // A vector orthogonal to the axis, thus parallel to the optics plane
    // and to be used as reference for spider structure. It defines angle = 0
    fReference = TVector3(0, fAxis.Z(), -fAxis.Y()).Unit();
}

///////////////////////////////////////////////
/// \brief It initializes the fMaxRingRadius and fMinRingRadius using the ring ring definitions
///
void TRestAxionOptics::SetMaxAndMinRingRadius() {
    if (fMinRingRadius == -1 && fRingsRadii.size() > 0) fMinRingRadius = fRingsRadii[0].first;
    if (fMaxRingRadius == -1 && fRingsRadii.size() > 0) fMaxRingRadius = fRingsRadii[0].second;

    for (const auto& ringRadius : fRingsRadii) {
        if (ringRadius.first < fMinRingRadius) fMinRingRadius = ringRadius.first;
        if (ringRadius.second > fMaxRingRadius) fMaxRingRadius = ringRadius.second;
    }
}

///////////////////////////////////////////////
/// \brief It initializes the fMaxRingRadius and fMinRingRadius using the ring ring definitions
///
void TRestAxionOptics::InitializeSpiderAngles() {
    std::pair<Double_t, Double_t> additional_negative = {-1, -1};

    Double_t angle = fSpiderOffsetAngle;
    do {
        Double_t angle_down = angle - fSpiderWidth / 2.;
        Double_t angle_up = angle + fSpiderWidth / 2.;

        if (angle_down < 0) {
            additional_negative = {2 * TMath::Pi() + angle_down, 2 * TMath::Pi()};
            fSpiderPositiveRanges.push_back({0, angle_up});

        } else if (angle_up > TMath::Pi() && angle_down < TMath::Pi()) {
            fSpiderPositiveRanges.push_back({angle_down, TMath::Pi()});
            fSpiderNegativeRanges.push_back({TMath::Pi(), angle_up});
        } else if (angle_up < TMath::Pi()) {
            fSpiderPositiveRanges.push_back({angle_down, angle_up});
        } else if (angle_down >= TMath::Pi()) {
            fSpiderNegativeRanges.push_back({angle_down, angle_up});
        }

        angle += fSpiderArmsSeparationAngle;

    } while (angle + 1.e-3 < 2 * TMath::Pi());

    if (additional_negative.first != -1 && additional_negative.second != -1)
        fSpiderNegativeRanges.push_back(additional_negative);

    debug << "Printing positive spider angles" << endl;
    debug << "-------------------------------" << endl;
    for (int n = 0; n < fSpiderPositiveRanges.size(); n++) {
        debug << "n : " << n << " from : " << 180 * fSpiderPositiveRanges[n].first / TMath::Pi() << " to "
              << 180 * fSpiderPositiveRanges[n].second / TMath::Pi() << endl;
    }

    debug << "Printing negative spider angles" << endl;
    debug << "-------------------------------" << endl;
    for (int n = 0; n < fSpiderNegativeRanges.size(); n++) {
        debug << "n : " << n << " from : " << 180 * fSpiderNegativeRanges[n].first / TMath::Pi() << " to "
              << 180 * fSpiderNegativeRanges[n].second / TMath::Pi() << endl;
    }

    for (int n = 0; n < fSpiderNegativeRanges.size(); n++) {
        fSpiderNegativeRanges[n].first = TMath::Cos(fSpiderNegativeRanges[n].first);
        fSpiderNegativeRanges[n].second = TMath::Cos(fSpiderNegativeRanges[n].second);
    }

    for (int n = 0; n < fSpiderPositiveRanges.size(); n++) {
        fSpiderPositiveRanges[n].first = TMath::Cos(fSpiderPositiveRanges[n].first);
        fSpiderPositiveRanges[n].second = TMath::Cos(fSpiderPositiveRanges[n].second);
    }

    debug << "Printing positive spider angles" << endl;
    debug << "-------------------------------" << endl;
    for (int n = 0; n < fSpiderPositiveRanges.size(); n++) {
        debug << "n : " << n << " from : " << fSpiderPositiveRanges[n].first << " to "
              << fSpiderPositiveRanges[n].second << endl;
    }

    debug << "Printing negative spider cosines" << endl;
    debug << "--------------------------------" << endl;
    for (int n = 0; n < fSpiderNegativeRanges.size(); n++) {
        debug << "n : " << n << " from : " << fSpiderNegativeRanges[n].first << " to "
              << fSpiderNegativeRanges[n].second << endl;
    }
}

///////////////////////////////////////////////
/// \brief It moves a given particle with position `pos` and direction `dir` to the entrance of the optics
///
TVector3 TRestAxionOptics::GetPositionAtEntrance(const TVector3& pos, const TVector3& dir) {
    return REST_Physics::MoveToPlane(pos, dir, fAxis, fEntrance);
}

///////////////////////////////////////////////
/// \brief It determines if a given particle with position `pos` (that **should be** already placed
/// at the entrance plane of the optics) will be inside a ring with radius between `Rin` and `Rout`.
///
/// This method is private since it is just a helper method used by GetRing, and it will only
/// work under the assumption that the particle position is already at the entrance optics plane.
///
Bool_t TRestAxionOptics::IsInsideRing(const TVector3& pos, Double_t Rout, Double_t Rin) {
    Double_t d = REST_Physics::DistanceToAxis(fEntrance, fAxis, pos);

    if (d < Rout && d >= Rin) return true;

    return false;
}

///////////////////////////////////////////////
/// \brief It returns the optics ring index for a given particle with position `pos` and direction `dir`.
///
/// If the particle is not inside any ring it will return -1
///
/// This method is protected since it should be used only by derived classes
///
Int_t TRestAxionOptics::GetEntranceRing(const TVector3& pos, const TVector3& dir) {
    TVector3 posEntrance = GetPositionAtEntrance(pos, dir);
    if (HitsSpider(posEntrance)) return -1;
    if (!IsInsideRing(posEntrance, fMaxRingRadius, fMinRingRadius)) return -1;

    int n = 0;
    for (const auto& ringRadius : fRingsRadii) {
        if (IsInsideRing(posEntrance, ringRadius.second, ringRadius.first)) return n;
        n++;
    }

    return -1;
}

///////////////////////////////////////////////
/// \brief It returns `true` if the particle with position `pos` and direction `dir`
/// encounters the spider net.
///
/// This method is private since it is just a helper method used by GetEntranceRing, and it will only
/// work under the assumption that the particle position is already at the entrance optics plane.
///
Bool_t TRestAxionOptics::HitsSpider(const TVector3& pos) {
    Double_t d = REST_Physics::DistanceToAxis(fEntrance, fAxis, pos);
    if (fSpiderArmsSeparationAngle == 0 || d < fSpiderStartRadius) return false;

    TVector3 posUnit = (pos - fEntrance).Unit();
    Double_t cos_angle = posUnit.Dot(fReference);

    if (posUnit.X() >= 0) {
        for (const auto& ang : fSpiderPositiveRanges) {
            if (cos_angle < ang.first && cos_angle > ang.second) return true;
        }
    } else {
        for (const auto& ang : fSpiderNegativeRanges) {
            if (cos_angle > ang.first && cos_angle < ang.second) return true;
        }
    }

    return false;
}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionOptics field members through a RML file
///
void TRestAxionOptics::InitFromConfigFile() {
    TRestMetadata::InitFromConfigFile();

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

    // If we recover the metadata class from ROOT file we will need to call Initialize ourselves
    this->Initialize();
}

///////////////////////////////////////////////
/// \brief Prints on screen the information about the metadata members of TRestAxionOptics
///
void TRestAxionOptics::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    metadata << "Optics length: " << fLength << " mm" << endl;
    metadata << "Optics entrance: (" << fEntrance.X() << ", " << fEntrance.Y() << ", " << fEntrance.Z()
             << ") mm" << endl;
    metadata << "Optics center: (" << fCenter.X() << ", " << fCenter.Y() << ", " << fCenter.Z() << ") mm"
             << endl;
    metadata << "Optics exit: (" << fExit.X() << ", " << fExit.Y() << ", " << fExit.Z() << ") mm" << endl;
    metadata << "Optics axis: (" << fAxis.X() << ", " << fAxis.Y() << ", " << fAxis.Z() << ")" << endl;
    metadata << " " << endl;
    metadata << "Relation of mirror rings integrated in the optics:" << endl;
    metadata << "---------" << endl;
    int n = 0;
    for (const auto& ringRadius : fRingsRadii) {
        metadata << "Ring " << n << ": Rmin = " << ringRadius.first << "mm , Rmax = " << ringRadius.second
                 << "mm" << endl;
        n++;
    }
    if (fSpiderArmsSeparationAngle != 0) {
        metadata << " " << endl;
        metadata << "Spider net structure parameters:" << endl;
        metadata << "--------------------------------" << endl;
        metadata << " - Arms separation angle : " << 180. * fSpiderArmsSeparationAngle / TMath::Pi()
                 << " degrees" << endl;
        metadata << " - First arm offset angle : " << 180. * fSpiderOffsetAngle / TMath::Pi() << " degrees"
                 << endl;
        metadata << " - Arm angular width : " << 180. * fSpiderWidth / TMath::Pi() << " degrees" << endl;
        metadata << " - Spider start radius : " << fSpiderStartRadius * units("cm") << " cm" << endl;
    }
    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}
