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
/// TRestAxionWolterOptics is a class that inherits from TRestAxionOptics.
///
/// ToDO: Write what happens here
///
///
/// ### RML definition
///
/// Example 1:
/// \code
/// <TRestAxionWolterOptics name="dummy">
///   	<parameter name="center" value="(0,0,200)mm" />
///		<parameter name="axis" value="(0,0.02,0.98)" />
///		<parameter name="length" value="22cm" />
///
///		<!-- We build mirror shells with 0.1mm thickness -->
///		<parameter name="shellMinRadii" value="5,10,15,20,25" />
///		<parameter name="shellMaxRadii" value="9.9,14.9,19.9,24.9,29.9" />
/// <TRestAxionWolterOptics/>
/// \endcode
///
/// Example 2:
/// \code
/// <TRestAxionWolterOptics center="(0,0,950)mm" axis="(0,0,1)" />
/// \endcode
///
/// Example 3:
/// \code
/// <TRestAxionWolterOptics center="(0,0,95)" units="cm" />
/// \endcode
///
/// \htmlonly <style>div.image img[src="Wolter.png"]{width:800px;}</style> \endhtmlonly
///
/// ![Wolter optics schematic figure](Wolter.png)
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2022-February: First concept and implementation of TRestAxionWolterOptics class.
///            	  Johanna von Oy
///
/// 2022-May: Final integration
///            	  Javier Galan
///
/// \class      TRestAxionWolterOptics
/// \author     Johanna von Oy <vonoy@physik.uni-bonn.de>
/// \author     Javier Galan <javier.galan@unizar.es>
///
/// <hr>
///

#include "TRestAxionWolterOptics.h"

using namespace std;
#include <cmath>
#include "TRestPhysics.h"
using namespace REST_Physics;
ClassImp(TRestAxionWolterOptics);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionWolterOptics::TRestAxionWolterOptics() : TRestAxionOptics() { Initialize(); }

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionWolterOptics::~TRestAxionWolterOptics() {}

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
/// corresponding TRestAxionWolterOptics section inside the RML.
///
TRestAxionWolterOptics::TRestAxionWolterOptics(const char* cfgFileName, string name)
    : TRestAxionOptics(cfgFileName) {
    RESTDebug << "Entering TRestAxionWolterOptics constructor( cfgFileName, name )" << RESTendl;

    Initialize();

    LoadConfigFromFile(fConfigFileName, name);

    if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Info) PrintMetadata();
}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionWolterOptics members
///
void TRestAxionWolterOptics::Initialize() {
    TRestAxionOptics::Initialize();

    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    if (fAlpha.size() == 0) return;

    fCosAlpha.clear();
    for (const auto& a : fAlpha) fCosAlpha.push_back(TMath::Cos(a));

    fCosAlpha_3.clear();
    for (const auto& a : fAlpha) fCosAlpha_3.push_back(TMath::Cos(3 * a));

    fFrontVertex.clear();
    for (unsigned int n = 0; n < fAlpha.size(); n++) fFrontVertex.push_back(fR3[n] / TMath::Tan(fAlpha[n]));

    fBackVertex.clear();
    for (unsigned int n = 0; n < fAlpha.size(); n++)
        fBackVertex.push_back(fR3[n] / TMath::Tan(3. * fAlpha[n]));

    /// Initializing Entrance mask
    if (fEntranceRingsMask) {
        delete fEntranceRingsMask;
        fEntranceRingsMask = nullptr;
    }
    fEntranceRingsMask = new TRestRingsMask();

    std::vector<Double_t> inner, outer;
    for (unsigned int n = 0; n < fR1.size() - 1; n++) {
        std::cout << "n : " << n << " R1: " << fR1[n] << " Thickness : " << fThickness[n] << std::endl;
        inner.push_back(fR1[n] + fThickness[n]);
        outer.push_back(fR1[n + 1]);
    }
    fEntranceRingsMask->SetRadii(inner, outer);

    fEntranceMask->AddMask(fSpiderMask);
    fEntranceMask->AddMask(fEntranceRingsMask);

    /// Initializing Middle mask
    if (fMiddleRingsMask) {
        delete fMiddleRingsMask;
        fMiddleRingsMask = nullptr;
    }
    fMiddleRingsMask = new TRestRingsMask();

    inner.clear();
    outer.clear();
    for (unsigned int n = 0; n < fR3.size() - 1; n++) {
        inner.push_back(fR3[n] + fThickness[n]);
        outer.push_back(fR3[n + 1]);
    }
    fMiddleRingsMask->SetRadii(inner, outer);

    fMiddleMask->AddMask(fSpiderMask);
    fMiddleMask->AddMask(fMiddleRingsMask);

    /// Initializing Exit mask
    if (fExitRingsMask) {
        delete fExitRingsMask;
        fExitRingsMask = nullptr;
    }
    fExitRingsMask = new TRestRingsMask();

    inner.clear();
    outer.clear();
    for (unsigned int n = 0; n < fR5.size() - 1; n++) {
        inner.push_back(fR5[n] + fThickness[n]);
        outer.push_back(fR5[n + 1]);
    }
    fExitRingsMask->SetRadii(inner, outer);

    fExitMask->AddMask(fSpiderMask);
    fExitMask->AddMask(fExitRingsMask);

    fEntranceMask->Print();
    fMiddleMask->Print();
    fExitMask->Print();
}

void TRestAxionWolterOptics::FirstMirrorReflection(TVector3& pos, TVector3& dir) {
    Int_t mirror = GetMirror();
    if (mirror < 0) {
        RESTError << "TRestAxionWolterOptics::FirstMirrorReflection. Mirror index cannot be negative!"
                  << RESTendl;
        return 0;
    }

    if (mirror >= fFrontVertex.size()) {
        RESTError << "TRestAxionWolterOptics::FirstMirrorReflection. Mirror index above number of mirrors!"
                  << RESTendl;
        return 0;
    }

    TVector3 vertex(0, 0, fFrontVertex[mirror]);
    Double_t cosA = fCosAlpha[mirror];

    //// Reflection on first mirror
    pos = fEntrancePosition +
          fEntranceDirection * REST_Physics::GetConeVectorIntersection(fEntrancePosition, fEntranceDirection,
                                                                       TVector3(0, 0, 1), vertex, cosA);

    TVector3 coneNormal = REST_Physics::GetConeNormal(pos, fAlpha[mirror]);

    dir = GetVectorReflection(fEntranceDirection, coneNormal);
}

void TRestAxionWolterOptics::SecondMirrorReflection(TVector3& pos, TVector3& dir) {
    Int_t mirror = GetMirror();
    if (mirror < 0) {
        RESTError << "TRestAxionWolterOptics::FirstMirrorReflection. Mirror index cannot be negative!"
                  << RESTendl;
        return 0;
    }

    if (mirror >= fFrontVertex.size()) {
        RESTError << "TRestAxionWolterOptics::FirstMirrorReflection. Mirror index above number of mirrors!"
                  << RESTendl;
        return 0;
    }

    TVector3 vertex(0, 0, fBackVertex[mirror]);
    Double_t cosA = fCosAlpha_3[mirror];

    //// Reflection on first mirror
    pos = fMiddlePosition +
          fMiddleDirection * REST_Physics::GetConeVectorIntersection(fMiddlePosition, fMiddleDirection,
                                                                     TVector3(0, 0, 1), vertex, cosA);

    TVector3 coneNormal = REST_Physics::GetConeNormal(pos, 3 * fAlpha[mirror]);

    dir = GetVectorReflection(fMiddleDirection, coneNormal);
}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionWolterOptics field members through a RML file
///
void TRestAxionWolterOptics::InitFromConfigFile() {
    TRestAxionOptics::InitFromConfigFile();

    if (fOpticsFile != "") {
        std::string fullPathFileName = SearchFile(fOpticsFile);

        std::vector<std::vector<Double_t>> opticsData;
        TRestTools::ReadASCIITable(fullPathFileName, opticsData, 3);

        TRestTools::PrintTable(opticsData);

        // The relevant parameters will be fR3 and fAlpha
        // TODO We should check that the rings have been defined in increasing order
        fR1 = TRestTools::GetColumnFromTable(opticsData, 0);
        fR2 = TRestTools::GetColumnFromTable(opticsData, 1);
        fR3 = TRestTools::GetColumnFromTable(opticsData, 2);
        fR4 = TRestTools::GetColumnFromTable(opticsData, 3);
        fR5 = TRestTools::GetColumnFromTable(opticsData, 4);

        fAlpha = TRestTools::GetColumnFromTable(opticsData, 5);
        for (auto& x : fAlpha) x = x * units("rad") / units("deg");

        // For the moment we will only consider fixed length mirrors
        // fLength = TRestTools::GetColumnFromTable(opticsData, 6);

        fThickness = TRestTools::GetColumnFromTable(opticsData, 7);
    }

    if (fSpiderMask) {
        delete fSpiderMask;
        fSpiderMask = nullptr;
    }
    fSpiderMask = (TRestSpiderMask*)this->InstantiateChildMetadata("TRestSpiderMask");

    if (fSpiderMask == nullptr) {
        RESTWarning << "TRestAxionWolterOptics requires usually a TRestSpiderMask definition" << RESTendl;
    } else {
    }

    // If we recover the metadata class from ROOT file we will need to call Initialize ourselves
    this->Initialize();
}

void TRestAxionWolterOptics::PrintParameters() {
    if (fR3.size() > 0) {
        for (unsigned int n = 0; n < fR3.size(); n++) {
            Double_t dR1 = fR1[n] - fR3[n] - fMirrorLength * TMath::Sin(fAlpha[n]);
            Double_t dR5 = fR5[n] - fR3[n] + fMirrorLength * TMath::Sin(3 * fAlpha[n]);
            if (n % 10 == 0)
                std::cout << "## R1\tdelta R1\tR3\tR5\tdelta R5\talpha\tCosAlpha\tFrontVertex\tBackVertex"
                          << std::endl;

            std::cout << fR1[n] << "\t" << dR1 << "\t" << fR3[n] << "\t" << fR5[n] << "\t" << dR5 << "\t"
                      << fAlpha[n] << "\t" << fCosAlpha[n] << "\t" << fFrontVertex[n] << "\t"
                      << fBackVertex[n] << std::endl;
        }
    }
}

void TRestAxionWolterOptics::PrintSpider() {
    if (fSpiderMask) fSpiderMask->PrintMetadata();
}

///////////////////////////////////////////////
/// \brief Prints on screen the information about the metadata members of TRestAxionWolterOptics
///
void TRestAxionWolterOptics::PrintMetadata() {
    TRestAxionOptics::PrintMetadata();

    RESTMetadata << "---------" << RESTendl;
    RESTMetadata << " - Optics file : " << fOpticsFile << RESTendl;
    RESTMetadata << " " << RESTendl;
    RESTMetadata << " Use \"this->PrintMasks()\" to get masks info" << RESTendl;
    RESTMetadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << RESTendl;
}
