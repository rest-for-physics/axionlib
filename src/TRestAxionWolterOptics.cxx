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
/// This class defines an optics device using conical aproximation as it is
/// ilustrated in the following figure.
///
/// \htmlonly <style>div.image img[src="Wolter.png"]{width:800px;}</style> \endhtmlonly
///
/// ![Wolter optics schematic figure](Wolter.png)
///
/// The parameters shown there correspond to the geometrical description of each of
/// the mirrors that build the optics device. Those parameters will be retrieved by
/// TRestAxionOptics, and they will be placed inside TRestAxionOptics::fOpticsData.
/// This class will use those parameters to implement the ray-tracing interactions
/// inside the mirrors, and to define the entrance/middle/interface masks. The
/// format of the file must follow the same specifications as the file
/// [XMM.Wolter](https://github.com/rest-for-physics/axionlib-data/blob/master/optics/XMM.Wolter).
///
/// The optical parameters extracted from this file, in particular `R1`, `R3` and `R5`,
/// will be used to generate the TRestRingsMask at each of the optical interfaces,
/// entrance, middle, and exit. An additional mask, a TRestSpiderMask, common to all
/// the interfaces must be defined in the RML configuration. The TRestRingsMask and the
/// TRestSpiderMask will be used to build a TRestCombinedMask at each interface.
///
/// ### RML definition
///
/// Examples of RML optics definitions can be found at the [axion-lib data
/// repository](https://github.com/rest-for-physics/axionlib-data/tree/master/optics).
///
/// \code
///	<TRestAxionWolterOptics name="xmm" verboseLevel="warning" >
///		<parameter name="opticsFile" value="XMM.Wolter" />
///
///		<parameter name="mirrorLength" value="300" />
///
///		<TRestAxionOpticsMirror name="XMM">
///			<parameter name="mirrorType" value="Single" />
///
///			<parameter name="layerTop" value="Au" />
///			<parameter name="layerThicknessTop" value="250" />
///			<parameter name="sigmaTop" value="0.4" />
///
///			<parameter name="substrate" value="Ni" />
///		</TRestAxionOpticsMirror>
///
///		<TRestSpiderMask name="spider" verboseLevel="warning">
///			<parameter name="maskRadius" value="70cm"/>
///			<parameter name="offset" value="(0,0)cm"/>
///			<parameter name="rotationAngle" value="0"/>
///			<parameter name="armsWidth" value="2.5deg"/>
///			<parameter name="armsSeparationAngle" value="180./16degrees"/>
///			<parameter name="initialRadius" value="0cm"/>
///		</TRestSpiderMask>
///	</TRestAxionWolterOptics>
/// \endcode
///
/// The previous definition corresponds to the following image, that was
/// generated using the method TRestAxionOptics::DrawDensityMaps,
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
#include <TAxis.h>
#include <TGraph.h>
#include <TH1F.h>

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

    fR1 = GetR1();
    fR3 = GetR3();
    fR5 = GetR5();
    fAlpha = GetAlpha();
    fThickness = GetThickness();
    fXSep = 2 * (fR1 - fR3 - fMirrorLength * TMath::Sin(fAlpha)) / TMath::Tan(fAlpha);

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

    if (fRandom != nullptr) {
        delete fRandom;
        fRandom = nullptr;
    }
    fRandom = new TRandom3(0);
}

///////////////////////////////////////////////
/// \brief Implementation of first mirror interaction. It updates fFirstInteractionPosition and
/// fMiddleDirection making use of fEntrancePosition and fEntranceDirection
///
Int_t TRestAxionWolterOptics::FirstMirrorReflection(const TVector3& pos, const TVector3& dir) {
    Int_t mirror = GetMirror();
    RESTDebug << "--> Entering TRestAxionWolterOptics::FirstMirrorReflection" << RESTendl;
    RESTDebug << "Mirror: " << mirror << RESTendl;
    if (mirror < 0) {
        RESTError << "TRestAxionWolterOptics::FirstMirrorReflection. Mirror index cannot be negative!"
                  << RESTendl;
        return -1;
    }

    if (mirror >= 0 && (unsigned int)mirror >= fFrontVertex.size()) {
        RESTError << "TRestAxionWolterOptics::FirstMirrorReflection. Mirror index above number of mirrors!"
                  << RESTendl;
        return -1;
    }

    RESTDebug << "Vertex Z: " << fFrontVertex[mirror] << RESTendl;
    RESTDebug << "Cos Alpha: " << fCosAlpha[mirror] << RESTendl;
    TVector3 vertex(0, 0, fFrontVertex[mirror]);
    Double_t cosA = fCosAlpha[mirror];

    //// Reflection on first mirror
    fFirstInteractionPosition =
        pos + dir * REST_Physics::GetConeVectorIntersection(pos, dir, TVector3(0, 0, -1), vertex, cosA);

    if (fFirstInteractionPosition.Z() < GetEntrancePositionZ() || fFirstInteractionPosition.Z() > 0) {
        RESTDebug << "TRestAxionWolterOptics::FirstMirrorReflection. No interaction!" << RESTendl;
        fFirstInteractionPosition = REST_Physics::MoveByDistance(pos, dir, fMirrorLength / 2.);
        fMiddleDirection = fEntranceDirection;
        fFirstInteraction = false;
        return 0;
    }

    TVector3 coneNormal = REST_Physics::GetConeNormal(fFirstInteractionPosition, fAlpha[mirror]);
    RESTDebug << "Cone normal: (" << coneNormal.X() << ", " << coneNormal.Y() << ", " << coneNormal.Z() << ")"
              << RESTendl;

    fMiddleDirection = GetVectorReflection(fEntranceDirection, coneNormal);

    RESTDebug << "<-- Exiting TRestAxionWolterOptics::FirstMirrorReflection" << RESTendl;

    fFirstInteraction = true;
    return 1;
}

///////////////////////////////////////////////
/// \brief Implementation of first mirror interaction. It updates fSecondInteractionPosition and
/// fExitDirection making use of fMiddlePosition and fMiddleDirection
///
Int_t TRestAxionWolterOptics::SecondMirrorReflection(const TVector3& pos, const TVector3& dir) {
    Int_t mirror = GetMirror();
    if (mirror < 0) {
        RESTError << "TRestAxionWolterOptics::FirstMirrorReflection. Mirror index cannot be negative!"
                  << RESTendl;
        return 0;
    }

    if (mirror >= 0 && (unsigned int)mirror >= fFrontVertex.size()) {
        RESTError << "TRestAxionWolterOptics::FirstMirrorReflection. Mirror index above number of mirrors!"
                  << RESTendl;
        return 0;
    }

    TVector3 vertex(0, 0, fBackVertex[mirror]);
    Double_t cosA = fCosAlpha_3[mirror];

    //// Reflection on first mirror
    fSecondInteractionPosition =
        pos + dir * REST_Physics::GetConeVectorIntersection(pos, dir, TVector3(0, 0, -1), vertex, cosA);

    if (fSecondInteractionPosition.Z() > GetExitPositionZ() || fSecondInteractionPosition.Z() < 0) {
        RESTDebug << "TRestAxionWolterOptics::SecondMirrorReflection. No interaction!" << RESTendl;
        fSecondInteractionPosition = REST_Physics::MoveByDistance(pos, dir, fMirrorLength / 2.);
        fExitDirection = fMiddleDirection;
        fSecondInteraction = false;
        return 0;
    }

    TVector3 coneNormal = REST_Physics::GetConeNormal(fSecondInteractionPosition, 3 * fAlpha[mirror]);

    fExitDirection = GetVectorReflection(fMiddleDirection, coneNormal);

    fSecondInteraction = true;
    return 1;
}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionWolterOptics field members through a RML file
///
void TRestAxionWolterOptics::InitFromConfigFile() {
    if (fSpiderMask) {
        delete fSpiderMask;
        fSpiderMask = nullptr;
    }

    fSpiderMask = (TRestSpiderMask*)this->InstantiateChildMetadata("TRestSpiderMask");

    if (fSpiderMask == nullptr) {
        RESTWarning << "TRestAxionWolterOptics requires usually a TRestSpiderMask definition" << RESTendl;
    } else {
    }

    TRestAxionOptics::InitFromConfigFile();

    // If we recover the metadata class from ROOT file we will need to call Initialize ourselves
    this->Initialize();
}

///////////////////////////////////////////////
/// \brief It prints out the Wolter (relevant) parameters extracted from the optics data file,
/// and other parameters calculated after those input parameters, such as the vertex and angles.
///
/// It will also evaluate the precision loss due to angle to mirror raddius transformation.
///
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

///////////////////////////////////////////////
/// \brief It prints out the spider mask common to all the optical planes
///
void TRestAxionWolterOptics::PrintSpider() {
    if (fSpiderMask) fSpiderMask->PrintMetadata();
}

///////////////////////////////////////////////
/// \brief Prints on screen the information about the metadata members of TRestAxionWolterOptics
///
void TRestAxionWolterOptics::PrintMetadata() { TRestAxionOptics::PrintMetadata(); }

///////////////////////////////////////////////
/// \brief A method to to draw an optics schematic including the mirrors geometry.
///
TPad* TRestAxionWolterOptics::DrawMirrors() {
    TRestAxionOptics::CreatePad();

    fPad->cd();
    std::vector<TGraph*> graphCollection;
    for (unsigned int mirror = 0; mirror < fR3.size(); mirror++) {
        TGraph* gr = new TGraph();  //"Mirror" + IntegerToString(mirror + 1));

        Double_t lX = fMirrorLength * fCosAlpha[mirror];

        gr->SetPoint(0, -lX, fR1[mirror]);
        gr->SetPoint(1, 0, fR3[mirror]);
        gr->SetPoint(2, lX, fR5[mirror]);

        gr->GetXaxis()->SetLimits(-3.5 * lX, 3.5 * lX);
        gr->GetHistogram()->SetMaximum(fR1.back() * 1.15);
        gr->GetHistogram()->SetMinimum(fR1.front() * 0.8);

        gr->GetXaxis()->SetTitle("Z [mm]");
        gr->GetXaxis()->SetTitleSize(0.04);
        gr->GetXaxis()->SetLabelSize(0.04);
        gr->GetXaxis()->SetNdivisions(5);
        gr->GetYaxis()->SetTitle("R [mm]");
        gr->GetYaxis()->SetTitleOffset(1.4);
        gr->GetYaxis()->SetTitleSize(0.04);
        gr->GetYaxis()->SetLabelSize(0.04);
        gr->SetLineWidth(6 * fThickness[mirror]);
        gr->SetLineColor(20 + mirror % 20);
        if (mirror == 0)
            gr->Draw("AL");
        else
            gr->Draw("L");
    }

    return fPad;
}
