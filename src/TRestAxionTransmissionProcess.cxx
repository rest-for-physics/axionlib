/*************************************************************************
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
/// TRestAxionTransmissionProcess will serve to include photon transmission
/// efficiency (or response) through different media before reaching the
/// detection detector volume.
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2019-March:  First implementation
///             Javier Galan
///
/// \class      TRestAxionTransmissionProcess
/// \author
///
/// <hr>
///
#include "TRestAxionTransmissionProcess.h"
using namespace std;

ClassImp(TRestAxionTransmissionProcess);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionTransmissionProcess::TRestAxionTransmissionProcess() { Initialize(); }

///////////////////////////////////////////////
/// \brief Constructor loading data from a config file
///
/// If no configuration path is defined using TRestMetadata::SetConfigFilePath
/// the path to the config file must be specified using full path, absolute or relative.
///
/// The default behaviour is that the config file must be specified with
/// full path, absolute or relative.
///
/// \param cfgFileName A const char* giving the path to an RML file.
///
TRestAxionTransmissionProcess::TRestAxionTransmissionProcess(char* cfgFileName) {
    Initialize();

    LoadConfig(cfgFileName);
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionTransmissionProcess::~TRestAxionTransmissionProcess() { delete fAxionEvent; }

///////////////////////////////////////////////
/// \brief Function to load the default config in absence of RML input
///
void TRestAxionTransmissionProcess::LoadDefaultConfig() {
    SetName(this->ClassName());
    SetTitle("Default config");
}

///////////////////////////////////////////////
/// \brief Function to load the configuration from an external configuration file.
///
/// If no configuration path is defined in TRestMetadata::SetConfigFilePath
/// the path to the config file must be specified using full path, absolute or relative.
///
/// \param cfgFileName A const char* giving the path to an RML file.
/// \param name The name of the specific metadata. It will be used to find the
/// correspondig TRestGeant4AnalysisProcess section inside the RML.
///
void TRestAxionTransmissionProcess::LoadConfig(std::string cfgFilename, std::string name) {
    if (LoadConfigFromFile(cfgFilename, name)) LoadDefaultConfig();
}

///////////////////////////////////////////////
/// \brief Function to initialize input/output event members and define the section name
///
void TRestAxionTransmissionProcess::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    fAxionEvent = new TRestAxionEvent();
}

///////////////////////////////////////////////
/// \brief Process initialization. Data members that require initialization just before start processing
/// should be initialized here.
///
void TRestAxionTransmissionProcess::InitProcess() {
    RESTDebug << "Entering ... TRestAxionGeneratorProcess::InitProcess" << RESTendl;

    fXrayWindows.clear();
    for (const auto& wName : fWindowNames) {
        TRestAxionXrayWindow* w = (TRestAxionXrayWindow*)GetMetadata(wName);

        if (w == nullptr) {
            RESTError << "TRestAxionTransmissionProcess. Window definition with name : " << wName
                      << " not found!" << RESTendl;
        } else {
            fXrayWindows.push_back(w);
        }
    }

    // fRandom = new TRandom3(fSeed);
    // if (fSeed == 0) fSeed = fRandom->GetSeed();
}

///////////////////////////////////////////////
/// \brief The main processing event function
///
TRestEvent* TRestAxionTransmissionProcess::ProcessEvent(TRestEvent* evInput) {
    fAxionEvent = (TRestAxionEvent*)evInput;

    TVector3 inPos = fAxionEvent->GetPosition();
    TVector3 inDir = fAxionEvent->GetDirection();

    /// The component is placed at (0,0,0). It is TRestAxionEventProcess the responsible to translate the
    /// component (in reality the particle) according to fCenter.
    TVector3 newPos = REST_Physics::MoveToPlane(inPos, inDir, TVector3(0, 0, 1), TVector3(0, 0, 0));

    Double_t transmission = 1;
    Double_t x = newPos.X();
    Double_t y = newPos.Y();
    Double_t z = newPos.Z();
    Double_t en = fAxionEvent->GetEnergy();

    RESTDebug << "Particle position to evaluate window transmission. " << RESTendl;
    RESTDebug << "X : " << x << " Y: " << y << " Z: " << z << RESTendl;

    for (const auto& window : fXrayWindows) {
        transmission *= window->GetTransmission(en, x, y);
    }
    RESTDebug << "Transmission: " << transmission << RESTendl;

    SetObservableValue("transmission", transmission);

    if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) {
        fAxionEvent->PrintEvent();

        if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Extreme) GetChar();
    }

    return fAxionEvent;
}

///////////////////////////////////////////////
/// \brief Function reading input parameters from the RML TRestAxionTransmissionProcess metadata section
///
void TRestAxionTransmissionProcess::InitFromConfigFile() {
    TRestEventProcess::InitFromConfigFile();

    // This is the additional code required by the process to read window names
    TiXmlElement* windowDefinition = GetElement("window");
    while (windowDefinition) {
        fWindowNames.push_back(GetFieldValue("name", windowDefinition));

        windowDefinition = GetNextElement(windowDefinition);
    }
}
