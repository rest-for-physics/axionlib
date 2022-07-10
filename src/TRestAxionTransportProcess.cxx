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
/// TRestAxionTransportProcess simply translates the axion to a given
/// z-position without modifying the direction. This process will only
/// allow the movement on the sense of direction, particle cannot move
/// backwards.
///
/// It receives a single parameter that defines the Z-position where
/// the particle should be moved `zPosition`.
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2022-June:  Basic implementation of a photon transport process
///             Javier Galan
///
/// \class      TRestAxionTransportProcess
/// \author
///
/// <hr>
///
#include "TRestAxionTransportProcess.h"
using namespace std;

ClassImp(TRestAxionTransportProcess);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionTransportProcess::TRestAxionTransportProcess() { Initialize(); }

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
TRestAxionTransportProcess::TRestAxionTransportProcess(char* cfgFileName) {
    Initialize();

    LoadConfig(cfgFileName);
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionTransportProcess::~TRestAxionTransportProcess() { delete fAxionEvent; }

///////////////////////////////////////////////
/// \brief Function to load the default config in absence of RML input
///
void TRestAxionTransportProcess::LoadDefaultConfig() {
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
void TRestAxionTransportProcess::LoadConfig(std::string cfgFilename, std::string name) {
    if (LoadConfigFromFile(cfgFilename, name)) LoadDefaultConfig();
}

///////////////////////////////////////////////
/// \brief Function to initialize input/output event members and define the section name
///
void TRestAxionTransportProcess::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    fAxionEvent = new TRestAxionEvent();
}

///////////////////////////////////////////////
/// \brief Process initialization. Data members that require initialization just before start processing
/// should be initialized here.
///
void TRestAxionTransportProcess::InitProcess() {
    RESTDebug << "Entering ... TRestAxionGeneratorProcess::InitProcess" << RESTendl;
}

///////////////////////////////////////////////
/// \brief The main processing event function
///
TRestEvent* TRestAxionTransportProcess::ProcessEvent(TRestEvent* evInput) {
    fAxionEvent = (TRestAxionEvent*)evInput;

    TVector3 inPos = fAxionEvent->GetPosition();
    TVector3 inDir = fAxionEvent->GetDirection();

    if ((inPos.Z() - fZPosition) * inDir.Z() >= 0) {
        RESTWarning
            << "TRestAxionTransportProcess::ProcessEvent. Not appropiate particle direction to reach Z="
            << fZPosition << RESTendl;
        RESTWarning << "Axion position. Z:" << inPos.Z() << " direction Z" << inDir.Z() << RESTendl;
        return fAxionEvent;
    }

    TVector3 newPosition =
        REST_Physics::MoveToPlane(inPos, inDir, TVector3(0, 0, 1), TVector3(0, 0, fZPosition));

    fAxionEvent->SetPosition(newPosition);

    if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) {
        fAxionEvent->PrintEvent();

        if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Extreme) GetChar();
    }

    return fAxionEvent;
}
