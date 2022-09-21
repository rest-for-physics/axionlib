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
/// TRestAxionOpticsProcess TOBE documented
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2019-March:  First implementation of a dummy optics response
///             Javier Galan
///
/// \class      TRestAxionOpticsProcess
/// \author
///
/// <hr>
///
#include "TRestAxionOpticsProcess.h"
using namespace std;

ClassImp(TRestAxionOpticsProcess);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionOpticsProcess::TRestAxionOpticsProcess() { Initialize(); }

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
TRestAxionOpticsProcess::TRestAxionOpticsProcess(char* cfgFileName) {
    Initialize();

    LoadConfig(cfgFileName);
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionOpticsProcess::~TRestAxionOpticsProcess() { delete fAxionEvent; }

///////////////////////////////////////////////
/// \brief Function to load the default config in absence of RML input
///
void TRestAxionOpticsProcess::LoadDefaultConfig() {
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
void TRestAxionOpticsProcess::LoadConfig(std::string cfgFilename, std::string name) {
    if (LoadConfigFromFile(cfgFilename, name)) LoadDefaultConfig();
}

///////////////////////////////////////////////
/// \brief Function to initialize input/output event members and define the section name
///
void TRestAxionOpticsProcess::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    fAxionEvent = new TRestAxionEvent();
}

///////////////////////////////////////////////
/// \brief Process initialization. Data members that require initialization just before start processing
/// should be initialized here.
///
void TRestAxionOpticsProcess::InitProcess() {
    RESTDebug << "Entering ... TRestAxionGeneratorProcess::InitProcess" << RESTendl;

    fOptics = GetMetadata<TRestAxionOptics>();

    if (fOptics) {
        fOptics->PrintMetadata();
    } else {
        RESTError << "TRestAxionOptics::InitProcess. No sucess instantiating optics." << RESTendl;
    }
}

///////////////////////////////////////////////
/// \brief The main processing event function
///
TRestEvent* TRestAxionOpticsProcess::ProcessEvent(TRestEvent* evInput) {
    fAxionEvent = (TRestAxionEvent*)evInput;

    TVector3 inPos = fAxionEvent->GetPosition();
    TVector3 inDir = fAxionEvent->GetDirection();
    Double_t energy = fAxionEvent->GetEnergy();

    if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) fAxionEvent->PrintEvent();
    Double_t efficiency = fOptics->PropagatePhoton(inPos, inDir, energy);

    SetObservableValue("efficiency", efficiency);

    RESTDebug << "Optics efficiency: " << efficiency << RESTendl;

    // We register the event even if it is not properly reflected, i.e. efficiency = 0
    fAxionEvent->SetPosition(fOptics->GetLastGoodPosition());
    fAxionEvent->SetDirection(fOptics->GetLastGoodDirection());

    /// TODO set efficiency inside the event

    if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) {
        fAxionEvent->PrintEvent();

        if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Extreme) GetChar();
    }

    return fAxionEvent;
}

//////////////////////////////////////////////////////////////////////////
/// \brief End of event process.
///
void TRestAxionOpticsProcess::EndOfEventProcess(TRestEvent* evInput) {
    if (fOpticalAxis) {
        // The outgoing particle will be referenced to the optical axis
        TRestAxionEventProcess::EndOfEventProcess(evInput);
    } else {
        // The outgoing particle will be referenced to the universal axis
        TRestEventProcess::EndOfEventProcess(evInput);
    }
}
