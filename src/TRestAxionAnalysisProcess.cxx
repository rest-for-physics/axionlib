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
/// TRestAxionAnalysisProcess adds to the analysis tree information about
/// the particle state.
///
/// Observable values defined by this process:
/// - **energy**: The energy in keV.
/// - **posX**: The X-position in mm.
/// - **posY**: The Y-position in mm.
/// - **posZ**: The Z-position in mm.
/// - **thetaAngle**: The angle in radians respect to the propagation axis, defined in Z.
/// - **phiAngle**: The angle moving along the XY-plane.
/// - **R**: The distance to the Z-axis in mm.
///
/// The axion is generated with intensity proportional to g_ag = 1.0 x g10
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2019-March:  First implementation of shared memory buffer to rawsignal conversion.
///             Javier Galan
///
/// \class      TRestAxionAnalysisProcess
/// \author     Javier Galan
///
/// <hr>
///
#include "TRestAxionAnalysisProcess.h"
using namespace std;

ClassImp(TRestAxionAnalysisProcess);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionAnalysisProcess::TRestAxionAnalysisProcess() { Initialize(); }

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
TRestAxionAnalysisProcess::TRestAxionAnalysisProcess(char* cfgFileName) {
    Initialize();

    LoadConfig(cfgFileName);
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionAnalysisProcess::~TRestAxionAnalysisProcess() { delete fAxionEvent; }

///////////////////////////////////////////////
/// \brief Function to load the default config in absence of RML input
///
void TRestAxionAnalysisProcess::LoadDefaultConfig() {
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
void TRestAxionAnalysisProcess::LoadConfig(std::string cfgFilename, std::string name) {
    if (LoadConfigFromFile(cfgFilename, name)) LoadDefaultConfig();
}

///////////////////////////////////////////////
/// \brief Function to initialize input/output event members and define the section name
///
void TRestAxionAnalysisProcess::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    fAxionEvent = new TRestAxionEvent();

    fAxionEvent->Initialize();
}

///////////////////////////////////////////////
/// \brief The main processing event function
///
TRestEvent* TRestAxionAnalysisProcess::ProcessEvent(TRestEvent* evInput) {
    fAxionEvent = (TRestAxionEvent*)evInput;

    RESTDebug << "TRestAxionAnalysisProcess::ProcessEvent : " << fAxionEvent->GetID() << RESTendl;

    SetObservableValue("energy", fAxionEvent->GetEnergy());

    Double_t x = fAxionEvent->GetPosition().X();
    Double_t y = fAxionEvent->GetPosition().Y();
    SetObservableValue("posX", x);
    SetObservableValue("posY", y);
    SetObservableValue("posZ", fAxionEvent->GetPosition().Z());

    Double_t r = TMath::Sqrt(x * x + y * y);
    SetObservableValue("R", r);

    SetObservableValue("thetaAngle", fAxionEvent->GetDirection().Theta());
    SetObservableValue("phiAngle", fAxionEvent->GetDirection().Phi());

    if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) fAxionEvent->PrintEvent();

    return fAxionEvent;
}
