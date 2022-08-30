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
/// TRestAxionDeviationProcess simply produces a deviation given by an
/// angle, fDevAngle, that defines the cone directrix delimiting the
/// deviation along the main axion direction.
///
/// For the moment we simply produce a uniform distribution for the angle.
/// This will likely produce a non-homogeneous angular distribution, but
/// still good for small angles.
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2022-July:  Basic implementation of a photon transport process
///             Javier Galan
///
/// \class      TRestAxionDeviationProcess
/// \author
///
/// <hr>
///
#include "TRestAxionDeviationProcess.h"
using namespace std;

ClassImp(TRestAxionDeviationProcess);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionDeviationProcess::TRestAxionDeviationProcess() { Initialize(); }

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
TRestAxionDeviationProcess::TRestAxionDeviationProcess(char* cfgFileName) {
    Initialize();

    LoadConfig(cfgFileName);
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionDeviationProcess::~TRestAxionDeviationProcess() { delete fAxionEvent; }

///////////////////////////////////////////////
/// \brief Function to load the default config in absence of RML input
///
void TRestAxionDeviationProcess::LoadDefaultConfig() {
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
void TRestAxionDeviationProcess::LoadConfig(std::string cfgFilename, std::string name) {
    if (LoadConfigFromFile(cfgFilename, name)) LoadDefaultConfig();
}

///////////////////////////////////////////////
/// \brief Function to initialize input/output event members and define the section name
///
void TRestAxionDeviationProcess::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    fAxionEvent = new TRestAxionEvent();
}

///////////////////////////////////////////////
/// \brief Process initialization. Data members that require initialization just before start processing
/// should be initialized here.
///
void TRestAxionDeviationProcess::InitProcess() {
    RESTDebug << "Entering ... TRestAxionGeneratorProcess::InitProcess" << RESTendl;

    if (!fRandom) {
        delete fRandom;
        fRandom = nullptr;
    }

    fRandom = new TRandom3(fSeed);
    if (fSeed == 0) fSeed = fRandom->GetSeed();
}

///////////////////////////////////////////////
/// \brief The main processing event function
///
TRestEvent* TRestAxionDeviationProcess::ProcessEvent(TRestEvent* evInput) {
    fAxionEvent = (TRestAxionEvent*)evInput;

    TVector3 inPos = fAxionEvent->GetPosition();
    TVector3 inDir = fAxionEvent->GetDirection();

    Double_t theta = fDevAngle * fRandom->Rndm();
    Double_t phi = TMath::Pi() * fRandom->Rndm();

    TVector3 ortho = inDir.Orthogonal();
    TVector3 newDir = inDir;

    newDir.Rotate(theta, ortho);
    newDir.Rotate(phi, inDir);

    fAxionEvent->SetDirection(newDir);

    if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) {
        fAxionEvent->PrintEvent();

        if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Extreme) GetChar();
    }

    return fAxionEvent;
}
