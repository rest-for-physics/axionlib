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
///
/// This generator will produce solar axion events with the Sun at position
/// (0,0,-AU) and detector at position (0,0,0).
///
/// The axion generation properties, such as coupling type and intensity, are
/// defined by TRestAxionSolarFlux. That class defines the method
/// TRestAxionSolarFlux::GetRandomEnergyAndRadius that is exploited by this
/// process to define initial conditions for each axion, such as angular
/// direction or energy. The flux intensity and its dependency with the solar
/// radius and spectral energy will be fully coded inside TRestAxionSolarFlux,
/// while this method will be able to play with the Sun-generator position,
/// and the size of the target detector (which is placed at (0,0,0)) to
/// define an initial axion position, direction and energy.
///
/// The following metadata members can be tuned in this process:
/// - `axionMass:` It defines the axion mass in eV required later on by axion-photon
/// conversion processes. Its default value is 0.
/// - `detectorRadius:` All generated axion events will end up in a circular region
/// placed at the XY plane. This parameter defines the radius of the circular region.
/// The default value is 10cm.
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2022-February:  New generator implementation using TRestAxionSolarFlux
///             Javier Galan
///
/// \class      TRestAxionGeneratorProcess
/// \author     Javier Galan
///
/// <hr>
///
#include "TRestAxionGeneratorProcess.h"
using namespace std;

ClassImp(TRestAxionGeneratorProcess);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionGeneratorProcess::TRestAxionGeneratorProcess() { Initialize(); }

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
TRestAxionGeneratorProcess::TRestAxionGeneratorProcess(char* cfgFileName) {
    Initialize();

    LoadConfig(cfgFileName);
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionGeneratorProcess::~TRestAxionGeneratorProcess() { delete fOutputAxionEvent; }

///////////////////////////////////////////////
/// \brief Function to load the default config in absence of RML input
///
void TRestAxionGeneratorProcess::LoadDefaultConfig() {
    SetName("axionGenerator-Default");
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
void TRestAxionGeneratorProcess::LoadConfig(std::string cfgFilename, std::string name) {
    if (LoadConfigFromFile(cfgFilename, name)) LoadDefaultConfig();
}

///////////////////////////////////////////////
/// \brief Function to initialize input/output event members and define the section name
///
void TRestAxionGeneratorProcess::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    fOutputAxionEvent = new TRestAxionEvent();

    fIsExternal = true;
    fSingleThreadOnly = true;
}

///////////////////////////////////////////////
/// \brief Process initialization. Data members that require initialization just before start processing
/// should be initialized here.
///
void TRestAxionGeneratorProcess::InitProcess() {
    debug << "Entering ... TRestAxionGeneratorProcess::InitProcess" << endl;

    fAxionFlux = GetMetadata<TRestAxionSolarFlux>();

    if (fAxionFlux == nullptr) {
        if (!this->GetError()) this->SetError("The solar flux definition was not found.");
    }

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
TRestEvent* TRestAxionGeneratorProcess::ProcessEvent(TRestEvent* evInput) {
    debug << "TRestAxionGeneratorProcess::ProcessEvent : " << fCounter << endl;
    fOutputAxionEvent->SetID(fCounter);
    fCounter++;

    std::pair<Double_t, Double_t> p = fAxionFlux->GetRandomEnergyAndRadius();
    Double_t energy = p.first;
    Double_t radius = p.second;

    // We avoid the use of expensive trigonometric functions
    Double_t x, y, r;
    do {
        x = fRandom->Rndm() - 0.5;
        y = fRandom->Rndm() - 0.5;
        r = x * x + y * y;
    } while (r > 1 || r == 0);

    r = TMath::Sqrt(r);

    TVector3 axionPosition(REST_Physics::solarRadius * radius * x / r,
                           REST_Physics::solarRadius * radius * y / r, -REST_Physics::AU);

    TVector3 axionDirection = -axionPosition.Unit();

    fOutputAxionEvent->SetEnergy(energy);
    fOutputAxionEvent->SetPosition(axionPosition);
    fOutputAxionEvent->SetDirection(axionDirection);
    fOutputAxionEvent->SetMass(fAxionMass);

    if (GetVerboseLevel() >= REST_Debug) fOutputAxionEvent->PrintEvent();

    return fOutputAxionEvent;
}

void TRestAxionGeneratorProcess::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    metadata << "Axion mass: " << fAxionMass * units("eV") << " eV" << endl;
    metadata << "Detector radius: " << fDetectorRadius * units("cm") << " cm" << endl;
    metadata << "Random seed: " << (UInt_t)fSeed << endl;

    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

