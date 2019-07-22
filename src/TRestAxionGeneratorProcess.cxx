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
/// TRestAxionGeneratorProcess TOBE documented
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

    fInputEvent = NULL;
    fOutputEvent = fOutputAxionEvent;

    fRandom = new TRandom3(0);
}

void TRestAxionGeneratorProcess::InitProcess() {
    debug << "Entering ... TRestAxionGeneratorProcess::InitProcess" << endl;

    fAxionSolarModel = (TRestAxionSolarModel*)this->GetMetadata("TRestAxionSolarModel");

    if (!fAxionSolarModel) {
        error << "TRestAxionGeneratorProcess. Axion model was not defined!" << endl;
        exit(0);
    }
}

///////////////////////////////////////////////
/// \brief This method gets a random energy relying on the solar axion model defined in TRestAxionSolarModel
///
Double_t TRestAxionGeneratorProcess::GenerateEnergy() {
    debug << "Entering TRestAxionGeneratorProcess::GenerateEnergy() ..." << endl;
    Double_t solarFlux =
        fAxionSolarModel->GetSolarAxionFlux(fEnergyRange.X(), fEnergyRange.Y(), 1., fEnergyStep);

    Double_t random = solarFlux * fRandom->Uniform(0, 1.0);

    Double_t sum = 0;
    for (double en = fEnergyRange.X(); en < fEnergyRange.Y(); en += fEnergyStep) {
        sum += fAxionSolarModel->GetDifferentialSolarAxionFlux(en, 1.0) * fEnergyStep;

        if (random < sum) {
            debug << "TRestAxionGeneratorProcess::GenerateEnergy. Energy = " << en << endl;
            return en + fRandom->Uniform(0, fEnergyStep);
        }
    }

    return fEnergyRange.Y();
}

TVector3 TRestAxionGeneratorProcess::GenerateDirection() {
    TVector3 direction;

    // For the moment "circleWall" just generates in the XY-plane
    if (fAngularDistribution == "flux") {
        return fAngularDirection;
    }

    warning << "Angular distribution : " << fAngularDistribution << " is not defined!" << endl;

    return fAngularDirection;
}

TVector3 TRestAxionGeneratorProcess::GeneratePosition() {
    TVector3 position;

    // For the moment "circleWall" just generates in the XY-plane
    if (fSpatialDistribution == "circleWall") {
        Double_t r, x, y;

        do {
            x = fRandom->Uniform(-fSpatialRadius, fSpatialRadius);
            y = fRandom->Uniform(-fSpatialRadius, fSpatialRadius);

            r = x * x + y * y;
        } while (r > fSpatialRadius * fSpatialRadius);

        position = TVector3(x, y, 0) + fSpatialOrigin;

        return position;
    }

    warning << "Spatial distribution : " << fSpatialDistribution << " is not defined!" << endl;

    return position;
}

///////////////////////////////////////////////
/// \brief The main processing event function
///
TRestEvent* TRestAxionGeneratorProcess::ProcessEvent(TRestEvent* evInput) {
    debug << "TRestAxionGeneratorProcess::ProcessEvent : " << fCounter << endl;
    fOutputAxionEvent->SetID(fCounter);
    fCounter++;

    fOutputAxionEvent->SetEnergy(GenerateEnergy());
    fOutputAxionEvent->SetPosition(GeneratePosition());
    fOutputAxionEvent->SetDirection(GenerateDirection());

    // cout << "anaTree : " << fAnalysisTree << endl;
    // fAnalysisTree->SetObservableValue( this, "energy", fOutputAxionEvent->GetEnergy() );
    // fAnalysisTree->PrintObservables();

    if (GetVerboseLevel() >= REST_Debug) fOutputAxionEvent->PrintEvent();

    return fOutputAxionEvent;
}

///////////////////////////////////////////////
/// \brief Function reading input parameters from the RML TRestAxionGeneratorProcess metadata section
///
void TRestAxionGeneratorProcess::InitFromConfigFile() {
    // These 2-params should be moved to TRestAxionSolarModel
    fEnergyRange = Get2DVectorParameterWithUnits("energyRange", TVector2(0, 10));
    fEnergyStep = GetDoubleParameterWithUnits("energyStep", 1.e-3);

    fAngularDistribution = GetParameter("angularDistribution", "flux");
    fAngularDirection = StringTo3DVector(GetParameter("angularDirection"));//, "(0,0,1)"));

    fSpatialDistribution = GetParameter("spatialDistribution", "circleWall");
    fSpatialRadius = GetDoubleParameterWithUnits("spatialRadius", 10.e3);
    fSpatialOrigin = Get3DVectorParameterWithUnits("spatialOrigin");//, TVector3(0, 0, -30000));
}
