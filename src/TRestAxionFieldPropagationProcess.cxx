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
/// TRestAxionFieldPropagationProcess TOBE documented
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
/// \class      TRestAxionFieldPropagationProcess
/// \author     Javier Galan
///
/// <hr>
///
#include "TRestAxionFieldPropagationProcess.h"
using namespace std;

ClassImp(TRestAxionFieldPropagationProcess);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionFieldPropagationProcess::TRestAxionFieldPropagationProcess() { Initialize(); }

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
TRestAxionFieldPropagationProcess::TRestAxionFieldPropagationProcess(char* cfgFileName) {
    Initialize();

    LoadConfig(cfgFileName);
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionFieldPropagationProcess::~TRestAxionFieldPropagationProcess() {
    delete fInputAxionEvent;
    delete fOutputAxionEvent;
}

///////////////////////////////////////////////
/// \brief Function to load the default config in absence of RML input
///
void TRestAxionFieldPropagationProcess::LoadDefaultConfig() {
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
void TRestAxionFieldPropagationProcess::LoadConfig(std::string cfgFilename, std::string name) {
    if (LoadConfigFromFile(cfgFilename, name)) LoadDefaultConfig();
}

///////////////////////////////////////////////
/// \brief Function to initialize input/output event members and define the section name
///
void TRestAxionFieldPropagationProcess::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    fCharSizeExp = 30000.;

    fInputAxionEvent = new TRestAxionEvent();
    fOutputAxionEvent = new TRestAxionEvent();

    fInputEvent = fInputAxionEvent;
    fOutputEvent = fOutputAxionEvent;
}

///////////////////////////////////////////////
/// \brief The main processing event function
///
void TRestAxionFieldPropagationProcess::InitProcess() {
    debug << "Entering ... TRestAxionGeneratorProcess::InitProcess" << endl;

    fAxionMagneticField = (TRestAxionMagneticField*)this->GetMetadata("TRestAxionMagneticField");

    if (!fAxionMagneticField) {
        error << "TRestAxionFieldPropagationprocess. Magnetic Field was not defined!" << endl;
        exit(0);
    }

    fAxionBufferGas = (TRestAxionBufferGas*)this->GetMetadata("TRestAxionBufferGas");

    if (!fAxionBufferGas) {
        error << "TRestAxionBufferGas. Cannot access the buffer gas" << endl;
        exit(0);
    }

    fAxionPhotonConversion = new TRestAxionPhotonConversion();
    fAxionPhotonConversion->SetBufferGas(fAxionBufferGas);
}

TVector3 TRestAxionFieldPropagationProcess::MoveOneStep(TVector3 pos, TVector3 dir, Double_t step) {
    Double_t t;
    Double_t f;

    Int_t i, j, k;

    if (dir[1] != 0) {
        j = 1;
        k = 2;
        i = 0;
    } else {
        if (dir[0] != 0) {
            j = 0;
            i = 1;
            k = 2;
        } else {
            j = 2;
            i = 0;
            k = 1;
        }
    }

    f = pos[j] + (dir[j] / abs(dir[j])) * step;
    t = (f - pos[j]) / dir[j];

    pos[j] = f;
    pos[i] = pos[i] + t * dir[i];
    pos[k] = pos[k] + t * dir[k];

    return pos;
}

TVector3 TRestAxionFieldPropagationProcess::MoveVirtualBox(TVector3 pos, TVector3 dir, Double_t size) {
    if (size == -1) size = fCharSizeExp;

    Double_t step;
    Int_t i;

    if (dir[1] != 0)
        i = 1;
    else {
        if (dir[0] != 0)
            i = 0;
        else
            i = 2;
    }

    if (pos[i] < size && pos[i] > -size) return pos;

    if (pos[i] < 0) size = -size;

    step = abs(pos[i] - size);
    pos = MoveOneStep(pos, dir, step);

    return pos;
}

bool TRestAxionFieldPropagationProcess::ConditionStop(TVector3 pos, TVector3 dir) {
    if (dir[1] != 0)
        return ((pos[1] < 0 && (dir[1] < 0)) || (pos[1] > fCharSizeExp && (dir[1] > 0)));
    else {
        if (dir[0] != 0)
            return (pos[0] < -fCharSizeExp || pos[0] > fCharSizeExp);

        else
            return (pos[2] < -fCharSizeExp || pos[2] > fCharSizeExp);
    }
}

std::vector<TVector3> TRestAxionFieldPropagationProcess::FindBoundariesVolume(TVector3 pos, TVector3 dir,
                                                                              Double_t minStep) {
    std::vector<TVector3> boundary;
    TVector3 boundaryIn;
    TVector3 boundaryOut;

    Double_t dr = 5.;

    pos = MoveOneStep(pos, dir, minStep);

    while (dr > minStep) {
        while (fAxionMagneticField->GetMagneticField(pos[0], pos[1], pos[2]) == TVector3(0, 0, 0) &&
               not(ConditionStop(pos, dir)))
            pos = MoveOneStep(pos, dir, dr);

        if (ConditionStop(pos, dir)) {
            boundaryIn = pos;
            boundaryOut = TVector3(0, 0, 0);
            boundary.push_back(boundaryIn);
            boundary.push_back(boundaryOut);
            return boundary;
        }

        else {
            pos = MoveOneStep(pos, dir, -dr);
            dr = dr / 2.0;
        }
    }

    pos = MoveOneStep(pos, dir, dr * 2.0);
    boundaryIn = pos;

    while (fAxionMagneticField->GetMagneticField(pos[0], pos[1], pos[2]) != TVector3(0, 0, 0))
        pos = MoveOneStep(pos, dir, dr);

    boundaryOut = pos;

    boundary.push_back(boundaryIn);
    boundary.push_back(boundaryOut);

    return boundary;
}

std::vector<std::vector<TVector3>> TRestAxionFieldPropagationProcess::FindFieldBoundaries(Double_t minStep) {
    if (minStep == -1) minStep = 0.01;

    std::vector<std::vector<TVector3>> boundaryCollection;

    TVector3 posInitial = *(fInputAxionEvent->GetPosition());
    TVector3 direction = *(fInputAxionEvent->GetDirection());
    direction = direction.Unit();

    if (direction == TVector3(0, 0, 0) || (direction[1] * posInitial[1] > 0)) return boundaryCollection;

    std::vector<TVector3> buffVect;

    posInitial = MoveVirtualBox(posInitial, direction);

    std::vector<TVector3> bInt;
    bInt = FindBoundariesVolume(posInitial, direction, minStep);

    if (ConditionStop(bInt[0], direction))
        return boundaryCollection;

    else {
        buffVect.push_back(bInt[0]);
        buffVect.push_back(bInt[1]);
        boundaryCollection.push_back(buffVect);
        buffVect.clear();

        bInt = FindBoundariesVolume(bInt[1], direction, minStep);

        while (not(ConditionStop(bInt[0], direction))) {
            buffVect.push_back(bInt[0]);
            buffVect.push_back(bInt[1]);
            boundaryCollection.push_back(buffVect);
            buffVect.clear();

            bInt = FindBoundariesVolume(bInt[1], direction, minStep);
        }

        bInt.clear();
        return boundaryCollection;
    }
}

TVectorD TRestAxionFieldPropagationProcess::GetFieldVector(TVector3 in, TVector3 out, Int_t N) {
    if (N == 0) N = TMath::Power(10, 4);

    TVectorD Bt(N);
    TVector3 B;
    TVector3 direction = *(fInputAxionEvent->GetDirection());

    TVector3 differential = out - in;

    B = fAxionMagneticField->GetMagneticField(in[0], in[1], in[2]);
    Bt[0] = abs(B.Perp(direction));

    for (Int_t i = 1; i < N; i++) {
        in = in + differential * (Double_t(i) / Double_t(N - 1));
        B = fAxionMagneticField->GetMagneticField(in[0], in[1], in[2]);
        Bt[i] = abs(B.Perp(direction));
    }

    return Bt;
}

TRestEvent* TRestAxionFieldPropagationProcess::ProcessEvent(TRestEvent* evInput) {
    fInputAxionEvent = (TRestAxionEvent*)evInput;
    fOutputAxionEvent = fInputAxionEvent;

    Double_t Ea = fInputAxionEvent->GetEnergy();
    Double_t ma = 1.0;  // TODO : Add the axion mass information in the AxionEvent later ?

    std::vector<std::vector<TVector3>> boundaries;
    boundaries = FindFieldBoundaries();
    Int_t NofVolumes = boundaries.size();

    Double_t probability = 0.;

    if ( NofVolumes !=0 )
    {
         TVectorD B;
         TVectorD probabilities(NofVolumes);

         TVector3 lengthVector ;
         Double_t length;

         for ( Int_t i = 0; i < NofVolumes; i++ ) {
               lengthVector = boundaries[i][0]-boundaries[i][1];
               length = sqrt( lengthVector.Mag2() );
               B = GetFieldVector(boundaries[i][0], boundaries[i][1], 0);
               probabilities[i] = fAxionPhotonConversion->GammaTransmissionProbability(Ea, B, ma, length);
         }

         for (Int_t i = 0; i < NofVolumes; i++) probability = probability + probabilities[i];
    }


    fOutputAxionEvent->SetGammaProbability(probability); // Sum of probabilities or difference with the initial probability ?
    fOutputAxionEvent->SetPosition(boundaries[NofVolumes][1]); // TODO : Set the position after propagation of the axion at a fixed diatnce 

    if (GetVerboseLevel() >= REST_Debug) fOutputAxionEvent->PrintEvent();

    return fOutputEvent;
}

///////////////////////////////////////////////
/// \brief Function reading input parameters from the RML TRestAxionFieldPropagationProcess metadata section
///
void TRestAxionFieldPropagationProcess::InitFromConfigFile() {}
