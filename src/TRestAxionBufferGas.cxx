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
/// TRestAxionBufferGas is a class used to ...
///
///
/// TODO. Create an appropriate documentation here.
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2019-March: First concept and implementation of TRestAxionBufferGas class.
///             Javier Galan
///
/// \class      TRestAxionBufferGas
/// \author     Javier Galan
///
/// <hr>
///

#include "TRestAxionBufferGas.h"
using namespace std;

ClassImp(TRestAxionBufferGas);
//______________________________________________________________________________
TRestAxionBufferGas::TRestAxionBufferGas() : TRestMetadata() {
    // TRestAxionBufferGas default constructor
    Initialize();
}

//______________________________________________________________________________
TRestAxionBufferGas::TRestAxionBufferGas(const char* cfgFileName, string name) : TRestMetadata(cfgFileName) {
    cout << "Entering TRestAxionBufferGas constructor( cfgFileName, name )" << endl;

    Initialize();

    LoadConfigFromFile(fConfigFileName, name);

    PrintMetadata();
}

//______________________________________________________________________________
TRestAxionBufferGas::~TRestAxionBufferGas() {
    // TRestAxionBufferGas destructor
}

void TRestAxionBufferGas::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);
}

//______________________________________________________________________________
void TRestAxionBufferGas::InitFromConfigFile() {
    this->Initialize();

    // Initialize the metadata members from a configfile

    // fClassMember = GetParameter( "paramName", "defaultValue" );

    PrintMetadata();
}

void TRestAxionBufferGas::SetGasDensity(TString gasName, Double_t density) {
    Int_t gasIndex = FindGasIndex(gasName);

    fBufferGasDensity[gasIndex] = density;
}

// Returns the gas density in g/cm3
Double_t TRestAxionBufferGas::GetGasDensity(TString gasName) {
    Int_t gasIndex = FindGasIndex(gasName);

    return fBufferGasDensity[gasIndex];
}

void TRestAxionBufferGas::ReadGasData(TString gasName) {
    TString factorFileName = (TString)getenv("REST_PATH") + "/data/bufferGas/" + gasName + ".nff";

    debug << "TRestAxionBufferGas::ReadGasData. Reading factor file : " << factorFileName << endl;

    if (!fileExists((string)factorFileName)) {
        error << "TRestAxionBufferGas::ReadGasData( " << gasName << " )" << endl;
        error << "Gas factor file not found : " << factorFileName << endl;
        exit(1);
    }

    /* {{{ Reading form factor values */
    FILE* fin = fopen(factorFileName.Data(), "rt");

    double en, value;

    std::vector<Double_t> energyFactor;
    std::vector<Double_t> factor;

    // keV -- e/atom
    while (fscanf(fin, "%lf\t%lf\n", &en, &value) != EOF) {
        debug << "Energy : " << en << "keV -- Factor : " << value << endl;

        energyFactor.push_back(en);
        factor.push_back(value);
    }

    debug << "Items read : " << energyFactor.size() << endl;

    fclose(fin);
    /* }}} */

    TString absFileName = (TString)getenv("REST_PATH") + "/data/bufferGas/" + gasName + ".abs";

    debug << "TRestAxionBufferGas::ReadGasData. Reading factor file : " << absFileName << endl;

    if (!fileExists((string)absFileName)) {
        error << "TRestAxionBufferGas::ReadGasData( " << gasName << " )" << endl;
        error << "Gas absorption file not found : " << absFileName << endl;
        exit(1);
    }

    /* {{{ Reading absorption values */
    fin = fopen(absFileName.Data(), "rt");

    std::vector<Double_t> energyAbs;
    std::vector<Double_t> absorption;

    // keV -- e/atom
    while (fscanf(fin, "%lf\t%lf\n", &en, &value) != EOF) {
        debug << "Energy : " << en << "keV -- Absorption : " << value << endl;

        energyAbs.push_back(en);
        absorption.push_back(value);
    }

    debug << "Items read : " << energyAbs.size() << endl;

    fclose(fin);
    /* }}} */

    fBufferGasName.push_back(gasName);

    fFactorEnergy.push_back(energyFactor);
    fGasFormFactor.push_back(factor);

    fAbsEnergy.push_back(energyAbs);
    fGasAbsCoefficient.push_back(absorption);

    fBufferGasDensity.push_back(0);
}

Double_t TRestAxionBufferGas::GetFormFactor(TString gasName, Double_t energy) {
    Int_t gasIndex = FindGasIndex(gasName);
    debug << "TRestAxionBufferGas::GetFormFactor. Gas index = " << gasIndex << endl;

    if (gasIndex == -1) {
        ReadGasData(gasName);
        gasIndex = FindGasIndex(gasName);
    }

    Int_t energyIndex = GetEnergyIndex(fFactorEnergy[gasIndex], energy);
    debug << "Energy index : " << energyIndex << endl;

    if (energyIndex == -1) {
        error << "TRestAxionBufferGas::GetAbsorptionCoefficient. Energy out of range" << endl;
        exit(1);
    }

    // Absorption coefficient
    double y2 = fGasFormFactor[gasIndex][energyIndex + 1];
    double y1 = fGasFormFactor[gasIndex][energyIndex];

    // Normalized field
    double x2 = fFactorEnergy[gasIndex][energyIndex + 1];
    double x1 = fFactorEnergy[gasIndex][energyIndex];

    double m = (y2 - y1) / (x2 - x1);
    double n = y1 - m * x1;

    if (m * energy + n < 0) {
        error << "TRestAxionBufferGas::GetAbsorptionCoefficient. Negative coeffient" << endl;
        cout << "y2 : " << y2 << " y1 : " << y1 << endl;
        cout << "x2 : " << x2 << " x1 : " << x1 << endl;
        cout << "m : " << m << " n : " << n << endl;
        cout << "E : " << energy << " bin : " << energyIndex << endl;
        GetChar();
    }

    return (m * energy + n);
}

// energy in keV ---> returns Gamma in cm-1
Double_t TRestAxionBufferGas::GetPhotonAbsorptionLength(Double_t energy) {
    Double_t attLength = 0;
    for (unsigned int n = 0; n < fBufferGasName.size(); n++)
        attLength += fBufferGasDensity[n] * GetAbsorptionCoefficient(fBufferGasName[n], energy);

    return attLength;
}

Double_t TRestAxionBufferGas::GetPhotonAbsorptionLengthIneV(Double_t energy) {
    return cmToeV(GetPhotonAbsorptionLength(energy));
}

// Transforms cm-1 to eV
Double_t TRestAxionBufferGas::cmToeV(double l_Inv)  // E in keV, P in bar ---> Gamma in cm-1
{
    return l_Inv / REST_Physics::PhMeterIneV / 0.01;
}

Double_t TRestAxionBufferGas::GetPhotonMass(double en)  // in eV
{
    Double_t photonMass = 0;
    for (unsigned int n = 0; n < fBufferGasName.size(); n++) {
        // const double Wa_Helium = 4.002; // g/mol
        // const double Wa_Neon = 20.179; // g/mol
        // const double Wa_Argon = 39.948; // g/mol
        // const double Wa_Xenon = 131.293; // g/mol

        Double_t W_value = 0;
        if (fBufferGasName[n] == "He") W_value = 4.002;
        if (fBufferGasName[n] == "Ne") W_value = 20.179;
        if (fBufferGasName[n] == "Ar") W_value = 39.948;
        if (fBufferGasName[n] == "Xe") W_value = 131.293;

        if (W_value == 0) {
            error << "Gas name : " << fBufferGasName[n] << " is not implemented in TRestBufferGas!!" << endl;
            error << "W value must be defined in TRestAxionBufferGas::GetPhotonMass" << endl;
            error << "This gas will not contribute to the calculation of the photon mass!" << endl;
        } else {
            photonMass += fBufferGasDensity[n] * GetFormFactor(fBufferGasName[n], en) / W_value;
        }
    }

    return 28.77 * TMath::Sqrt(photonMass);
}

Double_t TRestAxionBufferGas::GetAbsorptionCoefficient(TString gasName, Double_t energy) {
    Int_t gasIndex = FindGasIndex(gasName);
    debug << "TRestAxionBufferGas::GetAbsorptionCoefficient. Gas index = " << gasIndex << endl;

    Int_t energyIndex = GetEnergyIndex(fAbsEnergy[gasIndex], energy);
    debug << "Energy index : " << energyIndex << endl;

    if (energyIndex == -1) {
        error << "TRestAxionBufferGas::GetAbsorptionCoefficient. Energy out of range" << endl;
        exit(1);
    }

    // Absorption coefficient
    double y2 = fGasAbsCoefficient[gasIndex][energyIndex + 1];
    double y1 = fGasAbsCoefficient[gasIndex][energyIndex];

    // Normalized field
    double x2 = fAbsEnergy[gasIndex][energyIndex + 1];
    double x1 = fAbsEnergy[gasIndex][energyIndex];

    double m = (y2 - y1) / (x2 - x1);
    double n = y1 - m * x1;

    if (m * energy + n < 0) {
        error << "TRestAxionBufferGas::GetAbsorptionCoefficient. Negative coeffient" << endl;
        cout << "y2 : " << y2 << " y1 : " << y1 << endl;
        cout << "x2 : " << x2 << " x1 : " << x1 << endl;
        cout << "m : " << m << " n : " << n << endl;
        cout << "E : " << energy << " bin : " << energyIndex << endl;
        GetChar();
    }

    return (m * energy + n);
}

Int_t TRestAxionBufferGas::GetEnergyIndex(std::vector<Double_t> enVector, Double_t energy) {
    for (unsigned int n = 0; n < enVector.size(); n++)
        if (energy < enVector[n]) return n - 1;

    return -1;
}

Int_t TRestAxionBufferGas::FindGasIndex(TString gasName) {
    Int_t nGases = (Int_t)fBufferGasName.size();

    for (int n = 0; n < nGases; n++)
        if (fBufferGasName[n] == gasName) return n;

    // If the gas is not found in the list we just force reading it from its gas datafile.
    ReadGasData(gasName);

    return FindGasIndex(gasName);
}

void TRestAxionBufferGas::PrintAbsorptionGasData(TString gasName) {
    Int_t gasIndex = FindGasIndex(gasName);

    for (unsigned int n = 0; n < fAbsEnergy[gasIndex].size(); n++)
        cout << "energy : " << fAbsEnergy[gasIndex][n]
             << " -- abs coeff : " << fGasAbsCoefficient[gasIndex][n] << endl;
}

void TRestAxionBufferGas::PrintFormFactorGasData(TString gasName) {
    Int_t gasIndex = FindGasIndex(gasName);

    for (unsigned int n = 0; n < fAbsEnergy[gasIndex].size(); n++)
        cout << "Energy : " << fFactorEnergy[gasIndex][n] << " -- Abs coeff : " << fGasFormFactor[gasIndex][n]
             << endl;
}

void TRestAxionBufferGas::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    metadata << "Number of buffer gases defined : " << fBufferGasName.size() << endl;
    if (fBufferGasName.size() == 0) {
        metadata << "Buffer medium is vacuum" << endl;
    } else {
        metadata << "Buffer gases defined : " << endl;
        metadata << "---------------------------" << endl;
        for (unsigned int n = 0; n < fBufferGasName.size(); n++) {
            metadata << " Gas name : " << fBufferGasName[n] << endl;
            metadata << " Gas density : " << fBufferGasDensity[n] << " g/cm3" << endl;
            metadata << " Form factor energy range : ( " << fFactorEnergy[n][0] << ", "
                     << fFactorEnergy[n].back() << " ) keV" << endl;
            metadata << " Absorption energy range : ( " << fAbsEnergy[n][0] << ", " << fAbsEnergy[n].back()
                     << " ) keV" << endl;
            metadata << " " << endl;
        }
    }

    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}
