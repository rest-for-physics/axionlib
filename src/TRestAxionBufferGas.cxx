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
/// TRestAxionBufferGas is a class used to create an interface to the gas
/// properties used in axion search calculations, such as the photon density
/// and photon absorption.
///
/// This class will allow us to create an arbitray gas mixture by specifying
/// the atomic components, and its contribution to the density. The available
/// elements can be found at the `data/bufferGas/` directory, which are
/// currenty He, Ne, Ar and Xe.
///
/// The following examples show how to create a particular gas mixture. Both
/// examples should lead to the same mixture.
///
/// Example 1:
///
/// \code
/// TRestAxionBufferGas *gas = new TRestAxionBufferGas();
///
/// //Density units must be expressed here in the default REST units, `kg/mm3`.
/// gas->SetGasDensity( "He", 2.6e-9 );
/// gas->SetGasDensity( "Xe", 5.6e-9 );
/// \endcode
///
/// Example 2:
///
/// \code
/// TRestAxionBufferGas *gas = new TRestAxionBufferGas();
///
/// gas->SetGasMixture( "He+Xe", "2.6e-6g/dm^3+5.6mg/m^3" );
/// \endcode
///
/// The corresponding RML section for initialization through a configuration
/// file would be as follows.
///
/// \code
///	<TRestAxionBufferGas name="heliumAndXenon" verboseLevel="warning" >
///		<gas name="He" density="2.6e-9"/>
///		<gas name="Xe" density="5.6mg/cm^3"/>
///	</TRestAxionBufferGas>
/// \endcode
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

#include "TRestSystemOfUnits.h"
using namespace REST_Units;

ClassImp(TRestAxionBufferGas);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionBufferGas::TRestAxionBufferGas() : TRestMetadata() { Initialize(); }

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
/// corresponding TRestGeant4Metadata section inside the RML.
///
TRestAxionBufferGas::TRestAxionBufferGas(const char* cfgFileName, string name) : TRestMetadata(cfgFileName) {
    RESTDebug << "Entering TRestAxionBufferGas constructor( cfgFileName, name )" << RESTendl;

    Initialize();

    LoadConfigFromFile(fConfigFileName, name);
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionBufferGas::~TRestAxionBufferGas() {}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionBufferGas members. It removes all gases.
///
void TRestAxionBufferGas::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    fBufferGasName.clear();
    fBufferGasDensity.clear();

    fAbsEnergy.clear();
    fGasAbsCoefficient.clear();

    fFactorEnergy.clear();
    fGasFormFactor.clear();
}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionBufferGas field members through a RML file
///
void TRestAxionBufferGas::InitFromConfigFile() {
    this->Initialize();

    auto gasDefinition = GetElement("gas");
    while (gasDefinition) {
        TString gasName = GetFieldValue("name", gasDefinition);
        if (gasName.Contains("+")) {
            TString gasDensity = GetFieldValue("density", gasDefinition);
            SetGasMixture(gasName, gasDensity);
        } else {
            Double_t gasDensity = GetDblParameterWithUnits("density", gasDefinition);
            SetGasDensity(gasName, gasDensity);
        }
        gasDefinition = GetNextElement(gasDefinition);
    }

    if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Info) PrintMetadata();
}

///////////////////////////////////////////////
/// \brief It adds a new gas component to the mixture. If it already exists it will update its density.
///
/// Density must be given in standard REST units: kg/mm^3
///
void TRestAxionBufferGas::SetGasDensity(TString gasName, Double_t density) {
    Int_t gasIndex = FindGasIndex(gasName);

    fBufferGasDensity[gasIndex] = density;
}

///////////////////////////////////////////////
/// \brief It re-initializes the gas mixture to the one provided by argument.
///
/// The first argument must be the gas components separated by the sign +.
/// The second argument are the corresponding densities with units.
///
/// Example : SetGasMixture("Ne+Xe", "2e-3g/cm3+3mg/cm3" );
///
/// If the second argument with densities is not given, the buffer gas will
/// add the gas components with zero density.
///
void TRestAxionBufferGas::SetGasMixture(TString gasMixture, TString gasDensities) {
    Initialize();

    std::vector<string> names = Split((string)gasMixture, "+");
    std::vector<string> densities;
    if (gasDensities == "0")
        densities.clear();
    else
        densities = Split((string)gasDensities, "+");

    if (!densities.empty() && names.size() == densities.size()) {
        for (unsigned int n = 0; n < names.size(); n++) {
            Double_t density = GetValueInRESTUnits(densities[n]);
            SetGasDensity(names[n], density);
        }
    } else if (densities.empty()) {
        for (unsigned int n = 0; n < names.size(); n++) {
            SetGasDensity(names[n], 0);
        }
    } else {
        this->SetError("SetGasMixture. Number of gases does not match the densities!");
    }
}

///////////////////////////////////////////////
/// \brief It returns the gas density - from the chosen gasName component - in g/cm3.
///
Double_t TRestAxionBufferGas::GetGasDensity(TString gasName) {
    Int_t gasIndex = FindGasIndex(gasName);

    if (gasIndex < 0) {
        RESTError << "TRestAxionBufferGas::GetGasDensity. Gas name : " << gasName << " not found!"
                  << RESTendl;
        return 0;
    }

    return fBufferGasDensity[gasIndex];
}

///////////////////////////////////////////////
/// \brief It reads the data files from the corresponding gas component.
///
/// This method will need to find under `data/bufferGas/` the corresponding
/// scattering form factors and absorption files for component X, i.e. X.nff
/// and X.abs.
///
/// TODO. This method might better use SearchFile method that searches in
/// globals <searchPath definition.
///
void TRestAxionBufferGas::ReadGasData(TString gasName) {
    TString factorFileName = SearchFile((string)gasName + ".nff");

    RESTDebug << "TRestAxionBufferGas::ReadGasData. Reading factor file : " << factorFileName << RESTendl;

    if (!TRestTools::fileExists((string)factorFileName)) {
        RESTError << "TRestAxionBufferGas::ReadGasData( " << gasName << " )" << RESTendl;
        RESTError << "Gas factor file not found : " << factorFileName << RESTendl;
        exit(1);
    }

    /* {{{ Reading form factor values */
    FILE* fin = fopen(factorFileName.Data(), "rt");

    double en, value;

    std::vector<Double_t> energyFactor;
    std::vector<Double_t> factor;

    // keV -- e/atom
    while (fscanf(fin, "%lf\t%lf\n", &en, &value) != EOF) {
        RESTDebug << "Energy : " << en << "keV -- Factor : " << value << RESTendl;

        energyFactor.push_back(en);
        factor.push_back(value);
    }

    RESTDebug << "Items read : " << energyFactor.size() << RESTendl;

    fclose(fin);
    /* }}} */

    TString absFileName = SearchFile((string)gasName + ".abs");

    RESTDebug << "TRestAxionBufferGas::ReadGasData. Reading factor file : " << absFileName << RESTendl;

    if (!TRestTools::fileExists((string)absFileName)) {
        RESTError << "TRestAxionBufferGas::ReadGasData( " << gasName << " )" << RESTendl;
        RESTError << "Gas absorption file not found : " << absFileName << RESTendl;
        exit(1);
    }

    /* {{{ Reading absorption values */
    fin = fopen(absFileName.Data(), "rt");

    std::vector<Double_t> energyAbs;
    std::vector<Double_t> absorption;

    // keV -- e/atom
    while (fscanf(fin, "%lf\t%lf\n", &en, &value) != EOF) {
        RESTDebug << "Energy : " << en << "keV -- Absorption : " << value << RESTendl;

        energyAbs.push_back(en);
        absorption.push_back(value);
    }

    RESTDebug << "Items read : " << energyAbs.size() << RESTendl;

    fclose(fin);
    /* }}} */

    fBufferGasName.push_back(gasName);

    fFactorEnergy.push_back(energyFactor);
    fGasFormFactor.push_back(factor);

    fAbsEnergy.push_back(energyAbs);
    fGasAbsCoefficient.push_back(absorption);

    fBufferGasDensity.push_back(0);
}

///////////////////////////////////////////////
/// \brief It returns the atomic form factor of the `gasName` component at the given energy.
///
/// Energy input parameter must be given in keV
///
Double_t TRestAxionBufferGas::GetFormFactor(TString gasName, Double_t energy) {
    // In case we are in vacuum
    if (GetNumberOfGases() == 0) return 0;

    Int_t gasIndex = FindGasIndex(gasName);
    RESTDebug << "TRestAxionBufferGas::GetFormFactor. Gas index = " << gasIndex << RESTendl;

    if (gasIndex == -1) {
        ReadGasData(gasName);
        gasIndex = FindGasIndex(gasName);
    }

    if (gasIndex == -1) {
        RESTError << "TRestAxionBufferGas::GetFormFactor. Gas: " << gasName << " Not Found!" << RESTendl;
        exit(1);
    }

    Int_t energyIndex = GetEnergyIndex(fFactorEnergy[gasIndex], energy);
    RESTDebug << "Energy index : " << energyIndex << RESTendl;

    if (energyIndex == -1) {
        RESTError << "TRestAxionBufferGas::GetFormFactor. Energy (" << energy << " keV) out of range"
                  << RESTendl;
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
        RESTError << "TRestAxionBufferGas::GetAbsorptionCoefficient. Negative coefficient" << RESTendl;
        cout << "y2 : " << y2 << " y1 : " << y1 << endl;
        cout << "x2 : " << x2 << " x1 : " << x1 << endl;
        cout << "m : " << m << " n : " << n << endl;
        cout << "E : " << energy << " bin : " << energyIndex << endl;
        GetChar();
    }

    return (m * energy + n);
}

///////////////////////////////////////////////
/// \brief It returns the inverse of the absorption lenght, for the gas mixture, in cm-1, for the given
/// energy in keV.
///
Double_t TRestAxionBufferGas::GetPhotonAbsorptionLength(Double_t energy) {
    Double_t attLength = 0;
    for (unsigned int n = 0; n < fBufferGasName.size(); n++)
        attLength +=
            fBufferGasDensity[n] * units("g/cm^3") * GetAbsorptionCoefficient(fBufferGasName[n], energy);

    return attLength;
}

///////////////////////////////////////////////
/// \brief It returns the inverse of the absorption lenght, for the gas mixture, in eV, for the given
/// energy in keV.
///
Double_t TRestAxionBufferGas::GetPhotonAbsorptionLengthIneV(Double_t energy) {
    return cmToeV(GetPhotonAbsorptionLength(energy));
}

///////////////////////////////////////////////
/// \brief It transforms cm-1 to eV
///
Double_t TRestAxionBufferGas::cmToeV(double l_Inv)  // E in keV, P in bar ---> Gamma in cm-1
{
    return l_Inv / REST_Physics::PhMeterIneV / 0.01;
}

///////////////////////////////////////////////
/// \brief It returns the equivalent photon mass (in eV) for the gas mixture at the given input energy
/// expressed in keV.
///
Double_t TRestAxionBufferGas::GetPhotonMass(double en) {
    Double_t photonMass = 0;
    for (unsigned int n = 0; n < fBufferGasName.size(); n++) {
        Double_t W_value = 0;
        if (fBufferGasName[n] == "H") W_value = 1.00794;   // g/mol
        if (fBufferGasName[n] == "He") W_value = 4.002;    // g/mol
        if (fBufferGasName[n] == "Ne") W_value = 20.179;   // g/mol
        if (fBufferGasName[n] == "Ar") W_value = 39.948;   // g/mol
        if (fBufferGasName[n] == "Xe") W_value = 131.293;  // g/mol

        if (W_value == 0) {
            RESTError << "Gas name : " << fBufferGasName[n] << " is not implemented in TRestBufferGas!!"
                      << RESTendl;
            RESTError << "W value must be defined in TRestAxionBufferGas::GetPhotonMass" << RESTendl;
            RESTError << "This gas will not contribute to the calculation of the photon mass!" << RESTendl;
        } else {
            photonMass +=
                fBufferGasDensity[n] * units("g/cm^3") * GetFormFactor(fBufferGasName[n], en) / W_value;
        }
    }

    return 28.77 * TMath::Sqrt(photonMass);
}

////////////////////////////////////////////
/// \brief It returns the equivalent gas density for a given photon mass expressed in eV. You have to define
/// previously the gas type.
///	The resulting density will be expressed in kg/mm^3, which are the standard REST Units.
///
Double_t TRestAxionBufferGas::GetMassDensity(double m_gamma) {
    Double_t massDensity = 0;
    for (unsigned int n = 0; n < fBufferGasName.size(); n++) {
        Double_t W_value = 0;
        Double_t Z_value = 0;
        if (fBufferGasName[n] == "H") {
            W_value = 1.00794;  // g/mol
            Z_value = 1;
        } else if (fBufferGasName[n] == "He") {
            W_value = 4.002602;  // g/mol
            Z_value = 2;
        } else if (fBufferGasName[n] == "Ne") {
            W_value = 20.179;  // g/mol
            Z_value = 10;
        } else if (fBufferGasName[n] == "Ar") {
            W_value = 39.948;  // g/mol
            Z_value = 18;
        } else if (fBufferGasName[n] == "Xe") {
            W_value = 131.293;  // g/mol
            Z_value = 54;
        }
        if (W_value == 0) {
            RESTError << "Gas name : " << fBufferGasName[n] << " is not implemented in TRestAxionBufferGas!!"
                      << RESTendl;
            RESTError << "W value must be defined in TRestAxionBufferGas::GetMassDensity" << RESTendl;
            RESTError << "This gas will not contribute to the calculation of the photon mass!" << RESTendl;
        } else {
            massDensity += pow(m_gamma, 2) * W_value / (Z_value * pow(28.77, 2));
        }
    }
    return massDensity / units("g/cm^3");
}

///////////////////////////////////////////////
/// \brief It returns the absorption coefficient, in cm2/g, for the given gas component and
/// energy in keV.
///
Double_t TRestAxionBufferGas::GetAbsorptionCoefficient(TString gasName, Double_t energy) {
    Int_t gasIndex = FindGasIndex(gasName);
    RESTDebug << "TRestAxionBufferGas::GetAbsorptionCoefficient. Gas index = " << gasIndex << RESTendl;

    if (gasIndex == -1) {
        ReadGasData(gasName);
        gasIndex = FindGasIndex(gasName);
    }

    if (gasIndex == -1) {
        RESTError << "TRestAxionBufferGas::GetAbsorptionCoefficient. Gas: " << gasName << " Not Found!"
                  << RESTendl;
        exit(1);
    }

    Int_t energyIndex = GetEnergyIndex(fAbsEnergy[gasIndex], energy);
    RESTDebug << "Energy index : " << energyIndex << RESTendl;

    if (energyIndex == -1) {
        RESTError << "TRestAxionBufferGas::GetAbsorptionCoefficient. Energy out of range" << RESTendl;
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
        RESTError << "TRestAxionBufferGas::GetAbsorptionCoefficient. Negative coeffient" << RESTendl;
        cout << "y2 : " << y2 << " y1 : " << y1 << endl;
        cout << "x2 : " << x2 << " x1 : " << x1 << endl;
        cout << "m : " << m << " n : " << n << endl;
        cout << "E : " << energy << " bin : " << energyIndex << endl;
        GetChar();
    }

    return (m * energy + n);
}

///////////////////////////////////////////////
/// \brief It returns the vector element index, from `enVector`, that is just below the given input energy.
///
Int_t TRestAxionBufferGas::GetEnergyIndex(std::vector<Double_t> enVector, Double_t energy) {
    for (unsigned int n = 0; n < enVector.size(); n++)
        if (energy < enVector[n]) return n - 1;

    return -1;
}

///////////////////////////////////////////////
/// \brief It returns the internal index of the gas component given by `gasName`.
///
Int_t TRestAxionBufferGas::FindGasIndex(TString gasName) {
    Int_t nGases = (Int_t)fBufferGasName.size();

    for (int n = 0; n < nGases; n++)
        if (fBufferGasName[n] == gasName) return n;

    // If the gas is not found in the list we just force reading it from its gas datafile.
    ReadGasData(gasName);

    return FindGasIndex(gasName);
}

///////////////////////////////////////////////
/// \brief Prints out the absorption coefficients as function of the energy for the given gas component, for
/// debugging pourposes.
///
void TRestAxionBufferGas::PrintAbsorptionGasData(TString gasName) {
    Int_t gasIndex = FindGasIndex(gasName);

    for (unsigned int n = 0; n < fAbsEnergy[gasIndex].size(); n++)
        cout << "energy : " << fAbsEnergy[gasIndex][n]
             << " -- abs coeff : " << fGasAbsCoefficient[gasIndex][n] << endl;
}

///////////////////////////////////////////////
/// \brief Prints out the atomic form factors as function of the energy for the given gas component, for
/// debugging pourposes.
///
void TRestAxionBufferGas::PrintFormFactorGasData(TString gasName) {
    Int_t gasIndex = FindGasIndex(gasName);

    for (unsigned int n = 0; n < fAbsEnergy[gasIndex].size(); n++)
        cout << "Energy : " << fFactorEnergy[gasIndex][n] << " -- Abs coeff : " << fGasFormFactor[gasIndex][n]
             << endl;
}

///////////////////////////////////////////////
/// \brief Prints on screen the information about the metadata members of TRestAxionBufferGas
///
void TRestAxionBufferGas::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    RESTMetadata << "Number of buffer gases defined : " << fBufferGasName.size() << RESTendl;
    if (fBufferGasName.size() == 0) {
        RESTMetadata << "Buffer medium is vacuum" << RESTendl;
    } else {
        RESTMetadata << "Photon mass at 4keV : " << this->GetPhotonMass(4.) << " eV" << RESTendl;
        RESTMetadata << " " << RESTendl;
        RESTMetadata << "Buffer gases inside mixture : " << RESTendl;
        RESTMetadata << "---------------------------" << RESTendl;
        for (unsigned int n = 0; n < fBufferGasName.size(); n++) {
            RESTMetadata << " Gas name : " << fBufferGasName[n] << RESTendl;
            RESTMetadata << " Gas density : " << fBufferGasDensity[n] * units("g/cm^3") << " g/cm3"
                         << RESTendl;
            RESTMetadata << " Form factor energy range : ( " << fFactorEnergy[n][0] << ", "
                         << fFactorEnergy[n].back() << " ) keV" << RESTendl;
            RESTMetadata << " Absorption energy range : ( " << fAbsEnergy[n][0] << ", "
                         << fAbsEnergy[n].back() << " ) keV" << RESTendl;
            RESTMetadata << " " << RESTendl;
        }
    }

    RESTMetadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << RESTendl;
}
