/******************** REST disclaimer ***********************************
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
/// TRestAxionSolarFlux will use a file in ASCII or binary format to initialize
/// a solar flux table that will describe the solar flux spectrum as a function
/// of the solar radius. It will be also possible to generate the solar table
/// by other means.
/// See for example the method TRestAxionSolarFlux::InitializeSolarTable using
/// a solar model description by TRestAxionSolarModel.
///
/// \code
///     <TRestAxionSolarFlux name="sunPrimakoff" verboseLevel="debug" >
///			<parameter name="couplingType" value="g_ag"/>
///			<parameter name="couplingStrength" value="1.e-10"/>
///			<parameter name="fluxDataFile" value="Primakoff_Gianotti_201904.dat"/>
///			<parameter name="fluxSptFile" value="Dummy_Galan_202202.spt"/>
///     </TRestAxionSolarFlux>
/// \endcode
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2022-February: Recovered from original TRestAxionSolarModel implementation
///                Javier Galan
///
/// \class      TRestAxionSolarFlux
/// \author     Javier Galan
///
/// <hr>
///

#include "TRestAxionSolarFlux.h"
using namespace std;

ClassImp(TRestAxionSolarFlux);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionSolarFlux::TRestAxionSolarFlux() : TRestMetadata() { Initialize(); }

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
/// corresponding TRestAxionMagneticField section inside the RML.
///
TRestAxionSolarFlux::TRestAxionSolarFlux(const char* cfgFileName, string name) : TRestMetadata(cfgFileName) {
    cout << "Entering TRestAxionSolarFlux constructor( cfgFileName, name )" << endl;

    Initialize();

    LoadConfigFromFile(fConfigFileName, name);

    if (GetVerboseLevel() >= REST_Info) PrintMetadata();
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionSolarFlux::~TRestAxionSolarFlux() {}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionSolarFlux members
///
void TRestAxionSolarFlux::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    string fullPathName = SearchFile((string)fFluxDataFile);

    debug << "Loading table from file : " << endl;
    debug << "File : " << fullPathName << endl;

    if (fullPathName == "") {
        ferr << "File not found : " << fFluxDataFile << endl;
        ferr << "Solar model table will not be loaded!!" << endl;
    } else {
        TRestTools::ReadASCIITable(fullPathName, fFluxTable);
        fFluxPerRadius.clear();
        for (int n = 0; n < fFluxTable.size(); n++) {
            double sum = 0;
            for (int m = 0; m < fFluxTable[n].size(); m++) {
                cout << fFluxTable[n][m] << "\t";
                sum += fFluxTable[n][m];
            }
            fFluxPerRadius.push_back(sum);
            cout << endl;
            cout << endl;
        }
        cout << endl;
        for (int n = 0; n < fFluxPerRadius.size(); n++)
            cout << "n : " << n << " flux : " << fFluxPerRadius[n] << endl;
        cout << endl;
    }
}

/*
Double_t TRestAxionSolarFlux::GetDifferentialSolarAxionFlux(Double_t energy, Double_t g10) {
    Double_t yearToSeconds = 3600. * 24. * 365.25;
    Double_t m2Tocm2 = 1.e4;

    if (fSolarAxionModel == "arXiv_0702006_Primakoff")
        return 6.02e10 * g10 * g10 * TMath::Power(energy, 2.481) * TMath::Exp(-energy / 1.205);
    if (fSolarAxionModel == "arXiv_1302.6283_Primakoff") {
        Double_t factor = 2.0e22 / yearToSeconds / m2Tocm2;
        return factor * g10 * g10 * TMath::Power(energy, 2.450) * TMath::Exp(-0.829 * energy);
    }
    if (fSolarAxionModel == "arXiv_1302.6283_BC") {
        Double_t factor_1 = 4.2e24 / yearToSeconds / m2Tocm2;
        Double_t factor_2 = 8.3e26 / yearToSeconds / m2Tocm2;

        return g10 * g10 *
               (factor_1 * TMath::Power(energy, 2.987) * TMath::Exp(-0.776 * energy) +
                factor_2 * energy * TMath::Exp(-.77 * energy) / (1 + 0.667 * TMath::Power(energy, 1.278)));
    }

    warning << "Solar model not recognized" << endl;
    warning << "--------------------------" << endl;
    warning << "Solar axion model : " << fSolarAxionModel << endl;

    return 0;
}

Double_t TRestAxionSolarFlux::GetSolarAxionFlux(Double_t eMin, Double_t eMax, Double_t g10, Double_t step) {
    if (fMode == "analytical" && fSolarEnergyFlux > 0)
        if (fg10 == g10 && fStep == step && eMin == fEnergyRange.X() && eMax == fEnergyRange.Y())
            return fSolarEnergyFlux;

    info << "TRestAxionSolarFlux::GetSolarAxionFlux re-calculating solar "
            "axion flux"
         << endl;

    fg10 = g10;
    fStep = step;
    fEnergyRange = TVector2(eMin, eMax);

    fSolarEnergyFlux = 0;
    for (Double_t en = eMin; en < eMax; en += step)
        fSolarEnergyFlux += GetDifferentialSolarAxionFlux(en, g10);

    fSolarEnergyFlux = fSolarEnergyFlux * step;

    return fSolarEnergyFlux;
} */

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionSolarFlux metadata members through a RML file
///
void TRestAxionSolarFlux::InitFromConfigFile() {
    debug << "Entering TRestAxionSolarFlux::InitFromConfigFile" << endl;

    fFluxDataFile = GetParameter("fluxDataFile", "");
    fFluxSptFile = GetParameter("fluxSptFile", "");
    fCouplingType = GetParameter("couplingType", "g_ag");
    fCouplingStrength = StringToDouble(GetParameter("couplingStrength", "1.e-10"));

    this->Initialize();
}

///////////////////////////////////////////////
/// \brief Prints on screen the information about the metadata members of TRestAxionMagneticField
///
void TRestAxionSolarFlux::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    metadata << " - Solar axion flux datafile : " << fFluxDataFile << endl;
    metadata << " - Coupling type : " << fCouplingType << endl;
    metadata << " - Coupling strength : " << fCouplingStrength << endl;
}
