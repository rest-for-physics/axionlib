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
/// TRestAxionSolarFluxs is a class used to define methods which define
/// theoretical model results that can be directly used in the calculation
/// of experimental observables.
///
/// TODO. Create an appropriate documentation here.
///
/// There will be two modes of defining the solar mode in this class:
/// *analytical* and *table*.
///
/// The following piece of code shows how to define an analytical solar model.
///
/// \code
///     <TRestAxionSolarFlux name="sunPrimakoff" verboseLevel="debug" >
///    <parameter name="mode" value="analytical" />
///    <parameter name="solarAxionSolarModel" value="arXiv_0702006_Primakoff" />
/// </TRestAxionSolarFlux>
/// \endcode
///
/// The available analytical solar axion models, for different production
/// mechanisms, is given in the following list.
///
/// - arXiv_0702006_Primakoff (https://arxiv.org/abs/hep-ex/0702006)
/// - arXiv_1302.6283_Primakoff ( https://arxiv.org/pdf/1302.6283v2.pdf )
/// - arXiv_1302.6283_BC ( https://arxiv.org/pdf/1302.6283v2.pdf )
///
/// The second mode, *table*, will provide further detail on the solar axion
/// production as a function of the solar radius. The tables will be
/// available as a file inside the data/solarModel/ directory. The
/// *solarAxionSolarModel* parameter
///
/// \code
///     <TRestAxionSolarFlux name="sunPrimakoff" verboseLevel="debug" >
///         <parameter name="mode" value="table" />
///         <parameter name="solarAxionSolarModel" value="arXiv_0702006_Primakoff" />
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
///             Javier Galan
///
/// \class      TRestAxionSolarFlux
/// \author     Javier Galan
///
/// <hr>
///

#include "TRestAxionSolarFlux.h"
using namespace std;

ClassImp(TRestAxionSolarFlux);

TRestAxionSolarFlux::TRestAxionSolarFlux() : TRestMetadata() {
    // TRestAxionSolarFlux default constructor
    Initialize();
}

TRestAxionSolarFlux::TRestAxionSolarFlux(const char* cfgFileName, string name) : TRestMetadata(cfgFileName) {
    cout << "Entering TRestAxionSolarFlux constructor( cfgFileName, name )" << endl;

    Initialize();

    LoadConfigFromFile(fConfigFileName, name);

    PrintMetadata();
}

TRestAxionSolarFlux::~TRestAxionSolarFlux() {}

void TRestAxionSolarFlux::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);
}

// Returns the solar axion flux in cm-2 keV-1 s-1 (on earth)
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
}

void TRestAxionSolarFlux::InitFromConfigFile() {
    debug << "Entering TRestAxionSolarFlux::InitFromConfigFile" << endl;

    this->Initialize();

    fMode = GetParameter("mode", "analytical");
    fSolarAxionModel = GetParameter("solarAxionModel", "arXiv_0702006_Primakoff");

    if (fMode == "table") {
        fStep = 0.1;
        debug << "Loading table from file : " << endl;
        string fullPathName = SearchFile((string)fSolarAxionModel);

        debug << "File : " << fullPathName << endl;

        if (fullPathName == "") {
            ferr << "File not found : " << fSolarAxionModel << endl;
            ferr << "Solar model table will not be loaded!!" << endl;
        } else {
            TRestTools::ReadASCIITable(fullPathName, fSolarTable);
            for (int n = 0; n < fSolarTable.size(); n++) {
                for (int m = 0; m < fSolarTable[n].size(); m++) cout << fSolarTable[n][m] << "\t";
                cout << endl;
                cout << endl;
            }
        }
    }

    PrintMetadata();
}

void TRestAxionSolarFlux::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    metadata << " - Mode : " << fMode << endl;
    metadata << " - Solar axion model : " << fSolarAxionModel << endl;
    metadata << "-------------------------------------------------" << endl;
    metadata << " - Axion-photon couping : " << fg10 << " x 10^{-10} GeV^{-1}" << endl;
    metadata << " - Integration step : " << fStep << " keV" << endl;
    metadata << " - Integration range : ( " << fEnergyRange.X() << ", " << fEnergyRange.Y() << " ) keV"
             << endl;
    metadata << " - Calculated solar flux : " << fSolarEnergyFlux << " cm-2 s-1" << endl;
    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}
