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

#ifndef _TRestAxionSpectrum
#define _TRestAxionSpectrum

#include "solaxflux/solar_model.hpp"

#include <TRestMetadata.h>

//! A metadata class to define a solar axion spectrum and functions to evaluate it.
class TRestAxionSpectrum : public TRestMetadata {
private:
    void Initialize();

    void InitFromConfigFile();

    TString fProcessName;

    TString fMetaDataFromFileHeader;

    std::vector<std::vector<double>> fSpectrumTable;  //->

public:
    bool isSpectrumTableLoaded() { return fSpectrumTable.size() > 0; }

    TString GetProcessName() { return fProcessName; }
    double GetSolarAxionFlux(double erg_lo, double erg_hi, double er_step_size);
    double GetDifferentialSolarAxionFlux(double erg);

    void PrintMetadata();

    // Constructors
    TRestAxionSpectrum();
    TRestAxionSpectrum(const char *cfgFileName, std::string name = "");
    // Destructor
    ~TRestAxionSpectrum();

    ClassDef(TRestAxionSpectrum, 1);
};

#endif

/*


/// The available analytical solar axion models, for different production
/// mechanisms, is given in the following list.
///
/// - arXiv_0702006_Primakoff (https://arxiv.org/abs/hep-ex/0702006)
/// - arXiv_1302.6283_Primakoff ( https://arxiv.org/pdf/1302.6283v2.pdf )
/// - arXiv_1302.6283_BC ( https://arxiv.org/pdf/1302.6283v2.pdf )
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

*/
