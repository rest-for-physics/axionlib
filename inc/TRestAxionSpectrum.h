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

#include "solaxflux/spectral_flux.hpp"

#include <TRestMetadata.h>

//! A metadata class to define a solar axion spectrum and functions to evaluate it.
class TRestAxionSpectrum : public TRestMetadata {
private:
    // Generic initialization routines.
    void Initialize();
    // The mode of this class. There are 3 different modes available: table, analytic, and solar_model.
    // Energy in keV, flux in axions / cm^2 s keV
    // table      : provide a text file with up to 4 columns (and reference values of the associated couplings g1 and g2).
    //              - 2 columns: energy | flux for coupling g1
    //              - 3 columns: energy | flux for coupling g1 | flux for coupling g2
    //              - 4 columns: energy | radius on solar disc | flux for coupling g1 | flux for coupling g2
    // analytical : provide parameters a, b, norm, g1ref for the following ansatz OR a named set of parameters available in the data/ folder
    //              flux = norm * (g1/g1ref)^2 * energy^a * exp(-b * energy)
    // solar_model: provide a text file with a solar model
    std::string sMode = "none";
    // InitFromConfigFile() generates an (interpolated) axion spectrum from a config file.
    void InitFromConfigFile();
    // An instance of the axion spectrum class from the AxionSolarFlux lib.
    AxionSpectrum spectrum;
    // Solar model if this needs to be defined locally.
    SolarModel sol;

    double fDefaultG1 = NAN;
    double fDefaultG2 = 0;
    std::string sTableFileName;

    //TString fProcessName;
    //TString fMetaDataFromFileHeader;
    //std::vector<std::vector<double>> fSpectrumTable;

public:
    // Constructors
    TRestAxionSpectrum();
    TRestAxionSpectrum(const char *cfgFileName, std::string name = "");
    ~TRestAxionSpectrum();

    double GetSolarAxionFlux(double erg_lo, double erg_hi, double er_step_size);
    double GetDifferentialSolarAxionFlux(double erg);

    void PrintMetadata();
    // bool isSpectrumTableLoaded() { return fSpectrumTable.size() > 0; }
    // TString GetProcessName() { return fProcessName; }

    ClassDef(TRestAxionSpectrum, 1);
};

#endif
