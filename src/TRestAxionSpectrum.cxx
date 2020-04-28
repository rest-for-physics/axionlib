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

/***************** DOXYGEN DOCUMENTATION ********************************
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// \class      TRestAxionSpectrum
/// \author     Sebastian Hoof
///
/// <hr>
///
 *************************************************************************/

// See this header file for more info on the functions.
#include "TRestAxionSpectrum.h"

// Map for named approximations of the spectrum
const std::map<std::string, std::vector<double>> avail_approximations = {
  {"arXiv_0702006_Primakoff", {6.02e10, 1.0e-10, 2.481, 1.0/1.205}},
  {"arXiv_1302.6283_Primakoff", {2.0e22/(3600.0*24.0*365.25*1.0e4), 1.0e-10, 2.450, 0.829}},
  {"arXiv_1302.6283_Compton", {4.2e24/(3600.0*24.0*365.25*1.0e4), 1.0e-13, 2.987, 0.776}}
};

ClassImp(TRestAxionSpectrum);

// Add basic TRest info
void TRestAxionSpectrum::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);
}

void TRestAxionSpectrum::InitFromConfigFile() {
    Initialize();
    sMode = GetParameter("mode");
    if (sMode == "table") {
      sTableFileName = GetParameter("spectrumTableFileName");
      double g1ref = GetDblParameterWithUnits("g1ref", NAN);
      double g2ref = GetDblParameterWithUnits("g2ref", NAN);
      std::string fullPathName = SearchFile((std::string)sTableFileName);
      if (fullPathName == "") {
        ferr << "File not found : " << sTableFileName << endl;
      } else if (std::isnan(g1ref)) {
        ferr << "You need to supply at least one reference value 'g1ref' in 'table' mode of TRestAxionSpectrum." << endl;
      } else {
        fDefaultG1 = g1ref;
        if (not(std::isnan(g2ref))) {
          // N.B. If g2 > 0 has no effect, a warning will be issued by SolarAxionFluxLib.
          spectrum = AxionSpectrum(sTableFileName, g1ref, g2ref);
          fDefaultG2 = g2ref;
        } else {
          // N.B. If g2 is not supplied, we initialise with the default value from SolarAxionFluxLib and issue a warning if the table submode is wrong.
          auto table_params = spectrum.get_table_parameters();
          int table_submode = get<0>(table_params);
          if (table_submode > 1) {
            fDefaultG2 = get<2>(table_params);
            cout << "WARNING! Your table for TRestAxionSpectrum supports two separate spectra, but you did not specify a value for coupling 'g2'.\n\
                              TRestAxionSpectrum will assume a default value of " << fDefaultG2 << " (in appropriate units)." << endl; };
          spectrum = AxionSpectrum(sTableFileName, g1ref);
        };
      };
    } else if (sMode == "analytical") {
      std::string named_approx = GetParameter("named_approx", "n/a");
      double norm = GetDblParameterWithUnits("norm", NAN);
      double g1ref = GetDblParameterWithUnits("g1ref", NAN);
      double a = GetDblParameterWithUnits("a", NAN);
      double b = GetDblParameterWithUnits("b", NAN);
      bool check_named_approx = not(avail_approximations.find(named_approx) == avail_approximations.end());
      bool check_numbers = not(std::isnan(norm) || std::isnan(a) || std::isnan(b) || std::isnan(g1ref));
      if (check_named_approx) {
        std::vector<double> p = avail_approximations.at(named_approx);
        spectrum = AxionSpectrum(p[0], p[1], p[2], p[3]);
        fDefaultG1 = p[1];
      } else if (check_numbers) {
        spectrum = AxionSpectrum(norm, g1ref, a, b);
        fDefaultG1 = g1ref;
      } else {
        std::string avail_approximation_names = "";
        for (auto el : avail_approximations) { avail_approximation_names += " " + el.first; };
        ferr << "You want to use the 'analytical' mode of TRestAxionSpectrum but neither supplied a known 'named_approx' OR the four required parameters 'a', 'b', 'norm', "
                "and 'g1ref'.\n The available known approximations are:"+avail_approximation_names << endl;
      }
    } else if (sMode == "solar_model") {
      // TODO: Ideally use TRestAxionSolarModel here...
      std::string sSolarModelFile = GetParameter("solarAxionModel");
      std::string fullPathName = SearchFile((std::string)sSolarModelFile);
      sol = SolarModel(fullPathName, OP, false);
      fDefaultG1 = sol.get_gagg_ref_value_in_inverse_GeV();
      fDefaultG2 = sol.get_gaee_ref_value();
      spectrum = AxionSpectrum(&sol);
    } else {
      ferr << "Mode for TRestAxionSpectrum not known! Choose one of 'table', 'analytical', and 'solar_model'." << endl;
      sMode = "none";
    };
}

TRestAxionSpectrum::TRestAxionSpectrum() : TRestMetadata() { Initialize(); }

TRestAxionSpectrum::TRestAxionSpectrum(const char* cfgFileName, std::string name) : TRestMetadata(cfgFileName) {
  debug << "Creating instance of TRestAxionSpectrum from file "+fConfigFileName+"..." << endl;
  Initialize();
  LoadConfigFromFile(fConfigFileName, name);
  PrintMetadata();
}

TRestAxionSpectrum::~TRestAxionSpectrum() { }

void TRestAxionSpectrum::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    if (sMode != "none") {
      metadata << " Solar spectrum created in mode '" << sMode << "'." << endl;
    } else {
      metadata << " Solar spectrum has not yet been initialised, yet!" << endl;
    };
    if (sMode == "solar_model") {
      metadata << " - Opacity code used : " << sol.get_solaxlib_name_and_version() << endl;
      metadata << " - Solar model file : " << sol.get_solar_model_name() << endl;
      metadata << " - Opacity code used : " << sol.get_opacitycode_name() << endl;
    } else if (sMode == "table") {
      metadata << " - Tabulted spectrum file used : " << sTableFileName << endl;
    };
    metadata << "-------------------------------------------------" << endl;
    metadata << " - Units of the solar axion flux from this class : axions / cm^2 s keV" << endl;
    if (not(std::isnan(fDefaultG1))) { metadata << " - Numerical value of coupling g1 (in appropriate units): " << fDefaultG1 << endl; };
    if (fDefaultG2 > 0) {
      metadata << " - Numerical value of coupling g2 (in appropriate units): " << fDefaultG2 << endl;
    } else {
      metadata << " - A secong coupling, g2, is not available in this class instance." << endl;
    };
    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

double TRestAxionSpectrum::GetSolarAxionFlux(double erg_lo, double erg_hi, double erg_step_size) { return 0; }
double TRestAxionSpectrum::GetDifferentialSolarAxionFlux(double erg) { return spectrum.axion_flux(erg, fDefaultG1, fDefaultG2); }
