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
/// TRestAxionSolarModel is a class used to calculate all axion interaction
/// rates in the Sun based on a solar model and opacity code. The main
/// purpose of the class is to return the two types of solar axion spectra,
/// i.e. from axion-photon and axion eletron interactions.
///
/// TODO. Create an appropriate documentation here.
///
/// The following piece of code shows how to define an analytical solar model.
///
/// \code
/// <TRestAxionSolarModel name="sunPrimakoff" verboseLevel="debug" >
///    <parameter name="mode" value="lib" />
///    <parameter name="solarAxionModel" value="SolarModel_B16-AGSS09.dat" />
/// </TRestAxionSolarModel>
/// \endcode
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2019-March: First concept and implementation of TRestAxionSolarModel class.
///             Javier Galan
/// 2020-April: Re-wrote TRestAxionSolarModel class as a wrapper for the SolarAxionFlux external library
///             Sebastian Hoof
///
/// \class      TRestAxionSolarModel
/// \author     Javier Galan, Sebastian Hoof
///
/// <hr>
///
 *************************************************************************/

#include "TRestAxionSolarModel.h"

ClassImp(TRestAxionSolarModel);

// Add basic TRest info
void TRestAxionSolarModel::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);
}

void TRestAxionSolarModel::InitFromConfigFile() {
    Initialize();
    sSolarModelFile = GetParameter("solarAxionModel", "SolarModel_B16-AGSS09.dat");
    std::string fullPathName = SearchFile((std::string)sSolarModelFile);
    if (fullPathName == "") {
        RESTError << "File not found : " << sSolarModelFile << RESTendl;
    } else {
#ifdef USE_SolaxFlux
        sol = SolarModel(fullPathName, OP, false);
        sExternalLibraryName = sol.get_solaxlib_name_and_version();
        sOpacityCodeName = sol.get_opacitycode_name();
        fRefPhotonCoupling = sol.get_gagg_ref_value_in_inverse_GeV();
        fRefElectronCoupling = sol.get_gaee_ref_value();
        bSolarModelInitialized = sol.is_initialised();
        if (bSolarModelInitialized) {
            RESTDebug << "Solar model file " << sSolarModelFile << " successfully loaded!" << RESTendl;
        } else {
            RESTError << "Solar model initialization was not successful!" << RESTendl;
        };
#endif
    };
}

// Default constructor
TRestAxionSolarModel::TRestAxionSolarModel() : TRestMetadata() { Initialize(); }

// From-file contructor
TRestAxionSolarModel::TRestAxionSolarModel(const char* cfgFileName, std::string name)
    : TRestMetadata(cfgFileName) {
    RESTDebug << "Creating instance of TRestAxionSolarModel from file " + fConfigFileName + "..." << RESTendl;
    Initialize();
    LoadConfigFromFile(fConfigFileName, name);
    PrintMetadata();
}

TRestAxionSolarModel::~TRestAxionSolarModel() {}  // SolarModel memory in sol will get destroyed automatically

void TRestAxionSolarModel::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    RESTMetadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << RESTendl;
    RESTMetadata << " Solar model created with " << sExternalLibraryName << "." << RESTendl;
    RESTMetadata << " - Solar model file : " << sSolarModelFile << RESTendl;
    RESTMetadata << " - Opacity code used : " << sOpacityCodeName << RESTendl;
    RESTMetadata << "-------------------------------------------------" << RESTendl;
    RESTMetadata << " - Reference value of the axion-photon coupling : " << fRefPhotonCoupling / 1.0e-10
                 << " x 10^{-10} / GeV" << RESTendl;
    RESTMetadata << " - Reference value of the axion-electron coupling : " << fRefElectronCoupling / 1.0e-13
                 << " x 10^{-13}" << RESTendl;
    RESTMetadata << " - Units of the solar axion flux from this class : axions / cm^2 s keV" << RESTendl;
    RESTMetadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << RESTendl;
}

std::string TRestAxionSolarModel::GetSolarModelFileName() { return sSolarModelFile; }

std::vector<double> TRestAxionSolarModel::GetSolarAxionFluxGAGamma(std::vector<double> energies,
                                                                   double r_max) {
    std::vector<double> result;
    if (bSolarModelInitialized) {
        RESTWarning << "TRestAxionSolarModel::GetSolarAxionFluxGAGamma." << RESTendl;
        RESTWarning << "This code has been commented to allow compilation and needs to be reviewed!"
                    << RESTendl;
        RESTWarning << "The result vector will be empty!" << RESTendl;
        // result = sol.calculate_spectral_flux_Primakoff(energies, r_max);
    } else {
        RESTError << "TRestAxionSolarModel not properly initialised for "
                     "RestAxionSolarModel::GetSolarAxionFluxGAGamma(...)!"
                  << RESTendl;
    };
    return result;
}

std::vector<double> TRestAxionSolarModel::GetSolarAxionFluxGAGamma(std::vector<double> energies,
                                                                   double g_agamma, double r_max) {
    std::vector<double> result = GetSolarAxionFluxGAGamma(energies);
    for (auto flux = result.begin(); flux != result.end(); flux++) {
        *flux *= pow(g_agamma / fRefPhotonCoupling, 2);
    };
    return result;
}

std::vector<double> TRestAxionSolarModel::GetSolarAxionFluxGAE(std::vector<double> energies, double r_max) {
    std::vector<double> result;
    if (bSolarModelInitialized) {
        RESTWarning << "TRestAxionSolarModel::GetSolarAxionFluxGAE." << RESTendl;
        RESTWarning << "This code has been commented to allow compilation and needs to be reviewed!"
                    << RESTendl;
        RESTWarning << "The result vector will be empty!" << RESTendl;
        // result = sol.calculate_spectral_flux_all_electron(energies, r_max);
    } else {
        RESTError << "TRestAxionSolarModel not properly initialised for "
                     "RestAxionSolarModel::GetSolarAxionFluxGAE(...)!"
                  << RESTendl;
    };
    return result;
}

std::vector<double> TRestAxionSolarModel::GetSolarAxionFluxGAE(std::vector<double> energies, double g_agae,
                                                               double r_max) {
    std::vector<double> result = GetSolarAxionFluxGAE(energies, r_max);
    for (auto flux = result.begin(); flux != result.end(); flux++) {
        *flux *= pow(g_agae / fRefElectronCoupling, 2);
    };
    return result;
}
