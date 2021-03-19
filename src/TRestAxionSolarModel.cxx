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
/// The TRestAxionSolarModel class is for calculating the solar axion
/// spectrum from first principles i.e. based on a solar model and an
/// appropriate opacity code. For this task it currently wraps the
/// external SolarAxionFlux library available at the following Github
/// repository: https://github.com/sebhoof/SolarAxionFlux. The resulting
/// spectra, fundametally arising from either ALP-electron or ALP-photon
/// interactions can then be used by e.g. the TRestAxionSolarSpectrum class in
/// other parts of the analysis.
///
/// The following piece of code shows how to set up an instance of
/// of the class with a the B16-AGSS09 solar model and OP opacity code.
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
/// 2020-April: TRestAxionSolarModel class now wrapping the SolarAxionFlux external library.
///             Sebastian Hoof
/// 2021-March: Overhaul of the class according to new strategy for the workflow.
///             Sebastian Hoof
///
/// \class      TRestAxionSolarModel
/// \author     Javier Galan, Sebastian Hoof
///
/// <hr>
///
 *************************************************************************/

#include "TRestAxionSolarModel.h"

//// Some local auxilliary functions
std::string format_number(double number, int precision = 5) {
    std::stringstream output;
    output << std::scientific << std::setprecision(precision) << number;
    return output.str();
}

ClassImp(TRestAxionSolarModel);

// Add basic TRest info
void TRestAxionSolarModel::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);
}

void TRestAxionSolarModel::InitFromConfigFile() {
    Initialize();
    fSolarModelFile = GetParameter("solarAxionModel");
    std::string fullPathName = SearchFile((std::string)fSolarModelFile);
    if (fullPathName == "") {
        ferr << "File not found : " << fSolarModelFile << endl;
    } else {
        fSol = SolarModel(fullPathName, OP, false);
        fExternalLibraryNameAndVersion = fSol.get_solaxlib_name_and_version();
        fOpacityCodeName = fSol.get_opacitycode_name();
        fRefPhotonCoupling = fSol.get_gagg_ref_value_in_inverse_GeV();
        fRefElectronCoupling = fSol.get_gaee_ref_value();
        fSolarModelInitialized = fSol.is_initialised();
        if (fSolarModelInitialized) {
            debug << "Solar model file " << fSolarModelFile << " successfully loaded." << endl;
        } else {
            ferr << "Solar model initialization was not successful!" << endl;
        }
    }
}

// Default constructor
TRestAxionSolarModel::TRestAxionSolarModel() : TRestMetadata() { Initialize(); }

// From-file contructor
TRestAxionSolarModel::TRestAxionSolarModel(const char* cfgFileName, std::string name)
    : TRestMetadata(cfgFileName) {
    debug << "Creating instance of TRestAxionSolarModel from file " + fConfigFileName + "..." << endl;
    Initialize();
    LoadConfigFromFile(fConfigFileName, name);
    PrintMetadata();
}

TRestAxionSolarModel::~TRestAxionSolarModel() {}  // SolarModel memory in fSol will be cleaned up
                                                  // automatically when fSol leaves scope

void TRestAxionSolarModel::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    metadata << " Solar model created with " << fExternalLibraryNameAndVersion << "." << endl;
    metadata << " - Solar model file : " << fSolarModelFile << endl;
    metadata << " - Opacity code used : " << fOpacityCodeName << endl;
    metadata << "-------------------------------------------------" << endl;
    metadata << " - Reference value of the axion-photon coupling : " << fRefPhotonCoupling / 1.0e-10
             << " x 10^{-10} / GeV" << endl;
    metadata << " - Reference value of the axion-electron coupling : " << fRefElectronCoupling / 1.0e-13
             << " x 10^{-13}" << endl;
    metadata << " - Units of the solar axion flux from this class : axions / cm^2 s keV" << endl;
    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}


std::vector<double> TRestAxionSolarModel::CalcSpectralSolarAxionFluxPrimakoff(std::vector<double> energies, double r_max, std::string outputfile, double g_agamma) {
    // Rescale the flux by the value of the coupling chosen by the user
    double rescaling = g_agamma/fRefPhotonCoupling;
    rescaling *= rescaling;

    std::vector<double> result = calculate_spectral_flux_Primakoff(energies, fSol, r_max);
    for (auto flux = result.begin(); flux != result.end(); flux++) { *flux *= rescaling; }

    std::vector<std::vector<double> > buffer = { energies, result };
    std::string header = "Spectral solar axion Primakoff flux dPhi/dE for g_agamma = "+format_number(g_agamma)
                         +" GeV^-1, calculated by "+fExternalLibraryNameAndVersion+".\nColumns:"
                         "Axion energy [keV] | Primakoff flux [axions / cm^2 s keV]";
    // Function from SolarAxionFlux lib that saves to file (with header) only if outputfile != ""
    save_to_file(outputfile, buffer, header);

    return result;
}

std::vector<double> TRestAxionSolarModel::CalcSpectralSolarAxionFluxABC(std::vector<double> energies, double r_max, std::string outputfile, double g_ae) {
    // Rescale the flux by the value of the coupling chosen by the user
    double rescaling = g_ae/fRefElectronCoupling;
    rescaling *= rescaling;

    std::vector<double> result = calculate_spectral_flux_axionelectron(energies, fSol, r_max);
    for (auto flux = result.begin(); flux != result.end(); flux++) { *flux *= rescaling; }

    std::vector<std::vector<double> > buffer = { energies, result };
    std::string header = "Spectral solar axion ABC flux dPhi/dE for g_ae = "+format_number(g_ae)
                         +", calculated by "+fExternalLibraryNameAndVersion+".\nColumns: "
                         "Axion energy [keV] | ABC flux [axions / cm^2 s keV]";
    // Function from SolarAxionFlux lib that saves to file (with header) only if outputfile != ""
    save_to_file(outputfile, buffer, header);

    return result;
}

std::vector<std::vector<double> > TRestAxionSolarModel::CalcSpectralSolarAxionFluxAll(std::vector<double> energies, double r_max, std::string outputfile, double g_agamma, double g_ae) {
    // Calculate fluxes for all individual components
    std::vector<double> flux_gag = CalcSpectralSolarAxionFluxPrimakoff(energies, r_max, "", g_agamma);
    std::vector<double> flux_gae = CalcSpectralSolarAxionFluxABC(energies, r_max, "", g_ae);

    std::vector<std::vector<double> > result = { energies, flux_gag, flux_gae };
    std::string header = "Spectral solar axion fluxes dPhi/dE for for g_agamma = "+format_number(g_agamma)
                         +" GeV^-1 and g_ae = "+format_number(g_ae)+", calculated by "
                         +fExternalLibraryNameAndVersion+".\nColumns: Axion energy [keV] | "
                         "Primakoff flux [axions / cm^2 s keV] | ABC flux [axions / cm^2 s keV]";
    // Function from SolarAxionFlux lib that saves to file (with header) only if outputfile != ""
    save_to_file(outputfile, result, header);

    return result;
}

std::vector<std::vector<double> > TRestAxionSolarModel::CalcSpectralAndSpatialSolarAxionFluxPrimakoff(std::vector<double> energies, std::vector<double> radii, std::string outputfile, double g_agamma) {
    // Rescale the flux by the value of the coupling chosen by the user
    double rescaling = g_agamma/fRefPhotonCoupling;
    rescaling *= rescaling;

    std::vector<std::vector<double> > result = calculate_spectral_flux_Primakoff(energies, radii, fSol);
    for (auto flux = result.begin(); flux != result.end(); flux++) { (*flux)[2] *= rescaling; }

    std::string header = "Spectral and spatial solar axion Primakoff flux d^2Phi/dEdr on the solar "
                         "disc for g_agamma = "+format_number(g_agamma)
                         +" GeV^-1, calculated by "+fExternalLibraryNameAndVersion+".\nColumns: "
                         "Radius on the solar disc [R_sol] | Axion energy [keV] | Primakoff flux [axions / cm^2 s keV]";
    // Function from SolarAxionFlux lib that saves to file (with header) only if outputfile != ""
    save_to_file(outputfile, result, header);

    return result;
}

std::vector<std::vector<double> > TRestAxionSolarModel::CalcSpectralAndSpatialSolarAxionFluxABC(std::vector<double> energies, std::vector<double> radii, std::string outputfile, double g_ae) {
    // Rescale the flux by the value of the coupling chosen by the user
    double rescaling = g_ae/fRefElectronCoupling;
    rescaling *= rescaling;

    std::vector<std::vector<double> > result = calculate_spectral_flux_axionelectron(energies, radii, fSol);
    for (auto flux = result.begin(); flux != result.end(); flux++) { (*flux)[2] *= rescaling; }

    std::string header = "Spectral and spatial solar axion ABC flux d^2Phi/dEdr on the solar "
                         "disc for g_ae = "+format_number(g_ae)
                         +", calculated by "+fExternalLibraryNameAndVersion+".\nColumns: "
                         "Radius on the solar disc [R_sol] | Axion energy [keV] | ABC flux [axions / cm^2 s keV]";
    // Function from SolarAxionFlux lib that saves to file (with header) only if outputfile != ""
    save_to_file(outputfile, result, header);

    return result;
}

std::vector<std::vector<double> > TRestAxionSolarModel::CalcSpectralAndSpatialSolarAxionFluxAll(std::vector<double> energies, std::vector<double> radii, std::string outputfile, double g_agamma, double g_ae) {
    std::vector<std::vector<double> > flux_gag = CalcSpectralAndSpatialSolarAxionFluxPrimakoff(energies, radii, "", g_agamma);
    std::vector<std::vector<double> > flux_gae = CalcSpectralAndSpatialSolarAxionFluxABC(energies, radii, "", g_ae);

    std::vector<std::vector<double> > result = { flux_gag[0], flux_gag[1], flux_gag[2], flux_gae[2] };
    std::string header = "Spectral and spatial solar axion fluxes d^2Phi/dEdr on the solar "
                         "disc for g_agamma = "+format_number(g_agamma)+" GeV^-1, g_ae = "
                         +format_number(g_ae)+", calculated by "+fExternalLibraryNameAndVersion+"."
                         "\nColumns: Radius on the solar disc [R_sol] | Axion energy [keV] | "
                         "Primakoff flux [axions / cm^2 s keV] | ABC flux [axions / cm^2 s keV]";
    // Function from SolarAxionFlux lib that saves to file (with header) only if outputfile != ""
    save_to_file(outputfile, result, header);

    return result;
}
