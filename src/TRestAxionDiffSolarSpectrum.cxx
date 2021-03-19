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
/// TODO: Not all features fully implemented, yet!
/// The TRestAxionDiffSolarSpectrum serves as a way to read in an perform calculations
/// with solar axion fluxes calculated and/or tabulated elsewhere e.g. from the
/// external SolarAxionFlux library (https://github.com/sebhoof/SolarAxionFlux)
/// via the TRestAxionSolarModel class. The class can handle (i) analytical
/// expressions and (ii) tabulated spectra. Other than interpolation and
/// (exact) integration, it can perform Monte Carlo draws of axion events.
/// The default units for the differential fluxes are axions / cm^2 s keV.
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2021-March: Restructured the existing code as standalone code i.e. independent from SolarAxionFlux
///             library and TRestAxionSolarModel. Introduced the notion of differential and integrated
///             spectra.
///             Sebastian Hoof
///
/// \class      TRestAxionDiffSolarSpectrum
/// \author     Sebastian Hoof
///
/// <hr>
///
 *************************************************************************/

// See this header file for more info on the functions.
#include "TRestAxionSolarSpectrum.h"
#include "TRestAxionDiffSolarSpectrum.h"

ClassImp(TRestAxionDiffSolarSpectrum);

// Add basic TRest info
void TRestAxionDiffSolarSpectrum::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);
}

void TRestAxionDiffSolarSpectrum::Init2DSplines(const int n_grids) {
    fGSLAccel2D.resize(n_grids);
    fGSLSpline2D.resize(n_grids);
    fGSLMem2D.resize(n_grids);

    std::vector<double> x_unique = sort_unique_values(fData[0]);
    std::vector<double> y_unique = sort_unique_values(fData[1]);
    int nx = x_unique.size();
    int ny = y_unique.size();
    int pts = fData[0].size();
    if (nx*ny != pts) {
        ferr << "Number of data points ("+std::to_string(pts)+") in '"+fTableFileName+"' inconsistent "
                "with the number of unique 'x' and 'y' values ("+std::to_string(nx)+" and "
                +std::to_string(ny)+")! Check the formatting of the file." << endl;
    }
    // Determine first and last "x" and "y" values and grid step size.
    double x_lo = x_unique.front();
    double x_up = x_unique.back();
    double x_delta = (x_up-x_lo) / (nx-1);
    double y_lo = y_unique.front();
    double y_up = y_unique.back();
    double y_delta = (y_up-y_lo) / (ny-1);

    const double* x = &x_unique[0];
    const double* y = &y_unique[0];
    for (int i=0; i<n_grids; ++i) {
        fGSLMem2D[i] = (double*) malloc(nx * ny * sizeof(double));
        fGSLSpline2D[i] = gsl_spline2d_alloc(gsl_interp2d_bilinear, nx, ny);
        fGSLAccel2D[i].first = gsl_interp_accel_alloc();
        fGSLAccel2D[i].second = gsl_interp_accel_alloc();
    }

    // Intialise grids
    for (int j=0; j<pts; ++j) {
        // Determine appropriate indices for the grid points.
        double temp = (fData[0][j]-x_lo) / x_delta;
        int ind_x = (int) (temp+0.5);
        temp = (fData[1][j]-y_lo) / y_delta;
        int ind_y = (int) (temp+0.5);
        for (int i=0; i<n_grids; ++i) { gsl_spline2d_set(fGSLSpline2D[i], fGSLMem2D[i], ind_x, ind_y, fData[i+2][j]); }
    }

    for (int i=0; i<n_grids; ++i) { gsl_spline2d_init(fGSLSpline2D[i], x, y, fGSLMem2D[i], nx, ny); }
}

void TRestAxionDiffSolarSpectrum::InitFromConfigFile() {
    Initialize();

    fTableFileName = GetParameter("spectrumTableFileName");
    fRefPhotonCoupling = GetDblParameter("g_agamma", NAN);
    fRefElectronCoupling = GetDblParameter("g_ae", NAN);
    std::cout << GetDblParameter("g_agamma", NAN) << endl;
    std::cout << fRefPhotonCoupling << " " << fRefElectronCoupling << std::endl;
    bool g_agamma_nan = std::isnan(fRefPhotonCoupling) || not(fRefPhotonCoupling > 0);
    bool g_ae_nan = std::isnan(fRefElectronCoupling) || not(fRefElectronCoupling > 0);
    fTableSubMode = not(g_agamma_nan) + 2*not(g_ae_nan);

    fTableFileName = SearchFile((std::string)fTableFileName);
    if (fTableFileName == "") {
        ferr << "File not found : " << fTableFileName << endl;
    } else if (fTableSubMode == 0) {
        ferr << "You need to supply a reference value for either the axion-photon coupling "
                "'g_agamma' (in units of GeV^{-1}) or the axion-electron coupling 'g_ae' when "
                "initializing TRestAxionDiffSolarSpectrum."
             << endl;
    } else {
        // TRestTools::ReadASCIITable(fTableFileName, fData); -> unnatural rows vs columns
        read_natural_ascii_table(fTableFileName, fData);
        int n_cols = fData.size();
        if (n_cols < 3) {
            ferr << "Tabulated solar axion spectrum " << fTableFileName << " contains less than "
                    "three columns; this is incompatible with TRestAxionDiffSolarSpectrum." << endl;
        } else if (n_cols==3) {
            if (fTableSubMode < 3) {
                // Intialise one 2D array
                Init2DSplines(1);
            } else {
                ferr << "Tabulated solar axion spectrum " << fTableFileName << " contains three columns "
                        "but you supplied two non-NAN, non-zero values for the couplings 'g_agamma' "
                        "and 'g_ae'; this is incompatible with TRestAxionDiffSolarSpectrum."
                     << endl;
            }
        } else if (n_cols==4) {
            if (fTableSubMode==3) {
                // Initialise two 2D arrays
                Init2DSplines(2);
            } else if (fTableSubMode>=1) {
                // Initialise one 2D array
                Init2DSplines(1);
                debug << "Tabulated solar axion spectrum " << fTableFileName << " contains more than "
                         "three columns for one non-NAN, non-zero coupling. Excess columns will be "
                         "ignored." << endl;
            } else {
                ferr << "Tabulated solar axion spectrum " << fTableFileName << " contains four columns "
                        "but you supplied less than two non-NAN values for the couplings 'g_agamma' "
                        "and 'g_ae'; this is incompatible with TRestAxionDiffSolarSpectrum."
                     << endl;
            }
        } else {
            debug << "Tabulated solar axion spectrum " << fTableFileName << " contains more than "
                     "four columns. Excess columns will be ignored." << endl;
        }
    }
}

TRestAxionDiffSolarSpectrum::TRestAxionDiffSolarSpectrum() : TRestMetadata() { Initialize(); }

TRestAxionDiffSolarSpectrum::TRestAxionDiffSolarSpectrum(const char* cfgFileName, std::string name)
    : TRestMetadata(cfgFileName) {
    debug << "Creating instance of TRestAxionDiffSolarSpectrum from file " + fConfigFileName + "..." << endl;
    Initialize();
    LoadConfigFromFile(fConfigFileName, name);
    PrintMetadata();
}

// Need to free allocated memory by GSL routines!
TRestAxionDiffSolarSpectrum::~TRestAxionDiffSolarSpectrum() {
    for (auto s : fGSLSpline2D) { gsl_spline2d_free(s); }
    for (auto a : fGSLAccel2D) { gsl_interp_accel_free(a.first); gsl_interp_accel_free(a.second); }
    for (auto z : fGSLMem2D) { free(z); }
}

void TRestAxionDiffSolarSpectrum::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    metadata << " Instance of TRestAxionDiffSolarSpectrum." << endl;
    metadata << " - Tabulated spectrum file used : "
             << TRestTools::SeparatePathAndName(fTableFileName).second << endl;
    std::string sp_init = IsSpectrumReady() ? "yes" : "no";
    metadata << " - Class properly initialized : " << sp_init << endl;
    metadata << "-------------------------------------------------" << endl;
    metadata << " - Units of the solar axion flux from this class : axions / cm^2 s keV" << endl;
    if (not(std::isnan(fRefPhotonCoupling))) {
        metadata << " - Numerical value of coupling g_agamma : " << fRefPhotonCoupling*1.0e10 <<
        " x 10^{-10} GeV^{-1}."<< endl;
    } else {
      metadata << " - Axion-photon interactions (via g_agamma) are not considered in this spectrum." << endl;
    }
    if (not(std::isnan(fRefElectronCoupling))) {
        metadata << " - Numerical value of coupling g_ae : " << fRefElectronCoupling*1.0e13 <<
        " x 10^{-13}." << endl;
    } else {
        metadata << " - Axion-electron interactions (via g_ae) are not considered in this spectrum." << endl;
    }
    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

double TRestAxionDiffSolarSpectrum::GetDoubleDifferentialSolarAxionFlux(double r, double erg, double g_agamma, double g_ae) {
    double x_agamma2, x_ae2, result = 0;
    if (fTableSubMode % 2 == 1) { x_agamma2 = g_agamma/fRefPhotonCoupling; x_agamma2 *= x_agamma2; }
    if (fTableSubMode % 2 == 0) { x_ae2 = g_ae/fRefElectronCoupling; x_ae2 *= x_ae2; }
    if (fTableSubMode == 1) {
      std::cout << "Test" << std::endl;
        result = x_agamma2*gsl_spline2d_eval(fGSLSpline2D[0], r, erg, fGSLAccel2D[0].first, fGSLAccel2D[0].second);
    } else if (fTableSubMode == 2) {
        result = x_ae2*gsl_spline2d_eval(fGSLSpline2D[0], r, erg, fGSLAccel2D[0].first, fGSLAccel2D[0].second);
    } else if (fTableSubMode == 3) {
        result = x_agamma2*gsl_spline2d_eval(fGSLSpline2D[0], r, erg, fGSLAccel2D[0].first, fGSLAccel2D[0].second) +
                 x_ae2*gsl_spline2d_eval(fGSLSpline2D[1], r, erg, fGSLAccel2D[1].first, fGSLAccel2D[1].second);
    } else if (fTableSubMode == 0) {
        ferr << "TRestAxionDiffSolarSpectrum has not been properly initialized." << endl;
    } else {
        ferr << "Internal error in TRestAxionDiffSolarSpectrum. This is a serious bug; please report it." << endl;
    }
    std::cout << fData.size() << " " << fData[0].size() << " " << fData[0][300] << " | " << fTableSubMode << " | " << x_agamma2 << endl;
    return result;
}
//double TRestAxionDiffSolarSpectrum::GetSolarAxionFlux(double erg_lo, double erg_hi, double g_agamma, double g_ae);
//double TRestAxionDiffSolarSpectrum::GetSolarAxionFlux(double r, double erg_lo, double erg_hi, double g_agamma, double g_ae);
//double TRestAxionDiffSolarSpectrum::GetIntegratedSolarAxionFlux(double erg_lo, double erg_hi, double g_agamma, double g_ae);
//double TRestAxionDiffSolarSpectrum::GetIntegratedSolarAxionFlux(double r, double erg_lo, double erg_hi, double g_agamma, double g_ae);
//double TRestAxionDiffSolarSpectrum::GetMCSamples(std::vector<double> rand_vars);
//double TRestAxionDiffSolarSpectrum::GetMCSamplesGivenR(std::vector<double> rand_vars, double r);
