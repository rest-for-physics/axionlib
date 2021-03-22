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
/// The TRestAxionSolarSpectrum serves as a way to read in an perform calculations
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
/// \class      TRestAxionSolarSpectrum
/// \author     Sebastian Hoof
///
/// <hr>
///
 *************************************************************************/

// See this header file for more info on the functions.
#include "TRestAxionSolarSpectrum.h"

// Map for named approximations of the spectrum
// The parameters correspond to 'norm', 'g' (user-defined coupling),  and spectral parameters 'a'
// and 'b' (numbers). The ansatz for the flux is given by the following expression:
// flux = (norm / cm^2 s keV) * (g/gref)^2 * (energy / keV)^a * exp(-b * energy / keV)
const std::map<std::string, std::vector<double>> avail_approximations = {
    {"arXiv_0702006_Primakoff", {6.02e10, 1.0e-10, 2.481, 1.0 / 1.205}},
    {"arXiv_1302.6283_Primakoff", {2.0e22 / (3600.0 * 24.0 * 365.25 * 1.0e4), 1.0e-10, 2.450, 0.829}},
    {"arXiv_1302.6283_Compton", {4.2e24 / (3600.0 * 24.0 * 365.25 * 1.0e4), 1.0e-13, 2.987, 0.776}}
};

// Helper function to sort and get unique entries of a vector.
std::vector<double> sort_unique_values(std::vector<double> x) {
    std::vector<double> unique_x = x;
    std::sort(unique_x.begin(), unique_x.end());
    unique_x.erase(std::unique(unique_x.begin(), unique_x.end()), unique_x.end());
    return unique_x;
}

// Credit to Christoph Weniger
int read_natural_ascii_table(std::string filename, std::vector<std::vector<double> >& data) {
  data.clear();
  std::ifstream in (filename.c_str(), std::ios::binary);
  if (in.fail()) { ferr << "File "+filename+" does not exist or could not be opened." << endl; }
  std::string line;
  while (std::getline(in, line)) {
      if (line[0] == '#') continue;
      std::stringstream strstr (line);
      size_t i = 0;
      double tmp;
      while(strstr >> tmp) {
          if (i+1 > data.size()) data.resize(i+1);
          data[i].push_back(tmp);
          i++;
      }
  }
  in.close();
  return 0;
}

ClassImp(TRestAxionSolarSpectrum);

// Add basic TRest info
void TRestAxionSolarSpectrum::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);
}

void TRestAxionSolarSpectrum::Init1DSpline(const int index) {
    const double* x = &fData[0][0];
    const double* y = &fData[index+1][0];
    int pts = fData[0].size();
    fGSLAccel1D[index] = gsl_interp_accel_alloc();
    fGSLSpline1D[index] = gsl_spline_alloc(gsl_interp_linear, pts);
    gsl_spline_init(fGSLSpline1D[index], x, y, pts);
}

void TRestAxionSolarSpectrum::Init2DSplines(const int n_grids) {
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

void TRestAxionSolarSpectrum::InitFromConfigFile() {
    Initialize();
    fMode = GetParameter("mode");
    if (fMode == "table") {
        fTableFileName = GetParameter("spectrumTableFileName");
        fRefPhotonCoupling = GetDblParameterWithUnits("g_agamma", NAN);
        fRefElectronCoupling = GetDblParameterWithUnits("g_ae", NAN);
        bool g_agamma_nan = std::isnan(fRefPhotonCoupling) || not(fRefPhotonCoupling > 0);
        bool g_ae_nan = std::isnan(fRefElectronCoupling) || not(fRefElectronCoupling > 0);
        int temp_submode = not(g_agamma_nan) + 2*g_ae_nan;
        fTableFileName = SearchFile((std::string)fTableFileName);
        if (fTableFileName == "") {
            ferr << "File not found : " << fTableFileName << endl;
        } else if (g_agamma_nan && g_ae_nan) {
            ferr << "You need to supply a reference value for either the axion-photon coupling "
                    "'g_agamma' (in units of GeV^{-1}) or the axion-electron coupling 'g_ae' when "
                    "using the 'table' mode of TRestAxionSolarSpectrum."
                 << endl;
        } else {
            // TRestTools::ReadASCIITable(fTableFileName, fData); // -> unnatural rows vs columns
            read_natural_ascii_table(fTableFileName, fData);
            int n_cols = fData.size();
            if (n_cols < 2) {
                ferr << "Tabulated solar axion spectrum " << fTableFileName << "contains less than "
                        "two columns." << endl;
            } else if (n_cols==2) {
                if (temp_submode<3) {
                    // Intialise one 1D interpolation
                    fGSLAccel1D.resize(1);
                    fGSLSpline1D.resize(1);
                    Init1DSpline(0);
                } else {
                    ferr << "Tabulated solar axion spectrum " << fTableFileName << "contains two columns "
                            "but you supplied two non-NAN values for the couplings 'g_agamma' and 'g_ae'."
                         << endl;
                }
            } else if (n_cols==3) {
                if (temp_submode<3) {
                    temp_submode += 4;
                    // Initialise one 2D interpolation
                    Init2DSplines(1);
                } else {
                    // Initialise two 1D interpolations
                    fGSLAccel1D.resize(2);
                    fGSLSpline1D.resize(2);
                    Init1DSpline(0);
                    Init1DSpline(1);
                }
            } else if (n_cols==4) {
                if (temp_submode==3) {
                    // Initialise two 2D interpolations
                    temp_submode += 4;
                    Init2DSplines(2);
                } else {
                    ferr << "Tabulated solar axion spectrum " << fTableFileName << "contains four "
                            "columns; this is incompatible with TRestAxionSolarSpectrum. You need to "
                            "supply two non-NAN values for the couplings 'g_agamma' and 'g_ae'."
                        << endl;
                }
            } else {
                debug << "Tabulated solar axion spectrum " << fTableFileName << "contains more than "
                         "four columns. Excess columns will be ignored." << endl;
            }
            fSpectrumReady = true;
        }
    } else if (fMode == "analytical") {
        std::string named_approx = GetParameter("named_approx", "NN");
        fAnalyticalNorm = GetDblParameterWithUnits("norm", NAN);
        fAnalyticalRefG = GetDblParameterWithUnits("gref", NAN);
        fAnalyticalA = GetDblParameterWithUnits("a", NAN);
        fAnalyticalB = GetDblParameterWithUnits("b", NAN);
        bool check_named_approx = not(avail_approximations.find(named_approx) == avail_approximations.end());
        bool check_numbers = not( std::isnan(fAnalyticalNorm) || std::isnan(fAnalyticalRefG) ||
                                  std::isnan(fAnalyticalA) || std::isnan(fAnalyticalB) );
        if (check_named_approx && check_numbers) {
            ferr << "You want to use the 'analytical' mode of TRestAxionSolarSpectrum but supplied both "
                    "a 'named_approx' and explicit values for one or more of the parameters 'norm', "
                    "'g', 'a' or 'b'. You can only do one or the other." << endl;
        } else if (check_named_approx) {
            std::vector<double> p = avail_approximations.at(named_approx);
            fAnalyticalNorm = p[0];
            fAnalyticalRefG = p[1];
            fAnalyticalA = p[2];
            fAnalyticalB = p[3];
            fSpectrumReady = true;
        } else if (check_numbers) {
            fSpectrumReady = true;
        } else {
            std::string avail_approximation_names = "";
            for (auto el : avail_approximations) { avail_approximation_names += el.first + "\n"; }
            ferr << "You want to use the 'analytical' mode of TRestAxionSolarSpectrum but NEITHER supplied a "
                    "known 'named_approx' NOR ALL OF the four required parameters 'a', 'b', 'g', "
                    "and 'norm'.\nInformation: the approximations known to the code are :"
                    +avail_approximation_names
                 << endl;
        }
    } else {
        ferr << "Unknown mode for TRestAxionSolarSpectrum! Choose either 'table' or 'analytical'."
             << endl;
    }
}

TRestAxionSolarSpectrum::TRestAxionSolarSpectrum() : TRestMetadata() { Initialize(); }

TRestAxionSolarSpectrum::TRestAxionSolarSpectrum(const char* cfgFileName, std::string name)
    : TRestMetadata(cfgFileName) {
    debug << "Creating instance of TRestAxionSolarSpectrum from file " + fConfigFileName + "..." << endl;
    Initialize();
    LoadConfigFromFile(fConfigFileName, name);
    PrintMetadata();
}

// Need to free allocated memory by GSL routines!
TRestAxionSolarSpectrum::~TRestAxionSolarSpectrum() {
    for (auto s : fGSLSpline1D) { gsl_spline_free(s); }
    for (auto s : fGSLSpline2D) { gsl_spline2d_free(s); }
    for (auto a : fGSLAccel1D) { gsl_interp_accel_free(a); }
    for (auto a : fGSLAccel2D) { gsl_interp_accel_free(a.first); gsl_interp_accel_free(a.second); }
    for (auto z : fGSLMem2D) { free(z); }
}

void TRestAxionSolarSpectrum::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    if (fMode != "none") {
        metadata << " Solar spectrum created in mode '" << fMode << "'." << endl;
    } else {
        metadata << " Solar spectrum has not yet been initialised, yet!" << endl;
    }
    if (fMode == "table") {
        metadata << " - Tabulated spectrum file used : "
                 << TRestTools::SeparatePathAndName(fTableFileName).second << endl;
    }
    metadata << "-------------------------------------------------" << endl;
    metadata << " - Units of the solar axion flux from this class : axions / cm^2 s keV" << endl;
    if (fMode == "table") {
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
    }
    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

double TRestAxionSolarSpectrum::GetDifferentialSolarAxionFlux(double erg, double g_agamma, double g_ae) { return 0; }
double TRestAxionSolarSpectrum::GetDifferentialSolarAxionFlux(double r, double erg, double g_agamma, double g_ae) { return 0; }
double TRestAxionSolarSpectrum::GetSolarAxionFlux(double erg_lo, double erg_hi, double erg_delta, double g_agamma, double g_ae) { return 0; }
double TRestAxionSolarSpectrum::GetSolarAxionFlux(double r, double erg_lo, double erg_hi, double erg_delta, double g_agamma, double g_ae) { return 0; }

double TRestAxionSolarSpectrum::GetMCSampleFullSun(double rv) {
    double result = 0;
    if (fMode == "analytical") {
        static double ap1 = fAnalyticalA + 1.0;
        static double invb = 1.0/fAnalyticalB;
        // The spectral ansatz for the 'analytical' mode is actually just ~ gamma distribution. Can
        // get samples directly from the corresponding analytical inverse CDF.
        result = gsl_cdf_gamma_Pinv(rv, ap1, invb);
    //} else if (fMode == "table") {
    //    result = 0;
    } else {
        ferr << "Member function GetMCSampleFullSun cannot be used with mode " + fMode + "." << endl;
    }
    return result;
}
