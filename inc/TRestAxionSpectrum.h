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

#include <TRestTools.h>
#include <TRestMetadata.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

//! A metadata class to define a solar axion spectrum (dPhi/dE) and functions to evaluate it.
class TRestAxionSpectrum : public TRestMetadata {
private:
    // Generic initialization routines.
    void Initialize();
    // The mode of this class. There are 2 different modes available: table and analytic.
    // Energy in keV, flux in axions / cm^2 s keV
    // table      : Provide a text file with up to 4 columns (and reference values of the associated couplings).
    //              The numerical value of the sub-mode is computed via binary logic (true = 1, false = 0):
    //              fTableSubMode = 2^2 * (r values provided?) + 2 * (g_ae provided?) + (g_agamma provided?)
    // analytical : Provide EITHER the parameters norm, gref, a and b for the following ansatz
    //                      OR a named set of parameters available in the data/ folder
    //              Ansatz: flux = norm * (g/gref)^2 * (energy / keV)^a * exp(-b * energy / keV)
    void InitFromConfigFile();

    void Init1DSpline(const int index);
    void Init2DSplines(const int n_grids);

    std::vector<std::vector<double> > fData;
    std::vector<gsl_interp_accel*> fGSLAccel1D;
    std::vector<gsl_spline*> fGSLSpline1D;
    std::vector<std::pair<gsl_interp_accel*,gsl_interp_accel*>> fGSLAccel2D;
    std::vector<gsl_spline2d*> fGSLSpline2D;
    std::vector<double*> fGSLMem2D;

    double fRefPhotonCoupling = NAN;
    double fRefElectronCoupling = NAN;
    std::string fTableFileName;

    double fAnalyticalRefG = NAN;
    double fAnalyticalNorm = NAN;
    double fAnalyticalA = NAN;
    double fAnalyticalB = NAN;

    std::string fMode = "none";
    int fTableSubMode = 0;
    bool fSpectrumInitialized = false;

public:
    TRestAxionSpectrum();
    TRestAxionSpectrum(const char *cfgFileName, std::string name = "");
    ~TRestAxionSpectrum();

    double GetDifferentialSolarAxionFlux(double erg, double g_agamma = 1.0e-10, double g_ae = 0);
    double GetDifferentialSolarAxionFlux(double r, double erg, double g_agamma, double g_ae);
    double GetSolarAxionFlux(double erg_lo, double erg_hi, double erg_delta, double g_agamma  = 1.0e-10, double g_ae = 0);
    double GetSolarAxionFlux(double r, double erg_lo, double erg_hi, double erg_delta, double g_agamma, double g_ae);
    //double GetMCSamples(std::vector<double> rand_vars);
    //double GetMCSamplesGivenR(std::vector<double> rand_vars, double r);

    void PrintMetadata();
    bool IsSpectrumReady() { return fSpectrumInitialized; }
    std::string GetSpectrumMode() { return fMode; };

    ClassDef(TRestAxionSpectrum, 1);
};

//! A metadata class to define a differential solar axion spectrum (d^2Phi/dEdr) and functions to evaluate it.
//class TRestAxionDiffSpectrum : public TRestMetadata {
/*
private:
    void Initialize();
    void InitFromConfigFile();

    std::vector<std::vector<double> > fData;
    std::vector<std::pair<gsl_interp_accel*,gsl_interp_accel*>> fGSLAccel2D;
    std::vector<gsl_spline2d*> fGSLSpline2D;

    double fRefPhotonCoupling = NAN;
    double fRefElectronCoupling = NAN;
    std::string fTableFileName;

    int fMode = 0;
    bool fSpectrumInitialized = false;

public:
    TRestAxionSpectrum();
    TRestAxionSpectrum(const char *cfgFileName, std::string name = "");
    ~TRestAxionSpectrum();

    double GetDoubleDifferentialSolarAxionFlux(double r, double erg, double g_agamma = 1.0e-10, double g_ae = 0);
    //double GetSolarAxionFlux(double erg_lo, double erg_hi, double g_agamma  =1.0e-10, double g_ae = 0);
    //double GetSolarAxionFlux(double r, double erg_lo, double erg_hi, double g_agamma, double g_ae);
    //double GetIntegratedSolarAxionFlux(double erg_lo, double erg_hi, double g_agamma  =1.0e-10, double g_ae = 0);
    //double GetIntegratedSolarAxionFlux(double r, double erg_lo, double erg_hi, double g_agamma, double g_ae);
    //double GetMCSamples(std::vector<double> rand_vars);
    //double GetMCSamplesGivenR(std::vector<double> rand_vars, double r);

    void PrintMetadata();
    bool IsSpectrumReady() { return fSpectrumInitialized; }
    int GetDiffSpectrumMode() { return fMode; };

    ClassDef(TRestAxionDiffSpectrum, 1);
    */
//};

#endif
