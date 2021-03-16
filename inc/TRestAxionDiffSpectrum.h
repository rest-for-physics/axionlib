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

#ifndef _TRestAxionDiffSpectrum
#define _TRestAxionDiffSpectrum

// #include <TRestTools.h>
#include <TRestMetadata.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

//! A metadata class to define a differential solar axion spectrum (d^2Phi/dEdr) and functions to evaluate it.
class TRestAxionDiffSpectrum : public TRestMetadata {
private:
    void Initialize();
    void InitFromConfigFile();

    void Init2DSplines(const int n_grids);

    std::vector<std::vector<double> > fData;
    std::vector<std::pair<gsl_interp_accel*,gsl_interp_accel*>> fGSLAccel2D;
    std::vector<gsl_spline2d*> fGSLSpline2D;
    std::vector<double*> fGSLMem2D;

    double fRefPhotonCoupling = NAN;
    double fRefElectronCoupling = NAN;
    std::string fTableFileName;

    int fTableSubMode = 0;

public:
    TRestAxionDiffSpectrum();
    TRestAxionDiffSpectrum(const char *cfgFileName, std::string name = "");
    ~TRestAxionDiffSpectrum();

    double GetDoubleDifferentialSolarAxionFlux(double r, double erg, double g_agamma = 1.0e-10, double g_ae = 0);
    //double GetSolarAxionFlux(double erg_lo, double erg_hi, double g_agamma = 1.0e-10, double g_ae = 0);
    //double GetSolarAxionFlux(double r, double erg_lo, double erg_hi, double g_agamma, double g_ae);
    //double GetIntegratedSolarAxionFlux(double erg_lo, double erg_hi, double g_agamma = 1.0e-10, double g_ae = 0);
    //double GetIntegratedSolarAxionFlux(double r, double erg_lo, double erg_hi, double g_agamma, double g_ae);
    //double GetMCSamples(std::vector<double> rand_vars);
    //double GetMCSamplesGivenR(std::vector<double> rand_vars, double r);

    void PrintMetadata();
    bool IsSpectrumReady() { return (fTableSubMode > 0); }
    int GetDiffSpectrumMode() { return fTableSubMode; }

    ClassDef(TRestAxionDiffSpectrum, 1);
};

#endif
