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

#ifndef _TRestAxionSolarModel
#define _TRestAxionSolarModel

#include "solaxflux/solar_model.hpp"

#include <TRestMetadata.h>

//! A metadata class to define theoretical axion models and calculations related
class TRestAxionSolarModel : public TRestMetadata {
   private:
    void Initialize();
    void InitFromConfigFile();

    // The solar model object
    SolarModel sol;

    // The range integration step size for solar energy spectrum (in keV)
    //double fStep = 0.01;
    //std::vector<double> fEnergyRange = {0.1,10.0};
    //double r_max;

    // Variables for diagnostics and metadata
    bool bSolarModelInitialized = false;
    double fRefPhotonCoupling;
    double fRefElectronCoupling;
    std::string sExternalLibraryName;
    std::string sSolarModelFile;
    std::string sOpacityCodeName;

  public:
    // Constructors and destructors
    TRestAxionSolarModel();
    TRestAxionSolarModel(const char *cfgFileName, std::string name = "");
    ~TRestAxionSolarModel();

    // Solar axion flux calculators
    std::vector<double> GetSolarAxionFluxGAGamma(std::vector<double> energies, double g_agamma, double r_max);
    std::vector<double> GetSolarAxionFluxGAGamma(std::vector<double> energies, double r_max=1.0);
    //std::vector<double> GetSolarAxionFluxGAGamma();
    std::vector<double> GetSolarAxionFluxGAE(std::vector<double> energies, double g_agae, double r_max);
    std::vector<double> GetSolarAxionFluxGAE(std::vector<double> energies, double r_max=1.0);
    //std::vector<double> GetSolarAxionFluxGAE();

    // Diagnostics and metadata
    void PrintMetadata();
    bool isSolarModelClassReady() { return bSolarModelInitialized; }
    std::string GetSolarModelFileName();
    //std::string GetOpacityCodeName() { return fOpacityCodeName; }

    ClassDef(TRestAxionSolarModel, 1);
};

#endif
