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

#include <sstream>
#include <iomanip>

#include "solaxflux/utils.hpp"
#include "solaxflux/solar_model.hpp"
#include "solaxflux/spectral_flux.hpp"

#include <TRestMetadata.h>

//! A metadata class to define theoretical axion models and calculations related
class TRestAxionSolarModel : public TRestMetadata {
   private:
    void Initialize();
    void InitFromConfigFile();

    // The SolarModel class object from the external SolarAxionFlux library
    SolarModel fSol;

    // Variables for diagnostics and metadata
    bool fSolarModelInitialized = false;
    double fRefPhotonCoupling;
    double fRefElectronCoupling;
    std::string fExternalLibraryNameAndVersion;
    std::string fSolarModelFile;
    std::string fOpacityCodeName;

  public:
    // Constructors and destructors
    TRestAxionSolarModel();
    TRestAxionSolarModel(const char *cfgFileName, std::string name="");
    ~TRestAxionSolarModel();

    //// Member functions to calculate the solar axion flux

    // The following member functions calculate the spectral solar axion flux i.e. dPhi/dE after
    // after integrating up to(!) the radius 'r_max' on the solar disc for all values of 'energies'.
    // Returns a vector of flux values in units of axions/(cm^2 s keV). Optionally saves the result
    // as 'outputfile' if a path is provided.
    std::vector<double> CalcSpectralSolarAxionFluxPrimakoff(std::vector<double> energies,
                            double r_max=1.0, std::string outputfile="", double g_agamma=1.0e-10);
    // std::vector<double> GetSolarAxionFluxLP(std::vector<double> energies, double r_max=1.0, std::string outputfile="", double g_agamma=1.0e-10);
    // std::vector<double> GetSolarAxionFluxTP(std::vector<double> energies, double r_max=1.0, std::string outputfile="", double g_agamma=1.0e-10);
    std::vector<double> CalcSpectralSolarAxionFluxABC(std::vector<double> energies,
                            double r_max=1.0, std::string outputfile="", double g_ae=1.0e-13);
    // Different to the functions above, the next function returns three columns containing the energies, ALP-photon, and ALP-electron flux, respectively
    std::vector<std::vector<double> > CalcSpectralSolarAxionFluxAll(std::vector<double> energies,
                            double r_max=1.0, std::string outputfile="",
                            double g_agamma=1.0e-10, double g_ae=1.0e-13);

    // The following member functions calculate the spectral and spatial solar axion flux
    // i.e. d^2Phi/dEdr where r is the dimensionless(!) radius on the solar disc for all values of
    // 'energies' and 'radii'.
    // Returns a vector of vectors, containing the selected radii, energies, and (potentially
    // multiple) flux values in units of axions/(cm^2 s keV).
    std::vector<std::vector<double> > CalcSpectralAndSpatialSolarAxionFluxPrimakoff(std::vector<double> energies, std::vector<double> radii,
                                          std::string outputfile="", double g_agamma=1.0e-10);
    // std::vector<std::vector<double> > CalcSpectralAndSpatialSolarAxionFluxLP(std::vector<double> energies, std::vector<double> radii, std::string outputfile="", double g_agamma=1.0e-10);
    // std::vector<std::vector<double>> CalcSpectralAndSpatialSolarAxionFluxTP(std::vector<double> energies, std::vector<double> radii, std::string outputfile="", double g_agamma=1.0e-10);
    std::vector<std::vector<double> > CalcSpectralAndSpatialSolarAxionFluxABC(std::vector<double> energies, std::vector<double> radii,
                                          std::string outputfile="", double g_ae=1.0e-13);
    std::vector<std::vector<double> > CalcSpectralAndSpatialSolarAxionFluxAll(std::vector<double> energies, std::vector<double> radii,
                                          std::string outputfile="", double g_agamma=1.0e-10, double g_ae=1.0e-13);

    //// Diagnostics and metadata
    void PrintMetadata();
    bool IsSolarModelClassReady() { return fSolarModelInitialized; };
    std::string GetSolarModelFileName() { return fSolarModelFile; };
    std::string GetOpacityCodeName() { return fOpacityCodeName; };
    std::string GetExternalLibraryVersion() { return fExternalLibraryNameAndVersion; };

    ClassDef(TRestAxionSolarModel, 1);
};

#endif
