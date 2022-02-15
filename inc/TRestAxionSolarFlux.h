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

#ifndef _TRestAxionSolarFlux
#define _TRestAxionSolarFlux

#include <TRestMetadata.h>

#include <TRestAxionSolarModel.h>

//! A metadata class to load tabulated solar axion fluxes
class TRestAxionSolarFlux : public TRestMetadata {
   private:
    void Initialize();

    void InitFromConfigFile();

    // The filename containning the solar flux table with continuum spectrum
    std::string fFluxDataFile = "";  //<

    // The filename containning the solar flux spectra for monochromatic spectrum
    std::string fFluxSptFile = "";  //<

    // Axion coupling. Defines coupling type and strength.
    std::string fCouplingType;  //<

    // Axion coupling strength
    Double_t fCouplingStrength;  //<

    // Contains the tabulated solar flux in cm-2 s-1 keV-1 versus solar radius and energy
    std::vector<std::vector<Double_t>> fFluxTable;  //!

    // Contains the tabulated solar flux in cm-2 s-1 versus solar radius for any number of spectral energies
    std::map<Double_t, std::vector<Double_t>> fFluxLines;  //!

    // Contains the contribution to the flux for each solar ring (continuum + monochromatic)
    std::vector<Double_t> fFluxPerRadius;  //!

    // The total solar disk integrated flux in cm-2 s-1
    Double_t fTotalFlux = 0;  //!

    void LoadContinuumFluxTable();
    void LoadMonoChromaticFluxTable();
    void IntegrateSolarRingsFluxes();

   public:
    Bool_t isSolarTableLoaded() { return fFluxTable.size() > 0; }

    Bool_t isSolarSpectrumLoaded() { return fFluxLines.size() > 0; }

    std::pair<Double_t, Double_t> GetRandomEnergyAndRadius() {
        std::pair<Double_t, Double_t> result = {0, 0};
        return result;
    }

    /// Tables might be loaded using a solar model description by TRestAxionSolarModel
    void InitializeSolarTable(TRestAxionSolarModel* model) {
        // TOBE implemented
        // This method should initialize the tables fFluxTable and fFluxLines
    }

    void ExportTables(std::string fname) {
        // TOBE implemented. Creates fname.dat and fname.spt
        // If we have external methods to initialize solar flux tables this method
        // might be used to generate the tables that can be used later on directly
        //
        // Check data/solarFlux/README.md for data format and file naming conventions
    }

    void PrintMetadata();

    void PrintContinuumSolarTable();
    void PrintIntegratedRingFlux();
    void PrintMonoChromaticFlux();

    TRestAxionSolarFlux();
    TRestAxionSolarFlux(const char* cfgFileName, std::string name = "");
    ~TRestAxionSolarFlux();

    ClassDef(TRestAxionSolarFlux, 1);
};
#endif
