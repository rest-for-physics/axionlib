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

#include <TH1F.h>
#include <TH2F.h>
#include <TRandom3.h>

#include <TRestMetadata.h>

#include <TRestAxionSolarModel.h>

//! A metadata class to load tabulated solar axion fluxes
class TRestAxionSolarFlux : public TRestMetadata {
   private:
    void Initialize();

    /// The filename containning the solar flux table with continuum spectrum
    std::string fFluxDataFile = "";  //<

    /// The filename containning the solar flux spectra for monochromatic spectrum
    std::string fFluxSptFile = "";  //<

    /// Axion coupling. Defines coupling type and strength.
    std::string fCouplingType;  //<

    /// Axion coupling strength
    Double_t fCouplingStrength;  //<

    /// Seed used in random generator
    Int_t fSeed = 0;  //<

    /// It will be used when loading `.flux` files to define the input file energy binsize in eV.
    Double_t fBinSize = 0;  //<

    /// It will be used when loading `.flux` files to define the threshold for peak identification
    Double_t fPeakThreshold = 0;  //<

    /// The tabulated solar flux continuum spectra TH1F(100,0,20)keV in cm-2 s-1 keV-1 versus solar radius
    std::vector<TH1F*> fFluxTable;  //!

    /// The tabulated solar flux in cm-2 s-1 for a number of monochromatic energies versus solar radius
    std::map<Double_t, TH1F*> fFluxLines;  //!

    /// Accumulative integrated solar flux for each solar ring for continuum spectrum (renormalized to unity)
    std::vector<Double_t> fFluxTableIntegrals;  //!

    /// Accumulative integrated solar flux for each monochromatic energy (renormalized to unity)
    std::vector<Double_t> fFluxLineIntegrals;  //!

    /// Total solar flux for monochromatic contributions
    Double_t fTotalMonochromaticFlux = 0;  //!

    /// Total solar flux for monochromatic contributions
    Double_t fTotalContinuumFlux = 0;  //!

    /// The ratio between monochromatic and total flux
    Double_t fFluxRatio = 0;  //!

    /// Random number generator
    TRandom3* fRandom = nullptr;  //!

    void ReadFluxFile();
    void LoadContinuumFluxTable();
    void LoadMonoChromaticFluxTable();
    void IntegrateSolarFluxes();

   public:
    /// It returns true if continuum flux spectra was loaded
    Bool_t isSolarTableLoaded() { return fFluxTable.size() > 0; }

    /// It returns true if monochromatic flux spectra was loaded
    Bool_t isSolarSpectrumLoaded() { return fFluxLines.size() > 0; }

    /// It returns the total integrated flux at earth in cm-2 s-1
    Double_t GetTotalFlux() { return fTotalContinuumFlux + fTotalMonochromaticFlux; }

    std::pair<Double_t, Double_t> GetRandomEnergyAndRadius();

    void LoadTables();

    /// Tables might be loaded using a solar model description by TRestAxionSolarModel
    void InitializeSolarTable(TRestAxionSolarModel* model) {
        // TOBE implemented
        // This method should initialize the tables fFluxTable and fFluxLines
    }

    void ExportTables(std::string fname);

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
