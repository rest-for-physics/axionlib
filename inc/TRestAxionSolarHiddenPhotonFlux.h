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

#ifndef _TRestAxionSolarHiddenPhotonFlux
#define _TRestAxionSolarHiddenPhotonFlux

#include <TRestAxionSolarFlux.h>
#include <TRestAxionSolarModel.h>

//! A metadata class to load tabulated solar hidden photon fluxes. Kinetic mixing set to 1.
class TRestAxionSolarHiddenPhotonFlux : public TRestAxionSolarFlux {
   private:
    /// The filename containing the solar flux table with continuum spectrum
    std::string fFluxDataFile = "";  //<

    /// The filename containing the resonance width (wGamma)
    std::string fWidthDataFile = "";  //<

    /// The filename containing the plasma frequency (wp) table
    std::string fPlasmaFreqDataFile = "";  //<

    /// It will be used when loading `.flux` files to define the input file energy binsize in eV.
    Double_t fBinSize = 0;  //<

    /// The tabulated solar flux continuum spectra TH1F(200,0,20)keV in cm-2 s-1 keV-1 versus solar radius
    std::vector<TH1F*> fFluxTable;  //!

    /// The tabulated resonance width TH1F(200,0,20)keV in eV2 versus solar radius
    std::vector<TH1F*> fWGammaTable;  //!

    /// The solar plasma frequency vector in eV versus solar radius
    std::vector<TH1F*> fWpTable;  //!

    /// The total solar flux TH1F(200,0,20)keV in cm-2 s-1 keV-1 versus solar radius
    std::vector<TH1F*> fFullFluxTable;  //!

    /// Accumulative integrated solar flux for each solar ring for continuum spectrum (renormalized to unity)
    std::vector<Double_t> fFluxTableIntegrals;  //!

    /// Total solar flux for monochromatic contributions
    Double_t fTotalContinuumFlux = 0;  //!

    /// A pointer to the continuum spectrum histogram
    TH1F* fContinuumHist = nullptr;  //!

    /// A pointer to the superposed monochromatic and continuum spectrum histogram
    TH1F* fTotalHist = nullptr;  //!

    void ReadFluxFile();
    void LoadContinuumFluxTable();
    void LoadMonoChromaticFluxTable();
    void IntegrateSolarFluxes();

   public:
    /// It returns true if continuum flux spectra was loaded
    Bool_t isSolarTableLoaded() { return fFluxTable.size() > 0; }

    /// It returns true if width table was loaded
    Bool_t isWidthTableLoaded() { return fWGammaTable.size() > 0; }

    /// It returns true if plasma frequency table was loaded
    Bool_t isPlasmaFreqLoaded() { return fWpTable.size() > 0; }

    /// It returns the integrated flux at earth in cm-2 s-1 for the given energy range
    Double_t IntegrateFluxInRange(TVector2 eRange = TVector2(-1, -1), Double_t mass = 0) override;

    /// It defines how to generate Monte Carlo energy and radius values to reproduce the flux
    std::pair<Double_t, Double_t> GetRandomEnergyAndRadius(TVector2 eRange = TVector2(-1, -1),
                                                           Double_t mass = 0) override;

    /// It defines how to read the solar tables at the inhereted class for a given mass in eV
    Bool_t LoadTables(Double_t mass) override;
    
    void LoadContinuumTable();
    void LoadWidthTable();
    void LoadPlasmaFreqTable();
    
    // calculate solar HP flux from the 3 tables and mass
    void CalculateSolarFlux();
    /// It returns the total integrated flux at earth in cm-2 s-1
    Double_t GetTotalFlux(Double_t mass = 0) override { return fTotalContinuumFlux; }

    /// It returns an energy integrated spectrum in cm-2 s-1 keV-1
    TH1F* GetEnergySpectrum(Double_t m = 0) override { return GetTotalSpectrum(); }

    TH1F* GetContinuumSpectrum();
    TH1F* GetTotalSpectrum();

    virtual TCanvas* DrawSolarFlux() override;

    /// Tables might be loaded using a solar model description by TRestAxionSolarModel
    void InitializeSolarTable(TRestAxionSolarModel* model) {
        // TOBE implemented
        // This method should initialize the tables fFluxTable and fFluxLines
    }

    void ExportTables(Bool_t ascii = false) override;

    void PrintMetadata() override;

    void PrintContinuumSolarTable();
    void PrintIntegratedRingFlux();

    TRestAxionSolarHiddenPhotonFlux();
    TRestAxionSolarHiddenPhotonFlux(const char* cfgFileName, std::string name = "");
    ~TRestAxionSolarHiddenPhotonFlux();

    ClassDefOverride(TRestAxionSolarHiddenPhotonFlux, 1);
};
#endif
