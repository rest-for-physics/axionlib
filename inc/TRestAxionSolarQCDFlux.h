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

#ifndef _TRestAxionSolarQCDFlux
#define _TRestAxionSolarQCDFlux

#include <TRestAxionSolarFlux.h>
#include <TRestAxionSolarModel.h>

//! A metadata class to load tabulated solar axion fluxes. Mass independent.
class TRestAxionSolarQCDFlux : public TRestAxionSolarFlux {
   private:
    /// The filename containning the solar flux table with continuum spectrum
    std::string fFluxDataFile = "";  //<

    /// The filename containning the solar flux spectra for monochromatic spectrum
    std::string fFluxSptFile = "";  //<

    /// It will be used when loading `.flux` files to define the input file energy binsize in eV.
    Double_t fBinSize = 0;  //<

    /// It will be used when loading `.flux` files to define the threshold for peak identification
    Double_t fPeakSigma = 0;  //<

    /// The tabulated solar flux continuum spectra TH1F(200,0,20)keV in cm-2 s-1 keV-1 versus solar radius
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

    /// A pointer to the continuum spectrum histogram
    TH1F* fContinuumHist = nullptr;  //!

    /// A pointer to the monochromatic spectrum histogram
    TH1F* fMonoHist = nullptr;  //!

    /// A pointer to the superposed monochromatic and continuum spectrum histogram
    TH1F* fTotalHist = nullptr;  //!

    void ReadFluxFile();
    void LoadContinuumFluxTable();
    void LoadMonoChromaticFluxTable();
    void IntegrateSolarFluxes();

   public:
    /// It returns true if continuum flux spectra was loaded
    Bool_t isSolarTableLoaded() { return fFluxTable.size() > 0; }

    /// It returns true if monochromatic flux spectra was loaded
    Bool_t isSolarSpectrumLoaded() { return fFluxLines.size() > 0; }

    /// It returns the integrated flux at earth in cm-2 s-1 for the given energy range
    Double_t IntegrateFluxInRange(TVector2 eRange = TVector2(-1, -1), Double_t mass = 0) override;

    /// It defines how to generate Monte Carlo energy and radius values to reproduce the flux
    std::pair<Double_t, Double_t> GetRandomEnergyAndRadius(TVector2 eRange = TVector2(-1, -1),
                                                           Double_t mass = 0) override;

    /// It defines how to read the solar tables at the inhereted class
    Bool_t LoadTables() override;

    /// It returns the total integrated flux at earth in cm-2 s-1
    Double_t GetTotalFlux(Double_t mass = 0) override {
        return fTotalContinuumFlux + fTotalMonochromaticFlux;
    }

    /// It returns an energy integrated spectrum in cm-2 s-1 keV-1
    TH1F* GetEnergySpectrum(Double_t m = 0) override { return GetTotalSpectrum(); }

    TH1F* GetContinuumSpectrum();
    TH1F* GetMonochromaticSpectrum();
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
    void PrintMonoChromaticFlux();

    TRestAxionSolarQCDFlux();
    TRestAxionSolarQCDFlux(const char* cfgFileName, std::string name = "");
    ~TRestAxionSolarQCDFlux();

    ClassDefOverride(TRestAxionSolarQCDFlux, 1);
};
#endif
