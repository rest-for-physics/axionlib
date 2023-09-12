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

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TRestMetadata.h>

//! A metadata class to load tabulated solar axion fluxes
class TRestAxionSolarFlux : public TRestMetadata {
   private:
    /// Axion coupling. Defines coupling type and strength.
    std::string fCouplingType = "";  //<

    /// Axion coupling strength
    Double_t fCouplingStrength = 0;  //<

    /// Mass parameter
    Double_t fMass = 0;  //!

    /// Seed used in random generator
    Int_t fSeed = 0;  //<

    /// A metadata member to control if this class has been initialized
    Bool_t fTablesLoaded = false;  //!

   protected:
    /// A canvas pointer for drawing
    TCanvas* fCanvas = nullptr;  //!

    /// Random number generator
    TRandom3* fRandom = nullptr;  //!

    TRestAxionSolarFlux();
    TRestAxionSolarFlux(const char* cfgFileName, std::string name = "");

    /// It defines how to read the solar tables at the inhereted class
    virtual Bool_t LoadTables() = 0;

   public:
    /// It is required in order to load solar flux tables into memory
    void Initialize();

    /// It is required in order to load solar flux tables into memory for specific mass
    void InitializeMass( Double_t mass ) { SetMass(mass); Initialize(); }
    
    /// Set mass and reinitialise
	void SetMass( Double_t m ) { fMass = m; }	//Initialize(); }	

    /// It returns the integrated flux at earth in cm-2 s-1 for the given energy range
    virtual Double_t IntegrateFluxInRange(TVector2 eRange = TVector2(-1, -1)) = 0;

    /// It returns the total integrated flux at earth in cm-2 s-1
    virtual Double_t GetTotalFlux() = 0;

    /// It defines how to generate Monte Carlo energy and radius values to reproduce the flux
    virtual std::pair<Double_t, Double_t> GetRandomEnergyAndRadius(TVector2 eRange = TVector2(-1, -1)) = 0;

    /// It returns an energy integrated spectrum in cm-2 s-1 keV-1
    virtual TH1D* GetEnergySpectrum() = 0;

    virtual TCanvas* DrawSolarFlux();

    virtual void ExportTables(Bool_t ascii = false) {
        RESTWarning << "TRestAxionSolarFlux::ExportTables must be re-implemented in the inherited class"
                    << RESTendl;
    }

    Bool_t AreTablesLoaded() { return fTablesLoaded; }

    Double_t GetMass() { return fMass; }

    TH1D* GetFluxHistogram(std::string fname, Double_t binSize = 0.01);
    TCanvas* DrawFluxFile(std::string fname, Double_t binSize = 0.01);

    void PrintMetadata();

    ~TRestAxionSolarFlux();

    ClassDef(TRestAxionSolarFlux, 2);
};
#endif
