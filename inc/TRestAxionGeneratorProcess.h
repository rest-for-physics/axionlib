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

#ifndef RestCore_TRestAxionGeneratorProcess
#define RestCore_TRestAxionGeneratorProcess

#include "TRandom3.h"
#include "TRestAxionEvent.h"
#include "TRestAxionSolarFlux.h"
#include "TRestEventProcess.h"

//! A process to initialize the axion event (mainly through TRestAxionSolarFlux)
class TRestAxionGeneratorProcess : public TRestEventProcess {
   private:
    /// A pointer to the specific TRestAxionEvent output
    TRestAxionEvent* fOutputAxionEvent;  //!

    /// A pointer to the TRestAxionSolarFlux metadata description
    TRestAxionSolarFlux* fAxionFlux;  //!

    /// Used internally to define the event id
    Int_t fCounter = 0;  //!

    /// Internal process random generator
    TRandom3* fRandom = nullptr;  //!

    /// The axion mass range in keV
    TVector2 fAxionMassRange = TVector2(1.e-6, 1e-2);  //<

    /// The target size in mm (or generator source extension) for the generator.
    Double_t fTargetRadius = 800;  //<

    /// The generator type (solarFlux/flat)
    TString fGeneratorType = "solarFlux";  //<

    /// Seed used in random generator
    Int_t fSeed = 0;  //<

    /// It defines the energy range for the axion event generator. Default between 50eV and 10keV.
    TVector2 fEnergyRange = TVector2(0.05, 10);  //<

    /// Absolute solar flux (cm-2 s-1). Required for future calculations.
    Double_t fTotalFlux = 0;  //<

    void Initialize() override;

    void LoadDefaultConfig();

   public:
    RESTValue GetInputEvent() const override { return nullptr; }
    RESTValue GetOutputEvent() const override { return fOutputAxionEvent; }

    void InitProcess() override;
    TRestEvent* ProcessEvent(TRestEvent* eventInput) override;

    void PrintMetadata() override;

    /// Returns the name of this process
    const char* GetProcessName() const override { return "axionGenerator"; }

    TRestAxionGeneratorProcess();   // Constructor
    ~TRestAxionGeneratorProcess();  // Destructor

    ClassDefOverride(TRestAxionGeneratorProcess, 4);
};
#endif
