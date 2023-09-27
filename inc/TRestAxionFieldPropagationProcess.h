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

#ifndef RestCore_TRestAxionFieldPropagationProcess
#define RestCore_TRestAxionFieldPropagationProcess

#include "TRestAxionEvent.h"
#include "TRestAxionEventProcess.h"
#include "TRestAxionMagneticField.h"
#include "TRestAxionQCDField.h"
#include "TRestPhysics.h"
#include "TVector3.h"
#include "TVectorD.h"

//! A process to introduce the magnetic field profile integration along the track
class TRestAxionFieldPropagationProcess : public TRestAxionEventProcess {
   private:
    /// The differential length in mm used for the field integration
    Double_t fIntegrationStep = 50;  //<

    /// The additional length in mm that the converted photon propagates without magnetic field
    Double_t fBufferGasAdditionalLength = 0;  //<

    /// A pointer to the magnetic field description stored in TRestRun
    TRestAxionMagneticField* fMagneticField = nullptr;  //!

    /// A pointer to TRestAxionField that implements probability calculations
    TRestAxionQCDField* fAxionField = nullptr;  //!

    /// A pointer to TRestBufferGas given to TRestAxionField to perform calculations in a particular gas
    TRestAxionBufferGas* fBufferGas = nullptr;  //!

    void Initialize() override;

   public:
    void InitProcess() override;

    RESTValue GetInputEvent() const override { return fAxionEvent; }
    RESTValue GetOutputEvent() const override { return fAxionEvent; }

    TRestEvent* ProcessEvent(TRestEvent* eventInput) override;

    /// It prints out the process parameters stored in the metadata structure
    void PrintMetadata() override {
        BeginPrintProcess();

        RESTMetadata << "Integration step length : " << fIntegrationStep << " mm" << RESTendl;
        RESTMetadata << "Buffer gas additional length : " << fBufferGasAdditionalLength * units("m") << " m"
                     << RESTendl;

        EndPrintProcess();
    }

    /// Returns the name of this process
    const char* GetProcessName() const override { return "axionFieldPropagation"; }

    // Constructor
    TRestAxionFieldPropagationProcess();
    TRestAxionFieldPropagationProcess(char* cfgFileName);

    // Destructor
    ~TRestAxionFieldPropagationProcess();

    ClassDefOverride(TRestAxionFieldPropagationProcess, 2);
};
#endif
