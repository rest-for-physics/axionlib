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

#ifndef RestCore_TRestAxionTransportProcess
#define RestCore_TRestAxionTransportProcess

#include "TRestAxionEvent.h"
#include "TRestAxionEventProcess.h"

//! A process to transport the axion to a given z-position without changing direction
class TRestAxionTransportProcess : public TRestAxionEventProcess {
   private:
    /// The Z-position where we will move the particle
    Double_t fZPosition;  //<

    void LoadDefaultConfig();

   protected:
   public:
    void InitProcess() override;

    void Initialize() override;

    TRestEvent* ProcessEvent(TRestEvent* evInput) override;

    void LoadConfig(std::string cfgFilename, std::string name = "");

    /// It prints out the process parameters stored in the metadata structure
    void PrintMetadata() override {
        BeginPrintProcess();

        RESTMetadata << "Moving to Z: " << fZPosition << " mm" << RESTendl;

        EndPrintProcess();
    }

    /// Returns the name of this process
    const char* GetProcessName() const override { return "axionTransport"; }

    // Constructor
    TRestAxionTransportProcess();
    TRestAxionTransportProcess(char* cfgFileName);

    // Destructor
    ~TRestAxionTransportProcess();

    ClassDefOverride(TRestAxionTransportProcess, 1);
};
#endif
