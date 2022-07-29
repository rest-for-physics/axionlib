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

#ifndef RestCore_TRestAxionOpticsProcess
#define RestCore_TRestAxionOpticsProcess

#include "TRestAxionEvent.h"
#include "TRestAxionEventProcess.h"
#include "TRestAxionOptics.h"

//! A process to introduce the response from optics in the axion signal generation chain
class TRestAxionOpticsProcess : public TRestAxionEventProcess {
   private:
    /// A variable to determine if the new axis will be optical or universal axis
    Bool_t fOpticalAxis = false;  //<

    /// A pointer to the optics description defined inside TRestRun
    TRestAxionOptics* fOptics;  //!

    void Initialize() override;

    void LoadDefaultConfig();

   protected:
   public:
    void InitProcess() override;

    /// This InitFromConfigFile could be removed as soon as PR rest-for-physics/framework#275
    /// is solved
    void InitFromConfigFile() override {
        TRestEventProcess::InitFromConfigFile();

        if (ToUpper(GetParameter("opticalAxis", "false")) == "TRUE")
            fOpticalAxis = true;
        else
            fOpticalAxis = false;
    }

    TRestEvent* ProcessEvent(TRestEvent* evInput) override;

    /// End of event process. Called directly after ProcessEvent()
    void EndOfEventProcess(TRestEvent* evInput = nullptr) override;

    void LoadConfig(std::string cfgFilename, std::string name = "");

    /// It prints out the process parameters stored in the metadata structure
    void PrintMetadata() override {
        BeginPrintProcess();

        if (fOpticalAxis)
            RESTMetadata << "Output particle is described respect to the optical axis" << RESTendl;
        else
            RESTMetadata << "Output particle is described respect to the universal axis" << RESTendl;

        EndPrintProcess();
    }

    /// Returns the name of this process
    const char* GetProcessName() const override { return "axionOptics"; }

    // Constructor
    TRestAxionOpticsProcess();
    TRestAxionOpticsProcess(char* cfgFileName);

    // Destructor
    ~TRestAxionOpticsProcess();

    ClassDefOverride(TRestAxionOpticsProcess, 1);
};
#endif
