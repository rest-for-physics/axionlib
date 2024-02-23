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

#ifndef RestCore_TRestAxionAnalysisProcess
#define RestCore_TRestAxionAnalysisProcess

#include "TRestAxionEvent.h"
#include "TRestEventProcess.h"

//! An analyis process to add TRestAxionEvent observables to the analysis tree
class TRestAxionAnalysisProcess : public TRestEventProcess {
   private:
    /// A pointer to the specific TRestAxionEvent
    TRestAxionEvent* fAxionEvent;  //!

    /// The analysis position in mm with regards to the sun at (0,0,-AU).
    TVector3 fAnalysisPosition = TVector3(0, 0, 0);  //<

    void Initialize() override;

    void LoadDefaultConfig();

   protected:
   public:
    RESTValue GetInputEvent() const override { return fAxionEvent; }
    RESTValue GetOutputEvent() const override { return fAxionEvent; }

    TRestEvent* ProcessEvent(TRestEvent* evInput) override;

    void LoadConfig(std::string cfgFilename, std::string name = "");

    /// It prints out the process parameters stored in the metadata structure
    void PrintMetadata() override {
        BeginPrintProcess();

        EndPrintProcess();
    }

    /// Returns the name of this process
    const char* GetProcessName() const override { return "axionAnalysis"; }

    // Constructor
    TRestAxionAnalysisProcess();
    TRestAxionAnalysisProcess(char* cfgFileName);

    // Destructor
    ~TRestAxionAnalysisProcess();

    ClassDefOverride(TRestAxionAnalysisProcess, 1);
};
#endif
