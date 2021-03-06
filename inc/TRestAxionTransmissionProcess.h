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

#ifndef RestCore_TRestAxionTransmissionProcess
#define RestCore_TRestAxionTransmissionProcess

#include "TRestAxionEvent.h"
#include "TRestEventProcess.h"

//! A process to include photon transmission from different interfaces found till reaching the detector. E.g.
//! differential vaccuum windows
class TRestAxionTransmissionProcess : public TRestEventProcess {
   private:
    /// A pointer to the specific TRestAxionEvent
    TRestAxionEvent* fAxionEvent;  //!

    void InitFromConfigFile();

    void Initialize();

    void LoadDefaultConfig();

   protected:
   public:
    any GetInputEvent() { return fAxionEvent; }
    any GetOutputEvent() { return fAxionEvent; }

    TRestEvent* ProcessEvent(TRestEvent* evInput);

    void LoadConfig(std::string cfgFilename, std::string name = "");

    /// It prints out the process parameters stored in the metadata structure
    void PrintMetadata() {
        BeginPrintProcess();

        EndPrintProcess();
    }

    /// Returns a new instance of this class
    TRestEventProcess* Maker() { return new TRestAxionTransmissionProcess; }

    /// Returns the name of this process
    TString GetProcessName() { return (TString) "axionTransmission"; }

    // Constructor
    TRestAxionTransmissionProcess();
    TRestAxionTransmissionProcess(char* cfgFileName);

    // Destructor
    ~TRestAxionTransmissionProcess();

    ClassDef(TRestAxionTransmissionProcess, 1);
};
#endif
