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

#include "TVector3.h"
#include "TVectorD.h"

#include "TRestAxionEvent.h"
#include "TRestAxionEventProcess.h"
#include "TRestAxionMagneticField.h"
#include "TRestPhysics.h"

//! A process to introduce the magnetic field profile integration along the track
class TRestAxionFieldPropagationProcess : public TRestAxionEventProcess {
   private:
    void Initialize();

    /// A pointer to the magnetic field stored in TRestRun
    TRestAxionMagneticField* fAxionMagneticField;  //!

   protected:
   public:
    void InitProcess();

    any GetInputEvent() { return fAxionEvent; }
    any GetOutputEvent() { return fAxionEvent; }

    TRestEvent* ProcessEvent(TRestEvent* evInput);

    /// It prints out the process parameters stored in the metadata structure
    void PrintMetadata() {
        BeginPrintProcess();

        // To be implemented

        EndPrintProcess();
    }

    /// Returns a new instance of this class
    TRestEventProcess* Maker() { return new TRestAxionFieldPropagationProcess; }

    /// Returns the name of this process
    TString GetProcessName() { return (TString) "axionFieldPropagation"; }

    // Constructor
    TRestAxionFieldPropagationProcess();
    TRestAxionFieldPropagationProcess(char* cfgFileName);

    // Destructor
    ~TRestAxionFieldPropagationProcess();

    ClassDef(TRestAxionFieldPropagationProcess, 1);
};
#endif
