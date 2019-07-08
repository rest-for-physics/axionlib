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
#include "TRestEventProcess.h"
#include "TRestAxionMagneticField.h"
#include "TRestAxionPhotonConversion.h"

//! A process to introduce the axion-photon conversion probability in the signal generation chain
class TRestAxionFieldPropagationProcess : public TRestEventProcess {
private:

    void InitFromConfigFile();

    void Initialize();

    void LoadDefaultConfig();

/// A pointer to the specific TRestAxionEvent
    TRestAxionEvent *fInputAxionEvent; //!
    TRestAxionEvent *fOutputAxionEvent; //!

/// A pointer to the magnetic field stored in TRestRun
    TRestAxionMagneticField *fAxionMagneticField; //!

/// A pointer for the gamma conversion probability 
    TRestAxionPhotonConversion *fAxionPhotonConversion; //!



protected:

public:
    void InitProcess();

    TRestEvent *ProcessEvent( TRestEvent *evInput );

    void LoadConfig( std::string cfgFilename, std::string name = "" );

    /// It prints out the process parameters stored in the metadata structure
    void PrintMetadata()
        {
            BeginPrintProcess();

            EndPrintProcess();
        }

    std::vector <TVector3> FindOneVolume( TVector3 pos, TVector3 dir, Double_t minStep );
    std::vector <std::vector <TVector3> > FindFieldBoundaries( Double_t minStep = -1); 

    TVectorD GetFieldVector( TVector3 in, TVector3 out, Int_t N = 0 ); 

    /// Returns a new instance of this class
    TRestEventProcess *Maker() { return new TRestAxionFieldPropagationProcess; }

    /// Returns the name of this process
    TString GetProcessName() { return (TString) "axionFieldPropagation"; }

    //Constructor
    TRestAxionFieldPropagationProcess();
    TRestAxionFieldPropagationProcess( char *cfgFileName );

    //Destructor
    ~TRestAxionFieldPropagationProcess();

    ClassDef(TRestAxionFieldPropagationProcess, 1);
};
#endif
