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

#include "TRestAxionEvent.h"
#include "TRestEventProcess.h"

#include "TRestAxionSpectrum.h"

#include "TRandom3.h"

//! A process to generate axions following a particular solar axion model
class TRestAxionGeneratorProcess : public TRestEventProcess {
   private:
    /// A pointer to the specific TRestAxionEvent output
    TRestAxionEvent* fOutputAxionEvent;  //!

    /// Used internally to define the event id
    Int_t fCounter = 0;  //!

    /// A pointer to the axion model stored in TRestRun
    TRestAxionSpectrum *fAxionSpectrum; //!

    /// Random number generator
    TRandom3* fRandom;  //!

    /// The solar energy range
    TVector2 fEnergyRange;  //->

    /// Energy step
    Double_t fEnergyStep;  //->

    /// Axion mass
    Double_t fAxionMass;  //->

    /// The angular distribution generator type
    TString fAngularDistribution;  //->

    /// The main direction of the angular distribution
    TVector3 fAngularDirection;  //->

    /// The spatial distribution generator type
    TString fSpatialDistribution;  //->

    /// The spatial distribution generator type
    Double_t fSpatialRadius;  //->

    /// The spatial origin from the spatial distribution
    TVector3 fSpatialOrigin;  //->

    /// Rotation vector (according to x,y,z) for the circle wall defined from the plan (X,Y) - insead of
    /// rotating magnet Volumes
    TVector3 fRotation;  //->

    /// The normal vector of the wall (other way to do a rotation)
    TVector3 fNormalPlan;  //->

    /// Mode of rotated circle Wall construction
    TString fMode;  //->

    void InitFromConfigFile();

    void Initialize();

    void LoadDefaultConfig();

    Double_t GenerateEnergy();
    TVector3 GeneratePosition();
    TVector3 GenerateDirection();

   protected:
   public:
    void InitProcess();

    any GetInputEvent() { return (TRestEvent*)NULL; }
    any GetOutputEvent() { return fOutputAxionEvent; }

    TRestEvent* ProcessEvent(TRestEvent* eventInput);

    void LoadConfig(std::string cfgFilename, std::string name = "");

    /// It prints out the process parameters stored in the metadata structure
    void PrintMetadata() {
        BeginPrintProcess();

        metadata << "Energy distribution" << endl;
        metadata << "---------------------" << endl;
        metadata << "Energy range : (" << fEnergyRange.X() << ", " << fEnergyRange.Y() << ") keV" << endl;
        metadata << "Energy step : " << fEnergyStep << " keV" << endl;
        metadata << " " << endl;

        metadata << "Angular distribution" << endl;
        metadata << "----------------------" << endl;
        metadata << "Type : " << fAngularDistribution << endl;
        metadata << "Main direction : (" << fAngularDirection.X() << "," << fAngularDirection.Y() << ","
                 << fAngularDirection.Z() << ")" << endl;
        metadata << " " << endl;
        metadata << "Spatial distribution" << endl;
        metadata << "----------------------" << endl;
        metadata << "Type : " << fSpatialDistribution << endl;
        metadata << "Radius : " << fSpatialRadius << " mm" << endl;
        metadata << "Origin : (" << fSpatialOrigin.X() << "," << fSpatialOrigin.Y() << ","
                 << fSpatialOrigin.Z() << ")" << endl;

        EndPrintProcess();
    }

    /// Returns a new instance of this class
    TRestEventProcess* Maker() { return new TRestAxionGeneratorProcess; }

    /// Returns the name of this process
    TString GetProcessName() { return (TString) "axionGenerator"; }

    // Constructor
    TRestAxionGeneratorProcess();
    TRestAxionGeneratorProcess(char* cfgFileName);

    // Destructor
    ~TRestAxionGeneratorProcess();

    ClassDef(TRestAxionGeneratorProcess, 1);
};
#endif
