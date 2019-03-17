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

#ifndef _TRestAxionModel
#define _TRestAxionModel

#include <TRestMetadata.h>

//! A metadata class to define theoretical axion models and calculations related
class TRestAxionModel:public TRestMetadata {
    private:
        void Initialize();

        void InitFromConfigFile();

        TString fSolarAxionFluxModel; //->
		TString fSolarProductionMechanism; //->

		// Integrated solar axion spectrum in cm-2 s-1
		Double_t fSolarEnergyFlux = 0; //->

		// The axion-photon g10 coupling 
		Double_t fg10 = 1.; //->

		// The integration step for solar energy spectrum
		Double_t fStep = 1.e-3; //->

		TVector2 fEnergyRange; //->

    public:
		
        TString GetSolarAxionModel() { return fSolarAxionFluxModel; }
        TString GetSolarProductionMechanism() { return fSolarProductionMechanism; }

		void SetSolarAxionModel( TString modelName ) { fSolarAxionFluxModel = modelName; }
		void SetSolarProductionMechanism( TString mName ) { fSolarProductionMechanism = mName; }

		void ResetSolarEnergyFlux() { fSolarEnergyFlux = 0; }

		Double_t GetSolarAxionFlux( Double_t eMin = 0., Double_t eMax = 10., Double_t g10 = 1., Double_t step = 0.001 );
		Double_t GetDifferentialSolarAxionFlux( Double_t energy, Double_t g10 = 1. );

        void PrintMetadata( );

        //Constructors
        TRestAxionModel();
        TRestAxionModel( const char *cfgFileName, std::string name = "");
        //Destructor
        ~TRestAxionModel();

        ClassDef(TRestAxionModel, 1); 
};
#endif
