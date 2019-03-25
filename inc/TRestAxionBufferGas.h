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

#ifndef _TRestAxionBufferGas
#define _TRestAxionBufferGas

#include <TRestMetadata.h>

//! A metadata class to define theoretical axion models and calculations related
class TRestAxionBufferGas:public TRestMetadata {
    private:
        void Initialize();

        void InitFromConfigFile();

		// Name of the buffer gas (He, Ne, Ar, Xe, ..., etc )
		std::vector <TString> fBufferGasName; //->

		// Gas density of the corresponding gasName in g/cm3
		std::vector <Double_t> fBufferGasDensity; //->

		// Energy values for gas absorption coefficient
		std::vector <std::vector <Double_t> > fAbsEnergy; //->

		// Gas absorption coefficient in cm2/g
		std::vector <std::vector <Double_t> > fGasAbsCoefficient; //->

		// Energy values for gas form factor
		std::vector <std::vector <Double_t> > fFactorEnergy; //->
		
		// Gas form factor
		std::vector <std::vector <Double_t> > fGasFormFactor; //->

		void ReadGasData( TString gasName );

		Int_t FindGasIndex( TString gName );
		Int_t GetEnergyIndex( std::vector <Double_t> enVector, Double_t energy );

    public:

		void SetGasDensity( TString gasName, Double_t density );
		Double_t GetGasDensity( TString gasName );

		// Returns value in cm2/g
		Double_t GetAbsorptionCoefficient( TString gasName, Double_t energy );
		Double_t GetFormFactor( TString gasName, Double_t energy );

		// Returns absorption length in cm-1
		Double_t GetPhotonAbsorptionLength( Double_t energy );

		// Returns absorption length in eV
		Double_t GetPhotonAbsorptionLengthIneV( Double_t energy );

		// Transforms cm-1 to eV (This method might be transferred to REST_Physics)
		Double_t cmToeV( double l_Inv );

		// in eV
		Double_t GetPhotonMass( double en );

		void PrintAbsorptionGasData( TString gasName );
		void PrintFormFactorGasData( TString gasName );

        void PrintMetadata( );

        //Constructors
        TRestAxionBufferGas();
        TRestAxionBufferGas( const char *cfgFileName, std::string name = "");
        //Destructor
        ~TRestAxionBufferGas();

        ClassDef(TRestAxionBufferGas, 1); 
};
#endif
