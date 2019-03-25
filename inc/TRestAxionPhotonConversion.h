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

#ifndef _TRestAxionPhotonConversion
#define _TRestAxionPhotonConversion

#include <TRestMetadata.h>

#include "TRestAxionBufferGas.h"

//! A metadata class to define analytical axion-photon conversion probabilities for axion helioscopes
class TRestAxionPhotonConversion:public TRestMetadata {
    private:
        void Initialize();

        void InitFromConfigFile();

		// Axion mass in eV
		Double_t fAxionMass = 0; //->

		// Coherence length in mm [REST default units]
		Double_t fCohLength = 0; //->

		// Magnet field intensity in T
		Double_t fBMag = 0; //->

		// The axion-photon g10 coupling 
		Double_t fg10 = 1.; //->

		// A pointer to the buffer gas definition
		TRestAxionBufferGas *fBufferGas = NULL; //!

    public:

		void AssignBufferGas( TRestAxionBufferGas *buffGas ) { fBufferGas = buffGas; }

		void SetAxionMass( Double_t m ) { fAxionMass = m; }
		void SetCoherenceLength( Double_t l ) { fCohLength = l; }
		void SetMagneticField( Double_t B ) { fBMag = B; }

		Double_t GetAxionMass( ) { return fAxionMass; }
		Double_t GetCoherenceLength( ) { return fCohLength; }
		Double_t GetMagneticField( ) { return fBMag; }

		// (BL/2)**2
		Double_t BLFactor(Double_t Lcoh = -1, Double_t Bmag = -1 );

		/// ma in eV, Ea in keV, Length in cm
		Double_t GammaTransmissionProbability( Double_t Ea, Double_t ma = -1, Double_t Lcoh = -1, Double_t Bmag = -1 );
		
        void PrintMetadata( );

        //Constructors
        TRestAxionPhotonConversion();
        TRestAxionPhotonConversion( const char *cfgFileName, std::string name = "");
        //Destructor
        ~TRestAxionPhotonConversion();

        ClassDef(TRestAxionPhotonConversion, 1); 
};
#endif
