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

//////////////////////////////////////////////////////////////////////////
/// TRestAxionPhotonConversion is a class used to ...
///
/// ... for the moment we implement here the equations from van Bibber paper
///
/// TODO. Create an appropriate documentation here.
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2019-March: First concept and implementation of TRestAxionPhotonConversion class.
///             Javier Galan
///
/// \class      TRestAxionPhotonConversion
/// \author     Javier Galan
///
/// <hr>
///

#include "TRestAxionPhotonConversion.h"
using namespace std;

using namespace REST_Physics;

ClassImp(TRestAxionPhotonConversion)
//______________________________________________________________________________
TRestAxionPhotonConversion::TRestAxionPhotonConversion() : TRestMetadata()
{
   // TRestAxionPhotonConversion default constructor
   Initialize();
}


//______________________________________________________________________________
TRestAxionPhotonConversion::TRestAxionPhotonConversion( const char *cfgFileName, string name ) : TRestMetadata (cfgFileName)
{
	cout << "Entering TRestAxionPhotonConversion constructor( cfgFileName, name )" << endl;

    Initialize();

    LoadConfigFromFile( fConfigFileName, name );

    PrintMetadata();
}


//______________________________________________________________________________
TRestAxionPhotonConversion::~TRestAxionPhotonConversion()
{
    // TRestAxionPhotonConversion destructor
}

void TRestAxionPhotonConversion::Initialize()
{
	SetSectionName( this->ClassName() );
}

double TRestAxionPhotonConversion::BLFactor( Double_t Lcoh, Double_t Bmag ) // (BL/2)**2
{
	// If we provide the values as argument we update the class members,
	// If not, we use the existing values in the class member
	if( Lcoh == -1 ) Lcoh = fCohLength; else fCohLength = Lcoh;
	if( Bmag == -1 ) Bmag = fBMag; else fBMag = Bmag;

	Double_t tm = lightSpeed / naturalElectron * 1.0e-9; // gev
	Double_t sol = Lcoh * Bmag * tm / 2;
	sol = sol* sol * 1.0e-20;

	return sol;
}

/// mgamma will be obtainned from buffer gas definition
///
/// If ma, Lcoh, Bmag = -1. Value will be taken from the class members definition.
/// Otherwise, members of the class will be updated
///
/// Ea in keV, ma in eV, mgamma in eV, Lcoh in m, Bmag in T
Double_t TRestAxionPhotonConversion::GammaTransmissionProbability( Double_t Ea, Double_t ma, Double_t Lcoh, Double_t Bmag )
{
	// If we provide the values as argument we update the class members,
	// If not, we use the existing values in the class member
	if( ma == -1 ) ma = fAxionMass; else fAxionMass = ma;
	if( Lcoh == -1 ) Lcoh = fCohLength; else fCohLength = Lcoh;
	if( Bmag == -1 ) Bmag = fBMag; else fBMag = Bmag;

	if( Bmag == 0 )
		warning << "TRestAxionPhotonConversion::GammaTransmissionProbability. Bmag = 0" << endl;
	if( Lcoh == 0 )
		warning << "TRestAxionPhotonConversion::GammaTransmissionProbability. Lcoh = 0" << endl;

	Double_t photonMass = 0.;
	if( !fBufferGas )
	{
		warning << "TRestAxionPhotonConversion::GammaTransmissionProbability. No buffer gas definition found!" << endl;
		warning << "Assuming vacuum medium. mgamma = 0" << endl;
	}
	else
	{
		photonMass = fBufferGas->GetPhotonMass( Ea );
	}

	debug << "+--------------------------------------------------------------------------+" << endl;
	debug << " TRestAxionPhotonConversion::GammaTransmissionProbability. Parameter summary" << endl;
    debug << " Photon mass : " << photonMass << " eV" << endl;
	debug << " Axion mass : " << ma << " eV" << endl;
   	debug << " Axion energy : " << Ea << " keV" << endl;
    debug << " Lcoh : " << Lcoh << " m" << endl;
    debug << " Bmag : " << Bmag << " T" << endl;
	debug << "+--------------------------------------------------------------------------+" << endl;

	if( ma == 0.0 && photonMass == 0.0 ) return BLFactor( );

	double q = (ma*ma - photonMass*photonMass)/2./Ea/1000.0;
	double l = Lcoh * PhMeterIneV;
	double phi = q * l;


	Double_t Gamma = 0.;
	if( fBufferGas )
		Gamma = fBufferGas->GetPhotonAbsorptionLength( Ea ); //cm-1

	Double_t GammaL = Gamma * Lcoh * 100;

	debug << "+------------------------+" << endl;
	debug << " Intermediate calculations" << endl;
    debug << " q : " << q << " eV" << endl;
    debug << " l : " << l << " eV-1" << endl;
	debug << " phi : " << phi << endl;
    debug << "Gamma : " << Gamma << endl;
    debug << "GammaL : " << GammaL << endl;
	debug << "+------------------------+" << endl;


	double MFactor = phi*phi + GammaL*GammaL/4.0;
	MFactor = 1.0/MFactor;

	double sol = MFactor * BLFactor( ) * ( 1 + exp(-GammaL) - 2*exp(-GammaL/2)*cos(phi));

	debug << "Axion-photon transmission probability : " << sol << endl;

	return sol;
}

//______________________________________________________________________________
void TRestAxionPhotonConversion::InitFromConfigFile()
{
    this->Initialize();

    // Initialize the metadata members from a configfile
	
    // fClassMember = GetParameter( "paramName", "defaultValue" );
	
    PrintMetadata();
}

void TRestAxionPhotonConversion::PrintMetadata( )
{
	TRestMetadata::PrintMetadata();

	metadata << " - Magnetic field : " << fBMag << " T" << endl;
	metadata << " - Coherence length : " << fCohLength * REST_Units::cm << " cm" << endl;
	metadata << " - Axion mass : " << fAxionMass << " eV" << endl;
	metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

	if( fBufferGas )
		fBufferGas->PrintMetadata();
}

