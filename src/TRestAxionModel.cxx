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
/// TRestAxionModels is a class used to define methods which define
/// theoretical model results that can be directly used in the calculation
/// of experimental observables.
///
/// TODO. Create an appropriate documentation here.
///
/// List of solar axion flux on Earth models available :
/// - arXiv_0702006
///
/// List of solar production mechanism models available :
/// - Primakoff
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2019-March: First concept and implementation of TRestAxionModel class.
///             Javier Galan
///
/// \class      TRestAxionModel
/// \author     Javier Galan
///
/// <hr>
///

#include "TRestAxionModel.h"
using namespace std;

ClassImp(TRestAxionModel)
//______________________________________________________________________________
TRestAxionModel::TRestAxionModel() : TRestMetadata()
{
   // TRestAxionModel default constructor
   Initialize();
}


//______________________________________________________________________________
TRestAxionModel::TRestAxionModel( const char *cfgFileName, string name ) : TRestMetadata (cfgFileName)
{
	cout << "Entering TRestAxionModel constructor( cfgFileName, name )" << endl;

    Initialize();

    LoadConfigFromFile( fConfigFileName, name );

    PrintMetadata();
}


//______________________________________________________________________________
TRestAxionModel::~TRestAxionModel()
{
    // TRestAxionModel destructor
}

void TRestAxionModel::Initialize()
{
	SetSectionName( this->ClassName() );
}

// Returns the solar axion flux in cm-2 keV-1 s-1 (on earth)
Double_t TRestAxionModel::GetDifferentialSolarAxionFlux( Double_t energy, Double_t g10 )
{
	// https://arxiv.org/abs/hep-ex/0702006
	if( fSolarAxionFluxModel == "arXiv_0702006" && fSolarProductionMechanism == "Primakoff" )
		return 6.02e10 * g10 * g10 * TMath::Power( energy, 2.481 ) * TMath::Exp( -energy/1.205 );

	warning << "Solar model not recognized" << endl;
	warning << "--------------------------" << endl;
	warning << "Solar axion model : " << fSolarAxionFluxModel << endl;
	warning << "Axion production : " << fSolarProductionMechanism << endl;

	return 0;

}

Double_t TRestAxionModel::GetSolarAxionFlux( Double_t eMin, Double_t eMax, Double_t g10, Double_t step )
{
	if( fSolarEnergyFlux > 0 )
		if( fg10 == g10 && fStep == step && eMin == fEnergyRange.X() && eMax == fEnergyRange.Y() )
			return fSolarEnergyFlux;

	info << "TRestAxionModel::GetSolarAxionFlux re-calculating solar axion flux" << endl;

	fg10 = g10;
	fStep = step;
	fEnergyRange = TVector2( eMin, eMax );

	fSolarEnergyFlux = 0;
	for( Double_t en = eMin; en < eMax; en += step )
		fSolarEnergyFlux += GetDifferentialSolarAxionFlux( en, g10 );

	fSolarEnergyFlux = fSolarEnergyFlux * step;

	return fSolarEnergyFlux;
}


//______________________________________________________________________________
void TRestAxionModel::InitFromConfigFile()
{
    this->Initialize();

    // Initialize the metadata members from a configfile
	
    // fClassMember = GetParameter( "paramName", "defaultValue" );
	
    fSolarAxionFluxModel = GetParameter( "solarAxionModel", "arXiv_0702006" );
	fSolarProductionMechanism = GetParameter( "solarProductionMechanism", "Primakoff" );

    PrintMetadata();
}

void TRestAxionModel::PrintMetadata( )
{
	TRestMetadata::PrintMetadata();

	metadata << " - Solar axion flux model : " << fSolarAxionFluxModel << endl;
	metadata << " - Solar axion production mechanism : " << fSolarProductionMechanism << endl;
	metadata << "-------------------------------------------------" << endl;
	metadata << " - Axion-photon couping : " << fg10 << " x 10^{-10} GeV^{-1}" << endl;
	metadata << " - Integration step : " << fStep << " keV" << endl;
	metadata << " - Integration range : ( " << fEnergyRange.X() << ", " << fEnergyRange.Y() << " ) keV" << endl;
	metadata << " - Calculated solar flux : " << fSolarEnergyFlux << " cm-2 s-1" << endl;
	metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

}

