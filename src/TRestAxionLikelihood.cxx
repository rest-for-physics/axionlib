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
/// TRestAxionLikelihood just a class to serve as example of a specific
/// implementation of a metadata class.
///
/// Just copy .cxx and .h files to a new class name,
/// cp TRestAxionLikelihood.cxx TRestAxionSpecificMetadata.cxx
/// cp TRestAxionLikelihood.h TRestAxionSpecificMetadata.h
///
/// Then, inside each file replace TRestAxionLikelihood by
/// TRestAxionSpecificMetadata.
///
/// Re-run "cmake ../" at the build directory to include the new
/// class to the compilation.
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2019-March: First concept and implementation of TRestAxionLikelihood class.
///             Javier Galan
///
/// \class      TRestAxionLikelihood
/// \author     Javier Galan
///
/// <hr>
///

#include "TRestAxionLikelihood.h"
using namespace std;

ClassImp(TRestAxionLikelihood)
//______________________________________________________________________________
TRestAxionLikelihood::TRestAxionLikelihood() : TRestMetadata()
{
   // TRestAxionLikelihood default constructor
   Initialize();
}


//______________________________________________________________________________
TRestAxionLikelihood::TRestAxionLikelihood( const char *cfgFileName, string name ) : TRestMetadata (cfgFileName)
{
	cout << "Entering TRestAxionLikelihood constructor( cfgFileName, name )" << endl;

    Initialize();

    LoadConfigFromFile( fConfigFileName, name );

    PrintMetadata();
}

//______________________________________________________________________________
TRestAxionLikelihood::~TRestAxionLikelihood()
{
    // TRestAxionLikelihood destructor
}

void TRestAxionLikelihood::Initialize()
{
	SetSectionName( this->ClassName() );
}

//______________________________________________________________________________
void TRestAxionLikelihood::InitFromConfigFile()
{
	this->Initialize();

	// Initialize the metadata members from a configfile
	fBmag = GetDblParameterWithUnits( "Bmag", 0. );
	fRmag = GetDblParameterWithUnits( "Rmag", 0. );
	fLmag = GetDblParameterWithUnits( "Lmag", 0. );

	fEfficiency = StringToDouble( GetParameter( "efficiency", "1" ) );

	fBackgroundLevel = StringToDouble( GetParameter( "bckLevel", "1.e-7" ) ); // cpd keV-1 s-1

	fTExpVacuum = StringToDouble( GetParameter( "expTimeVacuum", "1000" ) ); // In hours (vacuum phase)

	fTExpPerStep = StringToDouble( GetParameter( "expTimePerStep", "10" ) ); // In hours (gas phase, time at each step)

	fNSteps = StringToInteger( GetParameter( "pressureSteps", "100" ) );
}

void TRestAxionLikelihood::PrintMetadata( )
{
	TRestMetadata::PrintMetadata();

	metadata << " Magnetic field : " << fBmag << " T" << endl;
	metadata << " Magnet radius : " << fRmag << " mm" << endl;
	metadata << " Magnet length : " << fLmag << " mm" << endl;
	metadata << " Signal overall efficiency : " << fEfficiency << endl;
	metadata << " Background level : " << fBackgroundLevel << " cpd keV-1 s-1" << endl;
	metadata << " Vacuum phase exposure time : " << fTExpVacuum << " hours" << endl;
	metadata << " Gas phase exposure time per step : " << fTExpPerStep << endl;
	metadata << " Total number of steps : " << fNSteps << endl;

	metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

