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
/// TRestAxionOpticsResponseProcess TOBE documented
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2019-March:  First implementation of a dummy optics response
///             Javier Galan
///
/// \class      TRestAxionOpticsResponseProcess
/// \author     
///
/// <hr>
///
#include "TRestAxionOpticsResponseProcess.h"
using namespace std;

ClassImp( TRestAxionOpticsResponseProcess )

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionOpticsResponseProcess::TRestAxionOpticsResponseProcess()
{
    Initialize();
}

///////////////////////////////////////////////
/// \brief Constructor loading data from a config file
/// 
/// If no configuration path is defined using TRestMetadata::SetConfigFilePath
/// the path to the config file must be specified using full path, absolute or relative.
///
/// The default behaviour is that the config file must be specified with 
/// full path, absolute or relative.
///
/// \param cfgFileName A const char* giving the path to an RML file.
///
TRestAxionOpticsResponseProcess::TRestAxionOpticsResponseProcess( char *cfgFileName )
{
    Initialize();

    LoadConfig( cfgFileName );
}

///////////////////////////////////////////////
/// \brief Default destructor 
/// 
TRestAxionOpticsResponseProcess::~TRestAxionOpticsResponseProcess()
{
    delete fInputAxionEvent;
    delete fOutputAxionEvent;
}


///////////////////////////////////////////////
/// \brief Function to load the default config in absence of RML input
/// 
void TRestAxionOpticsResponseProcess::LoadDefaultConfig( )
{
    SetName( this->ClassName() );
    SetTitle( "Default config" );
}

///////////////////////////////////////////////
/// \brief Function to load the configuration from an external configuration file.
/// 
/// If no configuration path is defined in TRestMetadata::SetConfigFilePath
/// the path to the config file must be specified using full path, absolute or relative.
///
/// \param cfgFileName A const char* giving the path to an RML file.
/// \param name The name of the specific metadata. It will be used to find the 
/// correspondig TRestGeant4AnalysisProcess section inside the RML.
///
void TRestAxionOpticsResponseProcess::LoadConfig( std::string cfgFilename, std::string name )
{
    if( LoadConfigFromFile( cfgFilename, name ) ) LoadDefaultConfig( );
}

///////////////////////////////////////////////
/// \brief Function to initialize input/output event members and define the section name
/// 
void TRestAxionOpticsResponseProcess::Initialize()
{
    SetSectionName( this->ClassName() );

    fInputAxionEvent = new TRestAxionEvent();
    fOutputAxionEvent = new TRestAxionEvent();

    fInputEvent = fInputAxionEvent;
    fOutputEvent = fOutputAxionEvent;
}

///////////////////////////////////////////////
/// \brief The main processing event function
/// 
TRestEvent* TRestAxionOpticsResponseProcess::ProcessEvent( TRestEvent *evInput )
{
	fInputAxionEvent = (TRestAxionEvent *) evInput;

	*fOutputAxionEvent = *fInputAxionEvent;

	if( GetVerboseLevel() >= REST_Debug ) 
	{
		fOutputAxionEvent->PrintEvent();

		if ( GetVerboseLevel() >= REST_Extreme )
			GetChar();
	}

    return fOutputEvent;
}

///////////////////////////////////////////////
/// \brief Function reading input parameters from the RML TRestAxionOpticsResponseProcess metadata section
/// 
void TRestAxionOpticsResponseProcess::InitFromConfigFile( )
{
}

