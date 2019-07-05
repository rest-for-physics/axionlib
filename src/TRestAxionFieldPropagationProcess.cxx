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
/// TRestAxionFieldPropagationProcess TOBE documented
///
/// The axion is generated with intensity proportional to g_ag = 1.0 x g10
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2019-March:  First implementation of shared memory buffer to rawsignal conversion.
///             Javier Galan
///
/// \class      TRestAxionFieldPropagationProcess
/// \author     Javier Galan
///
/// <hr>
///
#include "TRestAxionFieldPropagationProcess.h"
using namespace std;

ClassImp(TRestAxionFieldPropagationProcess);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionFieldPropagationProcess::TRestAxionFieldPropagationProcess() { Initialize(); }

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
TRestAxionFieldPropagationProcess::TRestAxionFieldPropagationProcess(char* cfgFileName) {
    Initialize();

    LoadConfig(cfgFileName);
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionFieldPropagationProcess::~TRestAxionFieldPropagationProcess() {
    delete fInputAxionEvent;
    delete fOutputAxionEvent;
}

///////////////////////////////////////////////
/// \brief Function to load the default config in absence of RML input
///
void TRestAxionFieldPropagationProcess::LoadDefaultConfig() {
    SetName(this->ClassName());
    SetTitle("Default config");
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
void TRestAxionFieldPropagationProcess::LoadConfig(std::string cfgFilename, std::string name) {
    if (LoadConfigFromFile(cfgFilename, name)) LoadDefaultConfig();
}

///////////////////////////////////////////////
/// \brief Function to initialize input/output event members and define the section name
///
void TRestAxionFieldPropagationProcess::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    fInputAxionEvent = new TRestAxionEvent();
    fOutputAxionEvent = new TRestAxionEvent();

    fInputEvent = fInputAxionEvent;
    fOutputEvent = fOutputAxionEvent;
}

///////////////////////////////////////////////
/// \brief The main processing event function
///

void TRestAxionFieldPropagationProcess::InitProcess() {

    debug << "Entering ... TRestAxionGeneratorProcess::InitProcess" << endl;

    fAxionMagneticField = (TRestAxionMagneticField*)this->GetMetadata("TRestAxionMagneticField");

    if (!fAxionMagneticField) {
        error << "TRestAxionFieldPropagationprocess. Magnetic Field was not defined!" << endl;
        exit(0);
    }

}

std::vector <TVector3> TRestAxionFieldPropagationProcess::FindOneVolume( TVector3 pos, TVector3 dir, Double_t y, Double_t minStep )
{

    if ( dir[1] > 0 ) 
         error << " y is ascendant, and you entered the direction vector : " << "("<<dir[0]<<","<<dir[1]<<","<<dir[2]<< ")" <<endl;

    std::vector <TVector3> boundary;
    TVector3 boundaryIn;
    TVector3 boundaryOut;

    double dr = 10.; // Can be between 1 and 10

    double t = (y-pos[1])/dir[1];
    pos[1] = y; 
    pos[0] = pos[0]+t*dir[0];
    pos[2] = pos[2]+t*dir[2];

  
    while ( dr > minStep )
    {
    	while ( fAxionMagneticField->GetMagneticField(pos[0],pos[1],pos[2]) == TVector3(0,0,0) && pos[1]>0 ) 
        {
                y = pos[1]-dr;
	        t = (y-pos[1])/dir[1]; 
                pos[1] = y; 
                pos[0] = pos[0]+t*dir[0];
                pos[2] = pos[2]+t*dir[2];
        }

	if ( pos[1] <= 0 ) { 

	     boundaryIn = pos;
 	     boundaryOut = TVector3(0,0,0);
             boundary.push_back(boundaryIn);
             boundary.push_back(boundaryOut);
             return boundary; }

	else {

	     y = pos[1]+dr;
             t = (y-pos[1])/dir[1]; 
             pos[1] = y; 
             pos[0] = pos[0]+t*dir[0];
             pos[2] = pos[2]+t*dir[2];
             dr = dr/2.0;
             }
    }
    
    boundaryIn = pos;

    while ( fAxionMagneticField->GetMagneticField(pos[0],pos[1],pos[2]) != TVector3(0,0,0) ) 
    {
                y = pos[1]-dr;
	        t = (y-pos[1])/dir[1]; 
                pos[1] = y; 
                pos[0] = pos[0]+t*dir[0];
                pos[2] = pos[2]+t*dir[2];
    }
    
    boundaryOut = pos;

    boundary.push_back(boundaryIn);
    boundary.push_back(boundaryOut);
               
    return boundary;

} 

std::vector <std::vector <TVector3>> TRestAxionFieldPropagationProcess::FindFieldBoundaries( Double_t minStep )
{

    if ( minStep == -1 )
         minStep = 0.01 ; 

    std::vector <std::vector <TVector3>> boundaryCollection;

    TVector3 posInitial = *(fInputAxionEvent->GetPosition());
    TVector3 direction = *(fInputAxionEvent->GetDirection());

    Double_t t = (10000.-posInitial[1])/direction[1]; // Translation of the axion at the plane y=10 m, it could be 20m, or 5m etc. ...
    posInitial[1] = 10000.; 
    posInitial[0] = posInitial[0]+t*direction[0];
    posInitial[2] = posInitial[2]+t*direction[2];

    std::vector <TVector3> bInt;
    bInt = FindOneVolume(posInitial,direction,posInitial[1],minStep);

    if ( bInt[0][1] <= 0. ) 
         return boundaryCollection ;

    else 
    {
         boundaryCollection.push_back( bInt );
         bInt = FindOneVolume(bInt[1],direction,bInt[1][1],minStep);

         while ( bInt[0][1] > 0 ) 
         {
                 boundaryCollection.push_back( bInt );
                 bInt = FindOneVolume(bInt[1],direction,bInt[1][1],minStep);
         }

         bInt.clear();  
         return boundaryCollection;
    }
   
}

TRestEvent* TRestAxionFieldPropagationProcess::ProcessEvent(TRestEvent* evInput) {
    fInputAxionEvent = (TRestAxionEvent*)evInput;

    *fOutputAxionEvent = *fInputAxionEvent;
     

    // fOutputAxionEvent->SetGammaProbability( 0.1 );

    /*if (GetVerboseLevel() >= REST_Debug) {
        fOutputAxionEvent->PrintEvent();
        GetChar();
    }*/

    return fOutputEvent;
}

///////////////////////////////////////////////
/// \brief Function reading input parameters from the RML TRestAxionFieldPropagationProcess metadata section
///
void TRestAxionFieldPropagationProcess::InitFromConfigFile() {




}
