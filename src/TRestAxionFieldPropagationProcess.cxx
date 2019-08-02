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

    fFinalNormalPlan = TVector3();
    fFinalPositionPlan = TVector3();
    fDistance = 0.0;
}

///////////////////////////////////////////////
/// \brief The main processing event function
///
void TRestAxionFieldPropagationProcess::InitProcess() {
    debug   << "Entering ... TRestAxionGeneratorProcess::InitProcess" << endl;

    fAxionMagneticField = (TRestAxionMagneticField*)this->GetMetadata("TRestAxionMagneticField");

    if (!fAxionMagneticField) {
        error << "TRestAxionFieldPropagationprocess. Magnetic Field was not defined!" << endl;
        exit(0);
    }

    fAxionBufferGas = (TRestAxionBufferGas*)this->GetMetadata("TRestAxionBufferGas");

    if (!fAxionBufferGas) {
        error << "TRestAxionBufferGas. Cannot access the buffer gas" << endl;
        exit(0);
    }

    fAxionPhotonConversion = new TRestAxionPhotonConversion();
    fAxionPhotonConversion->SetBufferGas(fAxionBufferGas);
}

TVector3 TRestAxionFieldPropagationProcess::MoveToPlan(TVector3 pos, TVector3 dir, Double_t f, Int_t i) 
{
    if ( dir[i] == 0 ) error << "The component of the direction you chose is equal to 0 " << endl;
    
    Double_t t = (f-pos[i])/dir[i];
    pos[i] = f;
    pos[(i+1)%3] = pos[(i+1)%3] + t * dir[(i+1)%3];    
    pos[(i+2)%3] = pos[(i+2)%3] + t * dir[(i+2)%3];

    return pos;
}

TVector3 TRestAxionFieldPropagationProcess::FinalPositionInPlan(TVector3 pos, TVector3 dir,TVector3 normalPlan,TVector3 pointPlan)
{
   if ( normalPlan.Dot(dir)==0 ) return pos;

   Double_t t, d ;
   d = normalPlan.Dot(pos);
   t = - ( normalPlan.Dot(pos) + d ) / ( normalPlan.Dot(dir) ) ;

   pos = pos + t * dir;
   
   return pos;
}

TVector3 TRestAxionFieldPropagationProcess::MoveToFinalDistance(TVector3 pos, TVector3 dir, Double_t distance)
{
   Double_t t = distance / dir.Mag();
   pos = pos + t * dir;

   return pos;
}

bool TRestAxionFieldPropagationProcess::IsInBoundedPlan(TVector3 pos, Int_t i, Int_t p) 
{
    Double_t minCond1, maxCond1, minCond2, maxCond2;
    Int_t j, k;

    if ( i==0 ) 
    {
         minCond1 = ( fAxionMagneticField->GetYmin() )[p]; maxCond1 = ( fAxionMagneticField->GetYmax() )[p];
         minCond2 = ( fAxionMagneticField->GetZmin() )[p]; maxCond2 = ( fAxionMagneticField->GetZmax() )[p];
         j = 1; k = 2;
    }
    if ( i==1 ) 
    {
         minCond1 = ( fAxionMagneticField->GetXmin() )[p]; maxCond1 = ( fAxionMagneticField->GetXmax() )[p];
         minCond2 = ( fAxionMagneticField->GetZmin() )[p]; maxCond2 = ( fAxionMagneticField->GetZmax() )[p];
         j = 0; k = 2;
    }
    if ( i==2 ) 
    {
         minCond1 = ( fAxionMagneticField->GetXmin() )[p]; maxCond1 = ( fAxionMagneticField->GetXmax() )[p];
         minCond2 = ( fAxionMagneticField->GetYmin() )[p]; maxCond2 = ( fAxionMagneticField->GetYmax() )[p];
         j = 0 ; k = 1;  
    }
  
    bool cond1 = ( pos[j] <= maxCond1 && pos[j] >= minCond1 );
    bool cond2 = ( pos[k] <= maxCond2 && pos[k] >= minCond2 );

    return cond1 && cond2;

}

std::vector <TVector3> TRestAxionFieldPropagationProcess::InOut(std::vector <TVector3> bounds, TVector3 dir ) 
{

    Int_t i = 0;

    TVector3 in = bounds[0];
    TVector3 out = bounds[1];

    while ( i < 3 && dir[i] * ( out[i]-in[i] ) >= 0 )
            i = i+1 ;
 
    if ( i < 3 )
    {
         bounds[0] = out;
         bounds[1] = in;
    }

    return bounds;
    
}

std::vector <TVector3 > TRestAxionFieldPropagationProcess::FindBoundariesOneVolume(TVector3 pos, TVector3 dir, Int_t p)
{
     std::vector <TVector3> boundaries;
     TVector3 in;
     TVector3 out;

     Int_t i = 0;
     Int_t j = 0;

     TVectorD f(6);
     f[0] = ( fAxionMagneticField -> GetXmin() )[p]; f[1] = ( fAxionMagneticField -> GetXmax() )[p];
     f[2] = ( fAxionMagneticField -> GetYmin() )[p]; f[3] = ( fAxionMagneticField -> GetYmax() )[p];
     f[4] = ( fAxionMagneticField -> GetZmin() )[p]; f[5] = ( fAxionMagneticField -> GetZmin() )[p];

     Int_t nFace = 0; 

     while ( nFace < 6 && i <= 0 )
     {
             if ( dir[Int_t(nFace/2)] != 0 )
             {
                  if ( (pos[Int_t(nFace/2)] - f[nFace]) * dir[Int_t(nFace/2)] <= 0 )
                  {  
                        in = MoveToPlan(pos,dir,f[nFace],Int_t(nFace/2)); 
                        if ( IsInBoundedPlan(in,Int_t(nFace/2),p) ) 
                             i = i+1; 
         
                  }
             }
             nFace = nFace + 1;
     
     }

     if ( i == 1 )   
     {
          while ( nFace < 6 && j <= 0 )
          {
                  if ( dir[Int_t(nFace/2)] != 0 )
                  {
                       if ( (pos[Int_t(nFace/2)] - f[nFace]) * dir[Int_t(nFace/2)] <= 0 )
                      {
                             out = MoveToPlan(pos,dir,f[nFace],Int_t(nFace/2)); 
                             if ( IsInBoundedPlan(out,Int_t(nFace/2),p) )  
                                  j = j+1;
                      }
                  }
                  nFace = nFace + 1;
          }
     }     
     
     if ( i+j == 2 && in != out)
     {
        boundaries.push_back(in);
        boundaries.push_back(out);
        return boundaries;
     }
    
     else return boundaries;
             
}
std::vector <TVector3> TRestAxionFieldPropagationProcess::FieldBoundary(std::vector <TVector3> boundaries, Double_t minStep)
{
    std::vector <TVector3> boundariesField;

    TVector3 in = boundaries[0];
    TVector3 out = boundaries[1];
    TVector3 diff = out-in;

    Int_t N=10;
    
    Int_t i = 0;

    while ( 1.0/Double_t(N) > minStep )
    {
            while ( ( fAxionMagneticField -> GetMagneticField(in[0], in[1], in[2]) ) == TVector3(0, 0, 0) && i<N)
            {
                    in = in + Double_t(i) / Double_t(N) * diff;
                    i = i+1;
            }
            
            if ( i==N ) return boundariesField;
            
            in = in - 1.0 / Double_t(N) * diff;
            N = 10 * N;
            i = 0;
     }

     in = in + 10.0 / Double_t (N) * diff;
     boundariesField.push_back(in);
   
     N = 10;
     i = 0;

    while ( 1.0/Double_t(N) > minStep )
    {
            while ( ( fAxionMagneticField -> GetMagneticField(in[0], in[1], in[2]) == TVector3(0,0,0) ) && i<N)
            {
                    out = out - Double_t(i) / Double_t(N) * diff;
                    i = i+1;
            }
            
            out = out + 1.0 / Double_t(N) * diff;
            N = 10 * N;
            i = 0;
     }

     out = out - 10.0 / Double_t (N) * diff;
     boundariesField.push_back(out);

     return boundariesField;        
}


std::vector<std::vector<TVector3>> TRestAxionFieldPropagationProcess::FindFieldBoundaries(Double_t minStep) 
{
    if (minStep == -1) minStep = 0.01;

    std::vector<std::vector<TVector3>> boundaryCollection;
    std::vector<std::vector<TVector3>> boundaryFinalCollection;
    std::vector<TVector3> buffVect;
    TVector3 boundaryIn;
    TVector3 boundaryOut;

    TVector3 posInitial = *(fInputAxionEvent->GetPosition());
    TVector3 direction = *(fInputAxionEvent->GetDirection());
    direction = direction.Unit();

    if (direction == TVector3(0, 0, 0) ) // No moves 
        return boundaryCollection;

    if ( ( fAxionMagneticField -> GetXmin() ).size() == 0 ) 
         error << " The magnetic field has not been loaded " << endl;
   
    /*** Find global boundaries volume ***/
    
    Int_t N = ( fAxionMagneticField -> GetXmin() ).size();


    for (Int_t p=0 ; p<N ; p++)
    {
          buffVect = FindBoundariesOneVolume( posInitial, direction, p);
          if ( buffVect.size() == 2 )
          {
              buffVect = InOut( buffVect,direction );
              boundaryCollection.push_back(buffVect);
              buffVect.clear();
          }
    }
   
    debug << "+------------------------+" << endl;
    debug   << " Number of volume boundaries : " << 2 * boundaryCollection.size() << endl;
    debug << "+------------------------+" << endl;

    debug << "+------------------------+" << endl;
    for (Int_t p=0 ; p < boundaryCollection.size() ; p++) 
    {
         debug   << "for" << p << " in : (" << boundaryCollection[p][0].X() << ","  << boundaryCollection[p][0].Y() << "," << boundaryCollection[p][0].Z() << ")" << endl;
         debug   << "for" << p << " out : (" << boundaryCollection[p][1].X() << ","  << boundaryCollection[p][1].Y() << "," << boundaryCollection[p][1].Z() << ")" << endl; 
    }
    debug << "+------------------------+" << endl;

   /*** Find precise boundaries field ***/ 
   
    for ( Int_t i = 0 ; i < boundaryCollection.size() ; i ++ )
    {
              buffVect = FieldBoundary( boundaryCollection[i],minStep) ;
              if ( buffVect.size() == 2)
              {
                   boundaryFinalCollection.push_back(buffVect);
                   buffVect.clear();
              }
    }

    boundaryCollection.clear();

    debug << "+------------------------+" << endl;
    debug   << " Number of field boundaries : " << 2 * boundaryFinalCollection.size() << endl;
    debug << "+------------------------+" << endl;

    debug << "+------------------------+" << endl;
    for (Int_t p=0 ; p < boundaryFinalCollection.size() ; p++) 
    {
         debug   << "for" << p << " in : (" << boundaryFinalCollection[p][0].X() << ","  << boundaryFinalCollection[p][0].Y() << "," << boundaryFinalCollection[p][0].Z() << ")" << endl;
         debug   << "for" << p << " out : (" << boundaryFinalCollection[p][1].X() << ","  << boundaryFinalCollection[p][1].Y() << "," << boundaryFinalCollection[p][1].Z() << ")" << endl; 
    }
    debug << "+------------------------+" << endl;
    
    return boundaryFinalCollection;  

}


TVectorD TRestAxionFieldPropagationProcess::GetFieldVector(TVector3 in, TVector3 out, Int_t N) {
    if (N == 0) N = TMath::Power(10, 4);

    TVectorD Bt(N);
    TVector3 B;
    TVector3 direction = *(fInputAxionEvent->GetDirection());

    TVector3 differential = out - in;

    B = fAxionMagneticField->GetMagneticField(in[0], in[1], in[2]);
    Bt[0] = abs(B.Perp(direction));

    for (Int_t i = 1; i < N; i++) {
        in = in + differential * (Double_t(i) / Double_t(N - 1));
        B = fAxionMagneticField->GetMagneticField(in[0], in[1], in[2]);
        Bt[i] = abs(B.Perp(direction));
    }

    return Bt;
}

TRestEvent* TRestAxionFieldPropagationProcess::ProcessEvent(TRestEvent* evInput) {
    fInputAxionEvent = (TRestAxionEvent*)evInput;
    fOutputAxionEvent = fInputAxionEvent;

    TVector3 position = *(fInputAxionEvent->GetPosition());
    TVector3 direction = *(fInputAxionEvent->GetDirection());

    debug << "+------------------------+" << endl;
    debug << "Initial position of the axion input : " << endl;
    debug << "(" << position.X() << ","<< position.Y()  << ","<< position.Z() << ")";
    debug << "Direction of the axion input : " << endl;
    debug << "(" << direction.X() << ","<< direction.Y() << ","<< direction.Z() << ")";
    debug << "+------------------------+" << endl;
     
    Double_t Ea = fInputAxionEvent->GetEnergy();
    Double_t ma = fInputAxionEvent->GetMass();
    
    std::vector<std::vector<TVector3>> boundaries;
    boundaries = FindFieldBoundaries();
    Int_t NofVolumes = boundaries.size();

    debug << "+------------------------+" << endl;
    debug  << "Number of magnetic field through which the axion passes : " << NofVolumes << endl; 
    debug << "+------------------------+" << endl;

    Double_t probability = 0.;

    Int_t N = TMath::Power(10,4);

    if ( NofVolumes !=0 )
    {
         TVectorD B(N);
         TVectorD probabilities(NofVolumes);

         TVector3 lengthVector ;
         Double_t length;

         for ( Int_t i = 0; i < NofVolumes; i++ ) {
               lengthVector = boundaries[i][0]-boundaries[i][1];
               length = sqrt( lengthVector.Mag2() );
               B = GetFieldVector(boundaries[i][0], boundaries[i][1], 0);
               probabilities[i] = fAxionPhotonConversion->GammaTransmissionProbability(Ea, B, ma, length);
         }

         for (Int_t i = 0; i < NofVolumes; i++) probability = probability + probabilities[i];
    }

    fOutputAxionEvent->SetGammaProbability(probability);

    if (fMode == "plan") fOutputAxionEvent->SetPosition( FinalPositionInPlan(position,direction,fFinalNormalPlan,fFinalPositionPlan) );
    if (fMode == "distance") fOutputAxionEvent->SetPosition( MoveToFinalDistance(position,direction,fDistance) );

    debug << "+------------------------+" << endl;
    debug << "Final position of the axion input : " << endl;
    debug << "(" << fOutputAxionEvent->GetPositionX() << ","<< fOutputAxionEvent->GetPositionY()<< ","<< fOutputAxionEvent->GetPositionZ() << ")";
    debug << "+------------------------+" << endl;

    if (GetVerboseLevel() >= REST_Debug   ) fOutputAxionEvent->PrintEvent();
    
    boundaries.clear();
    
    return fOutputAxionEvent;
}

///////////////////////////////////////////////
/// \brief Function reading input parameters from the RML TRestAxionFieldPropagationProcess metadata section
///
void TRestAxionFieldPropagationProcess::InitFromConfigFile() {
    this->Initialize();

    fMode = GetParameter("mode");
    fFinalNormalPlan = Get3DVectorParameterWithUnits("finalNPlan");
    fFinalPositionPlan = Get3DVectorParameterWithUnits("finalPositionPlan");
    fDistance = GetDoubleParameterWithUnits("distance");

    PrintMetadata();
}
