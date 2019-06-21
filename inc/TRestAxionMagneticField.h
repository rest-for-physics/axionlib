///______________________________________________________________________________
///______________________________________________________________________________
///______________________________________________________________________________
///
///
///             RESTSoft : Software for Rare Event Searches with TPCs
///
///             myMetadata.h
///
///             Metadata class to be used to store basic info
///             inherited from TRestMetadata
///
///             jun 2016:   First concept. Javier Galan
//
///_______________________________________________________________________________

#ifndef _TRestAxionMagneticField
#define _TRestAxionMagneticField

#include <iostream>
#include <TRestMetadata.h>

#include "TVector3.h"
#include "TVectorD.h"

#if defined USE_Garfield
#include <ComponentVoxel.hh>
#include <ComponentBase.hh>
#include <Sensor.hh>
using namespace Garfield;
#else
class Sensor;
class ComponentVoxel;
class ComponentBase;
#endif

using namespace std;

class TRestAxionMagneticField : public TRestMetadata {
private:
    void Initialize();

    void InitFromConfigFile();

    Sensor* fSetOfField;  //!
    ComponentVoxel* fMesh; //!

    Double_t fXmin; //<
    Double_t fXmax; //<
    Double_t fYmin; //<
    Double_t fYmax; //<
    Double_t fZmin; //<
    Double_t fZmax; //<
    Double_t sizeMesh; //<<
	
    Int_t fNofVolumes; //

    TString fFileName;
    TVector3 fPos;
    vector<TVector3> fPositions;

public:
     void LoadMagneticVolume(TVector3 Pos);
     void AddFieldComponent(); 
     void PrintMetadata();

    // Constructors
    TRestAxionMagneticField();
    TRestAxionMagneticField(const char* cfgFileName, std::string name = "");
    // Destructor
    ~TRestAxionMagneticField();

    ClassDef(TRestAxionMagneticField, 1);
};
#endif
