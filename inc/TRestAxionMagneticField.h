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

#include <TRestMetadata.h>
#include <iostream>

#include "TVector3.h"
#include "TVectorD.h"

#if defined USE_Garfield
#include <ComponentBase.hh>
#include <ComponentVoxel.hh>
#include <Sensor.hh>
class Garfield::Sensor;
class Garfield::ComponentVoxel;
class Garfield::ComponentBase;
#endif

class TRestAxionMagneticField : public TRestMetadata {
   private:
    void Initialize();

    void InitFromConfigFile();

#ifdef USE_Garfield
    Garfield::Sensor* fSetOfField;    //!
#endif

    Int_t fNofVolumes;  //<

    std::vector<TVector3> fPositions;
    std::vector<TString> fFileNames;

    std::vector<Double_t> fXmax;
    std::vector<Double_t> fXmin;
    std::vector<Double_t> fYmax;
    std::vector<Double_t> fYmin;
    std::vector<Double_t> fZmin;
    std::vector<Double_t> fZmax;
    std::vector<Double_t> fSizeMesh;

   public:
    void LoadMagneticVolumes();
    void PrintMetadata();

    // Constructors
    TRestAxionMagneticField();
    TRestAxionMagneticField(const char* cfgFileName, std::string name = "");
    // Destructor
    ~TRestAxionMagneticField();

    ClassDef(TRestAxionMagneticField, 1);
};
#endif
