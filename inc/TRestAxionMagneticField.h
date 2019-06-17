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

#include "TRestAxionMagneticField.h"
#include "TVectorD.h"
#include <iostream>

using namespace std;

const double fXmin=-0.24;
const double fXmax=0.24;
const double fYmin=-0.24;
const double fYmax=0.24;
const double fZmin=-5.0;
const double fZmax=5.0;

class TRestAxionMagneticField : public TRestMetadata {
private:
    void Initialize();

    void InitFromConfigFile();

    string fType;
    string fFileName;
    Double_t fPosX;
    Double_t fPosY;
    Double_t fPosZ;
    Double_t fRotX;
    Double_t fRotY;
    Double_t fRotZ;

public:
     void LoadMagneticVolumeRegPar();
     void PrintMetadata();

    // Constructors
    TRestAxionMagneticField();
    TRestAxionMagneticField(const char* cfgFileName, std::string name = "");
    // Destructor
    ~TRestAxionMagneticField();

    ClassDef(TRestAxionMagneticField, 1);
};
#endif
