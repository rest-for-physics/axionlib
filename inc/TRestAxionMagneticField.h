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

#include "TCanvas.h"
#include "TH2D.h"
#include "TVector3.h"
#include "TVectorD.h"

#if defined USE_Garfield
#include <ComponentBase.hh>
#include <ComponentVoxel.hh>
#include <Sensor.hh>
#else
class Sensor;
class ComponentVoxel;
class ComponentBase;
#endif

class TRestAxionMagneticField : public TRestMetadata {
   private:
    void Initialize();

    void InitFromConfigFile();

#ifdef USE_Garfield
    Garfield::Sensor* fSetOfField;  //!
#endif

    Int_t fNofVolumes;  //<

    std::vector<TVector3> fPositions;  //<
    std::vector<TString> fFileNames;   //<

    std::vector<Double_t> fXmax;      //<
    std::vector<Double_t> fXmin;      //<
    std::vector<Double_t> fYmax;      //<
    std::vector<Double_t> fYmin;      //<
    std::vector<Double_t> fZmin;      //<
    std::vector<Double_t> fZmax;      //<
    std::vector<Double_t> fSizeMesh;  //<

    TH2D* fHisto;      //!
    TCanvas* fCanvas;  //!

   public:
    void LoadMagneticVolumes();
    TCanvas* DrawHistogram(TString projection, TString Bcomp, Int_t VIndex = -1, Double_t step = -1);
    void PrintMetadata();

    std::vector<Double_t> GetXmin() { return fXmin; }
    std::vector<Double_t> GetYmin() { return fYmin; }
    std::vector<Double_t> GetZmin() { return fZmin; }
    std::vector<Double_t> GetXmax() { return fXmax; }
    std::vector<Double_t> GetYmax() { return fYmax; }
    std::vector<Double_t> GetZmax() { return fZmax; }

    TVector3 GetMagneticField(Double_t x, Double_t y, Double_t z);

    // Constructors
    TRestAxionMagneticField();
    TRestAxionMagneticField(const char* cfgFileName, std::string name = "");
    // Destructor
    ~TRestAxionMagneticField();

    ClassDef(TRestAxionMagneticField, 1);
};
#endif
