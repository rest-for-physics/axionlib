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

#include "TRestMesh.h"

#if defined USE_Garfield
#include <ComponentBase.hh>
#include <ComponentVoxel.hh>
#include <Sensor.hh>
using namespace Garfield;
#else
class Sensor;
#endif

struct MagneticFieldVolume {
    // Offset position applied to the volume
    TVector3 position;

    // A description of the grid, mesh size, and number of nodes
    TRestMesh mesh;

    // The field data connected to the grid defined by the mesh
    std::vector<std::vector<std::vector<TVector3>>> field;
};

/// A class to load magnetic field maps and provide an easy access to it.
class TRestAxionMagneticField : public TRestMetadata {
   private:
    /// The name of the filenames containing the field data
    std::vector<string> fFileNames;  //<

    /// The absolute position of each of the magnetic volumes defined in this class
    std::vector<TVector3> fPositions;  //<

    /// A magnetic field volume structure to store field data and mesh.
    std::vector<MagneticFieldVolume> fMagneticFieldVolumes;  //!

    /// A helper histogram to plot the field
    TH2D* fHisto;  //!

    /// A canvas to insert the histogram drawing
    TCanvas* fCanvas;  //!

    void Initialize();

    void InitFromConfigFile();

    void LoadMagneticVolumes();

    void LoadMagneticFieldData(MagneticFieldVolume& mVol, std::vector<std::vector<Double_t>> data);

    TVector3 GetMagneticVolumeNode(MagneticFieldVolume mVol, TVector3 pos);

   public:
    /// The number of magnetic volumes loaded into the object
    Int_t GetNumberOfVolumes() { return fMagneticFieldVolumes.size(); }

    TVector3 GetMagneticField(Double_t x, Double_t y, Double_t z);

    TVector3 GetMagneticField(TVector3 pos);

    TCanvas* DrawHistogram(TString projection, TString Bcomp, Int_t VIndex = -1, Double_t step = -1);

    void PrintMetadata();

    TRestAxionMagneticField();
    TRestAxionMagneticField(const char* cfgFileName, std::string name = "");
    ~TRestAxionMagneticField();

    ClassDef(TRestAxionMagneticField, 1);
};
#endif
