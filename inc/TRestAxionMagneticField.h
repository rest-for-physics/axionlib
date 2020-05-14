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

#ifndef _TRestAxionMagneticField
#define _TRestAxionMagneticField

#include <TRestMetadata.h>
#include <iostream>

#include "TCanvas.h"
#include "TH2D.h"
#include "TVector3.h"
#include "TVectorD.h"

#include "TRestAxionBufferGas.h"
#include "TRestMesh.h"

/// A structure to define the properties and store the field data of a single magnetic volume inside
/// TRestAxionMagneticField
struct MagneticFieldVolume {
    /// Offset position applied to the volume
    TVector3 position;

    /// A description of the grid, mesh size, and number of nodes
    TRestMesh mesh;

    /// The field data connected to the grid defined by the mesh
    std::vector<std::vector<std::vector<TVector3>>> field;

    /// A pointer to the gas properties
    TRestAxionBufferGas* bGas = NULL;
};

/// A class to load magnetic field maps and evaluate the field on those maps including interpolation.
class TRestAxionMagneticField : public TRestMetadata {
   private:
    /// The name of the filenames containing the field data
    std::vector<string> fFileNames;  //<

    /// The absolute position of each of the magnetic volumes defined in this class
    std::vector<TVector3> fPositions;  //<

    /// A constant field component that will be added to the field map
    std::vector<TVector3> fConstantField;  //<

    /// The size of a grid element from the mesh in mm
    std::vector<TVector3> fMeshSize;  //<

    /// The type of the mesh used (default is cylindrical)
    std::vector<TString> fMeshType;  //<

    /// The gas mixture components that define the medium in each magnetic volume
    std::vector<TString> fGasMixtures;  //<

    /// The gas components densities corresponding to the gas mixture defined for each volume
    std::vector<TString> fGasDensities;  //<

    /// A vector to store the maximum bounding box values
    std::vector<TVector3> fBoundMax;  //<

    /// A magnetic field volume structure to store field data and mesh.
    std::vector<MagneticFieldVolume> fMagneticFieldVolumes;  //!

    /// A helper histogram to plot the field
    TH2D* fHisto;  //!

    /// A canvas to insert the histogram drawing
    TCanvas* fCanvas;  //!

    void Initialize();

    void InitFromConfigFile();

    void LoadMagneticFieldData(MagneticFieldVolume& mVol, std::vector<std::vector<Float_t>> data);

    TVector3 GetMagneticVolumeNode(MagneticFieldVolume mVol, TVector3 pos);

    /// \brief This private method returns true if the magnetic field volumes loaded are the same as
    /// the volumes defined.
    Bool_t FieldLoaded() { return GetNumberOfVolumes() == fMagneticFieldVolumes.size(); }

    /// It returns a  pointer to the corresponding magnetic volume id
    MagneticFieldVolume* GetMagneticVolume(Int_t id) {
        if (!FieldLoaded()) LoadMagneticVolumes();
        if (fMagneticFieldVolumes.size() > id)
            return &fMagneticFieldVolumes[id];
        else {
            ferr << "TRestAxionMagneticField::GetMagneticVolume. Id outside limits!" << endl;
            return NULL;
        }
    }

   public:
    void LoadMagneticVolumes();

    /// It returns true if no magnetic field map was loaded for that volume
    Bool_t IsFieldConstant(Int_t id) {
        if (GetMagneticVolume(id)) return GetMagneticVolume(id)->field.size() == 0;
        return true;
    }

    /// The number of magnetic volumes loaded into the object
    Int_t GetNumberOfVolumes() { return fPositions.size(); }

    Bool_t CheckOverlaps();

    std::vector<TVector3> GetVolumeBoundaries(Int_t id, TVector3 pos, TVector3 dir);
    std::vector<TVector3> GetFieldBoundaries(Int_t id, TVector3 pos, TVector3 dir, Double_t precision = 0);

    TVector3 GetMagneticField(Double_t x, Double_t y, Double_t z);
    TVector3 GetMagneticField(TVector3 pos);

    Double_t GetPhotonMass(Double_t x, Double_t y, Double_t z, Double_t en);
    Double_t GetPhotonMass(TVector3 pos, Double_t en);
    Double_t GetPhotonMass(Int_t id, Double_t en);

    Double_t GetPhotonAbsorptionLength(Double_t x, Double_t y, Double_t z, Double_t en);
    Double_t GetPhotonAbsorptionLength(TVector3 pos, Double_t en);
    Double_t GetPhotonAbsorptionLength(Int_t id, Double_t en);

    Int_t GetVolumeIndex(TVector3 pos);

    TVector3 GetVolumePosition(Int_t id);
    TVector3 GetVolumeCenter(Int_t id);

    Double_t GetTransversalComponent(TVector3 position, TVector3 direction);

    std::vector<Double_t> GetTransversalComponentAlongPath(TVector3 from, TVector3 to, Double_t dl = 1.,
                                                           Int_t Nmax = 0);

    Double_t GetTransversalFieldAverage(TVector3 from, TVector3 to, Double_t dl = 1., Int_t Nmax = 0);

    TCanvas* DrawHistogram(TString projection, TString Bcomp, Int_t volIndex = -1, Double_t step = -1,
                           TString style = "COLZ0", Double_t depth = -100010.0);

    void PrintMetadata();

    TRestAxionMagneticField();
    TRestAxionMagneticField(const char* cfgFileName, std::string name = "");
    ~TRestAxionMagneticField();

    ClassDef(TRestAxionMagneticField, 3);
};
#endif
