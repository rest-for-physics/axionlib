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

#include "TRestMesh.h"

/// A class to load magnetic field maps and provide an evaluate the field including interpolation.
class TRestAxionMagneticField : public TRestMetadata {
   private:
    struct MagneticFieldVolume {
        // Offset position applied to the volume
        TVector3 position;

        // A description of the grid, mesh size, and number of nodes
        TRestMesh mesh;

        // The field data connected to the grid defined by the mesh
        std::vector<std::vector<std::vector<TVector3>>> field;
    };

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

    Double_t GetTransversalComponent(TVector3 position, TVector3 direction);

    std::vector<Double_t> GetTransversalComponentAlongPath(TVector3 from, TVector3 to, Double_t dl = 1.,
                                                           Int_t Nmax = 0);

    Double_t GetTransversalFieldAverage(TVector3 from, TVector3 to, Double_t dl = 1., Int_t Nmax = 0);

    TCanvas* DrawHistogram(TString projection, TString Bcomp, Int_t VIndex = -1, Double_t step = -1);

    void PrintMetadata();

    TRestAxionMagneticField();
    TRestAxionMagneticField(const char* cfgFileName, std::string name = "");
    ~TRestAxionMagneticField();

    ClassDef(TRestAxionMagneticField, 1);
};
#endif
