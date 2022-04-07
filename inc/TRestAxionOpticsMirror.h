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

#ifndef _TRestAxionOpticsMirror
#define _TRestAxionOpticsMirror

#include <TRestMetadata.h>
#include <iostream>

#include <TCanvas.h>

/// A metadata class accessing the Henke database to load reflectivity data
class TRestAxionOpticsMirror : public TRestMetadata {
   private:
    /// The mirror type (Thick, Single, Bilayer, Multilayer). Only `Single` is supported now.
    std::string fMirrorType = "Single";  //<

    /// The mirror layer material (chemical forumula).
    std::string fLayer = "C";  //<

    /// The layer thickness in nm
    std::string fLayerThickness = "30";  //<

    /// The substrate material
    std::string fSubstrate = "SiO2";  //<

    /// Top surface roughness in nm
    std::string fSigma1 = "0";  //<

    /// A set of key-value pairs sent to the Henke website for data request
    std::map<std::string, std::string> fHenkeKeys;  //!

    /// The reflectivity loaded as a table with angle versus energy
    std::vector<std::vector<Float_t> > fReflectivityTable;  //!

    /// The transmission loaded as a table with angle versus energy
    std::vector<std::vector<Float_t> > fTransmissionTable;  //!

    /// A canvas to insert the optic properties drawing
    TCanvas* fCanvas = nullptr;  //!

    std::string DownloadHenkeFile();

    Int_t ExportTables();

    std::string GetReflectivityFilename();
    std::string GetTransmissionFilename();

   public:
    void Initialize();

    void SetMirrorType(const std::string& type) { fMirrorType = type; }
    void SetLayerMaterial(const std::string& layer) { fLayer = layer; }
    void SetLayerThickness(const std::string& thickness) { fLayerThickness = thickness; }
    void SetSubstrateMaterial(const std::string& substrate) { fSubstrate = substrate; }
    void SetRoughness(const std::string& roughness) { fSigma1 = roughness; }

    void LoadTables();

    Double_t GetReflectivity(const Double_t angle, const Double_t energy);

    Double_t GetTransmission(const Double_t angle, const Double_t energy);

    void PrintMetadata();

    TCanvas* DrawOpticsProperties(std::string options = "", Double_t lowRange = 1.e-5,
                                  Double_t lowRange2 = 1.e-3);

    TRestAxionOpticsMirror();
    TRestAxionOpticsMirror(const char* cfgFileName, std::string name = "");
    ~TRestAxionOpticsMirror();

    ClassDef(TRestAxionOpticsMirror, 1);
};
#endif
