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

#ifndef _TRestAxionXrayWindow
#define _TRestAxionXrayWindow

#include <TRestMetadata.h>
#include <TRestPatternMask.h>

//! A metadata class to create x-ray transmission window definitions
class TRestAxionXrayWindow : public TRestMetadata {
   private:
    /// Position of the center of the window in mm
    TVector3 fCenter = TVector3(0, 0, 0);  //<

    /// Thicknesss of window material in mm
    Double_t fThickness = 0.0;  //<

    /// Window material name
    std::string fMaterial = "Vacuum";  //<

    /// A mask defining a pattern where the transmission will be effective
    TRestPatternMask* fMask = nullptr;  //<

    /// A vector with the energies loaded from the material file. Not stored in disk.
    std::vector<Double_t> fEnergy;  //!

    /// A vector with the transmission already renormalized using the material thickness. Not stored in disk.
    std::vector<Double_t> fTransmission;  //!

    void Initialize();

    void ReadMaterial();

    Bool_t HitsPattern(Double_t x, Double_t y);

    Int_t GetEnergyIndex(Double_t energy);

   public:
    Double_t GetWindowRadius() {
        if (!fMask) return 0;
        return fMask->GetMaskRadius();
    }

    TRestPatternMask* GetMask() const { return fMask; }

    Double_t GetTransmission(Double_t energy, Double_t x, Double_t y);

    void PrintTransmissionData();

    void InitFromConfigFile();
    void PrintMetadata();

    TRestAxionXrayWindow();
    TRestAxionXrayWindow(const char* cfgFileName, std::string name = "");

    ~TRestAxionXrayWindow();

    ClassDef(TRestAxionXrayWindow, 1);
};
#endif
