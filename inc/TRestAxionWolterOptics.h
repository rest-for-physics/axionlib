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

#ifndef _TRestAxionWolterOptics
#define _TRestAxionWolterOptics
#include <TRestAxionOptics.h>
#include <TRestRingsMask.h>
#include <TRestSpiderMask.h>
#include <iostream>

/// A class calculate the reflection path and probability of X-rays through a Wolter 1 telescope
class TRestAxionWolterOptics : public TRestAxionOptics {
   private:
    /// An optics file that contains all the specific Wolter optics parameters
    std::string fOpticsFile = "";

    /// Entrance radius R1 in mm. See schematic figure.
    std::vector<Double_t> fR1;  //<

    /// Radius R2 in mm. See schematic figure.
    std::vector<Double_t> fR2;  //<

    /// Radius R3 in mm. See schematic figure.
    std::vector<Double_t> fR3;  //<

    /// Radius R4 in mm. See schematic figure.
    std::vector<Double_t> fR4;  //<

    /// Radius R5 in mm. See schematic figure.
    std::vector<Double_t> fR5;  //<

    /// Mirror angle (alpha) in radians. See schematic figure.
    std::vector<Double_t> fAlpha;  //<

    /// Mirror thickness in mm. See schematic figure.
    std::vector<Double_t> fThickness;  //<

    /// The spider structure to be used as an optical opaque mask (common to all planes)
    TRestSpiderMask* fSpiderMask = nullptr;  //<

    /// The rings structure to be used at entrance as an optical opaque mask
    TRestRingsMask* fEntranceRingsMask = nullptr;

    /// The rings structure to be used at entrance as an optical opaque mask
    TRestRingsMask* fMiddleRingsMask = nullptr;

    /// The rings structure to be used at entrance as an optical opaque mask
    TRestRingsMask* fExitRingsMask = nullptr;

    /// Mirror pre-calculated cosine angle (alpha). See schematic figure.
    std::vector<Double_t> fCosAlpha;  //!

    /// Mirror pre-calculated cosine angle (alpha). See schematic figure.
    std::vector<Double_t> fCosAlpha_3;  //!

    /// The Z-position of the cone vertex defined by the front mirrors
    std::vector<Double_t> fFrontVertex;  //!

    /// The Z-position of the cone vertex defined by the back mirrors
    std::vector<Double_t> fBackVertex;  //!

    Int_t FirstMirrorReflection(const TVector3& pos, const TVector3& dir) override;
    Int_t SecondMirrorReflection(const TVector3& pos, const TVector3& dir) override;

   public:
    void Initialize() override;

    /// It returns the entrance Z-position defined by the optical axis.
    Double_t GetEntranceZPosition() override {
        if (fCosAlpha.size() > 0) return -fMirrorLength * fCosAlpha[0];
        return 0;
    }

    /// It returns the exit Z-position defined by the optical axis
    Double_t GetExitZPosition() override {
        if (fCosAlpha_3.size() > 0) return fMirrorLength * fCosAlpha_3[0];
        return 0;
    }

    void SetMirror() override {
        fCurrentMirror = fEntranceRingsMask->GetRegion(fEntrancePosition.X(), fEntrancePosition.Y());
    }

    TPad* DrawMirrors() override { return fPad; }

    void PrintParameters();
    void PrintSpider();
    void PrintMetadata() override;
    void InitFromConfigFile() override;

    TRestAxionWolterOptics();
    TRestAxionWolterOptics(const char* cfgFileName, std::string name = "");
    ~TRestAxionWolterOptics();

    ClassDefOverride(TRestAxionWolterOptics, 1);
};
#endif
