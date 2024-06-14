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

#ifndef _TRestAxionTrueWolterOptics
#define _TRestAxionTrueWolterOptics
#include <TRestAxionOptics.h>
#include <TRestRingsMask.h>
#include <TRestSpiderMask.h>
#include <TRestTools.h>

#include <iostream>

/// A class that calculates the reflection path of X-rays through a Wolter 1 telescope
class TRestAxionTrueWolterOptics : public TRestAxionOptics {
   private:
    /// Entrance radius R1 in mm. See schematic figure.
    std::vector<Double_t> fR1;  //!

    /// Radius R2 in mm. See schematic figure.
    std::vector<Double_t> fR2;  //!

    /// Radius R3 in mm. See schematic figure.
    std::vector<Double_t> fR3;  //!

    /// Radius R4 in mm. See schematic figure.
    std::vector<Double_t> fR4;  //!

    /// Radius R5 in mm. See schematic figure.
    std::vector<Double_t> fR5;  //!

    /// Mirror angle (alpha) in radians. See schematic figure.
    std::vector<Double_t> fAlpha;  //!

    /// Mirror thickness in mm. See schematic figure.
    std::vector<Double_t> fThickness;  //!

    /// Distance between mirror stacks in mm. See schematic figure.
    /// This is here calculated using the functions from
    /// https://backend.orbit.dtu.dk/ws/portalfiles/portal/122353510/phdthesis_for_DTU_orbit.pdf.
    std::vector<Double_t> fXSep;  //!

    /// The spider structure to be used as an optical opaque mask (common to all planes)
    TRestSpiderMask* fSpiderMask = nullptr;  //<

    /// The rings structure to be used at entrance as an optical opaque mask
    TRestRingsMask* fEntranceRingsMask = nullptr;  //<

    /// The rings structure to be used at entrance as an optical opaque mask
    TRestRingsMask* fMiddleRingsMask = nullptr;  //<

    /// The rings structure to be used at entrance as an optical opaque mask
    TRestRingsMask* fExitRingsMask = nullptr;  //<

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

    /// It returns a vector with the values of R1
    std::vector<Double_t> GetR1() {
        std::vector<Double_t> r = TRestTools::GetColumnFromTable(fOpticsData, 0);
        return r;
    }

    /// It returns a vector with the values of R2
    std::vector<Double_t> GetR2() {
        std::vector<Double_t> r = TRestTools::GetColumnFromTable(fOpticsData, 1);
        return r;
    }

    /// It returns a vector with the values of R3
    std::vector<Double_t> GetR3() {
        std::vector<Double_t> r = TRestTools::GetColumnFromTable(fOpticsData, 2);
        return r;
    }

    /// It returns a vector with the values of R4
    std::vector<Double_t> GetR4() {
        std::vector<Double_t> r = TRestTools::GetColumnFromTable(fOpticsData, 3);
        return r;
    }

    /// It returns a vector with the values of R5
    std::vector<Double_t> GetR5() {
        std::vector<Double_t> r = TRestTools::GetColumnFromTable(fOpticsData, 4);
        return r;
    }

    /// It returns a vector with the values of alpha
    std::vector<Double_t> GetAlpha() {
        std::vector<Double_t> alpha = TRestTools::GetColumnFromTable(fOpticsData, 5);
        for (auto& x : alpha) x = x * units("rad") / units("deg");
        return alpha;
    }

    /// It returns a vector with the values of mirror thickness
    std::vector<Double_t> GetThickness() {
        std::vector<Double_t> t = TRestTools::GetColumnFromTable(fOpticsData, 7);
        return t;
    }

    /// It returns the entrance Z-position defined by the optical axis.
    Double_t GetEntrancePositionZ() override {
        if (fCosAlpha.size() > 0) return -fMirrorLength * fCosAlpha[0] - 0.5 * fXSep[0];
        return 0;
    }

    /// It returns the exit Z-position defined by the optical axis
    Double_t GetExitPositionZ() override {
        if (fCosAlpha_3.size() > 0) return fMirrorLength * fCosAlpha_3[0] + 0.5 * fXSep[0];
        return 0;
    }

    /// It returns the value of max entrance radius in mm
	Double_t GetMaxEntranceRadius() { return GetR1().back(); }
	
    /// It returns the value of min entrance radius in mm
	Double_t GetMinEntranceRadius() { return GetR1().front(); }

    /// It returns the value min/max entrance radius in mm as a std::pair
    std::pair<Double_t, Double_t> GetRadialLimits() override {
        std::pair<Double_t, Double_t> result(0, 0);
        if (!fR1.empty()) {
            result = {fR1.front(), fR1.back()};
        }
        return result;
    }

    TRestSpiderMask* const& GetSpiderMask () const { return fSpiderMask; }

    void SetMirror() override {
        Double_t x = fEntrancePosition.X();
        Double_t y = fEntrancePosition.Y();
        fCurrentMirror = fEntranceRingsMask->GetRegion(x, y);
    }

    TPad* DrawMirrors() override;

    void PrintParameters();
    void PrintSpider();
    void PrintMetadata() override;
    void InitFromConfigFile() override;

    TRestAxionTrueWolterOptics();
    TRestAxionTrueWolterOptics(const char* cfgFileName, std::string name = "");
    ~TRestAxionTrueWolterOptics();

    ClassDefOverride(TRestAxionTrueWolterOptics, 1);
};
#endif
