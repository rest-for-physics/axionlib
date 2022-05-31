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

    /// Mirror length in mm. See schematic figure.
    std::vector<Double_t> fLength;  //<

    /// Mirror thickness in mm. See schematic figure.
    std::vector<Double_t> fThickness;  //<

   public:
    void Initialize() override;

    void PrintMetadata() override;
    void InitFromConfigFile() override;

    /// It returns the position at the optics exit plane for the incoming particle
    TVector3 GetPositionAtExit(const TVector3& pos, const TVector3& dir) override {
        return TVector3(0, 0, 0);
    }

    /// It returns the direction at the optics exit plane for the incoming particle
    TVector3 GetDirectionAtExit(const TVector3& pos, const TVector3& dir) override {
        return TVector3(0, 0, 0);
    }

    /// It returns the efficiency for particle with position `pos` and direction `dir`.
    Double_t GetEfficiency(const TVector3& pos, const TVector3& dir) override { return 0.0; }

    TRestAxionWolterOptics();
    TRestAxionWolterOptics(const char* cfgFileName, std::string name = "");
    ~TRestAxionWolterOptics();

    /// get the interaction point of the photon with the mirror
    TVector3 GetInteractionPoint(const TVector3& pos, const TVector3& dir);
    ClassDefOverride(TRestAxionWolterOptics, 1);
};
#endif
