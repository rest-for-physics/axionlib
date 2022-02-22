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

#ifndef _TRestAxionOptics
#define _TRestAxionOptics

#include <TRestMetadata.h>
#include <iostream>

/// An abstract class to define common optics parameters and methods
class TRestAxionOptics : public TRestMetadata {
   private:
    /// It is the position of the center of the optics system.
    TVector3 fCenter = TVector3(0, 0, 0);  //<

    /// The axis of the optics system
    TVector3 fAxis = TVector3(0, 0, 1);  //<

    /// Optics physical mirror length in mm
    Double_t fLength = 250;  //<

    /// A vector containing the shells ring radius definitions. First element is the lower radius.
    std::vector<std::pair<Double_t, Double_t>> fShellsRadii;  //<

    /// The angle between two consecutive spider arms measured in radians.
    Double_t fSpiderArmsSeparationAngle = 0;  //<

    /// The position angle at which the spider arm will start
    Double_t fSpiderOffsetAngle = 0;  //<

    /// The width of each specific spider arm. Measured in radians. Default is 2.5 degrees.
    Double_t fSpiderWidth = TMath::Pi() / 18. / 4.;  //<

    /// The spider structure will be effective from this radius, in mm. Default is from 0 mm.
    Double_t fSpiderStartRadius = 0.;  //<

    /// An internal variable to register the maximum shell radius
    Double_t fMaxShellRadius = -1;  //!

    /// An internal variable to register the minimum shell radius
    Double_t fMinShellRadius = -1;  //!

    /// It is the calculated axis position at the entrance of the optics plane.
    TVector3 fEntrance = TVector3(0, 0, 0);  //!

    /// It is the calculated axis position at the exit of the optics plane.
    TVector3 fExit = TVector3(0, 0, 0);  //!

    /// A vector used to define a reference vector at the optics plane
    TVector3 fReference = TVector3(0, 0, 0);  //!

    /// It defines the forbidden (cosine) angular ranges imposed by the spider structure (0,Pi)
    std::vector<std::pair<Double_t, Double_t>> fSpiderPositiveRanges;  //!

    /// It defines the forbidden (cosine) angular ranges imposed by the spider structure (Pi,2Pi)
    std::vector<std::pair<Double_t, Double_t>> fSpiderNegativeRanges;  //!

    void SetMaxAndMinShellRadius();
    void InitializeSpiderAngles();

    Bool_t IsInsideRing(const TVector3& pos, Double_t Rout, Double_t Rin = 0);
    Bool_t HitsSpider(const TVector3& pos);

   public:
    void Initialize();

    /// It returns the center of the optics system
    TVector3 GetCenter() { return fCenter; }

    /// It returns the axis vector of the optics system
    TVector3 GetAxis() { return fAxis; }

    /// It returns the physical length of the optics system
    Double_t GetLength() { return fLength; }

    /// It returns the number of shells implemented in the optics system
    Int_t GetNumberOfShells() { return fShellsRadii.size(); }

    /// It returns the entrance position defined by the optical axis
    TVector3 GetEntrance() { return fEntrance; }

    /// It returns the exit position defined by the optical axis
    TVector3 GetExit() { return fExit; }

    TVector3 GetPositionAtEntrance(const TVector3& pos, const TVector3& dir);

    /// Pure abstract method to be implemented at inherited class
    virtual TVector3 GetPositionAtExit(const TVector3& pos, const TVector3& dir) { return TVector3(0, 0, 0); }

    /// Pure abstract method to be implemented at inherited class
    virtual TVector3 GetDirectionAtExit(const TVector3& pos, const TVector3& dir) {
        return TVector3(0, 0, 0);
    }

    /// Pure abstract method to be implemented at inherited class
    virtual Double_t GetEfficiency(const TVector3& pos, const TVector3& dir) { return 0.0; }

    Int_t GetEntranceShell(const TVector3& pos, const TVector3& dir);

    void PrintMetadata();

    void InitFromConfigFile();

    TRestAxionOptics();
    TRestAxionOptics(const char* cfgFileName, std::string name = "");
    ~TRestAxionOptics();

    ClassDef(TRestAxionOptics, 1);
};
#endif
