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

#ifndef _TRestAxionGenericOptics
#define _TRestAxionGenericOptics
#include <TRestAxionOptics.h>
#include <iostream>
/// A class calculate the reflection path and probability of X-rays through a Wolter 1 telescope
class TRestAxionGenericOptics : public TRestAxionOptics {
   private:
    /// A vector containing the shells seperations between the two stacks. First element is the lowest radius.
    std::vector<Double_t> fShellsSep;  //<

    /// A vector containing the shells angles. First element is the lowest radius. Note that the second stack
    /// has the tripple of this angle.
    std::vector<Double_t> fShellsAngle;  //<

    /// coating materials and surface roughness or better yet: direction to a file that gives the reflectivity
    std::string fReflectivityFileName;  //<

    void Initialize();

   public:
    void PrintMetadata();
    void InitFromConfigFile();

    /// It returns the position at the optics exit plane for the incoming particle
    TVector3 GetPositionAtExit(const TVector3& pos, const TVector3& dir) { return TVector3(0, 0, 0); }

    /// It returns the direction at the optics exit plane for the incoming particle
    TVector3 GetDirectionAtExit(const TVector3& pos, const TVector3& dir) { return TVector3(0, 0, 0); }

    /// It returns the efficiency for particle with position `pos` and direction `dir`.
    Double_t GetEfficiency(const TVector3& pos, const TVector3& dir) { return 0.0; }

    TRestAxionGenericOptics();
    TRestAxionGenericOptics(const char* cfgFileName, std::string name = "");
    ~TRestAxionGenericOptics();

    /// get the interaction point of the photon with the mirror
    // TVector3 GetInteractionPoint(const TVector3& pos, const TVector3& dir, ): TRestAxionOptics(){}
    ClassDef(TRestAxionGenericOptics, 1);
};
#endif
