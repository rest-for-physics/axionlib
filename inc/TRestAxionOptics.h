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
    /// It is the calculated axis position at the entrance of the optics plane.
    TVector3 fEntrance = TVector3(0, 0, 0);  //!

    /// It is the calculated axis position at the exit of the optics plane.
    TVector3 fExit = TVector3(0, 0, 0);  //!

    /// A vector used to define a reference vector at the optics plane
    TVector3 fReference = TVector3(0, 0, 0);  //!

   protected:
    TRestAxionOptics();
    TRestAxionOptics(const char* cfgFileName, std::string name = "");

   public:
    virtual void Initialize();

    /// It returns the physical length of one mirror stack; the whole optical system would be L=(fLength + 1/2
    /// * xSep) * (cos(angleRing) + cos(angleRing)) which doesn't work here because the angele hasn't been
    /// defined
    //   Double_t GetMirrLength() { return fLength; }

    /// It returns the physical length of the whole optics approximated
    // Double_t GetLength() { return fLength * 2; }

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

    Int_t GetEntranceRing(const TVector3& pos, const TVector3& dir);

    void PrintMetadata();

    void InitFromConfigFile();

    ~TRestAxionOptics();

    ClassDef(TRestAxionOptics, 1);
};
#endif
