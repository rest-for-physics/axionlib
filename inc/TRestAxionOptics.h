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

#include <TRestCombinedMask.h>
#include <TRestMetadata.h>
#include <iostream>

/// An abstract class to define common optics parameters and methods
class TRestAxionOptics : public TRestMetadata {
   protected:
    /// The mirror length. If all mirrors got the same length. Otherwise will be zero.
    Double_t fMirrorLength = 0;  //<

    /// The particle position at the origin
    TVector3 fOriginPosition;  //!

    /// The particle position at the optics plane entrance
    TVector3 fEntrancePosition;  //!

    /// The particle position at the front mirror interaction point
    TVector3 fFirstInteractionPosition;  //!

    /// The particle position at the optics plane middle
    TVector3 fMiddlePosition;  //!

    /// The particle position at the back mirror interaction point
    TVector3 fSecondInteractionPosition;  //!

    /// The particle position at the optics plane exit
    TVector3 fExitPosition;  //!

    /// The particle position at the optics plane entrance
    TVector3 fEntranceDirection;  //!

    /// The particle position at the optics plane middle
    TVector3 fMiddleDirection;  //!

    /// The particle position at the optics plane exit
    TVector3 fExitDirection;  //!

    /// The entrance optical mask that defines a pattern with regions identified with a number
    TRestCombinedMask* fEntranceMask = nullptr;  //!

    /// The middle optical mask that defines a pattern with regions identified with a number
    TRestCombinedMask* fMiddleMask = nullptr;  //!

    /// The exit optical mask that defines a pattern with regions identified with a number
    TRestCombinedMask* fExitMask = nullptr;  //!

    /// During the photon propagation it keeps track of the active mirror shell.
    Int_t fCurrentMirror = -1;  //!

    TRestAxionOptics();
    TRestAxionOptics(const char* cfgFileName, std::string name = "");

   private:
    void ResetPositions();

   protected:
    /// A pad pointer to be used by the drawing methods
    TPad* fPad = nullptr;

    Int_t TransportToEntrance(const TVector3& pos, const TVector3& dir);
    Int_t TransportToMiddle(const TVector3& pos, const TVector3& dir);
    Int_t TransportToExit(const TVector3& pos, const TVector3& dir);

    /// It returns the mirror index to be used in the photon reflection.
    Int_t GetMirror() { return fCurrentMirror; }

    /// It must be implemented at the inherited optics, making use of fEntrancePosition
    virtual void SetMirror() = 0;

   public:
    virtual void Initialize();

    /// It returns the entrance Z-position defined by the optical axis.
    virtual Double_t GetEntranceZPosition() = 0;

    /// It returns the exit Z-position defined by the optical axis
    virtual Double_t GetExitZPosition() = 0;

    /// It updates the values fFirstInteractionPosition and fMiddleDirection. Returns 0 if is not in region.
    virtual Int_t FirstMirrorReflection(const TVector3& pos, const TVector3& dir) = 0;

    /// It updates the values fSecondInteractionPosition and fExitDirection. Returns 0 if is not in region.
    virtual Int_t SecondMirrorReflection(const TVector3& pos, const TVector3& dir) = 0;

    /// It draws the mirrors using a TGraph
    virtual TPad* DrawMirrors() = 0;

    /// It updates the internal TRestAxionOptics particle positions/directions and returns efficiency
    Double_t PropagatePhoton(const TVector3& pos, const TVector3& dir, Double_t energy);

    /// Returns the entrance position from the latest propagated photon
    TVector3 GetEntrancePosition() { return fEntrancePosition; }

    /// Returns the middle position from the latest propagated photon
    TVector3 GetMiddlePosition() { return fEntrancePosition; }

    /// Returns the exit position from the latest propagated photon
    TVector3 GetExitPosition() { return fEntrancePosition; }

    /// Returns the entrance position from the latest propagated photon
    TVector3 GetEntranceDirection() { return fEntranceDirection; }

    /// Returns the middle position from the latest propagated photon
    TVector3 GetMiddleDirection() { return fMiddleDirection; }

    /// Returns the exit position from the latest propagated photon
    TVector3 GetExitDirection() { return fExitDirection; }

    void PrintEntranceMask() {
        if (fEntranceMask)
            fEntranceMask->PrintMetadata();
        else
            RESTWarning << "TRestAxionOptics::PrintEntranceMask. Not available" << RESTendl;
    }

    void PrintMetadata();

    void PrintMasks();

    void InitFromConfigFile();

    ~TRestAxionOptics();

    ClassDef(TRestAxionOptics, 1);
};
#endif
