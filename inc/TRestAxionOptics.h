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

#include <TRestAxionOpticsMirror.h>
#include <TRestCombinedMask.h>
#include <TRestMetadata.h>
#include <iostream>

#include "TRandom3.h"

/// An abstract class to define common optics parameters and methods
class TRestAxionOptics : public TRestMetadata {
   protected:
    /// An optics file that contains all the specific optics parameters
    std::string fOpticsFile = "";

    /// The mirror length. If all mirrors got the same length. Otherwise will be zero.
    Double_t fMirrorLength = 0;  //<

    /// The optics data table extracted from fOpticsFile
    std::vector<std::vector<Double_t>> fOpticsData;  //<

    /// The mirror properties
    TRestAxionOpticsMirror* fMirror = nullptr;  //<

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

    /// The calculated focal in mm. It is updated by the method TRestAxionOptics::FindFocal.
    Double_t fFocal = -1;  //!

    /// The entrance optical mask that defines a pattern with regions identified with a number
    TRestCombinedMask* fEntranceMask = nullptr;  //!

    /// The middle optical mask that defines a pattern with regions identified with a number
    TRestCombinedMask* fMiddleMask = nullptr;  //!

    /// The exit optical mask that defines a pattern with regions identified with a number
    TRestCombinedMask* fExitMask = nullptr;  //!

    /// During the photon propagation it keeps track of the active mirror shell.
    Int_t fCurrentMirror = -1;  //!

    /// During the photon propagation it tells us if the photon interacted in the first mirror
    Bool_t fFirstInteraction = false;  //!

    /// During the photon propagation it tells us if the photon interacted in the second mirror
    Bool_t fSecondInteraction = false;  //!

    /// Random number generator
    TRandom3* fRandom = nullptr;  //!

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

    /// It returns the lower/higher radius range where photons are allowed
    virtual std::pair<Double_t, Double_t> GetRadialLimits() = 0;

    /// It returns the entrance Z-position defined by the optical axis.
    virtual Double_t GetEntranceZPosition() = 0;

    /// It returns the exit Z-position defined by the optical axis
    virtual Double_t GetExitZPosition() = 0;

    /// It updates the values fFirstInteractionPosition and fMiddleDirection. Returns 0 if is not in region.
    virtual Int_t FirstMirrorReflection(const TVector3& pos, const TVector3& dir) = 0;

    /// It updates the values fSecondInteractionPosition and fExitDirection. Returns 0 if is not in region.
    virtual Int_t SecondMirrorReflection(const TVector3& pos, const TVector3& dir) = 0;

    /// It draws the mirrors using a TGraph. To be implemented at the inherited class.
    virtual TPad* DrawMirrors() = 0;

    virtual Double_t FindFocal(Double_t from, Double_t to, Double_t energy, Double_t precision = 1,
                               Bool_t recalculate = false, Int_t particles = 5000);

    Double_t CalculateSpotSize(Double_t energy, Double_t z, Int_t particles = 15000);

    TPad* CreatePad(Int_t nx = 1, Int_t ny = 1);

    TPad* DrawParticleTracks(Double_t deviation = 0, Int_t particles = 10);

    TPad* DrawScatterMaps(Double_t z, Double_t energy = 0, Double_t deviation = 0, Int_t particles = 1000,
                          Double_t focalHint = 7500);

    TPad* DrawDensityMaps(Double_t z, Double_t energy = 0, Double_t deviation = 0, Int_t particles = 1000,
                          Double_t focalHint = 7500);

    Double_t PropagatePhoton(const TVector3& pos, const TVector3& dir, Double_t energy);

    Int_t PropagateMonteCarloPhoton(Double_t energy, Double_t deviation);

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

    /// It returns true if the photon got reflected in the first mirror
    Bool_t IsFirstMirrorReflection() { return fFirstInteraction; }

    /// It returns true if the photon got reflected in the second mirror
    Bool_t IsSecondMirrorReflection() { return fSecondInteraction; }

    Int_t GetNumberOfReflections();

    /// Returns a pointer to access directly the entrance mask information
    TRestCombinedMask* const& GetEntranceMask() const { return fEntranceMask; }

    /// Returns a pointer to access directly the middle mask information
    TRestCombinedMask* const& GetMiddleMask() const { return fMiddleMask; }

    /// Returns a pointer to access directly the exit mask information
    TRestCombinedMask* const& GetExitMask() const { return fExitMask; }

    void PrintMetadata();

    void PrintMirror();
    void PrintMasks();

    void PrintEntranceMask();
    void PrintMiddleMask();
    void PrintExitMask();

    void PrintPhotonTrackingSummary();

    void InitFromConfigFile();

    ~TRestAxionOptics();

    ClassDef(TRestAxionOptics, 1);
};
#endif
