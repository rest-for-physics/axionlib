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

#ifndef TRestSoft_TRestAxionEvent
#define TRestSoft_TRestAxionEvent

#include <iostream>

#include "TMath.h"
#include "TObject.h"
#include "TVector3.h"

#include "TRestEvent.h"
#include "TRestSystemOfUnits.h"

/// An event data class to define the parameters related to an axion particle
class TRestAxionEvent : public TRestEvent {
   private:
    /// Particle position
    TVector3 fPosition;  //<

    /// Particle momentum. Unitary direction.
    TVector3 fDirection;  //<

    /// Initial energy of the axion.
    Double_t fEnergy = 0;  //<

    /// Axion mass in eV
    Double_t fMass = 0.;  //<

    /// The effective magnetic field fixed by TRestAxionFieldPropagationProcess
    Double_t fBSquared = 0;  //<

    /// The effective conversion length fixed by TRestAxionFieldPropagationProcess
    Double_t fLcoherence = 0;  //<

    /// It keeps track of efficiency introduced at different helioscope components
    std::map<std::string, Double_t> fEfficiencies;

    /// We may use it to integrate a detector response inside each event
    std::vector<Double_t> fResponse;  //<

   protected:
   public:
    TVector3 GetPosition() { return fPosition; }

    Double_t GetPositionX() { return fPosition.X(); }  // returns value in mm
    Double_t GetPositionY() { return fPosition.Y(); }  // returns value in mm
    Double_t GetPositionZ() { return fPosition.Z(); }  // returns value in mm

    TVector3 GetDirection() { return fDirection; }

    Double_t GetDirectionX() { return fDirection.X(); }  // returns normalized vector x-component.
    Double_t GetDirectionY() { return fDirection.Y(); }  // returns normalized vector y-component
    Double_t GetDirectionZ() { return fDirection.Z(); }  // returns normalized vector z-component

    Double_t GetEnergy() { return fEnergy; }            // returns value in keV
    Double_t GetMass() { return fMass * units("eV"); }  // returns value in eV

    Double_t GetBSquared() { return fBSquared; }
    Double_t GetLConversion() { return fLcoherence; }

    Double_t GetEfficiency(std::string name) { return fEfficiencies[name]; }

    void SetPosition(TVector3 pos) { fPosition = pos; }
    void SetPosition(Double_t x, Double_t y, Double_t z) { SetPosition(TVector3(x, y, z)); }

    void Translate(const TVector3& delta);

    void SetDirection(TVector3 dir) { fDirection = dir; }
    void SetDirection(Double_t px, Double_t py, Double_t pz) { SetDirection(TVector3(px, py, pz)); }

    void RotateZX(const TVector3& center, Double_t phi, Double_t theta);
    void RotateXZ(const TVector3& center, Double_t theta, Double_t phi);

    void RotateXY(const TVector3& center, Double_t pitch, Double_t yaw);
    void RotateYX(const TVector3& center, Double_t yaw, Double_t pitch);

    void SetEnergy(Double_t en) { fEnergy = en; }
    void SetMass(Double_t m) { fMass = m; }

    void SetBSquared(Double_t b) { fBSquared = b; }
    void SetLConversion(Double_t conv) { fLcoherence = conv; }

    void AddEfficiency(std::string name, Double_t value) { fEfficiencies[name] = value; }

    virtual void Initialize();

    virtual void PrintEvent();

    TPad* DrawEvent(TString option = "");

    TRestAxionEvent();
    ~TRestAxionEvent();

    ClassDef(TRestAxionEvent, 1);
};
#endif
