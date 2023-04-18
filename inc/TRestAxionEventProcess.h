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

#ifndef RestCore_TRestAxionEventProcess
#define RestCore_TRestAxionEventProcess

#include "TRestAxionEvent.h"
#include "TRestEventProcess.h"

/// A base class for any axion event process. Defines position, rotation and component displacement.
class TRestAxionEventProcess : public TRestEventProcess {
   private:
    /// The position around which the rotation will be applied
    TVector3 fCenter = TVector3(0, 0, 0);

    /// The rotation angle with respect to the Y-axis
    Double_t fYaw = 0;

    /// The rotation angle with respect to X-axis
    Double_t fPitch = 0;
   
   /// The position (different than the center of the object) around which the rotation will be applied
    TVector3 fCenterSetup = TVector3(0, 0, 0);

    /// The rotation angle around CenterSetup with respect to the Y-axis
    Double_t fYawSetup = 0;

    /// The rotation angle around CenterSetup with respect to X-axis
    Double_t fPitchSetup = 0;

    /// If enabled it will skip the end rotation that recovers the original axion trajectory direction
    Bool_t fSkipEndProcessRotation = false;  //!

   protected:
    /// A pointer to the specific TRestAxionEvent
    TRestAxionEvent* fAxionEvent;  //!

    void BeginPrintProcess();
    void EndPrintProcess();

    TVector3 GetCenter() const { return fCenter; }

    void SkipEndProcessRotation(Bool_t value = true) { fSkipEndProcessRotation = value; }

   public:
    RESTValue GetInputEvent() const override { return fAxionEvent; }
    RESTValue GetOutputEvent() const override { return fAxionEvent; }

    virtual void InitProcess() override {}

    /// Begin of event process, preparation work. Called right before ProcessEvent()
    virtual void BeginOfEventProcess(TRestEvent* evInput = nullptr) override;

    /// End of event process. Called directly after ProcessEvent()
    virtual void EndOfEventProcess(TRestEvent* evInput = nullptr) override;

    TRestAxionEventProcess();
    ~TRestAxionEventProcess();

    ClassDefOverride(TRestAxionEventProcess, 1);
};
#endif
