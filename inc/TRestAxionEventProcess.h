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
    /// The position respect wich the rotation will be applied
    TVector3 fCenter = TVector3(0, 0, 0);

    /// The rotation angle respect to the Y-axis
    Double_t fTheta = 0;

    /// The rotation angle with respect to Z-axis (propagation axis)
    Double_t fPhi = 0;

    /// The displacement applied to the process component
    TVector2 fDisplacement = TVector2(0, 0);

   protected:
    /// A pointer to the specific TRestAxionEvent
    TRestAxionEvent* fAxionEvent;  //!

    void BeginPrintProcess();
    void EndPrintProcess();

   public:
    /// To be executed at the beginning of the run (outside event loop)
    virtual void InitProcess() {}

    /// Begin of event process, preparation work. Called right before ProcessEvent()
    void BeginOfEventProcess(TRestEvent* evInput = nullptr);

    /// End of event process. Called directly after ProcessEvent()
    void EndOfEventProcess(TRestEvent* evInput = nullptr);

    TRestAxionEventProcess();
    ~TRestAxionEventProcess();

    ClassDef(TRestAxionEventProcess, 1);
};
#endif
