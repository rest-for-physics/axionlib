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

//////////////////////////////////////////////////////////////////////////
///
/// TOBE written
///
/// \class TRestAxionEventProcess
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2022-March: First concept.
///        		Javier Galan
///
/// <hr>
//////////////////////////////////////////////////////////////////////////

#include "TRestAxionEventProcess.h"

using namespace std;

ClassImp(TRestAxionEventProcess);

//////////////////////////////////////////////////////////////////////////
/// TRestAxionEventProcess default constructor
///
TRestAxionEventProcess::TRestAxionEventProcess() { fSingleThreadOnly = false; }

//////////////////////////////////////////////////////////////////////////
/// TRestAxionEventProcess destructor
///
TRestAxionEventProcess::~TRestAxionEventProcess() {}

//////////////////////////////////////////////////////////////////////////
/// \brief Begin of event process, preparation work. Called right before ProcessEvent()
///
/// This method is called before calling ProcessEvent(). We initialize the process's output
/// event if not null and not same as input event. The event's basic info (ID, timestamp, etc.)
/// will also be set to the same as input event
///
void TRestAxionEventProcess::BeginOfEventProcess(TRestEvent* inEv) {
    TRestEventProcess::BeginOfEventProcess(inEv);

    fAxionEvent = (TRestAxionEvent*)inEv;
    debug << "BoEP: Initial Position. X: " << fAxionEvent->GetPosition().X()
          << " Y: " << fAxionEvent->GetPosition().Y() << " Z: " << fAxionEvent->GetPosition().Z() << endl;
    debug << "BoEP: Initial Direction. X: " << fAxionEvent->GetDirection().X()
          << " Y: " << fAxionEvent->GetDirection().Y() << " Z: " << fAxionEvent->GetDirection().Z() << endl;
    fAxionEvent->RotateZX(fCenter, -fPhi, -fTheta);
    fAxionEvent->Translate(TVector3(-fDisplacement.X(), -fDisplacement.Y(), 0));
    debug << " ---- " << endl;
    debug << "BoEP: Final Position. X: " << fAxionEvent->GetPosition().X()
          << " Y: " << fAxionEvent->GetPosition().Y() << " Z: " << fAxionEvent->GetPosition().Z() << endl;
    debug << "BoEP: Final Direction. X: " << fAxionEvent->GetDirection().X()
          << " Y: " << fAxionEvent->GetDirection().Y() << " Z: " << fAxionEvent->GetDirection().Z() << endl;
    debug << " ++++ " << endl;
}

//////////////////////////////////////////////////////////////////////////
/// \brief End of event process. Validate the updated observable number matches total defined observable
/// number
///
void TRestAxionEventProcess::EndOfEventProcess(TRestEvent* evInput) {
    TRestEventProcess::EndOfEventProcess(evInput);

    debug << "EoEP: Initial Position. X: " << fAxionEvent->GetPosition().X()
          << " Y: " << fAxionEvent->GetPosition().Y() << " Z: " << fAxionEvent->GetPosition().Z() << endl;
    debug << "EoEP: Initial Direction. X: " << fAxionEvent->GetDirection().X()
          << " Y: " << fAxionEvent->GetDirection().Y() << " Z: " << fAxionEvent->GetDirection().Z() << endl;
    debug << " ---- " << endl;
    fAxionEvent->Translate(TVector3(fDisplacement.X(), fDisplacement.Y(), 0));
    fAxionEvent->RotateXZ(fCenter, fTheta, fPhi);
    debug << "EoEP: Final Position. X: " << fAxionEvent->GetPosition().X()
          << " Y: " << fAxionEvent->GetPosition().Y() << " Z: " << fAxionEvent->GetPosition().Z() << endl;
    debug << "EoEP: Final Direction. X: " << fAxionEvent->GetDirection().X()
          << " Y: " << fAxionEvent->GetDirection().Y() << " Z: " << fAxionEvent->GetDirection().Z() << endl;
    debug << " ++++ " << endl;
}

//////////////////////////////////////////////////////////////////////////
/// \brief Pre-defined printer, can be used at the beginning in the
/// implementation of PrintMetadata()
///
/// Prints process type, name, title, verboselevel, outputlevel, input/output
/// event type, and several separators
void TRestAxionEventProcess::BeginPrintProcess() {
    TRestEventProcess::BeginPrintProcess();

    metadata << "Center: (" << fCenter.X() << ", " << fCenter.Y() << ", " << fCenter.Z() << ")" << endl;
    metadata << "Theta angle: " << fTheta * 180. / TMath::Pi() << " degrees" << endl;
    metadata << "Phi angle: " << fPhi * 180. / TMath::Pi() << " degrees" << endl;
    metadata << "X-displacement: " << fDisplacement.X() << " mm" << endl;
    metadata << "Y-displacement: " << fDisplacement.Y() << " mm" << endl;
    metadata << " --------------------------- " << endl;
    metadata << " " << endl;
}

//////////////////////////////////////////////////////////////////////////
/// \brief Adds the footer for PrintMetadata
///
void TRestAxionEventProcess::EndPrintProcess() { TRestEventProcess::EndPrintProcess(); }