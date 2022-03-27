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
/// TRestAxionEvent is an event class used to define the properties of an
/// axion particle.
///
/// TODO. Create an appropriate documentation here (if needed).
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2019-March: First concept and implementation of TRestAxionEvent class.
///             Javier Galan
///
/// \class      TRestAxionEvent
/// \author     Javier Galan
///
/// <hr>
///
#include "TRestAxionEvent.h"
#include "TRestTools.h"

using namespace std;
using namespace TMath;

ClassImp(TRestAxionEvent);

TRestAxionEvent::TRestAxionEvent() {
    Initialize();
    fPad = NULL;
}

TRestAxionEvent::~TRestAxionEvent() {}

void TRestAxionEvent::Initialize() { TRestEvent::Initialize(); }

TPad* TRestAxionEvent::DrawEvent(TString option) {
    vector<string> optList = TRestTools::GetOptions((string)option);

    for (unsigned int n = 0; n < optList.size(); n++) {
        if (optList[n] == "print") this->PrintEvent();
    }

    optList.erase(std::remove(optList.begin(), optList.end(), "print"), optList.end());

    if (optList.size() == 0) optList.push_back("TODO");

    if (fPad != NULL) {
        delete fPad;
        fPad = NULL;
    }

    fPad = new TPad(this->GetName(), " ", 0, 0, 1, 1);
    fPad->Divide(1, 1);
    fPad->Draw();

    // TOBE implemented if necessary

    return fPad;
}

///////////////////////////////////////////////
/// \brief This method will produce a rotation respect to a `center` given by argument.
/// First we rotate the particle and direction along the X-axis by an angle `theta`, then
/// we rotate the particle and direction along the Z-axis by an angle `phi`.
///
void TRestAxionEvent::RotateXZ(const TVector3& center, Double_t theta, Double_t phi) {
    TVector3 ref = fPosition - center;

    ref.RotateX(theta);
    ref.RotateZ(phi);

    fPosition = ref + center;

    fDirection.RotateX(theta);
    fDirection.RotateZ(phi);
}

///////////////////////////////////////////////
/// \brief This method will produce a rotation respect to a `center` given by argument.
/// First we rotate the particle and direction along the Z-axis by an angle `phi`, then
/// we rotate the particle and direction along the X-axis by an angle `theta`.
///
void TRestAxionEvent::RotateZX(const TVector3& center, Double_t phi, Double_t theta) {
    TVector3 ref = fPosition - center;

    ref.RotateZ(phi);
    ref.RotateX(theta);

    fPosition = ref + center;

    fDirection.RotateZ(phi);
    fDirection.RotateX(theta);
}

void TRestAxionEvent::PrintEvent() {
    TRestEvent::PrintEvent();

    cout << "Energy : " << GetEnergy() << endl;
    cout << "Position : ( " << fPosition.X() << ", " << fPosition.Y() << ", " << fPosition.Z() << " )"
         << endl;
    cout << "Direction : ( " << fDirection.X() << ", " << fDirection.Y() << ", " << fDirection.Z() << " )"
         << endl;
    cout << "Gamma state probability : " << fGammaProbability << endl;
    cout << endl;
}
