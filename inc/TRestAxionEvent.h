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

#include "TObject.h"
#include "TMath.h"
#include "TVector3.h"

#include "TRestEvent.h"

/// An event data class to define the parameters related to an axion particle
class TRestAxionEvent : public TRestEvent
{
    private:
		TVector3 fPosition; //->
		TVector3 fDirection; //->
		Double_t fEnergy = 0; //->

		Double_t fGammaProbability = 0; //-> Its square is the P_{ag} conversion probability

    protected:

    public:
        
		TVector3 *GetPosition( ) { return &fPosition; }

        Double_t GetPositionX( ) { return fPosition.X(); } // returns value in mm
        Double_t GetPositionY( ) { return fPosition.Y(); } // returns value in mm
        Double_t GetPositionZ( ) { return fPosition.Z(); } // returns value in mm
        
		TVector3 *GetDirection( ) { return &fDirection; }

        Double_t GetDirectionX( ) { return fDirection.X(); } // returns normalized vector x-component.
        Double_t GetDirectionY( ) { return fDirection.Y(); } // returns normalized vector y-component 
        Double_t GetDirectionZ( ) { return fDirection.Z(); } // returns normalized vector z-component

        Double_t GetEnergy( ) { return fEnergy; } //returns value in keV

		Double_t GetGammaProbability( ) { return fGammaProbability; }

		void SetPosition( TVector3 pos ) { fPosition = pos; }
		void SetPosition( Double_t x, Double_t y, Double_t z ) { SetPosition( TVector3( x, y, z) ); }

		void SetDirection( TVector3 dir ) { fDirection = dir; }
		void SetDirection( Double_t px, Double_t py, Double_t pz ) { SetDirection( TVector3( px, py, pz) ); }

		void SetEnergy( Double_t en ) { fEnergy = en; }

		void SetGammaProbability( Double_t p ) { fGammaProbability = p; }

        virtual void Initialize();

        virtual void PrintEvent( );

        TPad *DrawEvent( TString option = "" );

        //Construtor
        TRestAxionEvent();
        //Destructor
        ~TRestAxionEvent();

        ClassDef(TRestAxionEvent, 1);
};
#endif
