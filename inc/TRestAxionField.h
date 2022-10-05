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

#ifndef _TRestAxionField
#define _TRestAxionField

#include "TRestAxionBufferGas.h"
#include "TRestComplex.h"

/*
/// MOVED TO TRestAxionFieldPropagationProcess class
/// A structure to define the two components of a complex number using real precision.
/// To be used inside TRestAxionField.
struct ComplexReal {
    /// The real part of the number
    mpfr::mpreal real = 0;

    /// The imaginary part of the number
    mpfr::mpreal img = 0;
};
*/
//! A basic class to define analytical axion-photon conversion calculations for axion helioscopes
class TRestAxionField : public TObject {
   private:
    /// A two component vector to store the complex EM field amplitude.
    /// MOVED to TRestAxionFieldPropagationProcess
    /// ComplexReal fAem;  //!

    /// A two component vector to store the complex axion field amplitude.
    /// MOVED to TRestAxionFieldPropagationProcess
    /// ComplexReal faxion;  //!

    Bool_t fDebug = false;  //!

    void Initialize();

    /// A pointer to the buffer gas definition
    TRestAxionBufferGas* fBufferGas = NULL;  //!

   public:
    Double_t BL(Double_t Bmag, Double_t Lcoh);
    Double_t BLHalfSquared(Double_t Bmag, Double_t Lcoh);

    /// It enables/disables debug mode
    void SetDebug(Bool_t v) { fDebug = v; }

    /// It assigns a gas buffer medium to the calculation
    void AssignBufferGas(TRestAxionBufferGas* buffGas) { fBufferGas = buffGas; }

    /// It assigns a gas buffer medium to the calculation
    void SetBufferGas(TRestAxionBufferGas* buffGas) { fBufferGas = buffGas; }

    Double_t GammaTransmissionProbability(Double_t Bmag, Double_t Lcoh, Double_t Ea, Double_t ma,
                                          Double_t mg = 0, Double_t absLength = 0);

    Double_t GammaTransmissionProbability(std::vector<Double_t> Bmag, Double_t deltaL, Double_t Ea,
                                          Double_t ma, Double_t mg = 0, Double_t absLength = 0);

    Double_t AxionAbsorptionProbability(Double_t Bmag, Double_t Lcoh, Double_t Ea, Double_t ma,
                                        Double_t mg = 0, Double_t absLength = 0);

    /// Commented because it uses ComplexReal structure that is moved to TRestAxionFieldPropagationProcess
    /// class
    /// void PropagateAxion(Double_t Bmag, Double_t Lcoh, Double_t Ea, Double_t ma, Double_t mg = 0,
    ///                    Double_t absLength = 0);

    TRestAxionField();
    ~TRestAxionField();

    ClassDef(TRestAxionField, 2);
};
#endif
