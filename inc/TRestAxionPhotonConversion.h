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

#ifndef _TRestAxionPhotonConversion
#define _TRestAxionPhotonConversion

#include "TRestAxionBufferGas.h"

//! A basic class to define analytical axion-photon conversion calculations for axion helioscopes
class TRestAxionPhotonConversion : public TObject {
   private:
    void Initialize();

    /// Axion mass in eV
    Double_t fAxionMass = 0;  //->

    /// Coherence length in mm [REST default units]
    Double_t fCohLength = 0;  //->

    /// Magnet field intensity in T
    Double_t fBMag = 0;  //->

    /// The axion-photon g10 coupling
    Double_t fg10 = 1.;  //->

    /// A pointer to the buffer gas definition
    TRestAxionBufferGas* fBufferGas = NULL;  //!

    void UpdateParameters(Double_t& Bmag, Double_t& ma, Double_t& Lcoh);

    Double_t BL(Double_t Lcoh = -1, Double_t Bmag = -1);
    Double_t BLHalfSquared(Double_t Lcoh = -1, Double_t Bmag = -1);

   public:
    /// It assigns a gas buffer medium to the calculation
    void AssignBufferGas(TRestAxionBufferGas* buffGas) { fBufferGas = buffGas; }

    /// It assigns a gas buffer medium to the calculation
    void SetBufferGas(TRestAxionBufferGas* buffGas) { fBufferGas = buffGas; }

    /// Sets the value of `fAxionMass` metadata member to `m` in `eV`.
    void SetAxionMass(Double_t m) { fAxionMass = m; }

    /// Sets the value of `fCohLength` metadata member to `l` in `mm`.
    void SetCoherenceLength(Double_t l) { fCohLength = l; }

    /// Sets the value of `fBmag` metadata member to `B` in `T`.
    void SetMagneticField(Double_t B) { fBMag = B; }

    /// Returns the value stored in `fAxionMass` in `eV`.
    Double_t GetAxionMass() { return fAxionMass; }

    /// Returns the value stored in `fCohLength` in `mm`.
    Double_t GetCoherenceLength() { return fCohLength; }

    /// Returns the value stored in `fBMag` in `T`.
    Double_t GetMagneticField() { return fBMag; }

    Double_t GammaTransmissionProbability(Double_t Ea, Double_t Bmag = -1, Double_t ma = -1,
                                          Double_t Lcoh = -1);

    Double_t GammaTransmissionProbability(Double_t Ea, TVectorD B, Double_t ma = -1, Double_t Lcoh = -1);

    TRestAxionPhotonConversion();
    ~TRestAxionPhotonConversion();

    ClassDef(TRestAxionPhotonConversion, 1);
};
#endif
