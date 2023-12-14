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

//! A basic class to define analytical axion-photon conversion calculations for axion helioscopes
class TRestAxionField : public TObject {
   private:
    Bool_t fDebug = false;  //!

    /// The magnetic field in Teslas
    Double_t fBmag = 2.5;

    /// The coherence lenght (in mm) where the magnetic field is defined
    Double_t fLcoh = 10000;

    /// The energy of the axion in keV
    Double_t fEa = 4.2;

    void Initialize();

    /// A pointer to the buffer gas definition
    TRestAxionBufferGas* fBufferGas = NULL;  //!

   public:
    void SetMagneticField(Double_t b) { fBmag = b; }
    void SetCoherenceLength(Double_t l) { fLcoh = l; }
    void SetAxionEnergy(Double_t e) { fEa = e; }

    Double_t GetMagneticField() const { return fBmag; }
    Double_t GetCoherenceLength() const { return fLcoh; }
    Double_t GetAxionEnergy() const { return fEa; }

    Double_t BL(Double_t Bmag, Double_t Lcoh);
    Double_t BLHalfSquared(Double_t Bmag, Double_t Lcoh);

    /// It enables/disables debug mode
    void SetDebug(Bool_t v) { fDebug = v; }

    /// It assigns a gas buffer medium to the calculation
    void AssignBufferGas(TRestAxionBufferGas* buffGas) { fBufferGas = buffGas; }

    /// It assigns a gas buffer medium to the calculation
    void SetBufferGas(TRestAxionBufferGas* buffGas) { fBufferGas = buffGas; }

    Double_t GammaTransmissionProbability(Double_t ma, Double_t mg = 0, Double_t absLength = 0);

    Double_t GammaTransmissionProbability(Double_t Bmag, Double_t Lcoh, Double_t Ea, Double_t ma,
                                          Double_t mg = 0, Double_t absLength = 0);

    Double_t AxionAbsorptionProbability(Double_t ma, Double_t mg = 0, Double_t absLength = 0);

    Double_t AxionAbsorptionProbability(Double_t Bmag, Double_t Lcoh, Double_t Ea, Double_t ma,
                                        Double_t mg = 0, Double_t absLength = 0);

    Double_t GammaTransmissionProbability(std::vector<Double_t> Bmag, Double_t deltaL, Double_t Ea,
                                          Double_t ma, Double_t mg = 0, Double_t absLength = 0);

    Double_t GammaTransmissionFWHM(Double_t step = 0.00001);

    std::vector<std::pair<Double_t, Double_t>> GetMassDensityScanning(std::string gasName = "He",
                                                                      Double_t maMax = 1, Double_t scale = 1);

    TRestAxionField();
    ~TRestAxionField();

    ClassDef(TRestAxionField, 2);
};
#endif
