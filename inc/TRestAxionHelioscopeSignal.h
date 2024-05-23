/*************************************************************************
 * This file is part of the REST software framework.                     *
 *                                                                       *
 * Copyright (C) 2016 GIFNA/TREX (University of Zaragoza)                *
 * For more information see https://gifna.unizar.es/trex                 *
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
 * If not, see https://www.gnu.org/licenses/.                            *
 * For the list of contributors see $REST_PATH/CREDITS.                  *
 *************************************************************************/

#ifndef REST_TRestAxionHelioscopeSignal
#define REST_TRestAxionHelioscopeSignal

#include <TRestAxionBufferGas.h>
#include <TRestAxionField.h>
#include <TRestAxionSolarFlux.h>

#include "TRestComponent.h"

/// It allows to build a signal model using a given helioscope configuration and solar axion flux component
class TRestAxionHelioscopeSignal : public TRestComponent {
   private:
    /// It defines the helioscope type (IAXO/AMELIE)
    std::string fConversionType = "IAXO";

    /// It defines the number of bores or TPCs (number of magnetic volumes)
    Int_t fBores = 2;

    /// It defines the magnet aperture radius in standard REST units (mm)
    Double_t fMagnetRadius = 350;

    /// It defines the magnet length in standard REST units (mm)
    Double_t fMagnetLength = 10000;

    /// It defines the magnetic field strength in T
    Double_t fMagnetStrength = 2;

    /// If optics is present we may add an efficiency for Ngamma calculation
    Double_t fOpticsEfficiency = 1;

    /// If an x-ray window is present we may add an efficiency for Ngamma calculation
    Double_t fWindowEfficiency = 1;

    /// It defines the gas mixture we use inside our magnetic field. Vacuum if it is nullptr
    TRestAxionBufferGas* fGas = nullptr;

    /// It defines the solar flux we use inside our magnetic field. Must be defined.
    TRestAxionSolarFlux* fFlux = nullptr;

    /// It provides access to axion-photon conversion probabilites calculations
    TRestAxionField* fField = nullptr;

   protected:
    void InitFromConfigFile() override;

    void FillHistograms() override;

   public:
    Double_t GetSignalRate(std::vector<Double_t> point, Double_t mass = 0);

	TRestAxionBufferGas *GetGas() { return fGas; }

    void PrintMetadata() override;

    void SetNumberOfBores(const Int_t& bores) { fBores = bores; }
    void SetMagnetRadius(const Double_t& radius) { fMagnetRadius = radius; }
    void SetMagnetLength(const Double_t& length) { fMagnetLength = length; }
    void SetMagnetStrength(const Double_t& strength) { fMagnetStrength = strength; }
    void SetType(const std::string& type) { fConversionType = type; }

    Int_t GetNumberOfBores() const { return fBores; }
    Double_t GetMagnetRadius() const { return fMagnetRadius; }
    Double_t GetMagnetLength() const { return fMagnetLength; }
    Double_t GetMagnetStrength() const { return fMagnetStrength; }
    std::string GetType() const { return fConversionType; }

    void Initialize() override;
    TRestAxionHelioscopeSignal(const char* cfgFileName, const std::string& name);
    TRestAxionHelioscopeSignal();
    ~TRestAxionHelioscopeSignal();

    ClassDefOverride(TRestAxionHelioscopeSignal, 1);
};
#endif
