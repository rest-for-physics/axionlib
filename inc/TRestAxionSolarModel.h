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

#ifndef _TRestAxionSolarModel
#define _TRestAxionSolarModel

#include <TRestMetadata.h>

//! A metadata class to define theoretical axion models and calculations related
class TRestAxionSolarModel : public TRestMetadata {
   private:
    void Initialize();

    void InitFromConfigFile();

    TString fSolarAxionModel;  //->

    TString fMode;  //->

    // Integrated solar axion spectrum in cm-2 s-1
    Double_t fSolarEnergyFlux = 0;  //->

    // The axion-photon g10 coupling
    Double_t fg10 = 1.;  //->

    // The integration step for solar energy spectrum
    Double_t fStep = 1.e-3;  //->

    TVector2 fEnergyRange;  //->

    // Contains the tabulated solar disk in solar radius and energy. As
    // described in data/solarModel.
    std::vector<std::vector<Double_t>> fSolarTable;  //->

   public:
    Bool_t isSolarTableLoaded() { return fSolarTable.size() > 0; }

    TString GetSolarAxionSolarModel() { return fSolarAxionModel; }

    void SetSolarAxionSolarModel(TString modelName) { fSolarAxionModel = modelName; }

    void ResetSolarEnergyFlux() { fSolarEnergyFlux = 0; }

    Double_t GetSolarAxionFlux(Double_t eMin = 0., Double_t eMax = 10., Double_t g10 = 1.,
                               Double_t step = 0.001);

    Double_t GetDifferentialSolarAxionFlux(Double_t energy, Double_t g10 = 1.);

    void PrintMetadata();

    // Constructors
    TRestAxionSolarModel();
    TRestAxionSolarModel(const char* cfgFileName, std::string name = "");
    // Destructor
    ~TRestAxionSolarModel();

    ClassDef(TRestAxionSolarModel, 1);
};
#endif
