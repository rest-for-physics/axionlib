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

#ifndef _TRestAxionLikelihood
#define _TRestAxionLikelihood

#include <TRestMetadata.h>

#include "TRestAxionBufferGas.h"
#include "TRestAxionPhotonConversion.h"
#include "TRestAxionSolarModel.h"
#include "TRestAxionSpectrum.h"

#include "TRandom3.h"

//! A metadata class deninning a particular implementation of the likelihood to obtain the experimental
//! sensitivity
class TRestAxionLikelihood : public TRestMetadata {
   private:
    void Initialize();

    void InitFromConfigFile();

    // We consider all this values as constant parameters
    // In the future we will use the signal generator event chain
    // to determine a more accurate signal generation
    Double_t fBmag = 0;        //-> Manget field in T
    Double_t fRmag = 0;        //-> Magnet radius
    Double_t fLmag = 0;        //-> Manget length in mm
    Double_t fEfficiency = 1;  //->

    Int_t fNbores = 0;  //->

    Double_t fBackgroundLevel = 0.;  //->

    Double_t fSpotArea = 0.;  //->

    TVector2 fErange;  //->

    Double_t fTExpVacuum = 0;  //->

    Double_t fTExpPerStep = 0;  //->
    Int_t fNSteps = 0;          //->

    Double_t fLastStepDensity = 0.;  //->

    TRestAxionPhotonConversion* fPhotonConversion;  //!
    TRestAxionBufferGas* fBufferGas;                //!
    TRestAxionSpectrum* fAxionSpectrum;             //!

    /// Random number generator
    TRandom3* fRandom;  //!

    Int_t fMeasuredCountsVacuum;
    std::vector<Int_t> fMeasuredCountsPerStep;
    std::vector<Double_t> fExposureTimePerStep;
    std::vector<Double_t> fDensityInStep;

   public:
    void GenerateMonteCarlo();
    Double_t LogLikelihood(Double_t ma, Double_t g10, Double_t Nmeas, Double_t rho, Double_t tExp);

    Double_t GetSignal(Double_t ma, Double_t g10_4, Double_t rho, Double_t tExp);

    void LikelihoodTest(string fname);

    void PrintMetadata();

    // Constructors
    TRestAxionLikelihood();
    TRestAxionLikelihood(const char* cfgFileName, std::string name = "");
    // Destructor
    ~TRestAxionLikelihood();

    ClassDef(TRestAxionLikelihood, 1);
};
#endif
