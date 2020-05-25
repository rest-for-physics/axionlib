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
/// TRestAxionLikelihood just a class to serve as example of a specific
/// implementation of a metadata class.
///
/// Just copy .cxx and .h files to a new class name,
/// cp TRestAxionLikelihood.cxx TRestAxionSpecificMetadata.cxx
/// cp TRestAxionLikelihood.h TRestAxionSpecificMetadata.h
///
/// Then, inside each file replace TRestAxionLikelihood by
/// TRestAxionSpecificMetadata.
///
/// Re-run "cmake ../" at the build directory to include the new
/// class to the compilation.
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2019-March: First concept and implementation of TRestAxionLikelihood class.
///             Javier Galan
///
/// \class      TRestAxionLikelihood
/// \author     Javier Galan
///
/// <hr>
///

#include "TRestAxionLikelihood.h"
using namespace std;

ClassImp(TRestAxionLikelihood);
//______________________________________________________________________________
TRestAxionLikelihood::TRestAxionLikelihood() : TRestMetadata() {
    // TRestAxionLikelihood default constructor
    Initialize();
}

//______________________________________________________________________________
TRestAxionLikelihood::TRestAxionLikelihood(const char* cfgFileName, string name)
    : TRestMetadata(cfgFileName) {
    cout << "Entering TRestAxionLikelihood constructor( cfgFileName, name )" << endl;

    Initialize();

    LoadConfigFromFile(fConfigFileName, name);

    PrintMetadata();
}

//______________________________________________________________________________
TRestAxionLikelihood::~TRestAxionLikelihood() {
    // TRestAxionLikelihood destructor
}

void TRestAxionLikelihood::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    // Buffer gas properties definition (e.g. equivalent photon mass, absorption, etc)
    fBufferGas = new TRestAxionBufferGas();

    // Conversion probabilities as defined by van Bibber paper
    fPhotonConversion = new TRestAxionPhotonConversion();

    fPhotonConversion->AssignBufferGas(fBufferGas);

    // Solar axion flux on earth
    fAxionSpectrum = new TRestAxionSpectrum(fConfigFileName.c_str());

    fRandom = new TRandom3(0);
}

void TRestAxionLikelihood::GenerateMonteCarlo() {
    ////// Vacuum phase
    Double_t expTimeInSeconds = fTExpVacuum * 3600.;
    Double_t energyRange = fErange.Y() - fErange.X();

    debug << "Energy range : " << energyRange << endl;

    Double_t mean = fBackgroundLevel * fSpotArea * expTimeInSeconds * energyRange;

    fMeasuredCountsVacuum = fRandom->Poisson(mean);
    debug << "Vacuum phase. Mean counts : " << mean << endl;
    debug << "Vacuum phase. Measured counts : " << fMeasuredCountsVacuum << endl;

    ////// Gas phase

    fMeasuredCountsPerStep.clear();

    if (fTExpPerStep >= 0) {
        for (int n = 0; n < fNSteps; n++) {
            // Time
            Double_t exposureTime = fTExpPerStep;
            if (n < 3) exposureTime = 0;
            fExposureTimePerStep.push_back(exposureTime);

            // Density
            Double_t rho = fLastStepDensity * (n + 1) / fNSteps;
            fDensityInStep.push_back(rho);
        }
    } else if (fTExpPerStep == -2)  // TODO KSVZ exclusion
    {
        Double_t totalDays = 0;
        for (int n = 0; n < fNSteps; n++) {
            // Density
            Double_t rho = fLastStepDensity * (n + 1) / fNSteps;
            // Double_t rho = fLastStepDensity*(n+1)/fNSteps;

            fDensityInStep.push_back(rho);

            // Time
            fBufferGas->SetGasDensity("He", rho);

            Double_t mgamma = fBufferGas->GetPhotonMass(3.5);
            Double_t nGamma = GetSignal(mgamma, 1., rho, 1.);

            Double_t ksvzFactor = 3.75523 * mgamma;
            Double_t exposureTime = TMath::Log(20.) / TMath::Power(ksvzFactor, 4) / nGamma;  // in hours

            if (mgamma < 0.0475) exposureTime = 0;
            totalDays += exposureTime;

            fExposureTimePerStep.push_back(exposureTime);
        }

        for (int n = 0; n < fNSteps; n++) {
            cout << "Step : " << n << " Time : " << fExposureTimePerStep[n] / 12 << endl;
        }

        Double_t shareTime =
            (fExposureTimePerStep[4] + fExposureTimePerStep[5] + fExposureTimePerStep[6]) / 7.;
        cout << "Share time : " << shareTime / 12. << endl;
        for (int n = 0; n < 7; n++) fExposureTimePerStep[n] = shareTime;

        /////////

        Double_t totalTime = 0;
        for (int n = 0; n < fNSteps; n++) {
            cout << "Step : " << n << " Time : " << fExposureTimePerStep[n] / 12 << endl;
            totalTime += fExposureTimePerStep[n] / 12.;
        }

        GetChar();
    }

    /*
      debug << "Gas phase. Mean counts : " << mean << endl;
      GetChar();
    */

    for (int n = 0; n < fNSteps; n++) {
        mean = fBackgroundLevel * fSpotArea * fExposureTimePerStep[n] * 3600. * energyRange;
        Int_t counts = fRandom->Poisson(mean);
        fMeasuredCountsPerStep.push_back(counts);

        debug << "Step : " << n << " measured : " << counts << " :: " << fMeasuredCountsPerStep.back()
              << " counts" << endl;
        debug << "Time : " << fExposureTimePerStep[n] / 12 << " days" << endl;
    }

    debug << "Number of steps : " << fNSteps << " :: " << fMeasuredCountsPerStep.size() << endl;
}

void TRestAxionLikelihood::LikelihoodTest(string fname) {
    FILE* f = fopen(fname.c_str(), "wt");

    for (Double_t m = 0.008; m < 10; m = m * 1.04) {
        Double_t integral = 0;
        Double_t gBef = 0.;

        cout << "Calculating mass : " << m << " eV" << endl;
        cout << "-------------------------------_" << endl;
        for (Double_t g4 = 1.e-12; g4 < 1.e10; g4 = g4 * 1.02) {
            double l = LogLikelihood(m, g4, fMeasuredCountsVacuum, 0.0, fTExpVacuum);
            if (l == 0) break;

            // cout << "N steps : " << fNSteps << endl;
            for (int n = 0; n < fNSteps; n++) {
                Double_t rho = fDensityInStep[n];

                fBufferGas->SetGasDensity("He", rho);

                // We consider only neighbout steps
                if (fBufferGas->GetPhotonMass(3.5) - m > 0.004 || m - fBufferGas->GetPhotonMass(3.5) > 0.004)
                    continue;

                //                      cout << "Step n : " << n << endl;

                if (g4 == 1.e-12) cout << "Calculating Lhood for step " << n << endl;
                if (g4 == 1.e-12) cout << "step time " << fExposureTimePerStep[n] / 12. << " days" << endl;
                l = l * LogLikelihood(m, g4, fMeasuredCountsPerStep[n], rho, fExposureTimePerStep[n]);
                // cout << "n : " << n << " lH : " << l << endl;

                if (l == 0) break;
            }

            if (l == 0) break;
            //      GetChar();

            integral += l * (g4 - gBef);
            gBef = g4;

            // cout << "g4 : " << g4 << " l : " << l << endl;
            //      fprintf( f, "%e\t%e\n", g4, l );
        }
        //      GetChar();
        // fclose( f );

        Double_t int95 = 0;
        Double_t gLimit = 0;
        gBef = 0;
        for (Double_t g4 = 1.e-12; g4 < 1.e10; g4 = g4 * 1.02) {
            double l = LogLikelihood(m, g4, fMeasuredCountsVacuum, 0.0, fTExpVacuum);
            // double l = LogLikelihoodAverage( m, g4 , 0.0, fTExpVacuum );
            if (l == 0) break;

            // cout << "N steps : " << fNSteps << endl;
            for (int n = 0; n < fNSteps; n++) {
                Double_t rho = fDensityInStep[n];
                fBufferGas->SetGasDensity("He", rho);

                // We consider only neighbout steps
                if (fBufferGas->GetPhotonMass(3.5) - m > 0.004 || m - fBufferGas->GetPhotonMass(3.5) > 0.004)
                    continue;

                l = l * LogLikelihood(m, g4, fMeasuredCountsPerStep[n], rho, fExposureTimePerStep[n]);
                // l = l * LogLikelihoodAverage( m, g4, rho, fExposureTimePerStep[n] );
                // cout << "n : " << n << " lH : " << l << endl;

                if (l == 0) break;
            }

            if (l == 0) break;

            int95 += l * (g4 - gBef);
            gBef = g4;

            // cout << "g4 : " << g4 << " l : " << l << endl;
            // cout << "Int95 : " << int95 << " integral : " << integral << endl;
            if (int95 / integral > 0.95) {
                gLimit = g4;
                cout << "gLimit : " << g4 << endl;
                break;
            }

            //      fprintf( f, "%e\t%e\n", g4, l );
        }

        printf("ma : %e\t gL %e\n", m, sqrt(sqrt(gLimit)) * 1.e-10);
        // GetChar();
        fprintf(f, "%lf\t%e\n", m, sqrt(sqrt(gLimit)) * 1.e-10);
        fflush(f);
    }

    fclose(f);
}

Double_t TRestAxionLikelihood::GetSignal(Double_t ma, Double_t g10_4, Double_t rho, Double_t tExp) {
    fBufferGas->SetGasDensity("He", rho);

    Double_t magnetArea = fNbores * TMath::Pi() * TMath::Power(fRmag * REST_Units::cm, 2);

    Double_t signal = 0;
    Double_t dE = 0.01;
    for (Double_t en = fErange.X(); en < fErange.Y(); en = en + dE) {
        Double_t Phi_a = fAxionSpectrum->GetDifferentialSolarAxionFlux(en);
        // TODO Needs to be readapted to the new TRestAxionPhotonConversion implementation
        // Double_t Pa_g = fPhotonConversion->GammaTransmissionProbability(en, fBmag, ma);

        Double_t Pa_g = 0;  // Just dummy 0 probability!!
        Double_t nGamma = Pa_g * Phi_a;

        signal += nGamma;

        /*
          cout << "Energy : " << en << endl;
          cout << "Phi_a : " << Phi_a << endl;
          cout << "Pa_g : " << Pa_g << endl;
          cout << "Efficiency : " << fEfficiency << endl;
          cout << "Area magnet : " << magnetArea/10000 << endl;
          cout << "tExp : " << tExp << endl;
        */
    }

    return signal * dE * tExp * 3600. * magnetArea * fEfficiency * g10_4;
}

Double_t TRestAxionLikelihood::LogLikelihood(Double_t ma, Double_t g10_4, Double_t Nmeas, Double_t rho,
                                             Double_t tExp) {
    Double_t signal = GetSignal(ma, g10_4, rho, tExp);

    // cout << "Texp : " << tExp << endl;
    // cout << "g10_4 : " << g10_4 << " signal " << signal << " time : " <<  endl;

    // debug << "Ngamma vacuum : " << signal << endl;

    Double_t bckMean = fBackgroundLevel * (tExp * 3600.) * fSpotArea * (fErange.Y() - fErange.X());
    Double_t mu = signal + bckMean;

    Double_t lhood = TMath::Exp(-mu + Nmeas);
    if (Nmeas > 0) lhood += TMath::Exp(-mu + Nmeas) * pow(mu / Nmeas, Nmeas);

    // cout << "g10_4 : " << g10_4 << " signal " << signal << " Lh : " << lhood <<  endl;

    /*
      if( lhood == 0 )
      GetChar();

    */
    if (isinf(lhood) || isnan(lhood)) {
        cout << "IS INF" << endl;
        GetChar();
    }

    /*
      cout << "mu : " << mu << endl;
      cout << "Nmeas : " << Nmeas << endl;
      cout << "signal : " << signal << endl;
      cout << "----------------" << endl;
      cout << "lhood : " << lhood << endl;
    */

    // Double_t Ngamma = solarFluxG10 * magnetArea * (fTExpVacuum * 3600.) * fEfficiency;

    // Double_t nGammaVacuum =
    // fPhotonConversion->SetAxionMass( Double_t m ) { fAxionMass = m; }

    return lhood;
}

//______________________________________________________________________________
void TRestAxionLikelihood::InitFromConfigFile() {
    this->Initialize();

    // Initialize the metadata members from a configfile
    fBmag = GetDblParameterWithUnits("Bmag", 4.5);
    fRmag = GetDblParameterWithUnits("Rmag", 300);
    fLmag = GetDblParameterWithUnits("Lmag", 4500.);

    /// This needs to be reviewed. TRestAxionPhotonConversion does not store those parameter anymore.
    // fPhotonConversion->SetCoherenceLength(fLmag);
    // fPhotonConversion->SetMagneticField(fBmag);

    fEfficiency = StringToDouble(GetParameter("efficiency", "1"));

    fBackgroundLevel = StringToDouble(GetParameter("bckLevel", "1.e-8"));  // cpd keV-1 s-1

    fErange = StringTo2DVector(GetParameter("energyRange", "(1,8)"));  // keV

    fTExpVacuum = StringToDouble(GetParameter("expTimeVacuum", "12000"));  // In hours (vacuum phase)

    fTExpPerStep =
        StringToDouble(GetParameter("expTimePerStep", "54"));  // In hours (gas phase, time at each step)

    fSpotArea = StringToDouble(GetParameter("spotArea", "0.2"));  // in cm2

    fNSteps = StringToInteger(GetParameter("pressureSteps", "100"));

    fNbores = StringToInteger(GetParameter("bores", "2"));

    fLastStepDensity = StringToDouble(GetParameter("lastStepDensity", "0.1786e-3"));  // in g/cm3
}

void TRestAxionLikelihood::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    metadata << " Number of magnet bores : " << fNbores << endl;

    metadata << " Magnetic field : " << fBmag << " T" << endl;
    metadata << " Magnet length : " << fLmag << " mm" << endl;

    metadata << " Magnet radius : " << fRmag << " mm" << endl;
    metadata << " Signal overall efficiency : " << fEfficiency << endl;

    metadata << " Background level : " << fBackgroundLevel << " cpd keV-1 s-1" << endl;
    metadata << " Energy range : (" << fErange.X() << ", " << fErange.Y() << ")" << endl;

    metadata << " Vacuum phase exposure time : " << fTExpVacuum << " hours" << endl;
    metadata << " Gas phase exposure time per step : " << fTExpPerStep << endl;
    metadata << " Total number of steps : " << fNSteps << endl;
    metadata << " Last step Helium density : " << fLastStepDensity << " g/cm3" << endl;

    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}
