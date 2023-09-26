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
/// TRestAxionHiddenPhotonField is a class used to calculate the axion-photon mixing
/// and determine the probability of the particle being in an axion or photon
/// state.
///
/// A peculiarity from this class is that it encapsulates internally the high
/// precision calculations using the real precisions types using TRestComplex.
/// It is known that double precision is not good enough in some scenarios.
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2019-March: First concept and implementation of TRestAxionHiddenPhotonField class.
///             Javier Galan
/// 2023-September: Separation of QCD and HiddenPhoton classes.
///             Tom O'Shea
///
/// \class      TRestAxionHiddenPhotonField
/// \author     Javier Galan
///
/// <hr>
///
#include "TRestAxionHiddenPhotonField.h"

#include <TVectorD.h>

#include "TH1F.h"

#ifdef USE_MPFR
#include "TRestComplex.h"
#endif

#include <numeric>

using namespace std;

ClassImp(TRestAxionHiddenPhotonField);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionHiddenPhotonField::TRestAxionHiddenPhotonField() { Initialize(); }

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionHiddenPhotonField::~TRestAxionHiddenPhotonField() {}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionHiddenPhotonField class
///
/// It sets the default real precision to be used with mpfr types. Now it is 30 digits.
/// So that we can still calculate numbers such as : 1.0 - 1.e-30
///
void TRestAxionHiddenPhotonField::Initialize() {
#ifdef USE_MPFR
    TRestComplex::SetPrecision(30);
#endif

    fBufferGas = NULL;

    /// MOVED TO TRestAxionHiddenPhotonFieldPropagationProcess class
    /// faxion = SetComplexReal(1, 0);
    /// fAem = SetComplexReal(0, 0);
}


///////////////////////////////////////////////
///	vacuum conversion probability
/// ie. when photonMass = 0
Double_t TRestAxionHiddenPhotonField::VacuumConversion(Double_t Lcoh, Double_t Ea, Double_t ma) {
	Double_t q = momentumTransfer(Ea, ma, 0.);
	return pow(2*sin(q*Lcoh/2), 2);
}



///////////////////////////////////////////////
/// \brief Performs the calculation of hidden photon - photon conversion probability.
///
/// If m_gamma (mg) is not given as an argument, i.e. it is equal to zero, then m_gamma
/// will be obtainned from the buffer gas definition. If no buffer gas has been assigned
/// then the medium will be assumed to be vacuum.
///
/// Ea in keV, ma in eV, Lcoh in mm, mg in eV, absLength in cm-1
///
/// The returned value is given for chi = 1
///
Double_t TRestAxionHiddenPhotonField::GammaTransmissionProbability(Double_t Lcoh, Double_t Ea, Double_t ma,
                                                       Double_t mg, Double_t absLength) {
#ifndef USE_MPFR
    RESTWarning
        << "MPFR libraries not linked to REST libraries. Try adding -DREST_MPFR=ON to your REST compilation"
        << RESTendl;
    RESTWarning << "TRestAxionHiddenPhotonField::GammaTransmissionProbability will return 0" << RESTendl;
    return 0;
#else
    mpfr::mpreal hiddenPhotonMass = ma;
    mpfr::mpreal cohLength = Lcoh / 1000.;  // Default REST units are mm;

    mpfr::mpreal photonMass = mg;

    if (mg == 0 && fBufferGas) photonMass = fBufferGas->GetPhotonMass(Ea);

    RESTDebug << "+--------------------------------------------------------------------------+" << RESTendl;
    RESTDebug << " TRestAxionHiddenPhotonField::GammaTransmissionProbability. Parameter summary" << RESTendl;
    RESTDebug << " Photon mass : " << photonMass << " eV" << RESTendl;
    RESTDebug << " Hidden photon mass : " << ma << " eV" << RESTendl;
    RESTDebug << " Hidden photon energy : " << Ea << " keV" << RESTendl;
    RESTDebug << " Lcoh : " << Lcoh << " mm" << RESTendl;
    RESTDebug << " Bmag : " << Bmag << " T" << RESTendl;
    RESTDebug << "+--------------------------------------------------------------------------+" << RESTendl;

    if (photonMass == 0.0) return VacuumConversion(Lcoh, Ea, ma);

    mpfr::mpreal q = momentumTransfer(Ea, ma, mg)			// eV
    mpfr::mpreal l = cohLength * REST_Physics::PhMeterIneV;		// eV-1

    mpfr::mpreal Gamma = absLength;		// cm-1
    if (absLength == 0 && fBufferGas) Gamma = fBufferGas->GetPhotonAbsorptionLength(Ea);  // cm-1
    mpfr::mpreal GammaL = Gamma * cohLength * 100;
	Gamma *= 100 / REST_Physics::PhMeterIneV;	// eV

    if (fDebug) {
        RESTDebug << "+------------------------+" << RESTendl;
        RESTDebug << " Intermediate calculations" << RESTendl;
        RESTDebug << " q : " << q << " eV" << RESTendl;
        RESTDebug << " l : " << l << " eV-1" << RESTendl;
        RESTDebug << " ql : " << q*l << RESTendl;
        RESTDebug << "Gamma : " << Gamma << RESTendl;
        RESTDebug << "GammaL : " << GammaL << RESTendl;
        RESTDebug << "+------------------------+" << RESTendl;
    }

    mpfr::mpreal MFactor = pow(ma,4)/(pow(ma*ma - mg*mg, 2) + pow(Ea*Gamma, 2));	// all should be in eV


    if (fDebug) {
        RESTDebug << "Mfactor : " << MFactor << RESTendl;
        RESTDebug << "cos(ql) : " << cos(q*l) << RESTendl;
        RESTDebug << "Exp(-GammaL) : " << exp(-GammaL) << RESTendl;
    }

    double sol =
        (double)(MFactor * (1 + exp(-GammaL) - 2 * exp(-GammaL / 2) * cos(q*l)));

    RESTDebug << "Axion-photon transmission probability : " << sol << RESTendl;

    return sol;
#endif
}



///////////////////////////////////////////////
/// \brief Performs the calculation of the photon absorbtion probability in the buffer gas.
///
/// If not provided as argument, m_g it will be attempted to obtain the value from the buffer gas
/// definition. If no buffer gas has been assigned the medium will be assumed to be vacuum.
///
/// Ea in keV, ma in eV, Lcoh in mm, mg in eV, absLength in cm-1
///
/// The returned value is given for chi = 1
///
Double_t TRestAxionHiddenPhotonField::AxionAbsorptionProbability(Double_t Bmag, Double_t Lcoh, Double_t Ea, Double_t ma,
                                                     Double_t mg, Double_t absLength) {
#ifndef USE_MPFR
    RESTWarning
        << "MPFR libraries not linked to REST libraries. Try adding -DREST_MPFr=ON to your REST compilation"
        << RESTendl;
    RESTWarning << "TRestAxionHiddenPhotonField::GammaTransmissionProbability will return 0" << RESTendl;
    return 0;
#else
    mpfr::mpreal axionMass = ma;
    mpfr::mpreal cohLength = Lcoh / 1000.;  // Default REST units are mm;

    mpfr::mpreal photonMass = mg;
    if (mg == 0 && fBufferGas) photonMass = fBufferGas->GetPhotonMass(Ea);

    if (fDebug) {
        RESTDebug << "+--------------------------------------------------------------------------+"
                  << RESTendl;
        RESTDebug << " TRestAxionHiddenPhotonField::GammaTransmissionProbability. Parameter summary" << RESTendl;
        RESTDebug << " Photon mass : " << photonMass << " eV" << RESTendl;
        RESTDebug << " Axion mass : " << ma << " eV" << RESTendl;
        RESTDebug << " Axion energy : " << Ea << " keV" << RESTendl;
        RESTDebug << " Lcoh : " << Lcoh << " mm" << RESTendl;
        RESTDebug << " Bmag : " << Bmag << " T" << RESTendl;
        RESTDebug << "+--------------------------------------------------------------------------+"
                  << RESTendl;
    }

    if (ma == 0.0 && photonMass == 0.0) return BLHalfSquared(Bmag, Lcoh);

    mpfr::mpreal q = (ma * ma - photonMass * photonMass) / 2. / Ea / 1000.0;
    mpfr::mpreal l = cohLength * REST_Physics::PhMeterIneV;
    mpfr::mpreal phi = q * l;

    mpfr::mpreal Gamma = absLength;
    if (absLength == 0 && fBufferGas) Gamma = fBufferGas->GetPhotonAbsorptionLength(Ea);  // cm-1
    mpfr::mpreal GammaL = Gamma * cohLength * 100;

    if (fDebug) {
        RESTDebug << "+------------------------+" << RESTendl;
        RESTDebug << " Intermediate calculations" << RESTendl;
        RESTDebug << " q : " << q << " eV" << RESTendl;
        RESTDebug << " l : " << l << " eV-1" << RESTendl;
        RESTDebug << " phi : " << phi << RESTendl;
        RESTDebug << "Gamma : " << Gamma << RESTendl;
        RESTDebug << "GammaL : " << GammaL << RESTendl;
        RESTDebug << "+------------------------+" << RESTendl;
    }

    mpfr::mpreal MFactor = phi * phi + GammaL * GammaL / 4.0;
    MFactor = 1.0 / MFactor;

    if (fDebug) {
        RESTDebug << "Mfactor : " << MFactor << RESTendl;
        RESTDebug << "(BL/2)^2 : " << BLHalfSquared(Bmag, Lcoh) << RESTendl;
        RESTDebug << "cos(phi) : " << cos(phi) << RESTendl;
        RESTDebug << "Exp(-GammaL) : " << exp(-GammaL) << RESTendl;
    }

    double sol = (double)(MFactor * BLHalfSquared(Bmag, Lcoh) * GammaL);

    if (fDebug) RESTDebug << "Axion-photon absorption probability : " << sol << RESTendl;

    return sol;
#endif
}

/// Commented because it uses ComplexReal structure that is moved to TRestAxionHiddenPhotonFieldPropagationProcess class
/*
void TRestAxionHiddenPhotonField::PropagateAxion(Double_t Bmag, Double_t Lcoh, Double_t Ea, Double_t ma,
                                                Double_t mg, Double_t absLength) {
    mpfr::mpreal axionMass = ma;
    mpfr::mpreal cohLength = Lcoh / 1000.;  // Default REST units are mm;

    mpfr::mpreal photonMass = mg;
    if (mg == 0 && fBufferGas) photonMass = fBufferGas->GetPhotonMass(Ea);

    if (fDebug) {
        RESTDebug << "+--------------------------------------------------------------------------+" <<
RESTendl; RESTDebug << " TRestAxionHiddenPhotonField::GammaTransmissionProbability. Parameter summary" <<
RESTendl; RESTDebug << " Photon mass : " << photonMass << " eV" << RESTendl; RESTDebug << " Axion mass : " <<
ma << " eV" << RESTendl; RESTDebug << " Axion energy : " << Ea << " keV" << RESTendl; RESTDebug << " Lcoh : "
<< Lcoh << " mm" << RESTendl; RESTDebug << " Bmag : " << Bmag << " T" << RESTendl; RESTDebug <<
"+--------------------------------------------------------------------------+" << RESTendl;
    }

    mpfr::mpreal q = (ma * ma - photonMass * photonMass) / 2. / Ea / 1000.0;
    mpfr::mpreal l = cohLength * PhMeterIneV;
    mpfr::mpreal phi = q * l;

    mpfr::mpreal Gamma = absLength;
    if (absLength == 0 && fBufferGas) Gamma = fBufferGas->GetPhotonAbsorptionLength(Ea);  // cm-1
    mpfr::mpreal GammaL = Gamma * cohLength * 100;

    if (fDebug) {
        RESTDebug << "+------------------------+" << RESTendl;
        RESTDebug << " Intermediate calculations" << RESTendl;
        RESTDebug << " q : " << q << " eV" << RESTendl;
        RESTDebug << " l : " << l << " eV-1" << RESTendl;
        RESTDebug << " phi : " << phi << RESTendl;
        RESTDebug << "Gamma : " << Gamma << RESTendl;
        RESTDebug << "GammaL : " << GammaL << RESTendl;
        RESTDebug << "+------------------------+" << RESTendl;
    }

    mpfr::mpreal bl = BL(Bmag, Lcoh);

    /// We have now calculated the main quantities BL, QL, and GammaL

    ComplexReal I = SetComplexReal(0, 1);
    ComplexReal ExpPhi = SetComplexReal(cos(-phi), sin(-phi));

    ComplexReal Bcomplex = SetComplexReal(BL(Bmag, Lcoh), 0);
    ComplexReal Qcomplex = SetComplexReal(phi, -GammaL / 2);

    ComplexReal Bterm = ComplexCocient(Bcomplex, Qcomplex);
    Bterm = ComplexProduct(I, Bterm);

    mpfr::mpreal ExpGamma = exp(-GammaL / 2.);
    Double_t ExpGammaDouble = TMath::Exp((Double_t)-GammaL / 2.);

    cout.precision(30);

    if (fDebug) {
        cout << "ExpGamma : " << ExpGamma << RESTendl;
        cout << "ExpGammaDouble : " << ExpGammaDouble << RESTendl;
        RESTDebug << "(BL/2)^2 : " << BLHalfSquared(Bmag, Lcoh) << RESTendl;
        RESTDebug << "cos(phi) : " << cos(phi) << RESTendl;
        RESTDebug << "Exp(-GammaL) : " << exp(-GammaL) << RESTendl;
    }

    // if (fDebug) RESTDebug << "Axion-photon absorption probability : " << sol << RESTendl;
}
*/
