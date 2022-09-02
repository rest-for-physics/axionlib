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
/// TRestAxionPhotonConversion is a class used to calculate the axion-photon mixing
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
/// 2019-March: First concept and implementation of TRestAxionPhotonConversion class.
///             Javier Galan
///
/// \class      TRestAxionPhotonConversion
/// \author     Javier Galan
///
/// <hr>
///
#include "TRestAxionPhotonConversion.h"
#include <TVectorD.h>
#include "TH1F.h"

using namespace std;

// Better we keep specifying  mpfr:: explicitily
// using mpfr::mpreal;

ClassImp(TRestAxionPhotonConversion);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionPhotonConversion::TRestAxionPhotonConversion() { Initialize(); }

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionPhotonConversion::~TRestAxionPhotonConversion() {}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionPhotonConversion class
///
/// It sets the default real precision to be used with mpfr types. Now it is 30 digits.
/// So that we can still calculate numbers such as : 1.0 - 1.e-30
///
void TRestAxionPhotonConversion::Initialize() {
    TRestComplex::SetPrecision(30);

    fBufferGas = NULL;

    /// MOVED TO TRestAxionFieldPropagationProcess class
    /// faxion = SetComplexReal(1, 0);
    /// fAem = SetComplexReal(0, 0);
}

///////////////////////////////////////////////
/// \brief Performs the calculation of (BL) factor in natural units.
///
/// `Lcoh` should be expressed in `mm`, and `Bmag` in `T`.
/// The result will be given for an axion-photon coupling of 10^{-10} GeV^{-1}
///
double TRestAxionPhotonConversion::BL(Double_t Bmag, Double_t Lcoh) {
    Double_t lengthInMeters = Lcoh / 1000.;

    Double_t tm = REST_Physics::lightSpeed / REST_Physics::naturalElectron * 1.0e-9;  // GeV
    Double_t sol = lengthInMeters * Bmag * tm;
    sol = sol * 1.0e-10;

    return sol;
}

///////////////////////////////////////////////
/// \brief Performs the calculation of (BL/2)^2 factor in natural units.
///
/// `Lcoh` should be expressed in `mm`, and `Bmag` in `T`.
/// The result will be given for an axion-photon coupling of 10^{-10} GeV^{-1}
///
double TRestAxionPhotonConversion::BLHalfSquared(Double_t Bmag, Double_t Lcoh)  // (BL/2)**2
{
    Double_t lengthInMeters = Lcoh / 1000.;

    Double_t tm = REST_Physics::lightSpeed / REST_Physics::naturalElectron * 1.0e-9;  // gev
    Double_t sol = lengthInMeters * Bmag * tm / 2;
    sol = sol * sol * 1.0e-20;

    return sol;
}

///////////////////////////////////////////////
/// \brief Performs the calculation of axion-photon conversion probability using directly
/// equation (11) from van Bibber, Phys Rev D Part Fields. 1989.
///
/// m_gamma will be obtainned from buffer gas definition. If no buffer gas has been assigned the medium
/// will be assumed to be vacuum.
///
/// Ea in keV, ma in eV, mgamma in eV, Lcoh in mm, Bmag in T
///
/// mg in eV, absLength in cm-1
///
/// The returned value is given for g_ag = 10^-10 GeV-1
///
Double_t TRestAxionPhotonConversion::GammaTransmissionProbability(Double_t Bmag, Double_t Lcoh, Double_t Ea,
                                                                  Double_t ma, Double_t mg,
                                                                  Double_t absLength) {
    mpfr::mpreal axionMass = ma;
    mpfr::mpreal cohLength = Lcoh / 1000.;  // Default REST units are mm;

    mpfr::mpreal photonMass = mg;

    if (mg == 0 && fBufferGas) photonMass = fBufferGas->GetPhotonMass(Ea);

    RESTDebug << "+--------------------------------------------------------------------------+" << RESTendl;
    RESTDebug << " TRestAxionPhotonConversion::GammaTransmissionProbability. Parameter summary" << RESTendl;
    RESTDebug << " Photon mass : " << photonMass << " eV" << RESTendl;
    RESTDebug << " Axion mass : " << ma << " eV" << RESTendl;
    RESTDebug << " Axion energy : " << Ea << " keV" << RESTendl;
    RESTDebug << " Lcoh : " << Lcoh << " mm" << RESTendl;
    RESTDebug << " Bmag : " << Bmag << " T" << RESTendl;
    RESTDebug << "+--------------------------------------------------------------------------+" << RESTendl;

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

    double sol =
        (double)(MFactor * BLHalfSquared(Bmag, Lcoh) * (1 + exp(-GammaL) - 2 * exp(-GammaL / 2) * cos(phi)));

    RESTDebug << "Axion-photon transmission probability : " << sol << RESTendl;

    return sol;
}

///////////////////////////////////////////////
/// \brief Performs the calculation of axion-photon absorption probability using directly
/// equation (18) from van Bibber, Phys Rev D Part Fields. 1989.
///
/// If not provided as argument, m_gamma it will be attempted to obtain the value from the buffer gas
/// definition. If no buffer gas has been assigned the medium will be assumed to be vacuum.
///
/// Ea in keV, ma in eV, mgamma in eV, Lcoh in mm, Bmag in T
///
/// mg in eV, absLength in cm-1
///
/// The returned value is given for g_ag = 10^-10 GeV-1
///
Double_t TRestAxionPhotonConversion::AxionAbsorptionProbability(Double_t Bmag, Double_t Lcoh, Double_t Ea,
                                                                Double_t ma, Double_t mg,
                                                                Double_t absLength) {
    mpfr::mpreal axionMass = ma;
    mpfr::mpreal cohLength = Lcoh / 1000.;  // Default REST units are mm;

    mpfr::mpreal photonMass = mg;
    if (mg == 0 && fBufferGas) photonMass = fBufferGas->GetPhotonMass(Ea);

    if (fDebug) {
        RESTDebug << "+--------------------------------------------------------------------------+"
                  << RESTendl;
        RESTDebug << " TRestAxionPhotonConversion::GammaTransmissionProbability. Parameter summary"
                  << RESTendl;
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
}

/// Commented because it uses ComplexReal structure that is moved to TRestAxionFieldPropagationProcess class
/*
void TRestAxionPhotonConversion::PropagateAxion(Double_t Bmag, Double_t Lcoh, Double_t Ea, Double_t ma,
                                                Double_t mg, Double_t absLength) {
    mpfr::mpreal axionMass = ma;
    mpfr::mpreal cohLength = Lcoh / 1000.;  // Default REST units are mm;

    mpfr::mpreal photonMass = mg;
    if (mg == 0 && fBufferGas) photonMass = fBufferGas->GetPhotonMass(Ea);

    if (fDebug) {
        RESTDebug << "+--------------------------------------------------------------------------+" <<
RESTendl; RESTDebug << " TRestAxionPhotonConversion::GammaTransmissionProbability. Parameter summary" <<
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
