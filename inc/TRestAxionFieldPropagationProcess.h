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

#ifndef RestCore_TRestAxionFieldPropagationProcess
#define RestCore_TRestAxionFieldPropagationProcess

#include "TVector3.h"
#include "TVectorD.h"

#include "TRestAxionEvent.h"
#include "TRestAxionMagneticField.h"
#include "TRestAxionPhotonConversion.h"
#include "TRestEventProcess.h"
#include "TRestPhysics.h"
#include "mpreal.h"

/// A structure to define the two components of a complex number using real precision.
/// To be used inside TRestAxionFieldPropagationProcess. Copied from TRestAxionPhotonConversion.h
struct ComplexReal {
    /// The real part of the number
    mpfr::mpreal real = 0;

    /// The imaginary part of the number
    mpfr::mpreal img = 0;
};

//! A process to introduce the axion-photon conversion probability in the signal generation chain
class TRestAxionFieldPropagationProcess : public TRestEventProcess {
   private:
    /// complex axion field amplitude.
    ComplexReal faxionAmplitude;  //!

    /// complex amplitude of the photon field component parallel to the transverse magnetic field.
    ComplexReal fparallelPhotonAmplitude;  //!

    /// complex amplitude of the photon field component orthogonal to the transverse magnetic field.
    ComplexReal forthogonalPhotonAmplitude;  //!

    /////////////////////////////////////////////////////////////////////////
    // ----- Just a quick implementation of complex number operations ---- //
    // ------------ including mpfr real precision arithmetics ------------ //

    /// A function to calculate complex number addition with real precision
    ComplexReal ComplexAddition(const ComplexReal& a, const ComplexReal& b) {
        ComplexReal c;

        c.real = a.real + b.real;
        c.img = a.img + b.img;

        return c;
    }

    /// A function to calculate complex number substraction with real precision
    ComplexReal ComplexSubstraction(const ComplexReal& a, const ComplexReal& b) {
        ComplexReal c;

        c.real = a.real - b.real;
        c.img = a.img - b.img;

        return c;
    }

    /// A function to calculate complex number product with real precision
    ComplexReal ComplexProduct(const ComplexReal& a, const ComplexReal& b) {
        ComplexReal c;

        c.real = a.real * b.real - a.img * b.img;
        c.img = a.real * b.img + a.img * b.real;

        return c;
    }

    /// A function to calculate complex number product by a value with real precision
    ComplexReal ComplexProduct(const mpfr::mpreal& value, const ComplexReal& a) {
        ComplexReal c;

        c.real = value * a.real;
        c.img = value * a.img;

        return c;
    }

    /// A function to calculate complex number cocient with real precision
    ComplexReal ComplexCocient(const ComplexReal& a, const ComplexReal& b) {
        ComplexReal c = ComplexConjugate(b);
        c = ComplexProduct(a, c);

        mpfr::mpreal norm = 1. / Norm2(b);

        c = ComplexProduct(norm, c);

        return c;
    }

    /// A function to calculate complex conjugate with real precision
    ComplexReal ComplexConjugate(const ComplexReal& a) {
        ComplexReal c;

        c.real = a.real;
        c.img = -a.img;

        return c;
    }

    /// A function to calculate the norm squared from a complex number with real precision
    mpfr::mpreal Norm2(const ComplexReal& a) {
        mpfr::mpreal result = a.real * a.real + a.img * a.img;
        return result;
    }

    /// A function to calculate complex number product by a value with real precision
    ComplexReal SetComplexReal(const mpfr::mpreal& r, const mpfr::mpreal& i) {
        ComplexReal c;

        c.real = r;
        c.img = i;

        return c;
    }

    void PrintComplex(ComplexReal);

    void InitFromConfigFile() override;

    void Initialize() override;

    void LoadDefaultConfig();

    /// A pointer to the specific TRestAxionEvent
    TRestAxionEvent* fAxionEvent;  //!

    /// A pointer to the magnetic field stored in TRestRun
    TRestAxionMagneticField* fAxionMagneticField;  //!

    /// A pointer for the gamma conversion probability
    TRestAxionPhotonConversion* fAxionPhotonConversion;  //!

    /// A pointer for the gamma conversion probability
    TRestAxionBufferGas* fAxionBufferGas;  //!

    /// Variables for the output AxionEvent position
    TString fMode;
    TVector3 fFinalNormalPlan;
    TVector3 fFinalPositionPlan;
    Double_t fDistance;

   protected:
   public:
    void InitProcess() override;

    any GetInputEvent() const override { return fAxionEvent; }
    any GetOutputEvent() const override { return fAxionEvent; }

    TRestEvent* ProcessEvent(TRestEvent* evInput) override;

    void LoadConfig(std::string cfgFilename, std::string name = "");

    /// It prints out the process parameters stored in the metadata structure
    void PrintMetadata() override {
        BeginPrintProcess();

        metadata << "mode: " << fMode << endl;

        Double_t x = fFinalPositionPlan.X();
        Double_t y = fFinalPositionPlan.Y();
        Double_t z = fFinalPositionPlan.Z();
        metadata << "finalPositionPlan =  ( " << x << " ," << y << ", " << z << ") mm" << endl;

        x = fFinalNormalPlan.X();
        y = fFinalNormalPlan.Y();
        z = fFinalNormalPlan.Z();
        metadata << "finalNormalPlan = ( " << x << ", " << y << ", " << z << ")" << endl;

        metadata << "Distance: " << fDistance << " mm" << endl;

        EndPrintProcess();
    }

    bool IsInBoundedPlan(TVector3 pos, Int_t i, Int_t p);
    std::vector<TVector3> InOut(std::vector<TVector3> bounds, TVector3 dir);
    std::vector<TVector3> FindBoundariesOneVolume(TVector3 pos, TVector3 dir, Int_t p);
    std::vector<TVector3> FieldBoundary(std::vector<TVector3> boundaries, Double_t minStep);

    /// Returns the boundaries of the axion passed through magnetic fields
    std::vector<std::vector<TVector3>> FindFieldBoundaries(Double_t minStep = -1);

    // TVectorD GetFieldVector(TVector3 in, TVector3 out, Int_t N = 0);

    // Obsolete
    // TVector3 MoveToFinalDistance(TVector3 pos, TVector3 dir, Double_t distance);
    // TVector3 MoveToPlan(TVector3 pos, TVector3 dir, Double_t f, Int_t i);
    // TVector3 FinalPositionInPlan(TVector3 pos, TVector3 dir, TVector3 normalPlan, TVector3 pointPlan);

    // Calculates amplitudes of the axion field, parallel component of the photon field and orthogonal
    // component of the photon field (with respect to the transverse component of the magnetic field)
    // after passing one segment of the particle trajectory along which it is assumed
    // that the transverse component of the magnetic field with respect to the particle propagation direction
    // is not equal zero.
    void CalculateAmplitudesInSegment(ComplexReal& faxionAmplitude, ComplexReal& fparallelPhotonAmplitude,
                                      ComplexReal& forthogonalPhotonAmplitude, TVector3& averageBT,
                                      mpfr::mpreal axionMass, mpfr::mpreal photonMass, mpfr::mpreal Ea,
                                      TVector3 from, TVector3 to, mpfr::mpreal CommonPhase,
                                      mpfr::mpreal OrthogonalPhase);

    // Calculates amplitudes of the axion field, parallel component of the photon field and orthogonal
    // component
    // of the photon field after passing one subsegment of the particle trajectory along which it is assumed
    // that the magnetic field is constant
    void CalculateAmplitudesInSubsegment(ComplexReal& axionAmplitude, ComplexReal& parallelPhotonAmplitude,
                                         ComplexReal& orhogonalPhotonAmplitude, mpfr::mpreal theta,
                                         mpfr::mpreal lambda, Double_t length, mpfr::mpreal CommonPhase,
                                         mpfr::mpreal OrthogonalPhase);

    // Calculates amplitudes of the axion field, parallel component of the photon field and orthogonal
    // component
    // of the photon field after passing one segment of the particle trajectory along which the transverse
    // component of the magnetic field is zero
    void PropagateWithoutBField(ComplexReal& axionAmplitude, ComplexReal& parallelPhotonAmplitude,
                                ComplexReal& orthogonalPhotonAmplitude, mpfr::mpreal axionMass,
                                mpfr::mpreal photonMass, mpfr::mpreal Ea, TVector3 from, TVector3 to);

    /// Returns the name of this process
    const char* GetProcessName() const override { return "axionFieldPropagation"; }

    // Constructor
    TRestAxionFieldPropagationProcess();
    TRestAxionFieldPropagationProcess(char* cfgFileName);

    // Destructor
    ~TRestAxionFieldPropagationProcess();

    ClassDefOverride(TRestAxionFieldPropagationProcess, 1);
};
#endif
