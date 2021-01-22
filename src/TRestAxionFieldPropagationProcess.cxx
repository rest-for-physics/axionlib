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
///
/// This process is designed to work with the most general case of the magnetic field description using the
/// help of TRestAxionMagneticField metadata class. This means that the magnetic ﬁeld volume can be made of
/// several magnetic ﬁeld "regions" (as described in the documentation of TRestAxionMagneticField, where each
/// "region" is deﬁned by a separate `<addMagnetVolume...>` line in its RML conﬁguration ﬁle. One such example
/// of generic magnetic ﬁeld volume description is shown in the following figure, where the volume consists of
/// four regions labeled #1 . . . #4. In the following example the regions appear as connected regions
/// although in general they could be completely unconnected, or isolated regions. Anywhere where no region is
/// defined the field will be equal to zero.
///
/// \htmlonly <style>div.image img[src="AxionMagnetTrajectory.png"]{width:500px;}</style> \endhtmlonly
///
/// ![Schematic of the axion trajectory through the magnetic field volumes.](AxionMagnetTrajectory.png)
///
/// Here, the boundaries of each region are represented by blue lines. Also, in some (or all) regions the
/// magnetic ﬁeld can be zero in its outer parts in order to define a more complex field shape inside.
/// These parts where \f$B = 0\f$ are shown in blue, while the parts of the
/// regions with \f$B \neq 0\f$ are shown in green. In this case, it is possible that the particle, right
/// after entering the region, passes through the part where \f$B = 0\f$ before traversing the section with
/// \f$B \neq 0\f$. Also, just before exiting the region, the particle can again pass through the section
/// with \f$B = 0\f$.
///
/// This process will use the profile of the transversal magnetic field component at each of the volumes
/// along the path to calculate the probability the particle is in a photon state at the end of the
/// trajectory,
///
/// In a first approach this process will be only valid for the axion propagation inside a single magnetic
/// volume, until it is confirmed the process is valid for any number of volumes.
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2019-March:  First prototyping of the process.
///              Javier Galan
///
/// 2019-July:   Implementation of boundaries and magnetic field evaluation.
///              Eve Pachoud
///
/// 2020-March:  Review and validation of this process.
///              Javier Galan and Krešimir Jakovčić
///
///
/// \class      TRestAxionFieldPropagationProcess
/// \author     Javier Galan <javier.galan@unizar.es>
/// \author     Krešimir Jakovčić <kjakov@irb.hr>
///
/// <hr>
///
#include "TRestAxionFieldPropagationProcess.h"
#include <TVectorD.h>
#include "TComplex.h"
#include "TH1F.h"

using namespace std;
using namespace REST_Physics;

ClassImp(TRestAxionFieldPropagationProcess);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionFieldPropagationProcess::TRestAxionFieldPropagationProcess() { Initialize(); }

///////////////////////////////////////////////
/// \brief Constructor loading data from a config file
///
/// If no configuration path is defined using TRestMetadata::SetConfigFilePath
/// the path to the config file must be specified using full path, absolute or relative.
///
/// The default behaviour is that the config file must be specified with
/// full path, absolute or relative.
///
/// \param cfgFileName A const char* giving the path to an RML file.
///
TRestAxionFieldPropagationProcess::TRestAxionFieldPropagationProcess(char* cfgFileName) {
    Initialize();

    LoadConfig(cfgFileName);
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionFieldPropagationProcess::~TRestAxionFieldPropagationProcess() { delete fAxionEvent; }

///////////////////////////////////////////////
/// \brief Function to load the default config in absence of RML input
///
void TRestAxionFieldPropagationProcess::LoadDefaultConfig() {
    SetName(this->ClassName());
    SetTitle("Default config");
}

///////////////////////////////////////////////
/// \brief Function to load the configuration from an external configuration file.
///
/// If no configuration path is defined in TRestMetadata::SetConfigFilePath
/// the path to the config file must be specified using full path, absolute or relative.
///
/// \param cfgFileName A const char* giving the path to an RML file.
/// \param name The name of the specific metadata. It will be used to find the
/// corresponding TRestAxionFieldPropagationProcess section inside the RML.
///
void TRestAxionFieldPropagationProcess::LoadConfig(std::string cfgFilename, std::string name) {
    if (LoadConfigFromFile(cfgFilename, name)) LoadDefaultConfig();
}

///////////////////////////////////////////////
/// \brief Function to initialize input/output event members and define the section name
///
/// It sets the default real precision to be used with mpfr types. Now it is 30 digits.
/// So that we can still calculate numbers such as : 1.0 - 1.e-30
///
void TRestAxionFieldPropagationProcess::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    mpfr::mpreal::set_default_prec(mpfr::digits2bits(30));

    fAxionEvent = new TRestAxionEvent();

    fFinalNormalPlan = TVector3();
    fFinalPositionPlan = TVector3();
    fDistance = 0.0;
}

///////////////////////////////////////////////
/// \brief Process initialization. Data members that require initialization just before start processing
/// should be initialized here.
///
void TRestAxionFieldPropagationProcess::InitProcess() {
    debug << "Entering ... TRestAxionGeneratorProcess::InitProcess" << endl;

    fAxionMagneticField = (TRestAxionMagneticField*)this->GetMetadata("TRestAxionMagneticField");

    if (!fAxionMagneticField) {
        ferr << "TRestAxionFieldPropagationprocess. Magnetic Field was not defined!" << endl;
        exit(0);
    }

    fAxionBufferGas = (TRestAxionBufferGas*)this->GetMetadata("TRestAxionBufferGas");

    if (!fAxionBufferGas) {
        ferr << "TRestAxionBufferGas. Cannot access the buffer gas" << endl;
        exit(0);
    }

    fAxionPhotonConversion = new TRestAxionPhotonConversion();
    fAxionPhotonConversion->SetBufferGas(fAxionBufferGas);
}

///////////////////////////////////////////////
/// \brief This method will translate the vector with direction `dir` starting at position `pos` to the plane
/// defined by the normal vector plane, `n` that contains the point `a` in the plane.
///
/// This method has been migrated to TRestPhysics, and it will be accesible through REST_Physics
/*
TVector3 TRestAxionFieldPropagationProcess::MoveToPlane(TVector3 pos, TVector3 dir, TVector3 n, TVector3 a) {
    if (n * dir == 0) {
        ferr << "The vector is parallel to the plane!!" << endl;
        ferr << "Position will not be translated" << endl;
    } else {
        Double_t t = (n * a - n * pos) / (n * dir);

        return pos + t * dir;
    }
    return pos;
}*/

///////////////////////////////////////////////
/// \brief This method is OBSOLETE.
///
/// It is more intuitive to define the vector position and direction, and
/// plane vector and point. Then use this information to find the intersection. Instead of defining a
/// component and an impact factor.
/*
TVector3 TRestAxionFieldPropagationProcess::MoveToPlan(TVector3 pos, TVector3 dir, Double_t f, Int_t i) {
    if (dir[i] == 0) ferr << "The component of the direction you chose is equal to 0 " << endl;

    Double_t t = (f - pos[i]) / dir[i];
    pos[i] = f;
    pos[(i + 1) % 3] = pos[(i + 1) % 3] + t * dir[(i + 1) % 3];
    pos[(i + 2) % 3] = pos[(i + 2) % 3] + t * dir[(i + 2) % 3];

    return pos;
}
*/

///////////////////////////////////////////////
/// \brief This method is OBSOLETE.
///
/// Re-implementation of MoveToPlane
/*
TVector3 TRestAxionFieldPropagationProcess::FinalPositionInPlan(TVector3 pos, TVector3 dir,
                                                                TVector3 normalPlan, TVector3 pointPlan) {
    if (normalPlan.Dot(dir) == 0) return pos;

    Double_t t, d;
    d = -normalPlan.Dot(pointPlan);
    t = -(normalPlan.Dot(pos) + d) / (normalPlan.Dot(dir));

    pos = pos + t * dir;

    return pos;
}
*/

///////////////////////////////////////////////
/// \brief This method is OBSOLETE
///
/// Re-implementation in MoveByDistance. Just because the method name is more intuitive using BYDISTANCE
/*
TVector3 TRestAxionFieldPropagationProcess::MoveToFinalDistance(TVector3 pos, TVector3 dir,
                                                                Double_t distance) {
    Double_t t = distance / dir.Mag();
    pos = pos + t * dir;

    return pos;
}
*/

///////////////////////////////////////////////
/// \brief This method transports a position `pos` by a distance `d` in the direction defined by `dir`.
///
/// This method has been migrated to TRestPhysics, and it will be accesible through REST_Physics
/*
TVector3 TRestAxionFieldPropagationProcess::MoveByDistance(TVector3 pos, TVector3 dir, Double_t d) {
    return pos + d * dir.Unit();
} */

bool TRestAxionFieldPropagationProcess::IsInBoundedPlan(TVector3 pos, Int_t i, Int_t p) {
    Double_t minCond1, maxCond1, minCond2, maxCond2;
    Int_t j, k;

    /*
if (i == 0) {
    minCond1 = (fAxionMagneticField->GetYmin())[p];
    maxCond1 = (fAxionMagneticField->GetYmax())[p];
    minCond2 = (fAxionMagneticField->GetZmin())[p];
    maxCond2 = (fAxionMagneticField->GetZmax())[p];
    j = 1;
    k = 2;
}
if (i == 1) {
    minCond1 = (fAxionMagneticField->GetXmin())[p];
    maxCond1 = (fAxionMagneticField->GetXmax())[p];
    minCond2 = (fAxionMagneticField->GetZmin())[p];
    maxCond2 = (fAxionMagneticField->GetZmax())[p];
    j = 0;
    k = 2;
}
if (i == 2) {
    minCond1 = (fAxionMagneticField->GetXmin())[p];
    maxCond1 = (fAxionMagneticField->GetXmax())[p];
    minCond2 = (fAxionMagneticField->GetYmin())[p];
    maxCond2 = (fAxionMagneticField->GetYmax())[p];
    j = 0;
    k = 1;
}

    bool cond1 = (pos[j] <= maxCond1 && pos[j] >= minCond1);
    bool cond2 = (pos[k] <= maxCond2 && pos[k] >= minCond2);

    return cond1 && cond2;
    */
    return false;
}

std::vector<TVector3> TRestAxionFieldPropagationProcess::InOut(std::vector<TVector3> bounds, TVector3 dir) {
    Int_t i = 0;

    TVector3 in = bounds[0];
    TVector3 out = bounds[1];

    while (i < 3 && dir[i] * (out[i] - in[i]) >= 0) i = i + 1;

    if (i < 3) {
        bounds[0] = out;
        bounds[1] = in;
    }

    return bounds;
}

///////////////////////////////////////////////
/// \brief Finds the in/out particle trajectory boundaries for a particular magnetic volume.
///
/// This method checks if the particle (with the initial position `pos` and direction `dir`) passes through
/// the magnetic ﬁeld region speciﬁed by the input parameter p. It is done by searching for the points where
/// the particle trajectory intersects the boundary planes of that region. If two such points (entry point and
/// exit point) are found, their coordinates are stored in the vector boundaries. In the example shown in Fig.
/// 1 these points are: IN 1 and OUT 1 for the region #1 and IN2  and OUT 2 for the region #2.
/* This method has been moved to TRestAxionMagneticField::GetVolumeBoundaries
std::vector<TVector3> TRestAxionFieldPropagationProcess::FindBoundariesOneVolume(TVector3 pos, TVector3 dir,
                                                                                 Int_t p) {
    std::vector<TVector3> boundaries;
TVector3 in;
TVector3 out;

Int_t i = 0;
Int_t j = 0;

TVectorD f(6);
f[0] = (fAxionMagneticField->GetXmin())[p];
f[1] = (fAxionMagneticField->GetXmax())[p];
f[2] = (fAxionMagneticField->GetYmin())[p];
f[3] = (fAxionMagneticField->GetYmax())[p];
f[4] = (fAxionMagneticField->GetZmin())[p];
f[5] = (fAxionMagneticField->GetZmax())[p];

Int_t nFace = 0;

while (nFace < 6 && i <= 0) {
    if (dir[Int_t(nFace / 2)] != 0) {
        if ((pos[Int_t(nFace / 2)] - f[nFace]) * dir[Int_t(nFace / 2)] <= 0) {
            in = MoveToPlan(pos, dir, f[nFace], Int_t(nFace / 2));
            if (IsInBoundedPlan(in, Int_t(nFace / 2), p)) i = i + 1;
        }
    }
    nFace = nFace + 1;
}

if (i == 1) {
    while (nFace < 6 && j <= 0) {
        if (dir[Int_t(nFace / 2)] != 0) {
            if ((pos[Int_t(nFace / 2)] - f[nFace]) * dir[Int_t(nFace / 2)] <= 0) {
                out = MoveToPlan(pos, dir, f[nFace], Int_t(nFace / 2));
                if (IsInBoundedPlan(out, Int_t(nFace / 2), p)) j = j + 1;
            }
        }
        nFace = nFace + 1;
    }
}

if (i + j == 2 && in != out) {
    boundaries.push_back(in);
    boundaries.push_back(out);
    return boundaries;
}

else
    return boundaries;
    return boundaries;
}
*/

/*  This  method has been moved  to TRestAxionMagneticField::GetFieldBoundaries
std::vector<TVector3> TRestAxionFieldPropagationProcess::FieldBoundary(std::vector<TVector3> boundaries,
                                                                       Double_t minStep) {
    std::vector<TVector3> boundariesField;

    TVector3 in = boundaries[0];
    TVector3 out = boundaries[1];
    TVector3 diff = out - in;

    Int_t N = 10;

    Int_t i = 0;

    while (1.0 / Double_t(N) > minStep) {
        while ((fAxionMagneticField->GetMagneticField(in[0], in[1], in[2])) == TVector3(0, 0, 0) && i < N) {
            in = in + 1.0 / Double_t(N) * diff;
            i = i + 1;
        }

        if (i == N) return boundariesField;

        in = in - 1.0 / Double_t(N) * diff;
        N = 10 * N;
        i = 0;
    }

    in = in + 10.0 / Double_t(N) * diff;
    boundariesField.push_back(in);

    N = 10;
    i = 0;

    while (1.0 / Double_t(N) > minStep) {
        while ((fAxionMagneticField->GetMagneticField(out[0], out[1], out[2]) == TVector3(0, 0, 0)) &&
               i < N) {
            out = out - 1.0 / Double_t(N) * diff;
            i = i + 1;
        }

        out = out + 1.0 / Double_t(N) * diff;
        N = 10 * N;
        i = 0;
    }

    out = out - 10.0 / Double_t(N) * diff;
    boundariesField.push_back(out);

    return boundariesField;
}
*/

// This method is obsolete because it uses methods FindBoundariesOneVolume and FieldBoundary that are obsolete
// and replaced by the methods TRestAxionMagneticField::GetVolumeBoundaries and
// TRestAxionMagneticField::GetFieldBoundaries.
// The NEW version of this method is given below and uses TRestAxionMagneticField::GetFieldBoundaries
/*
std::vector<std::vector<TVector3>> TRestAxionFieldPropagationProcess::FindFieldBoundaries(Double_t minStep) {
    std::vector<std::vector<TVector3>> boundaryCollection;
    std::vector<std::vector<TVector3>> boundaryFinalCollection;

if (minStep == -1) minStep = 0.01;
std::vector<TVector3> buffVect;
TVector3 boundaryIn;
TVector3 boundaryOut;

TVector3 posInitial = *(fInputAxionEvent->GetPosition());
TVector3 direction = *(fInputAxionEvent->GetDirection());
direction = direction.Unit();

if (direction == TVector3(0, 0, 0))  // No moves
    return boundaryCollection;

if ((fAxionMagneticField->GetXmin()).size() == 0)
    ferr << " The magnetic field has not been loaded " << endl;

// Find global boundaries volume

Int_t N = (fAxionMagneticField->GetXmin()).size();

for (Int_t p = 0; p < N; p++) {
    buffVect = FindBoundariesOneVolume(posInitial, direction, p);
    if (buffVect.size() == 2) {
        buffVect = InOut(buffVect, direction);
        boundaryCollection.push_back(buffVect);
        buffVect.clear();
    }
}

debug << "+------------------------+" << endl;
debug << " Number of volume boundaries : " << 2 * boundaryCollection.size() << endl;
debug << "+------------------------+" << endl;

debug << "+------------------------+" << endl;
for (Int_t p = 0; p < boundaryCollection.size(); p++) {
    debug << "for" << p << " in : (" << boundaryCollection[p][0].X() << ","
          << boundaryCollection[p][0].Y() << "," << boundaryCollection[p][0].Z() << ")" << endl;
    debug << "for" << p << " out : (" << boundaryCollection[p][1].X() << ","
          << boundaryCollection[p][1].Y() << "," << boundaryCollection[p][1].Z() << ")" << endl;
}
debug << "+------------------------+" << endl;

// Find precise boundaries field

for (Int_t i = 0; i < boundaryCollection.size(); i++) {
    buffVect = FieldBoundary(boundaryCollection[i], minStep);
    if (buffVect.size() == 2) {
        boundaryFinalCollection.push_back(buffVect);
        buffVect.clear();
    }
}

boundaryCollection.clear();

debug << "+------------------------+" << endl;
debug << " Number of field boundaries : " << 2 * boundaryFinalCollection.size() << endl;
debug << "+------------------------+" << endl;

debug << "+------------------------+" << endl;
for (Int_t p = 0; p < boundaryFinalCollection.size(); p++) {
    debug << "for" << p << " in : (" << boundaryFinalCollection[p][0].X() << ","
          << boundaryFinalCollection[p][0].Y() << "," << boundaryFinalCollection[p][0].Z() << ")" << endl;
    debug << "for" << p << " out : (" << boundaryFinalCollection[p][1].X() << ","
          << boundaryFinalCollection[p][1].Y() << "," << boundaryFinalCollection[p][1].Z() << ")" << endl;
}
debug << "+------------------------+" << endl;

    return boundaryFinalCollection;
}
*/

// This is a NEW version of FindFieldBoundaries method and replaces the obsolete version of this method given
// above
///////////////////////////////////////////////
/// \brief Finds the boundary points for each segment along the particle trajectory where the **transversal**
/// component
/// of the magnetic field in not zero.
///
/// This method goes over all magnetic field regions and checks for each region if the particle (with the
/// initial
/// position `posInitial` and direction `direction`) passes through this region. If yes, it determines the
/// boundary
/// points of the trajectory segment through this region between which the **transversal** component of the
/// magnetic
/// field is not zero. The method returns the vector where each element contains a pair of boundary points,
/// i.e., the starting and ending points of one segment of the particle trajectory along which the
/// **transversal**
/// component of the magnetic field is not zero.

std::vector<std::vector<TVector3>> TRestAxionFieldPropagationProcess::FindFieldBoundaries(Double_t minStep) {
    std::vector<std::vector<TVector3>> boundaryFinalCollection;

    if (minStep == -1) minStep = 0.01;
    std::vector<TVector3> buffVect;

    TVector3 posInitial = *(fAxionEvent->GetPosition());
    TVector3 direction = *(fAxionEvent->GetDirection());
    direction = direction.Unit();

    if (direction == TVector3(0, 0, 0))  // No moves
        return boundaryFinalCollection;

    if (fAxionMagneticField->GetNumberOfVolumes() == 0)
        ferr << " The magnetic field has not been loaded " << endl;

    // Find field boundaries volume

    Int_t N = fAxionMagneticField->GetNumberOfVolumes();

    for (Int_t p = 0; p < N; p++) {
        buffVect.clear();
        buffVect = fAxionMagneticField->GetFieldBoundaries(p, posInitial, direction);
        if (buffVect.size() == 2) {
            boundaryFinalCollection.push_back(buffVect);
        }
    }

    debug << "+------------------------+" << endl;
    debug << " Number of field boundaries : " << 2 * boundaryFinalCollection.size() << endl;
    debug << "+------------------------+" << endl;

    debug << "+------------------------+" << endl;
    for (Int_t p = 0; p < boundaryFinalCollection.size(); p++) {
        debug << "for volume " << p << " in : (" << boundaryFinalCollection[p][0].X() << ","
              << boundaryFinalCollection[p][0].Y() << "," << boundaryFinalCollection[p][0].Z() << ")" << endl;
        debug << "for volume " << p << " out : (" << boundaryFinalCollection[p][1].X() << ","
              << boundaryFinalCollection[p][1].Y() << "," << boundaryFinalCollection[p][1].Z() << ")" << endl;
    }
    debug << "+------------------------+" << endl;

    return boundaryFinalCollection;
}

/* This method is now OBSOLETE. It seems to me there is a problem here, there are two directional vectors
defined.
The one defined by in/out coordinates, and the one defined by axion direction.

This method will be substituted by
std::vector<Double_t> TRestAxionMagneticField::GetTransversalComponentAlongPath(TVector3 from, TVector3 to,
Double_t dl, Int_t Nmax );

TVectorD TRestAxionFieldPropagationProcess::GetFieldVector(TVector3 in, TVector3 out, Int_t N) {
    if (N == 0) N = TMath::Power(10, 4);

    TVectorD Bt(N);
    TVector3 B;
    TVector3 direction = *(fInputAxionEvent->GetDirection());

    TVector3 differential = out - in;

    B = fAxionMagneticField->GetMagneticField(in[0], in[1], in[2]);
    Bt[0] = abs(B.Perp(direction));

    for (Int_t i = 1; i < N; i++) {
        in = in + differential * (1.0 / Double_t(N - 1));
        B = fAxionMagneticField->GetMagneticField(in[0], in[1], in[2]);
        Bt[i] = abs(B.Perp(direction));
    }

    return Bt;
} */

/// \brief Prints variables of the ComplexReal type, i.e, complex numbers
///
void TRestAxionFieldPropagationProcess::PrintComplex(ComplexReal p) {
    debug << p.real << " + " << p.img << "i" << endl;
}

/// \brief Calculates amplitudes of the axion field, parallel component of the photon field and orthogonal
/// component of the photon field (parallel and orthogonal are with respect to the transverse component of the
/// magnetic field)
/// after passing one segment of the particle trajectory along which it is assumed
/// that the transverse component of the magnetic field with respect to the particle propagation direction
/// is not equal zero. The segment is divided into a series of subsegments. It is assumed that the magnitude
/// and direction of the
/// transverse component of the magnetic field are constant along each subsegment, but they can change from
/// one subsegment
/// to another. The calculations are based on the procedure described in the Section 4 of the internal report
/// "Axion-photon conversion in the external magnetic fields" written by B. Lakic and K. Jakovcic.
/// For each subsegment, the values of certain parameters (e.g. theta, lambda) are calculated first and then
/// the method
/// `CalculateAmplitudesInSubsegment` is called to calculate the amplitudes after passing that subsegment,
/// i.e., to calculate the amplitudes
/// at the end of that subsegment. This procedure is repeated for each subsegment until
/// the end of segment is reached. The values of the amplitudes at the end of one subsegment are used to
/// calculate the initial
/// values of these amplitudes for the next subsegment by using equations (4.4)-(4.6).
///
/// NOTE: The amplitudes are calculated for the axion-photon coupling constant g_agg = 10^-10 GeV-1
/// The length of the subsegment is defined by variable `step`. It is currently set to 200 mm.

void TRestAxionFieldPropagationProcess::CalculateAmplitudesInSegment(
    ComplexReal& faxionAmplitude, ComplexReal& fparallelPhotonAmplitude,
    ComplexReal& forthogonalPhotonAmplitude, TVector3& averageBT, mpfr::mpreal axionMass,
    mpfr::mpreal photonMass, mpfr::mpreal Ea, TVector3 from, TVector3 to, mpfr::mpreal CommonPhase,
    mpfr::mpreal OrthogonalPhase) {
    mpfr::mpreal g_agg = 1.0e-10;                 // axion-photon coupling constant in GeV-1
    mpfr::mpreal TeslaineV = 195.3;               // conversion factor from Tesla to eV^2
    Double_t segment_length = (to - from).Mag();  // the length of the entire segment
    TVector3 subsegment_start,
        subsegment_end;          // coordinates of the starting and ending points of one subsegment
    Double_t subsegment_length;  // the length of one subsegment
    ComplexReal subsegment_A0_par, subsegment_A0_ort;  // values of the parallel and orthogonal photon
                                                       // amplitudes at the end of one subsegment (that is
                                                       // also the beginning of the next subsegment
    TVector3
        averageBT_0;  // transverse component of the average magnetic field vector in the previous subsegment

    Double_t step = 200.0;  // the length of the subsegment
    Double_t BTmag;  // average magnitude of the transverse component of the magnetic field in one subsegment
                     // (given in Tesla)
    Double_t BTangle;  // angle between the transverse component of the average magnetic field in one
                       // subsegment and the one in the previous subsegment

    subsegment_start = from;
    subsegment_A0_par = fparallelPhotonAmplitude;    // parallel photon amplitude at the beginning of the
                                                     // segment, i.e., at the point `from`
    subsegment_A0_ort = forthogonalPhotonAmplitude;  // orthogonal photon amplitude at the beginning of the
                                                     // segment, i.e., at the point `from`
    averageBT_0 = averageBT;
    TVector3 dir = (to - from).Unit();  // direction of the particle propagation

    while ((subsegment_start - from).Mag() < segment_length) {  // loop over subsegments in one segment
        subsegment_end = subsegment_start + step * dir;
        if ((subsegment_end - from).Mag() >= segment_length) subsegment_end = to;
        subsegment_length = (subsegment_end - subsegment_start).Mag();
        subsegment_length = subsegment_length / 1000.0;  // default REST units are mm

        // calculation of the average magnitude of the transverse magnetic field along the subsegment, i.e.,
        // between coordinates `subsegment_start` and `subsegment_end`
        BTmag =
            fAxionMagneticField->GetTransversalFieldAverage(subsegment_start, subsegment_end);  // in Tesla

        // calculation of the transverse component of the average magnetic field vector along the subsegment,
        // i.e., between coordinates `subsegment_start` and `subsegment_end`
        averageBT = fAxionMagneticField->GetFieldAverageTransverseVector(subsegment_start, subsegment_end);

        // calculation of the angle between the transverse component of the average magnetic field in one
        // subsegment and the one in the previous subsegment
        if ((averageBT_0.Mag() == 0.0) || (averageBT.Mag() == 0.0))
            BTangle = 0.0;
        else
            BTangle = averageBT.Angle(averageBT_0);
        if (averageBT.Mag() != 0.0) averageBT_0 = averageBT;

        // calculation of the initial values of the parallel and orthogonal photon amplitudes at the beginning
        // of the subsegment (at the point `subsegment_start`)
        // from the values at the end of the previous subsegment (also at the point `subsegment_start`)
        fparallelPhotonAmplitude = ComplexAddition(ComplexProduct(cos(BTangle), subsegment_A0_par),
                                                   ComplexProduct(sin(BTangle), subsegment_A0_ort));
        forthogonalPhotonAmplitude = ComplexAddition(ComplexProduct(-sin(BTangle), subsegment_A0_par),
                                                     ComplexProduct(cos(BTangle), subsegment_A0_ort));

        // calculation of the parameters theta and lambda for the subsegment
        mpfr::mpreal term_1 = 2 * (Ea * 1000.0) * (g_agg * 1.0e-9) * (BTmag * TeslaineV);  // in eV^2
        mpfr::mpreal term_2 = axionMass * axionMass - photonMass * photonMass;             // in ev^2
        mpfr::mpreal theta = 0.5 * atan(term_1 / term_2);
        mpfr::mpreal lambda = sqrt(term_1 * term_1 + term_2 * term_2) / (4. * Ea * 1000.0);  // in eV

        debug << "+--------------------------------------------------------------------------+" << endl;
        debug << " CalculateAmplitudesInSegment method: Parameter summary" << endl;
        debug << endl << "segment length = " << segment_length << " mm" << endl;
        debug << endl
              << "subsegment_start: (" << subsegment_start.x() << ", " << subsegment_start.y() << ", "
              << subsegment_start.z() << ") "
              << " mm" << endl;
        debug << "subsegment_end: ( " << subsegment_end.x() << ", " << subsegment_end.y() << ", "
              << subsegment_end.z() << ") "
              << " mm" << endl;
        debug << "subsegment length = " << subsegment_length << " m" << endl << endl;

        debug << " average magnitude of the transverse component of the magnetic field in the subsegment : "
              << BTmag << " T" << endl;
        debug << " angle of the transverse component of the average magnetic field in the subsegment with "
                 "respect to the previous subsegment : "
              << BTangle << " rad" << endl;
        debug << " g_agg : " << g_agg << " GeV-1" << endl;
        debug << " Theta : " << theta << endl;
        debug << " lambda : " << lambda << " eV" << endl;
        debug << " subsegment_A0_par : ";
        PrintComplex(subsegment_A0_par);
        debug << " subsegment_A0_ort : ";
        PrintComplex(subsegment_A0_ort);
        debug << " BEFORE calculating in subsegment: paralell photon component amplitude : ";
        PrintComplex(fparallelPhotonAmplitude);
        debug << " BEFORE calculating in subsegment: orthogonal photon component amplitude : ";
        PrintComplex(forthogonalPhotonAmplitude);
        debug << "+--------------------------------------------------------------------------+" << endl;

        CalculateAmplitudesInSubsegment(faxionAmplitude, fparallelPhotonAmplitude, forthogonalPhotonAmplitude,
                                        theta, lambda, subsegment_length, CommonPhase, OrthogonalPhase);
        subsegment_A0_par = fparallelPhotonAmplitude;    // parallel photon amplitude at the end of the
                                                         // subsegment, i.e. at the point `subsegment_end`
        subsegment_A0_ort = forthogonalPhotonAmplitude;  // orthogonal photon amplitude at the end of the
                                                         // subsegment, i.e. at the point `subsegment_end`
        debug << endl << " AFTER calculating in subsegment: subsegment_A0_par : ";
        PrintComplex(subsegment_A0_par);
        debug << " AFTER calculating in subsegment: subsegment_A0_ort : ";
        PrintComplex(subsegment_A0_ort);
        debug << " AFTER calculating in subsegment: paralell photon component amplitude : ";
        PrintComplex(fparallelPhotonAmplitude);
        debug << " AFTER calculating in subsegment: orthogonal photon component amplitude : ";
        PrintComplex(forthogonalPhotonAmplitude);
        subsegment_start = subsegment_end;
    }
}

/// \brief Calculates amplitudes of the axion field, parallel component of the photon field and orthogonal
/// component
/// of the photon field after passing one subsegment of the particle trajectory along which it is assumed
/// that the transverse component of the magnetic field is constant. It uses equations (3.15) and (3.38) given
/// in the internal report "Axion-photon conversion in the external magnetic fields" written by B. Lakic and
/// K. Jakovcic.
/// Also, amplitudes are of ComplexReal type, which stores complex numbers based on mpreal wrapper to allow
/// precise
/// calculation of small values.  At present, the  precision is set to 30 digits, so that we can still
/// calculate
/// numbers such as : 1.0 - 1.e-30
/// NOTE: The amplitudes are calculated for the axion-photon coupling constant g_agg = 10^-10 GeV-1
void TRestAxionFieldPropagationProcess::CalculateAmplitudesInSubsegment(
    ComplexReal& faxionAmplitude, ComplexReal& fparallelPhotonAmplitude,
    ComplexReal& forthogonalPhotonAmplitude, mpfr::mpreal theta, mpfr::mpreal lambda, Double_t length,
    mpfr::mpreal CommonPhase, mpfr::mpreal OrthogonalPhase) {
    // lambda is given in eV, theta has no dimension, length is in m, Common phase and Orthogonal phase are in
    // eV

    mpfr::mpreal::set_default_prec(mpfr::digits2bits(30));
    cout.precision(30);
    // setting initial parameters
    ComplexReal a0 =
        faxionAmplitude;  // axion field amplitude initial value at the beginning of the subsegment
    ComplexReal A0_par = fparallelPhotonAmplitude;  // photon field parallel component amplitude initial value
                                                    // at the beginning of the subsegment
    ComplexReal A0_ort = forthogonalPhotonAmplitude;  // photon field orthogonal component amplitude initial
                                                      // value at the beginning of the subsegment
    debug << "+--------------------------------------------------------------------------+" << endl;
    debug << " CalculateAmplitudesInSubsegment method: Parameter summary" << endl;
    debug << " Theta : " << theta << endl;
    debug << " lambda : " << lambda << " eV" << endl;
    debug << " subsegment length : " << length << " m" << endl;
    debug << " Common phase : " << CommonPhase << " eV" << endl;
    debug << " orthogonal photon component phase : " << OrthogonalPhase << " eV" << endl;
    debug << " a0 : ";
    PrintComplex(a0);
    debug << " A0_par : ";
    PrintComplex(A0_par);
    debug << " A0_ort : ";
    PrintComplex(A0_ort);
    debug << "+--------------------------------------------------------------------------+" << endl;

    // setting auxillary parameters used in calculations
    mpfr::mpreal cos2_theta = cos(theta) * cos(theta);
    mpfr::mpreal sin2_theta = sin(theta) * sin(theta);
    mpfr::mpreal sin_2theta = sin(2.0 * theta);

    mpfr::mpreal l = length * PhMeterIneV;  // length in eV-1

    mpfr::mpreal phi = lambda * l;
    ComplexReal exp_lambdaPlusZ = SetComplexReal(cos(-phi), sin(-phi));
    ComplexReal exp_lambdaMinusZ = SetComplexReal(cos(phi), sin(phi));

    mpfr::mpreal phi1 = CommonPhase * l;
    ComplexReal exp_CommonPhase = SetComplexReal(cos(phi1), sin(phi1));

    debug << "+--------------------------------------------------------------------------+" << endl;
    debug << " Intermediate calculations" << endl;
    debug << " cos^(2)_theta : " << cos2_theta << endl;
    debug << " sin^(2)_theta : " << sin2_theta << endl;
    debug << " sin(2*theta) : " << sin_2theta << endl;
    debug << " subsegment length : " << length << " m" << endl;
    debug << " l : " << l << " eV-1" << endl;
    debug << " phi : " << phi << endl;
    debug << " exp_lambdaPlusZ : ";
    PrintComplex(exp_lambdaPlusZ);
    debug << " exp_lambdaMinusZ : ";
    PrintComplex(exp_lambdaMinusZ);
    debug << " phi1 : " << phi1 << endl;
    debug << " exp_CommonPhase : ";
    PrintComplex(exp_CommonPhase);
    debug << "+--------------------------------------------------------------------------+" << endl;

    // calculation of the photon parallel component amplitude
    ComplexReal A_term_1_1 = ComplexProduct(cos2_theta, exp_lambdaPlusZ);
    ComplexReal A_term_1_2 = ComplexProduct(sin2_theta, exp_lambdaMinusZ);
    ComplexReal A_sum_1 = ComplexAddition(A_term_1_1, A_term_1_2);
    ComplexReal A_term_1 = ComplexProduct(A0_par, A_sum_1);
    ComplexReal A_sum_2 = ComplexSubstraction(exp_lambdaPlusZ, exp_lambdaMinusZ);
    ComplexReal temp = ComplexProduct(sin_2theta * 0.5, a0);
    ComplexReal A_term_2 = ComplexProduct(temp, A_sum_2);
    ComplexReal A_sum = ComplexAddition(A_term_1, A_term_2);
    fparallelPhotonAmplitude = ComplexProduct(exp_CommonPhase, A_sum);
    debug << "+--------------------------------------------------------------------------+" << endl;
    debug << " Intermediate calculations for the photon parallel component amplitude" << endl;
    debug << " A_term_1_1 : ";
    PrintComplex(A_term_1_1);
    debug << " A_term_1_2 : ";
    PrintComplex(A_term_1_2);
    debug << " A_sum_1 : ";
    PrintComplex(A_sum_1);
    debug << " A_term_1 : ";
    PrintComplex(A_term_1);
    debug << " A_sum_2 : ";
    PrintComplex(A_sum_2);
    debug << " A_term_2 : ";
    PrintComplex(A_term_2);
    debug << " A_sum : ";
    PrintComplex(A_sum);
    debug << " parallelPhotonAmplitude : ";
    PrintComplex(fparallelPhotonAmplitude);
    debug << "+--------------------------------------------------------------------------+" << endl;

    // calculation of the axion amplitude
    ComplexReal a_sum_1 = A_sum_2;
    temp = ComplexProduct(sin_2theta * 0.5, A0_par);
    ComplexReal a_term_1 = ComplexProduct(temp, a_sum_1);
    ComplexReal a_term_2_1 = ComplexProduct(sin2_theta, exp_lambdaPlusZ);
    ComplexReal a_term_2_2 = ComplexProduct(cos2_theta, exp_lambdaMinusZ);
    ComplexReal a_sum_2 = ComplexAddition(a_term_2_1, a_term_2_2);
    ComplexReal a_term_2 = ComplexProduct(a0, a_sum_2);
    ComplexReal a_sum = ComplexAddition(a_term_1, a_term_2);
    faxionAmplitude = ComplexProduct(exp_CommonPhase, a_sum);
    debug << "+--------------------------------------------------------------------------+" << endl;
    debug << " Intermediate calculations for the axion amplitude" << endl;
    debug << " a_sum_1 : ";
    PrintComplex(a_sum_1);
    debug << " a_term_1 : ";
    PrintComplex(a_term_1);
    debug << " a_term_2_1 : ";
    PrintComplex(a_term_2_1);
    debug << " a_term_2_2 : ";
    PrintComplex(a_term_2_2);
    debug << " a_sum_2 : ";
    PrintComplex(a_sum_2);
    debug << " a_term_2 : ";
    PrintComplex(a_term_2);
    debug << " a_sum : ";
    PrintComplex(a_sum);
    debug << " axionAmplitude : ";
    PrintComplex(faxionAmplitude);
    debug << "+--------------------------------------------------------------------------+" << endl;

    // calculation of the photon orthogonal component amplitude
    mpfr::mpreal phi2 = OrthogonalPhase * l;
    ComplexReal exp_OrthogonalPhase = SetComplexReal(cos(phi2), sin(phi2));
    forthogonalPhotonAmplitude = ComplexProduct(A0_ort, exp_OrthogonalPhase);
    debug << "+--------------------------------------------------------------------------+" << endl;
    debug << " Intermediate calculations for the photon orthogonal component amplitude" << endl;
    debug << " phi2 : " << phi2 << endl;
    debug << " exp_OrthogonalPhase : ";
    PrintComplex(exp_OrthogonalPhase);
    debug << " orthogonalPhotonAmplitude : ";
    PrintComplex(forthogonalPhotonAmplitude);
    debug << "+--------------------------------------------------------------------------+" << endl;
}

/// \brief Calculates amplitudes of the axion field, parallel component of the photon field and orthogonal
/// component
/// of the photon field after passing one segment of the particle trajectory along which the transverse
/// component
/// of the magnetic field is ZERO. It uses equation (3.15) given in the internal report "Axion-photon
/// conversion
/// in the external magnetic fields" written by B. Lakic and K. Jakovcic.
/// Also, amplitudes are of ComplexReal type, which stores complex numbers based on mpreal wrapper to allow
/// precise
/// calculation of small values.  At present, the  precision is set to 30 digits, so that we can still
/// calculate
/// numbers such as : 1.0 - 1.e-30
void TRestAxionFieldPropagationProcess::PropagateWithoutBField(ComplexReal& faxionAmplitude,
                                                               ComplexReal& fparallelPhotonAmplitude,
                                                               ComplexReal& forthogonalPhotonAmplitude,
                                                               mpfr::mpreal axionMass,
                                                               mpfr::mpreal photonMass, mpfr::mpreal Ea,
                                                               TVector3 from, TVector3 to) {
    mpfr::mpreal::set_default_prec(mpfr::digits2bits(30));
    cout.precision(30);
    mpfr::mpreal axionPhase = Ea * 1000.0 - (axionMass * axionMass) / (2. * Ea * 1000.0);     // in eV
    mpfr::mpreal photonPhase = Ea * 1000.0 - (photonMass * photonMass) / (2. * Ea * 1000.0);  // in eV
    Double_t length = (to - from).Mag();
    length = length / 1000.0;               // default REST units are mm
    mpfr::mpreal l = length * PhMeterIneV;  // length in eV-1

    debug << "+--------------------------------------------------------------------------+" << endl;
    debug << " Propagation without B field: " << endl;
    debug << " INITIAL VALUES " << endl;
    debug << " axionAmplitude : ";
    PrintComplex(faxionAmplitude);
    debug << " parallelPhotonAmplitude : ";
    PrintComplex(fparallelPhotonAmplitude);
    debug << " orthogonalPhotonAmplitude : ";
    PrintComplex(forthogonalPhotonAmplitude);
    debug << "+--------------------------------------------------------------------------+" << endl;

    mpfr::mpreal axionphi = axionPhase * l;
    ComplexReal exp_axionPhase = SetComplexReal(cos(axionphi), sin(axionphi));
    faxionAmplitude = ComplexProduct(faxionAmplitude, exp_axionPhase);

    mpfr::mpreal photonphi = photonPhase * l;
    ComplexReal exp_photonPhase = SetComplexReal(cos(photonphi), sin(photonphi));
    fparallelPhotonAmplitude = ComplexProduct(fparallelPhotonAmplitude, exp_photonPhase);
    forthogonalPhotonAmplitude = ComplexProduct(forthogonalPhotonAmplitude, exp_photonPhase);
    debug << "+--------------------------------------------------------------------------+" << endl;
    debug << " Intermediate calculations" << endl;
    debug << " axionPhase : " << axionPhase << endl;
    debug << " axionphi : " << axionphi << endl;
    debug << " exp_axionPhase : ";
    PrintComplex(exp_axionPhase);
    debug << "Norm2(exp_axionPhase) = " << Norm2(exp_axionPhase) << endl;
    debug << " photonPhase : " << photonPhase << endl;
    debug << " photonphi : " << photonphi << endl;
    debug << " exp_photonPhase : ";
    PrintComplex(exp_photonPhase);
    debug << "Norm2(exp_photonPhase) = " << Norm2(exp_photonPhase) << endl;
    debug << "+--------------------------------------------------------------------------+" << endl;
    debug << "+--------------------------------------------------------------------------+" << endl;
    debug << " Propagation without B field: " << endl;
    debug << " FINAL VALUES " << endl;
    debug << " axionAmplitude : ";
    PrintComplex(faxionAmplitude);
    debug << " parallelPhotonAmplitude : ";
    PrintComplex(fparallelPhotonAmplitude);
    debug << " orthogonalPhotonAmplitude : ";
    PrintComplex(forthogonalPhotonAmplitude);
    debug << "+--------------------------------------------------------------------------+" << endl;
}

TRestEvent* TRestAxionFieldPropagationProcess::ProcessEvent(TRestEvent* evInput) {
    mpfr::mpreal::set_default_prec(mpfr::digits2bits(30));
    cout.precision(30);
    fAxionEvent = (TRestAxionEvent*)evInput;

    TVector3 position = *(fAxionEvent->GetPosition());
    TVector3 direction = *(fAxionEvent->GetDirection());
    Double_t Ea = fAxionEvent->GetEnergy();
    Double_t ma = fAxionEvent->GetMass();

    faxionAmplitude = SetComplexReal(1.0, 0.0);
    fparallelPhotonAmplitude = SetComplexReal(0.0, 0.0);
    forthogonalPhotonAmplitude = SetComplexReal(0.0, 0.0);
    TVector3 averageBT = TVector3(0.0, 0.0, 0.0);  // initial value of the transverse component of the average magnetic field before entering the magnetic field region

    debug << "+------------------------+" << endl;
    debug << "Initial position of the axion input : " << endl;
    debug << "(" << position.X() << "," << position.Y() << "," << position.Z() << ")" << endl;
    debug << "Direction of the axion input : " << endl;
    debug << "(" << direction.X() << "," << direction.Y() << "," << direction.Z() << ")" << endl;
    debug << "Axion energy : " << Ea << endl;
    debug << "Axion mass : " << ma << endl;
    debug << "axion amplitude = " << faxionAmplitude.real << " + " << faxionAmplitude.img << "i" << endl;
    debug << "parallel photon amplitude = " << fparallelPhotonAmplitude.real << " + "
          << fparallelPhotonAmplitude.img << "i" << endl;
    debug << "orthogonal photon amplitude = " << forthogonalPhotonAmplitude.real << " + "
          << forthogonalPhotonAmplitude.img << "i" << endl;
    debug << "+------------------------+" << endl;

    std::vector<std::vector<TVector3>> boundaries;
    boundaries = FindFieldBoundaries();
    Int_t NofVolumes = boundaries.size();

    debug << "+------------------------+" << endl;
    debug << "Number of magnetic field regions through which the axion passes : " << NofVolumes << endl;
    debug << "+------------------------+" << endl;

    Double_t probability = 0.;  // initial value of the axion to photon conversion probability
    mpfr::mpreal axionMass = ma;  // in eV
    mpfr::mpreal photonMass = 0.0;

    for (Int_t i = 0; i < NofVolumes; i++) {    // loop over segments of trajectory where B is not equal to zero
        Int_t id = fAxionMagneticField->GetVolumeIndex(boundaries[i][0]);
        if (id < 0) {
            warning << "TRestAxionFieldPropagationProcess::ProcessEvent position is outside any volume. Setting photon mass to 0." << endl;
            photonMass = 0.0;  // in eV
        } else {
            Double_t mphoton = fAxionMagneticField->GetPhotonMass(id, Ea);  // It returns the effective photon mass in eV at the corresponding magnetic volume id.
            photonMass = mphoton;
        }
        debug << "Volume ID = " << id << "   photon mass = " << photonMass << endl;
        debug << "Volume ID = " << id << "   axion mass = " << axionMass << endl;
        debug << "Volume ID = " << id << "   axion energy = " << Ea << endl;

        // calculating common phase for the axion field and parallel component of the photon field (eqs. (4.1) and (4.2))
        // and phase for the orthogonal component of the photon field (eq. (4.3)) from the report "Axion-photon
        // conversion in the external magnetic fields"
        mpfr::mpreal CommonPhase = Ea * 1000.0 - (axionMass * axionMass + photonMass * photonMass) / (4. * Ea * 1000.0);   // in eV
        mpfr::mpreal OrthogonalPhase = Ea * 1000.0 - (photonMass * photonMass) / (2. * Ea * 1000.0);   // in eV

        debug << "Calculating amplitudes for one segment of trajectory where B is not zero. Segment boundaries are : (" << boundaries[i][0].X() << ","
              << boundaries[i][0].Y() << "," << boundaries[i][0].Z() << ") to (" << boundaries[i][1].X() << ","
              << boundaries[i][1].Y() << "," << boundaries[i][1].Z() << ")" << endl;
        CalculateAmplitudesInSegment(faxionAmplitude, fparallelPhotonAmplitude, forthogonalPhotonAmplitude, averageBT, axionMass, photonMass, Ea, boundaries[i][0], boundaries[i][1], CommonPhase, OrthogonalPhase);

    debug << endl << endl << endl << endl;
    debug << "----------------------------------------------------------------" << endl;
    debug << " Amplitude values after calculations in one segment: " << endl;
    debug << " axion amplitude : "; PrintComplex(faxionAmplitude);
    debug << " parallel photon component amplitude : "; PrintComplex(fparallelPhotonAmplitude);
    debug << " orthogonal photon component amplitude : "; PrintComplex(forthogonalPhotonAmplitude);
    debug << "+--------------------------------------------------------------------------+" << endl << endl;
        if ((i+1) < NofVolumes) {
            cout << "Calculating amplitudes along the part of trajectory where B = 0 between the two segments. Boundaries are : (" << boundaries[i][1].X() << ","
              << boundaries[i][1].Y() << "," << boundaries[i][1].Z() << ") to (" << boundaries[i+1][0].X() << ","
              << boundaries[i+1][0].Y() << "," << boundaries[i+1][0].Z() << ")" << endl;
            PropagateWithoutBField(faxionAmplitude, fparallelPhotonAmplitude, forthogonalPhotonAmplitude, axionMass, photonMass, Ea, boundaries[i][1], boundaries[i+1][0]);
        }
    }
    debug << endl << endl << endl << endl;
    debug << "----------------------------------------------------------------" << endl;
    debug << " FINAL AMPLITUDES: " << endl;
    debug << " axion amplitude : "; PrintComplex(faxionAmplitude);
    debug << " paralell photon component amplitude : "; PrintComplex(fparallelPhotonAmplitude);
    debug << " orthogonal photon component amplitude : "; PrintComplex(forthogonalPhotonAmplitude);
    debug << " PROBABILITY for axion to photon conversion (1-|a|^2): " << 1.0 - Norm2(faxionAmplitude) << endl;
    debug << "+--------------------------------------------------------------------------+" << endl;

    mpfr::mpreal probabilityHighPrecision = 1.0 - Norm2(faxionAmplitude);
    probability = probabilityHighPrecision.toDouble();
    fAxionEvent->SetGammaProbability(probability);
    debug << "+------------------------+" << endl;
    debug << "Conversion probability : " << endl;
    debug << "fAxionEvent->GetGammaProbability() = " << fAxionEvent->GetGammaProbability() << endl;
    debug << "+------------------------+" << endl;

    // if (fMode == "plan")
    //    fAxionEvent->SetPosition(MoveToPlane(position, direction, fFinalNormalPlan, fFinalPositionPlan));
    // if (fMode == "distance") fAxionEvent->SetPosition(MoveByDistance(position, direction, fDistance));

    debug << "+------------------------+" << endl;
    debug << "Final position of the axion : " << endl;
    debug << "(" << fAxionEvent->GetPositionX() << "," << fAxionEvent->GetPositionY() << ","
          << fAxionEvent->GetPositionZ() << ")" << endl;
    debug << "+------------------------+" << endl;

    if (GetVerboseLevel() >= REST_Debug) fAxionEvent->PrintEvent();

    boundaries.clear();

    return fAxionEvent;
}







///////////////////////////////////////////////
/// \brief Function reading input parameters from the RML TRestAxionFieldPropagationProcess metadata section
///
void TRestAxionFieldPropagationProcess::InitFromConfigFile() {
    this->Initialize();

    fMode = GetParameter("mode");
    fFinalNormalPlan = Get3DVectorParameterWithUnits("finalNPlan");
    fFinalPositionPlan = Get3DVectorParameterWithUnits("finalPositionPlan");
    fDistance = GetDblParameterWithUnits("distance");

    PrintMetadata();
}
