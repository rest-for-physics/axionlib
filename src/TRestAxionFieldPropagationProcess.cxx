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
/// help of TRestAxionMagneticField metadata class. This means that the magnetic field volume can be made of
/// several magnetic field "regions" (as described in the documentation of TRestAxionMagneticField, where each
/// "region" is defined by a separate `<addMagnetVolume...>` line in its RML configuration file. One such
/// example of generic magnetic field volume description is shown in the following figure, where the volume
/// consists of four regions labeled #1 . . . #4. In the following example the regions appear as connected
/// regions although in general they could be completely unconnected, or isolated regions. Anywhere where no
/// region is defined the field will be equal to zero.
///
/// \htmlonly <style>div.image img[src="AxionMagnetTrajectory.png"]{width:500px;}</style> \endhtmlonly
///
/// ![Schematic of the axion trajectory through the magnetic field volumes.](AxionMagnetTrajectory.png)
///
/// Here, the boundaries of each region are represented by blue lines. Also, in some (or all) regions the
/// magnetic field can be zero in its outer parts in order to define a more complex field shape inside.
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

using namespace std;

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
TRestAxionFieldPropagationProcess::TRestAxionFieldPropagationProcess(char* cfgFileName) { Initialize(); }

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionFieldPropagationProcess::~TRestAxionFieldPropagationProcess() { delete fAxionEvent; }

///////////////////////////////////////////////
/// \brief Function to initialize input/output event members and define the section name
///
void TRestAxionFieldPropagationProcess::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    fAxionEvent = new TRestAxionEvent();
}

///////////////////////////////////////////////
/// \brief Process initialization. Data members that require initialization just before start processing
/// should be initialized here.
///
void TRestAxionFieldPropagationProcess::InitProcess() {
    RESTDebug << "Entering ... TRestAxionGeneratorProcess::InitProcess" << RESTendl;

    fMagneticField = (TRestAxionMagneticField*)this->GetMetadata("TRestAxionMagneticField");

    if (!fMagneticField) {
        RESTError << "TRestAxionFieldPropagationprocess. Magnetic Field was not defined!" << RESTendl;
        exit(0);
    }

    if (!fAxionField) {
        fAxionField = new TRestAxionField();

        fBufferGas = (TRestAxionBufferGas*)this->GetMetadata("TRestAxionBufferGas");
        if (fBufferGas) fAxionField->AssignBufferGas(fBufferGas);
    }
}

TRestEvent* TRestAxionFieldPropagationProcess::ProcessEvent(TRestEvent* evInput) {
    // Already done by TRestAxionEventProcess
    fAxionEvent = (TRestAxionEvent*)evInput;
    RESTDebug << "TRestAxionFieldPropagationProcess::ProcessEvent : " << fAxionEvent->GetID() << RESTendl;

    std::vector<TVector3> trackBounds =
        fMagneticField->GetFieldBoundaries(fAxionEvent->GetPosition(), fAxionEvent->GetDirection());

    std::vector<Double_t> bProfile =
        fMagneticField->GetTransversalComponentAlongPath(trackBounds[0], trackBounds[1], fIntegrationStep);

    Double_t Ea = fAxionEvent->GetEnergy();
    Double_t ma = fAxionEvent->GetMass();

    Double_t prob = fAxionField->GammaTransmissionProbability(bProfile, fIntegrationStep, Ea, ma);
    SetObservableValue("probability", prob);

    Double_t lCoh = (bProfile.size() - 1) * fIntegrationStep;
    SetObservableValue("coherenceLength", lCoh);

    if (fBufferGas && fBufferGasAdditionalLength > 0) {
        Double_t Gamma = fBufferGas->GetPhotonAbsorptionLength(Ea);  // cm-1
        Double_t GammaL = Gamma * lCoh * units("cm");
        SetObservableValue("absorption", exp(-GammaL));
    }

    if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) fAxionEvent->PrintEvent();

    /// Missing to propagate the axion to the end of magnet bore?

    return fAxionEvent;
}
