/******************** REST disclaimer ***********************************
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
/// TRestAxionSolarFlux encapsulates common metadata members and methods that can
/// be later used as a common way to evaluate the axion solar flux from different
/// implementations.
///
/// At present, two different implemenations exist: TRestAxionSolarQCDFlux and
/// TRestAxionSolarHiddenPhotonFlux.
///
/// ### Common metadata members defined by this class
///
/// In order to trace the nature and intensity of the coupling in
/// future ray-tracking results we need to define the parameters `couplingType` and
/// `couplingStrength`. The ray-tracing processing will be done for different
/// coupling components in different event processing chains that will be combined
/// in a final sensitivity plot.
///
/// - *couplingType:* A string describing the coupling type, i.e. g_ag, g_ae, g_an, ...
/// - *couplingStrength:* The intensity of the coupling used to calculate the values
/// given in the solar flux tables.
///
/// ### Pre-generated solar axion flux tables
///
/// Some tables  will be available at the
/// [axionlib-data](https://github.com/rest-for-physics/axionlib-data/tree/master)
/// repository. The different RML flux definitions used to load those tables
/// will be found at the
/// [fluxes.rml](https://github.com/rest-for-physics/axionlib-data/blob/master/solarFlux/fluxes.rml)
/// file found at the axionlib-data repository.
///
/// Inside a local REST installation, the `fluxes.rml` file will be found at the REST
/// installation directory, and it will be located automatically by the
/// TRestMetadata::SearchFile method when initializing a TRestAxionSolarFlux class.
///
/// ### Re-implementing a specific TRestAxionSolarXYZFlux class
///
/// Different pure virtual methods must be re-implemented in the inherited class:
///
/// - **IntegratedFluxInRange**: It should return a double value with the total integrated
/// flux in cm-2 s-1 for a given energy range.
/// - **GetTotalFlux**: It should return the total flux in cm-2 s-1
/// - **GetRandomEnergyAndRadius**: It should return a random (energy,radius) sample based
/// on the flux distribution initialized.
/// - **LoadTables**: It is called by TRestAxionSolarFlux::Initialize to allow the inherited
/// class to load all the necessary tables in memory.
/// - **GetEnergySpectrum**: It should return a TH1D pointer with a energy spectrum histogram.
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2022-February: Recovered from original TRestAxionSolarModel implementation
///                Javier Galan
/// 2023-May:      This class has been transformed into a pure abstract class
///
/// \class      TRestAxionSolarFlux
/// \author     Javier Galan
///
/// <hr>
///

#include "TRestAxionSolarFlux.h"
using namespace std;

ClassImp(TRestAxionSolarFlux);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionSolarFlux::TRestAxionSolarFlux() : TRestMetadata() {}

///////////////////////////////////////////////
/// \brief Constructor loading data from a config file
///
/// If no configuration path is defined using TRestMetadata::SetConfigFilePath
/// the path to the config file must be specified using full path, absolute or
/// relative.
///
/// The default behaviour is that the config file must be specified with
/// full path, absolute or relative.
///
/// \param cfgFileName A const char* giving the path to an RML file.
/// \param name The name of the specific metadata. It will be used to find the
/// corresponding TRestAxionMagneticField section inside the RML.
///
TRestAxionSolarFlux::TRestAxionSolarFlux(const char* cfgFileName, string name) : TRestMetadata(cfgFileName) {
    RESTDebug << "Entering TRestAxionSolarFlux constructor( cfgFileName, name )" << RESTendl;
    RESTDebug << "File: " << cfgFileName << " Name: " << name << RESTendl;
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionSolarFlux::~TRestAxionSolarFlux() {}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionSolarFlux members
///
void TRestAxionSolarFlux::Initialize() {
    SetLibraryVersion(LIBRARY_VERSION);

    fTablesLoaded = false;
    if (LoadTables()) fTablesLoaded = true;

    if (fRandom) {
        delete fRandom;
        fRandom = nullptr;
    }

    fRandom = new TRandom3(fSeed);
    fSeed = fRandom->TRandom::GetSeed();
}


///////////////////////////////////////////////
/// \brief Initialization of TRestAxionSolarFlux members with specific mass
///
//void TRestAxionSolarFlux::InitializeMass( Double_t mass ) { SetMass(mass); RESTMetadata << GetMass() << RESTendl; }    // SetMass calls Initialize


///////////////////////////////////////////////
/// \brief It builds a histogram using the contents of the .flux file given
/// in the argument.
///
TH1D* TRestAxionSolarFlux::GetFluxHistogram(string fname, Double_t binSize) {
    string fullPathName = SearchFile(fname);

    std::vector<std::vector<Double_t>> fluxData;
    TRestTools::ReadASCIITable(fullPathName, fluxData, 3);

    TH2F* originalHist =
        new TH2F(Form("FluxTable_%s", GetName()), "", 100, 0., 1., (Int_t)(20. / binSize), 0., 20.);

    for (const auto& data : fluxData) {
        Double_t r = 0.005 + data[0];
        Double_t en = data[1] - 0.005;
        Double_t flux = data[2] * binSize;  // flux in cm-2 s-1 bin-1

        originalHist->Fill(r, en, flux);
    }

    return (TH1D*)originalHist->ProjectionY();
}

///////////////////////////////////////////////
/// \brief It draws the contents of a .flux file. This method just receives the
/// name of the .flux file and it works stand-alone.
///
TCanvas* TRestAxionSolarFlux::DrawFluxFile(string fname, Double_t binSize) {
    if (fCanvas != nullptr) {
        delete fCanvas;
        fCanvas = nullptr;
    }
    fCanvas = new TCanvas("canv", "This is the canvas title", 1400, 1200);
    fCanvas->Draw();

    TPad* pad1 = new TPad("pad1", "This is pad1", 0.01, 0.02, 0.99, 0.97);
    pad1->Draw();

    fCanvas->cd();
    pad1->cd();

    GetFluxHistogram(fname, binSize)->Draw("hist");

    return fCanvas;
}

///////////////////////////////////////////////
/// \brief It draws the contents of a .flux file. This method just receives the
///
TCanvas* TRestAxionSolarFlux::DrawSolarFlux() {
    if (fCanvas != nullptr) {
        delete fCanvas;
        fCanvas = nullptr;
    }
    fCanvas = new TCanvas("canv", "This is the canvas title", 1200, 500);
    fCanvas->Draw();

    TPad* pad1 = new TPad("pad1", "This is pad1", 0.01, 0.02, 0.99, 0.97);
    pad1->Draw();

    pad1->cd();
    pad1->SetLogy();
    pad1->SetRightMargin(0.09);
    pad1->SetLeftMargin(0.15);
    pad1->SetBottomMargin(0.15);

    TH1D* ht = GetEnergySpectrum();
    ht->SetLineColor(kBlack);
    ht->SetFillStyle(4050);
    ht->SetFillColor(kBlue - 10);

    ht->Draw("hist");

    return fCanvas;
}

///////////////////////////////////////////////
/// \brief Prints on screen the information about the metadata members of TRestAxionSolarFlux
///
void TRestAxionSolarFlux::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    RESTMetadata << " - Coupling type : " << fCouplingType << RESTendl;
    RESTMetadata << " - Coupling strength : " << fCouplingStrength << RESTendl;
    RESTMetadata << "--------" << RESTendl;
    RESTMetadata << " - Random seed : " << fSeed << RESTendl;
    RESTMetadata << "--------" << RESTendl;
}
