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
/// TRestAxionSolarFlux will use a file in ASCII or binary format to initialize
/// a solar flux table that will describe the solar flux spectrum as a function
/// of the solar radius. It will be also possible to generate the solar table
/// by other means.
///
/// Once the class has been initialized, the main use of this class will be provided
/// by the method TRestAxionSolarFlux::GetRandomEnergyAndRadius. This method will
/// return a random axion energy and position inside the solar radius following the
/// distributions given by the solar flux tables.
///
/// For the moment, in order to trace the nature and intensity of the coupling in
/// future ray-tracking results we need to define the parameters `couplingType` and
/// `couplingStrength`. The ray-tracing processing will be done for different
/// coupling components in different event processing chains.
///
/// Description of the parameters accepted by this metadata class.
/// - *couplingType:* A string describing the coupling type, i.e. g_ag, g_ae, g_an, ...
/// - *couplingStrength:* The intensity of the coupling used to calculate the values
/// given in the solar flux tables.
/// - *fluxDataFile:* A table with 100 rows representing the solar ring flux from the
/// center to the corona, and 200 columns representing the flux, measured in cm-2 s-1 keV-1,
/// for the range (0,20)keV in steps of 100eV. The table could be provided in ASCII format,
/// using `.dat` extension, or it might be a binary table using `.N200f` extension.
/// - *fluxSptFile:* A table where each column represents a monochromatic energy. The
/// column contains 101 rows, the first element is the energy of the monochromatic line
/// while the next 100 elements contain the flux, measured in cm-2 s-1, integrated to
/// each solar ring, being the second element the ring in the center of the sun.
///
/// Additionally this class will be able to read `.flux` files that are the original files
/// produced in 3-columns format (inner radius [solar units] / energy [keV] /
/// flux [cm-2 s-1 keV-1]). The `.flux` files may contain the full information,
/// continuum and spectral components. Those components will be splited into two independent
/// contributions by TRestAxionSolarFlux::ReadFluxFile to be managed internally. Two
/// additional parameters *will be required* to translate the `.flux` files into the tables
/// that are understood by this class.
/// - *binSize:* The energy binning used on the `.flux` file and inside the histogram used
/// for monochromatic lines identification.
/// - *peakSigma:* The ratio between the flux provided at the `.flux` file and the
/// average flux calculated in the peak surroundings. If the flux ratio is higher than
/// this value, the flux at that particular bin will be considered a peak.
///
/// Optionally, if we want to consider a different binning on the monochromatic/continuum
/// histogram used internally for the calculation we may specify optionally a new parameter.
/// In that case, `fBinSize` will be the binning of the internal histogram, while the new
/// parameter will be the binning given inside the `.flux` file.
/// - *fluxBinSize:* The bin size used on the `.flux` table.
///
/// Pre-generated solar axion flux tables will be available at the
/// [axionlib-data](https://github.com/rest-for-physics/axionlib-data/tree/master)
/// repository. The different RML flux definitions used to load those tables
/// will be found at the
/// [fluxes.rml](https://github.com/rest-for-physics/axionlib-data/blob/master/solarFlux/fluxes.rml)
/// file found at the axionlib-data repository.
///
/// Inside a local REST installation, the `fluxes.rml` file will be found at the REST
/// installation directory, and it will be located automatically by the
/// TRestMetadata::SearchFile method.
///
/// ### A basic RML definition
///
/// The following definition integrates an axion-photon component with a continuum
/// spectrum using a Primakoff production model, and a dummy spectrum file that
/// includes two monocrhomatic lines at different solar disk radius positions.
///
/// \code
///     <TRestAxionSolarFlux name="sunPrimakoff" verboseLevel="debug" >
///			<parameter name="couplingType" value="g_ag"/>
///			<parameter name="couplingStrength" value="1.e-10"/>
///			<parameter name="fluxDataFile" value="Primakoff_Gianotti_201904.dat"/>
///			<parameter name="fluxSptFile" value="Dummy_Galan_202202.spt"/>
///     </TRestAxionSolarFlux>
/// \endcode
///
/// \warning When the flux is loaded manually inside the `restRoot` interactive
/// shell, or inside a macro or script, after metadata initialization, it is necessary
/// to call the method TRestAxionSolarFlux::LoadTables to trigger the tables
/// initialization.
///
/// ### Performing MonteCarlo tests using pre-loaded tables
///
/// In order to test the response of different solar flux definitions we may use the script
/// `solarPlot.py` found at `pipeline/metadata/solarFlux/`. This script will generate a
/// number of particles and it will assign to each particle an energy and solar disk
/// location with the help of the method TRestAxionSolarFlux::GetRandomEnergyAndRadius.
///
/// \code
/// python3 solarPlot.py --fluxname LennertHoofABC --N 1000000
/// \endcode
///
/// By default, it will load the flux definition found at `fluxes.rml` from the
/// `axionlib-data` repository, and generate a `png` image with the resuts from the
/// Monte Carlo execution.
///
/// \htmlonly <style>div.image img[src="ABC_flux_MC.png"]{width:750px;}</style> \endhtmlonly
///
/// ![Solar flux distributions MC-generated with TRestAxionSolarFlux.](ABC_flux_MC.png)
///
/// ### Reading solar flux tables from `.flux` files
///
/// In a similar way we may initialize the class using a `.flux` file. The `.flux` described
/// previously will contain a high definition flux table measured in `cm-2 s-1 keV-1` as a
/// function of the energy (from 0 to 20keV) and the inner solar radius (from 0 to 1). The
/// binning of this table will be typically between 1eV and 10eV.
///
/// The class TRestAxionSolarFlux will be initialized directly using the `.flux` file
/// provided under the `fluxDataFile` parameter. During the initialization, the flux will be
/// splitted into two independent flux components. One smooth continuum component integrated
/// in 100 eV steps, and a monochromatic peak components.
///
/// In order to help with the identification of peaks we need to define the `binSize` used in
/// the `.flux` table and the `peakSigma` defining the number of sigmas over the average for
/// a bin to be considered a peak.
///
/// \code
///    <TRestAxionSolarFlux name="LennertHoofABC_Flux" verboseLevel="warning" >
///        <parameter name="couplingType" value="g_ae"/>
///        <parameter name="couplingStrength" value="1.e-13"/>
///        <parameter name="fluxDataFile" value="ABC_LennertHoof_202203.flux"/>
///
///        <parameter name="binSize" value="10eV" />
///        <parameter name="peakSigma" value="10" />
///
///        <parameter name="seed" value="137" />
///    </TRestAxionSolarFlux>
/// \endcode
///
/// We will be able to load this file as usual, using the following recipe inside
/// `restRoot`,
///
/// \code
///    TRestAxionSolarFlux *sFlux = new TRestAxionSolarFlux("fluxes.rml", "LennertHoofABC")
///    sFlux->LoadTables()
///    TCanvas *c = sFlux->DrawSolarFluxes()
///    c->Print("ABC_FluxTable.png" )
/// \endcode
///
/// will generate the following figure.
///
/// \htmlonly <style>div.image img[src="ABC_FluxTable.png"]{width:750px;}</style> \endhtmlonly
///
/// ![Solar flux distributions generated with TRestAxionSolarFlux::DrawSolarFlux.](ABC_FluxTable.png)
///
/// ### Exporting the solar flux tables
///
/// On top of that, we will be able to export those tables to the TRestAxionSolarFlux standard
/// format to be used in later occasions.
///
/// \code
///    TRestAxionSolarFlux *sFlux = new TRestAxionSolarFlux("fluxes.rml", "LennertHoofABC")
///    sFlux->LoadTables()
///    sFlux->ExportTables()
/// \endcode
///
/// which will produce two files, a binary table `.N200f` with the continuum flux, and an ASCII
/// table containning the `.spt` monochromatic lines. The filename root will be extracted from
/// the original `.flux` file. Optionally we may export the continuum flux to an ASCII file by
/// indicating it at the TRestAxionSolarFlux::ExportTables method call. The files will be placed
/// at the REST user space, at `$HOME/.rest/export/` directory.
///
/// TODO Implement the method TRestAxionSolarFlux::InitializeSolarTable using
/// a solar model description by TRestAxionSolarModel.
///
/// TODO Perhaps it would be interesting to replace fFluxTable for a TH2D
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2022-February: Recovered from original TRestAxionSolarModel implementation
///                Javier Galan
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

    if (!fRandom) {
        delete fRandom;
        fRandom = nullptr;
    }

    if (fRandom != nullptr) {
        delete fRandom;
        fRandom = nullptr;
    }

    fRandom = new TRandom3(fSeed);
    fSeed = fRandom->TRandom::GetSeed();
}

///////////////////////////////////////////////
/// \brief It builds a histogram using the contents of the .flux file given
/// in the argument.
///
TH1F* TRestAxionSolarFlux::GetFluxHistogram(string fname, Double_t binSize) {
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

    return (TH1F*)originalHist->ProjectionY();
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

    TH1F* ht = GetEnergySpectrum();
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
