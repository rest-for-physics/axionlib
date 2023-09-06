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
/// TRestAxionSolarQCDFlux will use a file in ASCII or binary format to initialize
/// a solar flux table that will describe the solar flux spectrum as a function
/// of the solar radius.
///
/// This class may serve to load any generic flux definition that is independent
/// from the axion mass. However, since the design of the class was motivated to
/// reproduce the standard QCD axion flux, therefore the name of the class.
///
/// Another scenario arises when the axion or an axion-like particle production
/// mechanism depends on its mass, and we may need to introduce a flux description
/// as the given in the class TRestAxionSolarHiddenPhotonFlux. Both classes are
/// prototyped by a pure base class TRestAxionSolarFlux that defines common methods
/// used to evaluate the flux, and generate Monte-Carlo events inside
/// TRestAxionGeneratorProcess.
///
/// ### Basic use
///
/// Once the class has been initialized, the main use of this class will be provided
/// by the method TRestAxionSolarQCDFlux::GetRandomEnergyAndRadius. This method will
/// return a random axion energy and position inside the solar radius following the
/// distributions given by the solar flux tables.
///
/// Description of the specific parameters accepted by this metadata class.
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
/// contributions by TRestAxionSolarQCDFlux::ReadFluxFile to be managed internally. Two
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
///     <TRestAxionSolarQCDFlux name="sunPrimakoff" verboseLevel="debug" >
///
///         <!-- Common parameters that belong to TRestAxionSolarFlux -->
///			<parameter name="couplingType" value="g_ag"/>
///			<parameter name="couplingStrength" value="1.e-10"/>
///
///         <!-- Specific parameters that belong to TRestAxionSolarQCDFlux -->
///			<parameter name="fluxDataFile" value="Primakoff_Gianotti_201904.dat"/>
///			<parameter name="fluxSptFile" value="Dummy_Galan_202202.spt"/>
///
///     </TRestAxionSolarQCDFlux>
/// \endcode
///
/// \warning When the flux is loaded manually inside the `restRoot` interactive
/// shell, or inside a macro or script, after metadata initialization, it is necessary
/// to call the method TRestAxionSolarQCDFlux::LoadTables to trigger the tables
/// initialization.
///
/// ### Performing MonteCarlo tests using pre-loaded tables
///
/// In order to test the response of different solar flux definitions we may use the script
/// `solarPlot.py` found at `pipeline/metadata/solarFlux/`. This script will generate a
/// number of particles and it will assign to each particle an energy and solar disk
/// location with the help of the method TRestAxionSolarQCDFlux::GetRandomEnergyAndRadius.
///
/// \code
/// python3 solarPlotQCD.py --fluxname LennertHoofABC --N 1000000
/// \endcode
///
/// By default, it will load the flux definition found at `fluxes.rml` from the
/// `axionlib-data` repository, and generate a `png` image with the resuts from the
/// Monte Carlo execution.
///
/// \htmlonly <style>div.image img[src="ABC_flux_MC.png"]{width:750px;}</style> \endhtmlonly
///
/// ![Solar flux distributions MC-generated with TRestAxionSolarQCDFlux.](ABC_flux_MC.png)
///
/// ### Reading solar flux tables from `.flux` files
///
/// In a similar way we may initialize the class using a `.flux` file. The `.flux` described
/// previously will contain a high definition flux table measured in `cm-2 s-1 keV-1` as a
/// function of the energy (from 0 to 20keV) and the inner solar radius (from 0 to 1). The
/// binning of this table will be typically between 1eV and 10eV.
///
/// The class TRestAxionSolarQCDFlux will be initialized directly using the `.flux` file
/// provided under the `fluxDataFile` parameter. During the initialization, the flux will be
/// splitted into two independent flux components. One smooth continuum component integrated
/// in 100 eV steps, and a monochromatic peak components.
///
/// In order to help with the identification of peaks we need to define the `binSize` used in
/// the `.flux` table and the `peakSigma` defining the number of sigmas over the average for
/// a bin to be considered a peak.
///
/// \code
///    <TRestAxionSolarQCDFlux name="LennertHoofABC_Flux" verboseLevel="warning" >
///        <parameter name="couplingType" value="g_ae"/>
///        <parameter name="couplingStrength" value="1.e-13"/>
///        <parameter name="fluxDataFile" value="ABC_LennertHoof_202203.flux"/>
///
///        <parameter name="binSize" value="10eV" />
///        <parameter name="peakSigma" value="10" />
///
///        <parameter name="seed" value="137" />
///    </TRestAxionSolarQCDFlux>
/// \endcode
///
/// We will be able to load this file as usual, using the following recipe inside
/// `restRoot`,
///
/// \code
///    TRestAxionSolarQCDFlux *sFlux = new TRestAxionSolarQCDFlux("fluxes.rml", "LennertHoofABC")
///    sFlux->Initialize()
///    TCanvas *c = sFlux->DrawSolarFluxes()
///    c->Print("ABC_FluxTable.png" )
/// \endcode
///
/// will generate the following figure.
///
/// \htmlonly <style>div.image img[src="ABC_FluxTable.png"]{width:750px;}</style> \endhtmlonly
///
/// ![Solar flux distributions generated with TRestAxionSolarQCDFlux::DrawSolarFlux.](ABC_FluxTable.png)
///
/// ### Exporting the solar flux tables
///
/// On top of that, we will be able to export those tables to the TRestAxionSolarQCDFlux standard
/// format to be used in later occasions.
///
/// \code
///    TRestAxionSolarQCDFlux *sFlux = new TRestAxionSolarQCDFlux("fluxes.rml", "LennertHoofABC")
///    sFlux->Initialize()
///    sFlux->ExportTables()
/// \endcode
///
/// which will produce two files, a binary table `.N200f` with the continuum flux, and an ASCII
/// table containning the `.spt` monochromatic lines. The filename root will be extracted from
/// the original `.flux` file. Optionally we may export the continuum flux to an ASCII file by
/// indicating it at the TRestAxionSolarQCDFlux::ExportTables method call. The files will be placed
/// at the REST user space, at `$HOME/.rest/export/` directory.
///
/// TODO Implement the method TRestAxionSolarQCDFlux::InitializeSolarTable using
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
/// 2023-May: Specific methods extracted from TRestAxionSolarFlux
///                Javier Galan
///
/// \class      TRestAxionSolarQCDFlux
/// \author     Javier Galan
///
/// <hr>
///

#include "TRestAxionSolarQCDFlux.h"
using namespace std;

ClassImp(TRestAxionSolarQCDFlux);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionSolarQCDFlux::TRestAxionSolarQCDFlux() : TRestAxionSolarFlux() {}

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
TRestAxionSolarQCDFlux::TRestAxionSolarQCDFlux(const char* cfgFileName, string name)
    : TRestAxionSolarFlux(cfgFileName) {
    LoadConfigFromFile(fConfigFileName, name);

    if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Info) PrintMetadata();
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionSolarQCDFlux::~TRestAxionSolarQCDFlux() {}

///////////////////////////////////////////////
/// \brief It will load the tables in memory by using the filename information provided
/// inside the metadata members.
///
Bool_t TRestAxionSolarQCDFlux::LoadTables(Double_t mass) {
    if (fFluxDataFile == "" && fFluxSptFile == "") return false;

    if (TRestTools::GetFileNameExtension(fFluxDataFile) == "flux") {
        ReadFluxFile();
    } else {
        LoadContinuumFluxTable();
        LoadMonoChromaticFluxTable();
    }

    IntegrateSolarFluxes();

    return true;
}

///////////////////////////////////////////////
/// \brief A helper method to load the data file containning continuum spectra as a
/// function of the solar radius. It will be called by TRestAxionSolarQCDFlux::Initialize.
///
void TRestAxionSolarQCDFlux::LoadContinuumFluxTable() {
    if (fFluxDataFile == "") {
        RESTDebug << "TRestAxionSolarflux::LoadContinuumFluxTable. No solar flux table was defined"
                  << RESTendl;
        return;
    }

    string fullPathName = SearchFile((string)fFluxDataFile);

    RESTDebug << "Loading table from file : " << RESTendl;
    RESTDebug << "File : " << fullPathName << RESTendl;

    std::vector<std::vector<Float_t>> fluxTable;
    if (TRestTools::GetFileNameExtension(fFluxDataFile) == "dat") {
        std::vector<std::vector<Double_t>> doubleTable;
        if (!TRestTools::ReadASCIITable(fullPathName, doubleTable)) {
            RESTError << "TRestAxionSolarQCDFlux::LoadContinuumFluxTable. " << RESTendl;
            RESTError << "Could not read solar flux table : " << fFluxDataFile << RESTendl;
            return;
        }
        for (const auto& row : doubleTable) {
            std::vector<Float_t> floatVec(row.begin(), row.end());
            fluxTable.push_back(floatVec);
        }
    } else if (TRestTools::IsBinaryFile(fFluxDataFile))
        TRestTools::ReadBinaryTable(fullPathName, fluxTable);
    else {
        fluxTable.clear();
        RESTError << "Filename extension was not recognized!" << RESTendl;
        RESTError << "Solar flux table will not be populated" << RESTendl;
        RESTError << "Filename extension: " << TRestTools::GetFileNameExtension(fFluxDataFile) << RESTendl;
    }

    if (fluxTable.size() != 100 && fluxTable[0].size() != 200) {
        fluxTable.clear();
        RESTError << "LoadContinuumFluxTable. The table does not contain the right number of rows or columns"
                  << RESTendl;
        RESTError << "Table will not be populated" << RESTendl;
    }

    for (unsigned int n = 0; n < fluxTable.size(); n++) {
        TH1F* h = new TH1F(Form("%s_ContinuumFluxAtRadius%d", GetName(), n), "", 200, 0, 20);
        for (unsigned int m = 0; m < fluxTable[n].size(); m++) h->SetBinContent(m + 1, fluxTable[n][m]);
        fFluxTable.push_back(h);
    }
}

///////////////////////////////////////////////
/// \brief A helper method to load the data file containning monochromatic spectral
/// lines as a function of the solar radius. It will be called by TRestAxionSolarQCDFlux::Initialize.
///
void TRestAxionSolarQCDFlux::LoadMonoChromaticFluxTable() {
    if (fFluxSptFile == "") {
        RESTDebug << "TRestAxionSolarflux::LoadMonoChromaticFluxTable. No solar flux monochromatic table was "
                     "defined"
                  << RESTendl;
        return;
    }

    string fullPathName = SearchFile((string)fFluxSptFile);

    RESTDebug << "Loading monochromatic lines from file : " << RESTendl;
    RESTDebug << "File : " << fullPathName << RESTendl;

    std::vector<std::vector<Double_t>> doubleTable;
    if (!TRestTools::ReadASCIITable(fullPathName, doubleTable)) {
        RESTError << "TRestAxionSolarQCDFlux::LoadMonoChromaticFluxTable." << RESTendl;
        RESTError << "Could not read solar flux table : " << fFluxSptFile << RESTendl;
        return;
    }

    std::vector<std::vector<Float_t>> asciiTable;
    for (const auto& row : doubleTable) {
        std::vector<Float_t> floatVec(row.begin(), row.end());
        asciiTable.push_back(floatVec);
    }

    fFluxLines.clear();

    if (asciiTable.size() != 101) {
        RESTError << "LoadMonoChromaticFluxTable. The table does not contain the right number of rows"
                  << RESTendl;
        RESTError << "Table will not be populated" << RESTendl;
        return;
    }

    for (unsigned int en = 0; en < asciiTable[0].size(); en++) {
        Float_t energy = asciiTable[0][en];
        TH1F* h = new TH1F(Form("%s_MonochromeFluxAtEnergy%6.4lf", GetName(), energy), "", 100, 0, 1);
        for (unsigned int r = 1; r < asciiTable.size(); r++) h->SetBinContent(r, asciiTable[r][en]);
        fFluxLines[energy] = h;
    }
}

///////////////////////////////////////////////
/// \brief It loads a .flux file. It will split continuum and monochromatic peaks, loading
/// both internal flux tables.
///
void TRestAxionSolarQCDFlux::ReadFluxFile() {
    if (fBinSize <= 0) {
        RESTError << "TRestAxionSolarflux::ReadFluxFile. Energy bin size of .flux file must be specified."
                  << RESTendl;
        RESTError << "Please, define binSize parameter in eV." << RESTendl;
        return;
    }

    if (fPeakSigma <= 0) {
        RESTWarning << "TRestAxionSolarflux::ReadFluxFile. Peak sigma must be specified to generate "
                       "monochromatic spectrum."
                    << RESTendl;
        RESTWarning << "Only continuum table will be generated. If this was intentional, please, ignore this "
                       "RESTWarning."
                    << RESTendl;
        return;
    }

    string fullPathName = SearchFile((string)fFluxDataFile);

    RESTDebug << "Loading flux table ...  " << RESTendl;
    RESTDebug << "File : " << fullPathName << RESTendl;
    std::vector<std::vector<Double_t>> fluxData;
    TRestTools::ReadASCIITable(fullPathName, fluxData, 3);

    RESTDebug << "Table loaded" << RESTendl;
    TH2F* originalHist = new TH2F("FullTable", "", 100, 0., 1., (Int_t)(20. / fBinSize), 0., 20.);
    TH2F* continuumHist = new TH2F("ContinuumTable", "", 100, 0., 1., (Int_t)(20. / fBinSize), 0., 20.);
    TH2F* spectrumHist = new TH2F("LinesTable", "", 100, 0., 1., (Int_t)(20. / fBinSize), 0., 20.);

    Double_t fluxBinSize = TRestTools::GetLowestIncreaseFromTable(fluxData, 1);

    for (const auto& data : fluxData) {
        Float_t r = 0.005 + data[0];
        Float_t en = data[1] - 0.005;
        Float_t flux = data[2] * fluxBinSize;  // flux in cm-2 s-1 bin-1

        originalHist->Fill(r, en, (Float_t)flux);
        continuumHist->Fill(r, en, (Float_t)flux);
    }
    RESTDebug << "Histograms filled" << RESTendl;

    Int_t peaks = 0;
    do {
        peaks = 0;
        // We just identify pronounced peaks, not smoothed gaussians!
        const int smearPoints = (Int_t)(5 / (fBinSize * 100));
        const int excludePoints = smearPoints / 5;
        for (const auto& data : fluxData) {
            Float_t r = 0.005 + data[0];
            Float_t en = data[1] - 0.005;

            Int_t binR = continuumHist->GetXaxis()->FindBin(r);
            Int_t binE = continuumHist->GetYaxis()->FindBin(en);

            Double_t avgFlux = 0;
            Int_t n = 0;
            for (int e = binE - smearPoints; e <= binE + smearPoints; e++) {
                if (e < 1 || (e < binE + excludePoints && e > binE - excludePoints)) continue;
                n++;
                avgFlux += continuumHist->GetBinContent(binR, e);
            }
            avgFlux /= n;

            Float_t targetBinFlux = continuumHist->GetBinContent(binR, binE);
            Float_t thrFlux = avgFlux + fPeakSigma * TMath::Sqrt(avgFlux);
            if (targetBinFlux > 0 && targetBinFlux > thrFlux) {
                continuumHist->SetBinContent(binR, binE, avgFlux);
                peaks++;
            }
        }
    } while (peaks > 0);

    for (int n = 0; n < originalHist->GetNbinsX(); n++)
        for (int m = 0; m < originalHist->GetNbinsY(); m++) {
            Float_t orig = originalHist->GetBinContent(n + 1, m + 1);
            Float_t cont = continuumHist->GetBinContent(n + 1, m + 1);

            spectrumHist->SetBinContent(n + 1, m + 1, orig - cont);
        }

    continuumHist->Rebin2D(1, (Int_t)(0.1 / fBinSize));  // cm-2 s-1 (100eV)-1
    continuumHist->Scale(10);                            // cm-2 s-1 keV-1
    // It could be over here if we would use directly a TH2D

    fFluxTable.clear();
    for (int n = 0; n < continuumHist->GetNbinsX(); n++) {
        TH1F* hc =
            (TH1F*)continuumHist->ProjectionY(Form("%s_ContinuumFluxAtRadius%d", GetName(), n), n + 1, n + 1);
        fFluxTable.push_back(hc);
    }

    fFluxLines.clear();
    for (int n = 0; n < spectrumHist->GetNbinsY(); n++) {
        if (spectrumHist->ProjectionX("", n + 1, n + 1)->Integral() > 0) {
            Double_t energy = spectrumHist->ProjectionY()->GetBinCenter(n + 1);
            TH1F* hm = (TH1F*)spectrumHist->ProjectionX(
                Form("%s_MonochromeFluxAtEnergy%5.3lf", GetName(), energy), n + 1, n + 1);
            fFluxLines[energy] = hm;
        }
    }

    cout << "Number of peaks identified: " << fFluxLines.size() << endl;
}

///////////////////////////////////////////////
/// \brief It builds a histogram with the continuum spectrum component.
/// The flux will be expressed in cm-2 s-1 keV-1. Binned in 100eV steps.
///
TH1F* TRestAxionSolarQCDFlux::GetContinuumSpectrum() {
    if (fContinuumHist != nullptr) {
        delete fContinuumHist;
        fContinuumHist = nullptr;
    }

    fContinuumHist = new TH1F("ContinuumHist", "", 200, 0, 20);
    for (const auto& x : fFluxTable) {
        fContinuumHist->Add(x);
    }

    fContinuumHist->SetStats(0);
    fContinuumHist->GetXaxis()->SetTitle("Energy [keV]");
    fContinuumHist->GetXaxis()->SetTitleSize(0.05);
    fContinuumHist->GetXaxis()->SetLabelSize(0.05);
    fContinuumHist->GetYaxis()->SetTitle("Flux [cm-2 s-1 keV-1]");
    fContinuumHist->GetYaxis()->SetTitleSize(0.05);
    fContinuumHist->GetYaxis()->SetLabelSize(0.05);

    return fContinuumHist;
}

///////////////////////////////////////////////
/// \brief It builds a histogram with the monochromatic spectrum component.
/// The flux will be expressed in cm-2 s-1 eV-1. Binned in 1eV steps.
///
TH1F* TRestAxionSolarQCDFlux::GetMonochromaticSpectrum() {
    if (fMonoHist != nullptr) {
        delete fMonoHist;
        fMonoHist = nullptr;
    }

    fMonoHist = new TH1F("MonochromaticHist", "", 20000, 0, 20);
    for (const auto& x : fFluxLines) {
        fMonoHist->Fill(x.first, x.second->Integral());  // cm-2 s-1 eV-1
    }

    fMonoHist->SetStats(0);
    fMonoHist->GetXaxis()->SetTitle("Energy [keV]");
    fMonoHist->GetXaxis()->SetTitleSize(0.05);
    fMonoHist->GetXaxis()->SetLabelSize(0.05);
    fMonoHist->GetYaxis()->SetTitle("Flux [cm-2 s-1 eV-1]");
    fMonoHist->GetYaxis()->SetTitleSize(0.05);
    fMonoHist->GetYaxis()->SetLabelSize(0.05);

    return fMonoHist;
}

///////////////////////////////////////////////
/// \brief It builds a histogram adding the continuum and the monochromatic
/// spectrum component. The flux will be expressed in cm-2 s-1 keV-1.
/// Binned in 1eV steps.
///
TH1F* TRestAxionSolarQCDFlux::GetTotalSpectrum() {
    TH1F* hm = GetMonochromaticSpectrum();
    TH1F* hc = GetContinuumSpectrum();

    if (fTotalHist != nullptr) {
        delete fTotalHist;
        fTotalHist = nullptr;
    }

    fTotalHist = new TH1F("fTotalHist", "", 20000, 0, 20);
    for (int n = 0; n < hc->GetNbinsX(); n++) {
        for (int m = 0; m < 100; m++) {
            fTotalHist->SetBinContent(n * 100 + 1 + m, hc->GetBinContent(n + 1));
        }
    }

    for (int n = 0; n < hm->GetNbinsX(); n++)
        // 1e-2 is the renormalization from 20000 bins to 200 bins
        fTotalHist->SetBinContent(n + 1, fTotalHist->GetBinContent(n + 1) + 100 * hm->GetBinContent(n + 1));

    fTotalHist->SetStats(0);
    fTotalHist->GetXaxis()->SetTitle("Energy [keV]");
    fTotalHist->GetXaxis()->SetTitleSize(0.05);
    fTotalHist->GetXaxis()->SetLabelSize(0.05);
    fTotalHist->GetYaxis()->SetTitle("Flux [cm-2 s-1 keV-1]");
    fTotalHist->GetYaxis()->SetTitleSize(0.05);
    fTotalHist->GetYaxis()->SetLabelSize(0.05);

    return fTotalHist;
}

///////////////////////////////////////////////
/// \brief It draws the contents of a .flux file. This method just receives the
///
TCanvas* TRestAxionSolarQCDFlux::DrawSolarFlux() {
    if (fCanvas != nullptr) {
        delete fCanvas;
        fCanvas = nullptr;
    }
    fCanvas = new TCanvas("canv", "This is the canvas title", 1200, 500);
    fCanvas->Draw();

    TPad* pad1 = new TPad("pad1", "This is pad1", 0.01, 0.02, 0.99, 0.97);
    pad1->Divide(2, 1);
    pad1->Draw();

    pad1->cd(1);
    pad1->cd(1)->SetLogy();
    pad1->cd(1)->SetRightMargin(0.09);
    pad1->cd(1)->SetLeftMargin(0.15);
    pad1->cd(1)->SetBottomMargin(0.15);

    TH1F* ht = GetTotalSpectrum();
    ht->SetLineColor(kBlack);
    ht->SetFillStyle(4050);
    ht->SetFillColor(kBlue - 10);

    TH1F* hm = GetMonochromaticSpectrum();
    hm->SetLineColor(kBlack);
    hm->Scale(100);  // renormalizing per 100eV-1

    ht->Draw("hist");
    hm->Draw("hist same");

    pad1->cd(2);
    pad1->cd(2)->SetRightMargin(0.09);
    pad1->cd(2)->SetLeftMargin(0.15);
    pad1->cd(2)->SetBottomMargin(0.15);

    ht->Draw("hist");
    hm->Draw("hist same");

    return fCanvas;
}

///////////////////////////////////////////////
/// \brief A helper method to initialize the internal class data members with the
/// integrated flux for each solar ring. It will be called by TRestAxionSolarQCDFlux::Initialize.
///
void TRestAxionSolarQCDFlux::IntegrateSolarFluxes() {
    fFluxLineIntegrals.clear();
    fTotalMonochromaticFlux = 0;

    for (const auto& line : fFluxLines) {
        fTotalMonochromaticFlux += line.second->Integral();
        fFluxLineIntegrals.push_back(fTotalMonochromaticFlux);
    }

    for (unsigned int n = 0; n < fFluxLineIntegrals.size(); n++)
        fFluxLineIntegrals[n] /= fTotalMonochromaticFlux;

    fTotalContinuumFlux = 0.0;
    for (unsigned int n = 0; n < fFluxTable.size(); n++) {
        fTotalContinuumFlux += fFluxTable[n]->Integral() * 0.1;  // We integrate in 100eV steps
        fFluxTableIntegrals.push_back(fTotalContinuumFlux);
    }

    for (unsigned int n = 0; n < fFluxTableIntegrals.size(); n++)
        fFluxTableIntegrals[n] /= fTotalContinuumFlux;

    fFluxRatio = fTotalMonochromaticFlux / (fTotalContinuumFlux + fTotalMonochromaticFlux);
}

///////////////////////////////////////////////
/// \brief It returns the integrated flux at earth in cm-2 s-1 for the given energy range
///
Double_t TRestAxionSolarQCDFlux::IntegrateFluxInRange(TVector2 eRange, Double_t mass) {
    if (eRange.X() == -1 && eRange.Y() == -1) {
        if (GetTotalFlux() == 0) IntegrateSolarFluxes();
        return GetTotalFlux();
    }

    Double_t flux = 0;
    for (const auto& line : fFluxLines)
        if (line.first > eRange.X() && line.first < eRange.Y()) flux += line.second->Integral();

    fTotalContinuumFlux = 0.0;
    for (unsigned int n = 0; n < fFluxTable.size(); n++) {
        flux += fFluxTable[n]->Integral(fFluxTable[n]->FindFixBin(eRange.X()),
                                        fFluxTable[n]->FindFixBin(eRange.Y())) *
                0.1;  // We integrate in 100eV steps
    }

    return flux;
}

///////////////////////////////////////////////
/// \brief It returns a random solar radius position and energy according to the
/// flux distributions defined inside the solar tables loaded in the class
///
std::pair<Double_t, Double_t> TRestAxionSolarQCDFlux::GetRandomEnergyAndRadius(TVector2 eRange,
                                                                               Double_t mass) {
    std::pair<Double_t, Double_t> result = {0, 0};
    if (!AreTablesLoaded()) return result;
    Double_t rnd = fRandom->Rndm();
    if (fTotalMonochromaticFlux == 0 || fRandom->Rndm() > fFluxRatio) {
        // Continuum
        for (unsigned int r = 0; r < fFluxTableIntegrals.size(); r++) {
            if (rnd < fFluxTableIntegrals[r]) {
                Double_t energy = fFluxTable[r]->GetRandom();
                if (eRange.X() != -1 && eRange.Y() != -1) {
                    if (energy < eRange.X() || energy > eRange.Y()) return GetRandomEnergyAndRadius(eRange);
                }
                Double_t radius = ((Double_t)r + fRandom->Rndm()) * 0.01;
                std::pair<Double_t, Double_t> p = {energy, radius};
                return p;
            }
        }
    } else {
        // Monochromatic
        int n = 0;
        for (const auto& line : fFluxLines) {
            if (rnd < fFluxLineIntegrals[n]) {
                std::pair<Double_t, Double_t> p = {line.first, line.second->GetRandom()};
                return p;
            }
            n++;
        }
    }
    return result;
}

///////////////////////////////////////////////
/// \brief It prints on screen the table that has been loaded in memory
///
void TRestAxionSolarQCDFlux::PrintContinuumSolarTable() {
    cout << "Continuum solar flux table: " << endl;
    cout << "--------------------------- " << endl;
    for (unsigned int n = 0; n < fFluxTable.size(); n++) {
        for (int m = 0; m < fFluxTable[n]->GetNbinsX(); m++)
            cout << fFluxTable[n]->GetBinContent(m + 1) << "\t";
        cout << endl;
        cout << endl;
    }
    cout << endl;
}

///////////////////////////////////////////////
/// \brief It prints on screen the integrated solar flux per solar ring
///
void TRestAxionSolarQCDFlux::PrintIntegratedRingFlux() {
    cout << "Integrated solar flux per solar ring: " << endl;
    cout << "--------------------------- " << endl;
    /*
    for (int n = 0; n < fFluxPerRadius.size(); n++)
    cout << "n : " << n << " flux : " << fFluxPerRadius[n] << endl;
    cout << endl;
    */
}

///////////////////////////////////////////////
/// \brief It prints on screen the spectral lines loaded in memory
///
void TRestAxionSolarQCDFlux::PrintMonoChromaticFlux() {
    //   cout << "Number of monochromatic lines: " << fFluxPerRadius.size() << endl;
    cout << "+++++++++++++++++++++++++++++++++++" << endl;
    for (auto const& line : fFluxLines) {
        cout << "Energy : " << line.first << " keV" << endl;
        cout << "-----------------" << endl;
        for (int n = 0; n < line.second->GetNbinsX(); n++)
            cout << "R : " << line.second->GetBinCenter(n + 1)
                 << " flux : " << line.second->GetBinContent(n + 1) << " cm-2 s-1" << endl;
    }
}

///////////////////////////////////////////////
/// \brief Prints on screen the information about the metadata members of TRestAxionSolarQCDFlux
///
void TRestAxionSolarQCDFlux::PrintMetadata() {
    TRestAxionSolarFlux::PrintMetadata();

    if (fFluxDataFile != "")
        RESTMetadata << " - Solar axion flux datafile (continuum) : " << fFluxDataFile << RESTendl;
    if (fFluxSptFile != "")
        RESTMetadata << " - Solar axion flux datafile (monochromatic) : " << fFluxSptFile << RESTendl;
    RESTMetadata << "-------" << RESTendl;
    RESTMetadata << " - Total monochromatic flux : " << fTotalMonochromaticFlux << " cm-2 s-1" << RESTendl;
    RESTMetadata << " - Total continuum flux : " << fTotalContinuumFlux << " cm-2 s-1" << RESTendl;
    RESTMetadata << "++++++++++++++++++" << RESTendl;

    if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Debug) {
        PrintContinuumSolarTable();
        PrintMonoChromaticFlux();
        PrintIntegratedRingFlux();
    }
}

///////////////////////////////////////////////
/// \brief It will create files with the continuum and spectral flux components to be used
/// in a later ocasion.
///
void TRestAxionSolarQCDFlux::ExportTables(Bool_t ascii) {
    string rootFilename = TRestTools::GetFileNameRoot(fFluxDataFile);

    string path = REST_USER_PATH + "/export/";

    if (!TRestTools::fileExists(path)) {
        std::cout << "Creating path: " << path << std::endl;
        int z = system(("mkdir -p " + path).c_str());
        if (z != 0) RESTError << "Could not create directory " << path << RESTendl;
    }

    if (fFluxTable.size() > 0) {
        std::vector<std::vector<Float_t>> table;
        for (const auto& x : fFluxTable) {
            std::vector<Float_t> row;
            for (int n = 0; n < x->GetNbinsX(); n++) row.push_back(x->GetBinContent(n + 1));

            table.push_back(row);
        }

        if (!ascii)
            TRestTools::ExportBinaryTable(path + "/" + rootFilename + ".N200f", table);
        else
            TRestTools::ExportASCIITable(path + "/" + rootFilename + ".dat", table);
    }

    if (fFluxLines.size() > 0) {
        std::vector<std::vector<Float_t>> table;
        for (const auto& x : fFluxLines) {
            std::vector<Float_t> row;
            row.push_back(x.first);
            for (int n = 0; n < x.second->GetNbinsX(); n++) row.push_back(x.second->GetBinContent(n + 1));

            table.push_back(row);
        }

        TRestTools::TransposeTable(table);

        TRestTools::ExportASCIITable(path + "/" + rootFilename + ".spt", table);
    }
}
