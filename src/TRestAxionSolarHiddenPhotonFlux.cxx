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
/// TRestAxionSolarHiddenPhotonFlux will use a file in binary format to initialize
/// a solar flux table that will describe the solar hidden photon flux spectrum as a function
/// of the solar radius.
///
/// This class loads the hidden photon flux that depends on the mass and kinetic mixing parameter.
/// For axion-like particle prodution independent of mass there is the class
/// TRestAxionSolarQCDFlux. Both classes are prototyped by a pure base class TRestAxionSolarFlux
/// that defines common methods used to evaluate the flux, and generate Monte-Carlo events inside
/// TRestAxionGeneratorProcess.
///
/// ### Basic use
///
/// Once the class has been initialized, the main use of this class will be provided
/// by the method TRestAxionSolarHiddenPhotonFlux::GetRandomEnergyAndRadius. This method will
/// return a random axion energy and position inside the solar radius following the
/// distributions given by the solar flux tables.
///
/// Description of the specific parameters accepted by this metadata class.
/// - *fluxDataFile:* A table with 1000 rows representing the solar ring flux from the
/// center to the corona, and 200 columns representing the flux, measured in cm-2 s-1 keV-1,
/// for the range (0,20)keV in steps of 100eV. The table is provided as a binary table using
/// `.N200f` extension.
/// - *widthDataFile:* A table with 1000 rows representing the width of the hidden photon
/// resonant production (wG) for each solar ring from the center to the corona, and 200 columns
/// representing the width, measured in eV2, for the range (0,20)keV in steps of 100eV. The
/// table is provided as a binary table using `.N200f` extension.
/// - *plasmaFreqDataFile:* A table with 1000 rows and only 1 column representing the solar
/// plasma frequency (wp) for each solar ring from the center to the corona, measured in eV. The
/// table is provided as a binary table using `.N1f` extension.
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
/// to call the method TRestAxionSolarHiddenPhotonFlux::LoadTables(mass) to trigger the tables
/// initialization.
///
/// ### Performing MonteCarlo tests using pre-loaded tables
///
/// In order to test the response of different solar flux definitions we may use the script
/// `solarPlot.py` found at `pipeline/metadata/solarFlux/`. This script will generate a
/// number of particles and it will assign to each particle an energy and solar disk
/// location with the help of the method TRestAxionSolarHiddenPhotonFlux::GetRandomEnergyAndRadius.
///
/// \code
/// python3 solarPlotHiddenPhoton.py --fluxname HiddenPhoton --N 1000000
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
/// ### Exporting the solar flux tables
///
/// On top of that, we will be able to export those tables to the TRestAxionSolarHiddenPhotonFlux
/// standard format to be used in later occasions.
///
/// \code
///    TRestAxionSolarHiddenPhotonFlux *sFlux = new TRestAxionSolarHiddenPhotonFlux("fluxes.rml",
///    "HiddenPhoton") sFlux->Initialize() sFlux->ExportTables()
/// \endcode
///
/// which will produce a binary table `.N200f` with the continuum flux. The filename root will be
/// extracted from the original `.flux` file. Optionally we may export the continuum flux to an
/// ASCII file by indicating it at the TRestAxionSolarHiddenPhotonFlux::ExportTables method call.
/// The files will be placed at the REST user space, at `$HOME/.rest/export/` directory.
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
/// 2023-June: TRestAxionSolarHiddenPhotonFlux created by editing TRestAxionSolarQCDFlux
///                Tomas O'Shea
///
/// \class      TRestAxionSolarHiddenPhotonFlux
/// \author     Javier Galan
///
/// <hr>
///

#include "TRestAxionSolarHiddenPhotonFlux.h"
using namespace std;

ClassImp(TRestAxionSolarHiddenPhotonFlux);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionSolarHiddenPhotonFlux::TRestAxionSolarHiddenPhotonFlux() : TRestAxionSolarFlux() {}

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
TRestAxionSolarHiddenPhotonFlux::TRestAxionSolarHiddenPhotonFlux(const char* cfgFileName, string name)
    : TRestAxionSolarFlux(cfgFileName) {
    LoadConfigFromFile(fConfigFileName, name);

    if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Info) PrintMetadata();
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionSolarHiddenPhotonFlux::~TRestAxionSolarHiddenPhotonFlux() {}

///////////////////////////////////////////////
/// \brief It will load the tables in memory by using the filename information provided
/// inside the metadata members, and calculate the solar flux for a given m.
///
Bool_t TRestAxionSolarHiddenPhotonFlux::LoadTables(Double_t mass) {
    if (fFluxDataFile == "" || fWidthDataFile == "" || fPlasmaFreqDataFile == "") return false;

    SetMass(mass);

    LoadContinuumFluxTable();
    LoadWidthTable();
    LoadPlasmaFreqTable();

    CalculateSolarFlux();

    IntegrateSolarFluxes();

    return true;
}

///////////////////////////////////////////////
/// \brief A helper method to load the data file containing continuum spectra as a
/// function of the solar radius. It will be called by TRestAxionSolarHiddenPhotonFlux::Initialize.
///
void TRestAxionSolarHiddenPhotonFlux::LoadContinuumTable() {
    if (fFluxDataFile == "") {
        RESTDebug
            << "TRestAxionSolarHiddenPhotonFlux::LoadContinuumFluxTable. No solar flux table was defined"
            << RESTendl;
        return;
    }

    string fullPathName = SearchFile((string)fFluxDataFile);

    RESTDebug << "Loading table from file : " << RESTendl;
    RESTDebug << "File : " << fullPathName << RESTendl;

    std::vector<std::vector<float>> fluxTable;

    if (!TRestTools::IsBinaryFile(fFluxDataFile)) {
        fluxTable.clear();
        RESTError << "File is not in binary format!" << RESTendl;
    }

    TRestTools::ReadBinaryTable(fullPathName, fluxTable);

    if (fluxTable.size() != 1000 && fluxTable[0].size() != 200) {
        fluxTable.clear();
        RESTError << "LoadContinuumFluxTable. The table does not contain the right number of rows or columns"
                  << RESTendl;
        RESTError << "Table will not be populated" << RESTendl;
    }

    for (unsigned int n = 0; n < fluxTable.size(); n++) {
        TH1F* h = new TH1F(Form("%s_ContinuumFluxAtRadius%d", GetName(), n), "", 200, 0, 20);
        for (unsigned int m = 0; m < fluxTable[n].size(); m++) h->SetBinContent(m + 1, fluxTable[n][m]);
        fContinuumTable.push_back(h);
    }
}

///////////////////////////////////////////////
/// \brief A helper method to load the data file containing resonance width as a
/// function of the solar radius. It will be called by TRestAxionSolarHiddenPhotonFlux::Initialize.
///
void TRestAxionSolarHiddenPhotonFlux::LoadWidthTable() {
    if (fFluxDataFile == "") {
        RESTDebug << "TRestAxionSolarHiddenPhotonFlux::LoadWidthTable. No width table was defined"
                  << RESTendl;
        return;
    }

    string fullPathName = SearchFile((string)fWidthDataFile);

    RESTDebug << "Loading table from file : " << RESTendl;
    RESTDebug << "File : " << fullPathName << RESTendl;

    std::vector<std::vector<float>> fluxTable;

    if (!TRestTools::IsBinaryFile(fWidthDataFile)) {
        fluxTable.clear();
        RESTError << "File is not in binary format!" << RESTendl;
    }

    TRestTools::ReadBinaryTable(fullPathName, fluxTable);

    if (fluxTable.size() != 1000 && fluxTable[0].size() != 200) {
        fluxTable.clear();
        RESTError << "LoadWidthTable. The table does not contain the right number of rows or columns"
                  << RESTendl;
        RESTError << "Table will not be populated" << RESTendl;
    }

    for (unsigned int n = 0; n < fluxTable.size(); n++) {
        TH1F* h = new TH1F(Form("%s_ResonanceWidthAtRadius%d", GetName(), n), "", 200, 0, 20);
        for (unsigned int m = 0; m < fluxTable[n].size(); m++) h->SetBinContent(m + 1, fluxTable[n][m]);
        fWidthTable.push_back(h);
    }
}

///////////////////////////////////////////////
/// \brief A helper method to load the data file containing resonance width as a
/// function of the solar radius. It will be called by TRestAxionSolarHiddenPhotonFlux::Initialize.
///
void TRestAxionSolarHiddenPhotonFlux::LoadPlasmaFreqTable() {
    if (fFluxDataFile == "") {
        RESTDebug
            << "TRestAxionSolarHiddenPhotonFlux::LoadPlasmaFreqTable. No plasma frequency table was defined"
            << RESTendl;
        return;
    }

    string fullPathName = SearchFile((string)fPlasmaFreqDataFile);

    RESTDebug << "Loading table from file : " << RESTendl;
    RESTDebug << "File : " << fullPathName << RESTendl;

    std::vector<std::vector<float>> fluxTable;

    if (!TRestTools::IsBinaryFile(fWidthDataFile)) {
        fluxTable.clear();
        RESTError << "File is not in binary format!" << RESTendl;
    }

    TRestTools::ReadBinaryTable(fullPathName, fluxTable);

    if (fluxTable.size() != 1000 && fluxTable[0].size() != 1) {
        fluxTable.clear();
        RESTError << "LoadWidthTable. The table does not contain the right number of rows or columns"
                  << RESTendl;
        RESTError << "Table will not be populated" << RESTendl;
    }

    for (unsigned int n = 0; n < fluxTable.size(); n++) {
        TH1F* h = new TH1F(Form("%s_PlasmaFreqAtRadius%d", GetName(), n), "", 1, 0, 20);
        for (unsigned int m = 0; m < fluxTable[n].size(); m++) h->SetBinContent(m + 1, fluxTable[n][m]);
        fPlasmaFreqTable.push_back(h);
    }
}

///////////////////////////////////////////////
/// \brief A helper method to calculate the real solar flux spectrum from the 3 tables, the
/// and the hidden photon mass for chi=1.
///
void TRestAxionSolarHiddenPhotonFlux::CalculateSolarFlux() {
    if (fMass == 0) {
        RESTError << "CalculateSolarFlux. The hidden photon mass is set to zero!" << RESTendl;
        return;
    }

    for (unsigned int n = 0; n < fluxTable.size(); n++) {
        // m4 * chi2 * wG * flux / ( (m2 - wp2)^2 + (w G)^2 )

        std::vector<float> mass2Vector(200, pow(fMass, 2));
        float wp = fPlasmaFreqTable[n].GetBinContent(1);
        std::vector<float> wp2Vector(200, pow(wp, 2));

        TH1F* hMass = new TH1D("hMass", "hMass", 200, 0, 20) TH1F* hWp =
            new TH1D("hWp", "hWp", 200, 0, 20) TH1F* hWg2 = (TH1F*)fWidthTable[n]->Clone();

        hMass->FillN(200, massVector);                 // m^2 hist
        hWp->FillN(200, wpVector);                     // wp^2 hist
        TH1F* hWg2 = fWidthTable[n] * fWidthTable[n];  // (w G)^2

        hMass->Add(hWp, -1);     // (m2 - wp2)
        hMass->Multiply(hMass);  // (m2 - wp2)^2
        hmass->Add(hWg2);        // (m2 - wp2)^2 + (w G)^2

        TH1F* h = fWidthTable[n] * fContinuumTable[n];
        h->Divide(hMass);
        h->Scale(pow(fMass, 4));

        fFluxTable.push_back(h);
    }
}

///////////////////////////////////////////////
/// \brief It builds a histogram with the continuum spectrum.
/// The flux will be expressed in cm-2 s-1 keV-1. Binned in 100eV steps.
///
TH1F* TRestAxionSolarHiddenPhotonFlux::GetContinuumSpectrum() {
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
/// \brief Same as GetContinuumSpectrum, the flux will be
/// expressed in cm-2 s-1 keV-1. Binned in 1eV steps.
///
TH1F* TRestAxionSolarHiddenPhotonFlux::GetTotalSpectrum() {
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
/// \brief A helper method to initialize the internal class data members with the
/// integrated flux for each solar ring. It will be called by TRestAxionSolarHiddenPhotonFlux::Initialize.
///
void TRestAxionSolarHiddenPhotonFlux::IntegrateSolarFluxes() {
    fTotalContinuumFlux = 0.0;
    for (unsigned int n = 0; n < fFluxTable.size(); n++) {
        fTotalContinuumFlux += fFluxTable[n]->Integral() * 0.1;  // We integrate in 100eV steps
        fFluxTableIntegrals.push_back(fTotalContinuumFlux);
    }

    for (unsigned int n = 0; n < fFluxTableIntegrals.size(); n++)
        fFluxTableIntegrals[n] /= fTotalContinuumFlux;
}

///////////////////////////////////////////////
/// \brief It returns the integrated flux at earth in cm-2 s-1 for the given energy range
///
Double_t TRestAxionSolarHiddenPhotonFlux::IntegrateFluxInRange(TVector2 eRange) {
    if (eRange.X() == -1 && eRange.Y() == -1) {
        if (GetTotalFlux() == 0) IntegrateSolarFluxes();
        return GetTotalFlux();
    }

    Double_t flux = 0;
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
std::pair<Double_t, Double_t> TRestAxionSolarHiddenPhotonFlux::GetRandomEnergyAndRadius(TVector2 eRange) {
   std::pair<Double_t, Double_t> result = {0, 0};
    if (!AreTablesLoaded()) return result;
    Double_t rnd = fRandom->Rndm();
    if (fRandom->Rndm() > fFluxRatio) {
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
    }
    return result;
}

///////////////////////////////////////////////
/// \brief It prints on screen the table that has been loaded in memory
///
void TRestAxionSolarHiddenPhotonFlux::PrintContinuumSolarTable() {
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
void TRestAxionSolarHiddenPhotonFlux::PrintIntegratedRingFlux() {
    cout << "Integrated solar flux per solar ring: " << endl;
    cout << "--------------------------- " << endl;
    /*
    for (int n = 0; n < fFluxPerRadius.size(); n++)
    cout << "n : " << n << " flux : " << fFluxPerRadius[n] << endl;
    cout << endl;
    */
}

///////////////////////////////////////////////
/// \brief Prints on screen the information about the metadata members of TRestAxionSolarHiddenPhotonFlux
///
void TRestAxionSolarHiddenPhotonFlux::PrintMetadata() {
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
        PrintIntegratedRingFlux();
    }
}

///////////////////////////////////////////////
/// \brief It will create files with spectra to be used
/// in a later ocasion.
///
void TRestAxionSolarHiddenPhotonFlux::ExportTables(Bool_t ascii) {
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
}
