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
/// for the range (0,20)keV in steps of 100eV.
/// - *fluxSptFile:* A table where each column represents a monochromatic energy. The
/// column contains 101 rows, the first element is the energy of the monochromatic line
/// while the next 100 elements contain the flux, measured in cm-2 s-1, integrated to
/// each solar ring, being the second element the ring in the center of the sun.
///
/// The following code shows how to define this class inside a RML file.
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
/// The previous definition was used to generate the following figure using the script
/// pipeline/solarFlux/solarFlux.py.
///
/// \htmlonly <style>div.image img[src="AxionSolarFlux.png"]{width:750px;}</style> \endhtmlonly
///
/// ![Solar flux distributions generated with TRestAxionSolarFlux.](AxionSolarFlux.png)
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
    cout << "Entering TRestAxionSolarFlux constructor( cfgFileName, name )" << endl;

    LoadConfigFromFile(fConfigFileName, name);

    if (GetVerboseLevel() >= REST_Info) PrintMetadata();
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionSolarFlux::~TRestAxionSolarFlux() {}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionSolarFlux members
///
void TRestAxionSolarFlux::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    if (fFluxDataFile != "" && TRestTools::GetFileNameExtension(fFluxDataFile) == ".flux") {
        ReadFluxFile();
    } else {
        LoadContinuumFluxTable();
        LoadMonoChromaticFluxTable();
    }

    IntegrateSolarFluxes();

    if (!fRandom) {
        delete fRandom;
        fRandom = nullptr;
    }

    fRandom = new TRandom3(fSeed);
    if (fSeed == 0) fSeed = fRandom->GetSeed();
}

///////////////////////////////////////////////
/// \brief A helper method to load the data file containning continuum spectra as a
/// function of the solar radius. It will be called by TRestAxionSolarFlux::Initialize.
///
void TRestAxionSolarFlux::LoadContinuumFluxTable() {
    if (fFluxDataFile == "") {
        debug << "TRestAxionSolarflux::LoadContinuumFluxTable. No solar flux table was defined" << endl;
        return;
    }

    string fullPathName = SearchFile((string)fFluxDataFile);

    debug << "Loading table from file : " << endl;
    debug << "File : " << fullPathName << endl;

    std::vector<std::vector<Float_t>> fluxTable;
    if (TRestTools::GetFileNameExtension(fFluxDataFile) == ".dat")
        TRestTools::ReadASCIITable(fullPathName, fluxTable);
    else if (TRestTools::IsBinaryFile(fFluxDataFile))
        TRestTools::ReadBinaryTable(fullPathName, fluxTable);
    else {
        fluxTable.clear();
        ferr << "Filename extension was not recognized!" << endl;
        ferr << "Solar flux table will not be populated" << endl;
    }

    if (fluxTable.size() != 100 && fluxTable[0].size() != 200) {
        fluxTable.clear();
        ferr << "LoadContinuumFluxTable. The table does not contain the right number of rows or columns"
             << endl;
        ferr << "Solar flux table will not be populated" << endl;
    }

    for (int n = 0; n < fluxTable.size(); n++) {
        TH1F* h = new TH1F(Form("%s_ContinuumFluxAtRadius%d", GetName(), n), "", 200, 0, 20);
        for (int m = 0; m < fluxTable[n].size(); m++) h->SetBinContent(m + 1, fluxTable[n][m]);
        fFluxTable.push_back(h);
    }
}

///////////////////////////////////////////////
/// \brief It loads a .flux file. It will split continuum and monochromatic peaks.
///
void TRestAxionSolarFlux::ReadFluxFile() {
    if (fFluxDataFile == "") {
        debug << "TRestAxionSolarflux::LoadContinuumFluxTable. No solar flux table was defined" << endl;
        return;
    }

    string fullPathName = SearchFile((string)fFluxDataFile);

    debug << "Loading table from file : " << endl;
    debug << "File : " << fullPathName << endl;

    /*
std::vector<std::vector<Float_t>> fluxTable;
std::vector<std::vector<Float_t>> sptTable;

TRestTools::ReadASCIITable(fname, fluxTable);
    */
}

///////////////////////////////////////////////
/// \brief A helper method to load the data file containning monochromatic spectral
/// lines as a function of the solar radius. It will be called by TRestAxionSolarFlux::Initialize.
///
void TRestAxionSolarFlux::LoadMonoChromaticFluxTable() {
    if (fFluxSptFile == "") {
        debug << "TRestAxionSolarflux::LoadMonoChromaticFluxTable. No solar flux monochromatic table was "
                 "defined"
              << endl;
        return;
    }

    string fullPathName = SearchFile((string)fFluxSptFile);

    debug << "Loading monochromatic lines from file : " << endl;
    debug << "File : " << fullPathName << endl;

    std::vector<std::vector<Float_t>> asciiTable;
    TRestTools::ReadASCIITable(fullPathName, asciiTable);

    fFluxLines.clear();

    if (asciiTable.size() != 101) {
        ferr << "LoadMonoChromaticFluxTable. The table does not contain the right number of rows" << endl;
        ferr << "Table will not be populated" << endl;
        return;
    }

    for (int en = 0; en < asciiTable[0].size(); en++) {
        Float_t energy = asciiTable[0][en];
        std::vector<Float_t> profile;
        TH1F* h = new TH1F(Form("%s_MonochromeFluxAtEnergy%4.2lf", GetName(), energy), "", 100, 0, 1);
        for (int r = 1; r < asciiTable.size(); r++) h->SetBinContent(r, asciiTable[r][en]);
        fFluxLines[energy] = h;
    }
}

///////////////////////////////////////////////
/// \brief A helper method to initialize the internal class data members with the
/// integrated flux for each solar ring. It will be called by TRestAxionSolarFlux::Initialize.
///
void TRestAxionSolarFlux::IntegrateSolarFluxes() {
    fFluxLineIntegrals.clear();
    fTotalMonochromaticFlux = 0;
    for (const auto& line : fFluxLines) {
        fTotalMonochromaticFlux += line.second->Integral();
        fFluxLineIntegrals.push_back(fTotalMonochromaticFlux);
    }

    for (int n = 0; n < fFluxLineIntegrals.size(); n++) fFluxLineIntegrals[n] /= fTotalMonochromaticFlux;

    fTotalContinuumFlux = 0.0;
    for (int n = 0; n < fFluxTable.size(); n++) {
        fTotalContinuumFlux += fFluxTable[n]->Integral() * 0.1;  // We integrate in 100eV steps
        fFluxTableIntegrals.push_back(fTotalContinuumFlux);
    }

    for (int n = 0; n < fFluxTableIntegrals.size(); n++) fFluxTableIntegrals[n] /= fTotalContinuumFlux;

    fFluxRatio = fTotalMonochromaticFlux / (fTotalContinuumFlux + fTotalMonochromaticFlux);
}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionSolarFlux metadata members through a RML file
///
void TRestAxionSolarFlux::InitFromConfigFile() {
    debug << "Entering TRestAxionSolarFlux::InitFromConfigFile" << endl;

    fFluxDataFile = GetParameter("fluxDataFile", "");
    fFluxSptFile = GetParameter("fluxSptFile", "");
    fCouplingType = GetParameter("couplingType", "g_ag");
    fCouplingStrength = StringToDouble(GetParameter("couplingStrength", "1.e-10"));
    fSeed = StringToInteger(GetParameter("seed", "0"));

    this->Initialize();
}

///////////////////////////////////////////////
/// \brief It returns a random solar radius position and energy according to the
/// flux distributions defined inside the solar tables loaded in the class
///
std::pair<Double_t, Double_t> TRestAxionSolarFlux::GetRandomEnergyAndRadius() {
    std::pair<Double_t, Double_t> result = {0, 0};
    Double_t rnd = fRandom->Rndm();
    if (fTotalMonochromaticFlux == 0 || fRandom->Rndm() > fFluxRatio) {
        // Continuum
        for (int r = 0; r < fFluxTableIntegrals.size(); r++) {
            if (rnd < fFluxTableIntegrals[r]) {
                Double_t energy = fFluxTable[r]->GetRandom();
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
void TRestAxionSolarFlux::PrintContinuumSolarTable() {
    cout << "Continuum solar flux table: " << endl;
    cout << "--------------------------- " << endl;
    for (int n = 0; n < fFluxTable.size(); n++) {
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
void TRestAxionSolarFlux::PrintIntegratedRingFlux() {
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
void TRestAxionSolarFlux::PrintMonoChromaticFlux() {
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
/// \brief Prints on screen the information about the metadata members of TRestAxionSolarFlux
///
void TRestAxionSolarFlux::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    if (fFluxDataFile != "")
        metadata << " - Solar axion flux datafile (continuum) : " << fFluxDataFile << endl;
    if (fFluxSptFile != "")
        metadata << " - Solar axion flux datafile (monochromatic) : " << fFluxSptFile << endl;
    metadata << " - Coupling type : " << fCouplingType << endl;
    metadata << " - Coupling strength : " << fCouplingStrength << endl;
    metadata << "-------" << endl;
    metadata << " - Total monochromatic flux : " << fTotalMonochromaticFlux << " cm-2 s-1" << endl;
    metadata << " - Total continuum flux : " << fTotalContinuumFlux << " cm-2 s-1" << endl;
    metadata << "--------" << endl;
    metadata << " - Random seed : " << fSeed << endl;
    metadata << "++++++++++++++++++" << endl;

    if (GetVerboseLevel() >= REST_Debug) {
        PrintContinuumSolarTable();
        PrintMonoChromaticFlux();
        PrintIntegratedRingFlux();
    }
}
