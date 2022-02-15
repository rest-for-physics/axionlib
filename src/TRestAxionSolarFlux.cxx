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
/// See for example the method TRestAxionSolarFlux::InitializeSolarTable using
/// a solar model description by TRestAxionSolarModel.
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
TRestAxionSolarFlux::TRestAxionSolarFlux() : TRestMetadata() { Initialize(); }

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

    Initialize();

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

    LoadContinuumFluxTable();

    LoadMonoChromaticFluxTable();

    IntegrateSolarRingsFluxes();
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

    TRestTools::ReadASCIITable(fullPathName, fFluxTable);

    if (fFluxTable.size() != 100 && fFluxTable[0].size() != 200) {
        fFluxTable.clear();
        ferr << "LoadContinuumFluxTable. The table does not contain the right number of rows or columns"
             << endl;
        ferr << "Table will not be populated" << endl;
    }
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

    std::vector<std::vector<Double_t>> asciiTable;
    TRestTools::ReadASCIITable(fullPathName, asciiTable);

    fFluxLines.clear();

    if (asciiTable.size() != 101) {
        ferr << "LoadMonoChromaticFluxTable. The table does not contain the right number of rows" << endl;
        ferr << "Table will not be populated" << endl;
        return;
    }

    for (int en = 0; en < asciiTable[0].size(); en++) {
        Double_t energy = asciiTable[0][en];
        std::vector<Double_t> profile;
        for (int r = 1; r < asciiTable.size(); r++) profile.push_back(asciiTable[r][en]);
        fFluxLines[energy] = profile;
    }
}

///////////////////////////////////////////////
/// \brief A helper method to initialize the internal class data members with the
/// integrated flux for each solar ring. It will be called by TRestAxionSolarFlux::Initialize.
///
void TRestAxionSolarFlux::IntegrateSolarRingsFluxes() {
    fFluxPerRadius.clear();
    for (int n = 0; n < fFluxTable.size(); n++) {
        double sum = 0;
        for (int m = 0; m < fFluxTable[n].size(); m++) sum += fFluxTable[n][m];
        for (const auto& line : fFluxLines) sum += line.second[n];
        fFluxPerRadius.push_back(sum);
    }
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

    this->Initialize();
}

///////////////////////////////////////////////
/// \brief It prints on screen the table that has been loaded in memory
///
void TRestAxionSolarFlux::PrintContinuumSolarTable() {
    cout << "Continuum solar flux table: " << endl;
    cout << "--------------------------- " << endl;
    for (int n = 0; n < fFluxTable.size(); n++) {
        for (int m = 0; m < fFluxTable[n].size(); m++) cout << fFluxTable[n][m] << "\t";
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
    for (int n = 0; n < fFluxPerRadius.size(); n++)
        cout << "n : " << n << " flux : " << fFluxPerRadius[n] << endl;
    cout << endl;
}

///////////////////////////////////////////////
/// \brief It prints on screen the spectral lines loaded in memory
///
void TRestAxionSolarFlux::PrintMonoChromaticFlux() {
    cout << "Number of monochromatic lines: " << fFluxPerRadius.size() << endl;
    cout << "+++++++++++++++++++++++++++++++++++" << endl;
    for (auto const& line : fFluxLines) {
        cout << "Energy : " << line.first << " keV" << endl;
        cout << "-----------------" << endl;
        for (int n = 0; n < line.second.size(); n++)
            cout << "R : " << n << " flux : " << line.second[n] << " cm-2 s-1" << endl;
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
    metadata << "++++++++++++++++++" << endl;

    if (GetVerboseLevel() >= REST_Debug) {
        PrintContinuumSolarTable();
        PrintMonoChromaticFlux();
        PrintIntegratedRingFlux();
    }
}
