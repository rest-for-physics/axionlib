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
/// TRestAxionOpticsMirror is a class that allows to load externally
/// defined optics response files. This metadata class will be a generic,
/// abstract, class that will be inherited by other more specific metadata
/// classes. This class will define few common metadata members helping to
/// describe the optics alignment, position, and basic geometry specifications,
/// such as number of mirror rings, or additional entrance masks, such as
/// spider mask.
///
/// The following metadata parameters define the optics position, size and
/// alignment:
/// * **center**: It defines the center of the optics, entrance and exit
/// optics planes will be defined using the half lenght and the `center`
/// position.
/// * **axis**: It defines the optical axis direction.
/// * **length**: It defines the size of the optics, used to calculate
/// the optics plane entrance, and the optics plane exit.
///
/// The following image is generated as a validation or a way to visualize the
/// TRestAxionOpticsMirror::GetEntranceRing method. Each color represents a particle
/// hitting in a different ring. The position is drawn at both, the generator
/// plane and the optics entrance plane. This creates an effect of diffusion at
/// the generator plane since the generator random direction is slightly tilted
/// respect to the optical axis.
///
/// \htmlonly <style>div.image img[src="xyz.png"]{width:750px;}</style> \endhtmlonly
///
/// ![Image description](xyz.png)
///
/// This image was generated using the xyz script.
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2022-February: First concept and implementation of TRestAxionOpticsMirror class.
///            	  Javier Galan
///
/// \class      TRestAxionOpticsMirror
/// \author     Javier Galan <javier.galan@unizar.es>
///
/// <hr>
///

#include "TRestAxionOpticsMirror.h"

using namespace std;

#include "TRestPhysics.h"
using namespace REST_Physics;

ClassImp(TRestAxionOpticsMirror);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionOpticsMirror::TRestAxionOpticsMirror() : TRestMetadata() { Initialize(); }

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
TRestAxionOpticsMirror::TRestAxionOpticsMirror(const char* cfgFileName, string name)
    : TRestMetadata(cfgFileName) {
    debug << "Entering TRestAxionOpticsMirror constructor( cfgFileName, name )" << endl;

    Initialize();

    LoadConfigFromFile(fConfigFileName, name);

    if (GetVerboseLevel() >= REST_Info) PrintMetadata();
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionOpticsMirror::~TRestAxionOpticsMirror() {}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionOpticsMirror members
///
void TRestAxionOpticsMirror::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    fHenkeKeys.clear();
    fHenkeKeys["Layer"] = fLayer;
    fHenkeKeys["Ldensity"] = "-1";
    fHenkeKeys["Thick"] = fLayerThickness;
    fHenkeKeys["Sigma1"] = fSigma1;
    fHenkeKeys["Substrate"] = fSubstrate;
    fHenkeKeys["Sdensity"] = "-1";
    fHenkeKeys["Sigma2"] = "0";
    fHenkeKeys["Pol"] = "0";
    fHenkeKeys["Scan"] = "Angle";
    fHenkeKeys["Min"] = "0";
    fHenkeKeys["Max"] = "90";
    fHenkeKeys["Npts"] = "90";
    fHenkeKeys["temp"] = "Energy+%28eV%29";
    fHenkeKeys["Fixed"] = "200";
    fHenkeKeys["Plot"] = "LinLog";
    fHenkeKeys["Output"] = "Plot";
}

///////////////////////////////////////////////
/// \brief Loads the reflectivity table
///
void TRestAxionOpticsMirror::LoadTables() {
    if (fHenkeKeys.size() == 0) Initialize();

    string mirrorFile = SearchFile(fMirrorType + "_" + fLayer + "_" + fLayerThickness + "_" + fSubstrate +
                                   "_" + fSigma1 + ".reflectivity");
    string windowFile = SearchFile(fMirrorType + "_" + fLayer + "_" + fLayerThickness + "_" + fSubstrate +
                                   "_" + fSigma1 + ".transmission");

    if (mirrorFile != "" && windowFile != "") {
        TRestTools::ReadASCIITable(mirrorFile, fReflectivityTable);
        TRestTools::ReadASCIITable(windowFile, fTransmissionTable);
        return;
    }

    cout << "--------------------------------- INFO ------------------------------" << endl;
    cout << "The optics material properties were not found in the database" << endl;
    cout << "Downloading from Henke database. This process will take a long time." << endl;
    cout << "After producing these files, please, contribute them to the " << endl;
    cout << "rest-for-physics/axionlib-data inside the optics directory ..." << endl;
    cout << "---------------------------------- INFO ------------------------------" << endl;

    fReflectivityTable.clear();
    fTransmissionTable.clear();
    map<Int_t, std::vector<Float_t>> reflectivity;
    map<Int_t, std::vector<Float_t>> transmission;
    for (int n = 0; n < 90; n += 9) {
        fHenkeKeys["Min"] = IntegerToString(n);
        fHenkeKeys["Max"] = IntegerToString(n + 9);

        cout.flush();
        vector<Float_t> reflect;
        vector<Float_t> transm;
        reflect.clear();
        transm.clear();

        cout << "Reading angles between " << n << " and " << n + 9 << ".";
        for (int e = 30; e <= 15000; e += 30) {
            fHenkeKeys["Fixed"] = IntegerToString(e);
            cout << ".";
            cout.flush();
            string fname = DownloadHenkeFile();

            std::vector<std::vector<Float_t>> data;
            TRestTools::ReadASCIITable(fname, data, 2);

            // we skip the last point if we are not at the latest angles file
            Int_t N = data.size() - 1;
            if (n == 81) N = N + 1;

            for (int m = 0; m < N; m++) reflectivity[e].push_back(data[m][1]);
            for (int m = 0; m < N; m++) transmission[e].push_back(data[m][2]);
        }
        cout << endl;
    }
    cout << "Number of angles stored : " << reflectivity[30].size() << endl;

    for (const auto& x : reflectivity) {
        fReflectivityTable.push_back(x.second);
    }
    for (const auto& x : transmission) {
        fTransmissionTable.push_back(x.second);
    }

    ExportTables();
}

Int_t TRestAxionOpticsMirror::ExportTables() {
    if (fReflectivityTable.size() == 0) {
        ferr << "Nothing to export!" << endl;
        return 1;
    }

    string path = REST_USER_PATH + "/export/";

    if (!TRestTools::fileExists(path)) {
        cout << "Creating path: " << path << endl;
        system(("mkdir -p " + path).c_str());
    }

    string fnameR = fMirrorType + "_" + fLayer + "_" + fLayerThickness + "_" + fSubstrate + "_" + fSigma1 +
                    ".reflectivity";
    TRestTools::ExportASCIITable(path + fnameR, fReflectivityTable);

    info << "Reflectivity table generated at: " << path + fnameR << endl;

    string fnameT = fMirrorType + "_" + fLayer + "_" + fLayerThickness + "_" + fSubstrate + "_" + fSigma1 +
                    ".transmission";
    TRestTools::ExportASCIITable(path + fnameT, fTransmissionTable);

    info << "Transmission table generated at: " << path + fnameT << endl;

    return 0;
}

///////////////////////////////////////////////
/// \brief It downloads the reflectivity file for the present mirror properties
/// defined at the metadata members.
///
/// \return It returns the location and filename of the downloaded file.
///
std::string TRestAxionOpticsMirror::DownloadHenkeFile() {
    string url = "https://henke.lbl.gov/cgi-bin/laymir.pl";
    string result = TRestTools::POSTRequest(url, fHenkeKeys);

    size_t start = result.find("HREF=\"") + 6;
    size_t length = result.find(".dat\">") + 4 - start;
    result = result.substr(start, length);

    return TRestTools::DownloadRemoteFile("https://henke.lbl.gov/" + result);
}

///////////////////////////////////////////////
/// \brief It returns the interpolated reflectivity for a given angle (in degrees)
/// and a given energy (in keV).
///
Double_t GetReflectivity(const Double_t angle, const Double_t energy) {
    Double_t en = energy;
    if (en < 0.030) {
        warning << "Energy is below 30eV! It should be between 30eV and 15keV" << endl;
        warning << "Setting energy to 30eV" << endl;
        en = 0.030;
    }

    if (en > 15) {
        warning << "Energy is above 15keV! It should be between 30eV and 15keV" << endl;
        warning << "Setting energy to 15keV" << endl;
        en = 15;
    }

    Double_t ang = angle;
    if (ang < 0.0) {
        warning << "Angle is below 0 degrees! It should be between 0 and 90 degrees" << endl;
        warning << "Setting angle to 0 degrees" << endl;
        ang = 0.0;
    }

    if (ang > 90) {
        warning << "Angle is above 90 degrees! It should be between 0 and 90 degrees" << endl;
        warning << "Setting angle to 90 degrees" << endl;
        ang = 90;
    }

    return energy;
}

///////////////////////////////////////////////
/// \brief
///
Double_t GetTransmission(const Double_t angle, const Double_t energy) { return 0.; }

///////////////////////////////////////////////
/// \brief Prints on screen the information about the metadata members of TRestAxionOpticsMirror
///
void TRestAxionOpticsMirror::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    metadata << "Mirror type: " << fMirrorType << endl;
    metadata << "Layer material: " << fLayer << endl;
    metadata << "Layer thickness: " << fLayerThickness << " nm" << endl;
    metadata << "Substrate material: " << fSubstrate << endl;
    metadata << "Roughness: " << fSigma1 << "nm" << endl;
    metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}
