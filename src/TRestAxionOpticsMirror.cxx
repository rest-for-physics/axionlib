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
/// TRestAxionOpticsMirror is a class that allows to define specific
/// mirror properties, such as the layer material and thickness, the
/// substrate material, or the surface roughness through different
/// class metadata members.
///
/// The following metadata parameters are used to define the mirror
/// properties:
/// * **mirrorType**: It defines the mirror type (Single, Thick, Bilayer,
/// Multilayer). For the moment only `Single` layer type has been implemented.
/// * **layer**: It defines the layer material. Chemical formula.
/// * **layerThickness**: It defines the layer thickness in nm.
/// * **substrate**: It defines the substrate material. Chemical formula.
/// * **sigmal1**: It defines the layer roughness in nm.
///
/// The metadata members defined will generate a request to the Henke
/// database to retrieve reflectivity and transmission values as a
/// function of the angle and energy of the incident photon. Once the
/// metadata members have been defined, it is necessary to call the
/// method TRestAxionOpticsMirror::LoadTables in order to initialize
/// the reflectivity and transmission tables.
///
/// The download process will take a reasonable amount of time, and after
/// all the data files have been downloaded, this class will generate a
/// unique table containning the reflectivity as a function of angle and
/// energy. The tables will be exported to the REST user path so that
/// in future calls these tables will be directly loaded without requiring
/// to download them again.
///
/// The following code shows how to define the mirror properties and launch
/// the drawing of few plots showing the reflectivity as a function of
/// incident angle and  energy.
///
/// \code
///     TRestAxionOpticsMirror *mirror = new TRestAxionOpticMirror();
///		mirror->SetMirrorType("Single");
///		mirror->SetLayer("SiO2");
///		mirror->SetLayerThickness("30");
///		mirror->SetSubstrateMaterial("C");
///		mirror->SetRoughness("C");
///
///     mirror->LoadTables();
///     mirror->DrawOpticsProperties();
/// \endcode
///
/// Alternatively we may use a RML definition and pass some options to the
/// TRestAxionOpticMirror::DrawOpticsProperties to define what it will be
/// drawn.
///
/// The following code will draw the reflectivity as a function of the angle
/// in the 0 to 5 degrees range at four different energies (1,4,7,10)keV, and
/// as a function of the energy in the 0 to 15 keV range at three different
/// angles (0.5,1,1.5) degrees. The lower y-axis value will be fixed to 1.e-4.
/// See the method documentation for more details.
///
/// \code
///     TRestAxionOpticsMirror *mirror = new TRestAxionOpticMirror("mirror.rml", "default");
///     mirror->DrawOpticsProperties("[1,4,7,10](0,15):[0.5,1,1.5](0,5)", 1.e-4);
/// \endcode
///
///
/// The following plots have been generated using the DrawOpticsProperties using those
/// options.
//
/// \htmlonly <style>div.image img[src="reflectivity.png"]{width:750px;}</style> \endhtmlonly
///
/// ![The reflectivity as a function of the incidence angle and the energy in keV.](reflectivity.png)
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2022-April: First concept and implementation of TRestAxionOpticsMirror class.
///            	  Javier Galan
///
/// \class      TRestAxionOpticsMirror
/// \author     Javier Galan <javier.galan@unizar.es>
///
/// <hr>
///

#include "TRestAxionOpticsMirror.h"
#include <TAxis.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TLegend.h>

using namespace std;

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
        TRestTools::ReadBinaryTable(mirrorFile, fReflectivityTable, 901);
        TRestTools::ReadBinaryTable(windowFile, fTransmissionTable, 901);
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

///////////////////////////////////////////////
/// \brief It is a private method to export the tables to a binary file once the tables
/// have been downloaded from Henke database
///
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
    TRestTools::ExportBinaryTable(path + fnameR, fReflectivityTable);

    info << "Reflectivity table generated at: " << path + fnameR << endl;

    string fnameT = fMirrorType + "_" + fLayer + "_" + fLayerThickness + "_" + fSubstrate + "_" + fSigma1 +
                    ".transmission";
    TRestTools::ExportBinaryTable(path + fnameT, fTransmissionTable);

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
Double_t TRestAxionOpticsMirror::GetReflectivity(const Double_t angle, const Double_t energy) {
    if (fReflectivityTable.size() == 0) LoadTables();

    Double_t en = energy;
    if (en < 0.030) {
        warning << "Energy is below 30eV! It should be between 30eV and 15keV" << endl;
        warning << "Setting energy to 30eV" << endl;
        en = 0.030;
    }

    if (en >= 15) {
        warning << "Energy is above 15keV! It should be between 30eV and 15keV" << endl;
        warning << "Setting energy to 15keV" << endl;
        en = 14.9999;
    }

    Int_t lowEnBin = (Int_t)((en - 0.03) / 0.03);
    Double_t deltaE = (en - (Double_t)(lowEnBin + 1) * 0.03) / 0.03;  // between 0 and 1

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

    Int_t lowAngBin = (Int_t)((ang) / 0.1);
    Double_t deltaAng = (ang - (Double_t)(lowAngBin)*0.1) / 0.1;  // between 0 and 1

    Double_t REnLowAngLow = fReflectivityTable[lowEnBin][lowAngBin];
    Double_t REnLowAngHi = fReflectivityTable[lowEnBin][lowAngBin + 1];
    Double_t REnHiAngLow = fReflectivityTable[lowEnBin + 1][lowAngBin];
    Double_t REnHiAngHi = fReflectivityTable[lowEnBin + 1][lowAngBin + 1];

    // We have renormalized the grid to unity and now we apply the equation z = f(x,y)
    // where x is associated with the energy, y is associated with the angle
    // z = f(x,y) = (1−x)(1−y) * v00+x(1−y) * v10+(1−x)y * v01+xy * v11
    // So that, for example, when x=1 and v=1 we get v11=REnHiAngHi
    return (1 - deltaE) * (1 - deltaAng) * REnLowAngLow + deltaE * (1 - deltaAng) * REnHiAngLow +
           (1 - deltaE) * deltaAng * REnLowAngHi + deltaE * deltaAng * REnHiAngHi;
}

///////////////////////////////////////////////
/// \brief It returns the interpolated transmission for a given angle (in degrees)
/// and a given energy (in keV).
///
Double_t TRestAxionOpticsMirror::GetTransmission(const Double_t angle, const Double_t energy) {
    if (fTransmissionTable.size() == 0) LoadTables();

    Double_t en = energy;
    if (en < 0.030) {
        warning << "Energy is below 30eV! It should be between 30eV and 15keV" << endl;
        warning << "Setting energy to 30eV" << endl;
        en = 0.030;
    }

    if (en >= 15) {
        warning << "Energy is above 15keV! It should be between 30eV and 15keV" << endl;
        warning << "Setting energy to 15keV" << endl;
        en = 14.9999;
    }

    Int_t lowEnBin = (Int_t)((en - 0.03) / 0.03);
    Double_t deltaE = (en - (Double_t)(lowEnBin + 1) * 0.03) / 0.03;  // between 0 and 1

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

    Int_t lowAngBin = (Int_t)((ang) / 0.01);
    Double_t deltaAng = (ang - (Double_t)(lowAngBin)*0.1) / 0.1;  // between 0 and 1

    Double_t REnLowAngLow = fTransmissionTable[lowEnBin][lowAngBin];
    Double_t REnLowAngHi = fTransmissionTable[lowEnBin][lowAngBin + 1];
    Double_t REnHiAngLow = fTransmissionTable[lowEnBin + 1][lowAngBin];
    Double_t REnHiAngHi = fTransmissionTable[lowEnBin + 1][lowAngBin + 1];

    // We have renormalized the grid to unity and now we apply the equation z = f(x,y)
    // where x is associated with the energy, y is associated with the angle
    // z = f(x,y) = (1−x)(1−y) * v00+𝑥(1−y) * v10+(1−𝑥)y * v01+xy * v11
    // So that, for example, when x=1 and v=1 we get v11=REnHiAngHi

    return (1 - deltaE) * (1 - deltaAng) * REnLowAngLow + deltaE * (1 - deltaAng) * REnHiAngLow +
           (1 - deltaE) * deltaAng * REnLowAngHi + deltaE * deltaAng * REnHiAngHi;
}

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

///////////////////////////////////////////////
/// \brief A method that creates a canvas where the mirror optics properties are drawn.
/// It generates two plots, on the left the reflectivity as a function of the angle, and
/// on the right the reflectivity as a function of the energy.
///
/// The first argument is a string where we may specify the energies and angles to be plotted
/// against the angle and energy. Between square brackets [ ] we define the values that will
/// be plotted, while between parenthesis we define the range of the x-axis to be plotted.
/// First, the energy curves to be plotted together with the energy range are given. Then,
/// the angle curves and the angular range is given.
///
/// For example the following definition will produce a plot with 3 energies (2,5,10) keV
/// as a function of the angle, in the range 0 to 45 degrees, and a second plot with 4
/// angles (1,2,5,15) degrees in the energy range 0 to 15keV.
///
/// ```
/// [2,5,10](0,15):[1,2,5,15](0,45)
/// ```
///
TCanvas* TRestAxionOpticsMirror::DrawOpticsProperties(std::string options, Double_t lowRange) {
    if (fReflectivityTable.size() == 0) LoadTables();

    std::vector<string> optList = TRestTools::GetOptions(options);

    if (optList.size() == 0) optList = TRestTools::GetOptions("[2,4,8](0,15):[0.25,0.5,1](0,5)");

    if (optList.size() != 2) {
        ferr << "TRestAxionOpticsMirror::DrawOpticsProperties. Wrong arguments!" << endl;
        return fCanvas;
    }

    std::vector<double> energies = StringToElements(optList[0], "[", ",", "]");
    std::vector<double> eRange = StringToElements(optList[0], "(", ",", ")");

    std::vector<double> angles = StringToElements(optList[1], "[", ",", "]");
    std::vector<double> aRange = StringToElements(optList[1], "(", ",", ")");

    if (eRange[0] < 0.03) eRange[0] = 0.03;
    if (eRange[1] > 15) eRange[1] = 15;
    if (aRange[0] < 0.0) aRange[0] = 0.0;
    if (aRange[1] > 90) aRange[1] = 90;

    //   Double_t lowReflec = TRestTools::GetMinValueFromTable(fReflectivityTable);
    //   Double_t highReflec = TRestTools::GetMaxValueFromTable(fReflectivityTable);

    if (fCanvas != NULL) {
        delete fCanvas;
        fCanvas = NULL;
    }
    fCanvas = new TCanvas("canv", "", 1400, 600);

    TPad* pad1 = new TPad("pad1", "This is pad1", 0.01, 0.02, 0.99, 0.97);
    pad1->Divide(2, 1);
    pad1->Draw();

    ////// Drawing reflectivity versus angle
    pad1->cd(1);
    pad1->cd(1)->SetLogy();
    pad1->cd(1)->SetRightMargin(0.09);
    pad1->cd(1)->SetLeftMargin(0.15);
    pad1->cd(1)->SetBottomMargin(0.15);

    std::vector<TGraph*> ref_vs_ang_graph;

    for (int n = 0; n < energies.size(); n++) {
        string grname = "gr" + IntegerToString(n);
        TGraph* gr = new TGraph();
        gr->SetName(grname.c_str());
        for (double a = aRange[0]; a <= aRange[1]; a += (aRange[1] - aRange[0]) / 100.) {
            gr->AddPoint(a, GetReflectivity(a, energies[n]));
        }
        gr->SetLineColor(49 - n * 3);
        gr->SetLineWidth(3);
        ref_vs_ang_graph.push_back(gr);
    }

    ref_vs_ang_graph[0]->GetXaxis()->SetLimits(aRange[0], aRange[1]);
    ref_vs_ang_graph[0]->GetHistogram()->SetMaximum(1);
    ref_vs_ang_graph[0]->GetHistogram()->SetMinimum(lowRange);

    ref_vs_ang_graph[0]->GetXaxis()->SetTitle("Angle [degrees]");
    ref_vs_ang_graph[0]->GetXaxis()->SetTitleSize(0.05);
    ref_vs_ang_graph[0]->GetXaxis()->SetLabelSize(0.05);
    ref_vs_ang_graph[0]->GetYaxis()->SetTitle("Reflectivity");
    ref_vs_ang_graph[0]->GetYaxis()->SetTitleOffset(1.5);
    ref_vs_ang_graph[0]->GetYaxis()->SetTitleSize(0.05);
    ref_vs_ang_graph[0]->GetYaxis()->SetLabelSize(0.05);
    pad1->cd(1)->SetLogy();
    ref_vs_ang_graph[0]->Draw("AL");
    for (int n = 1; n < energies.size(); n++) ref_vs_ang_graph[n]->Draw("L");

    TLegend* legend = new TLegend(0.6, 0.75, 0.9, 0.95);
    legend->SetTextSize(0.03);
    legend->SetHeader("Energies", "C");  // option "C" allows to center the header
    for (int n = 0; n < energies.size(); n++) {
        std::string lname = "gr" + IntegerToString(n);
        std::string ltitle = DoubleToString(energies[n]) + " keV";

        legend->AddEntry(lname.c_str(), ltitle.c_str(), "l");
    }
    legend->Draw();

    ////// Drawing reflectivity versus energy
    pad1->cd(2);
    pad1->cd(2)->SetLogy();
    pad1->cd(2)->SetRightMargin(0.09);
    pad1->cd(2)->SetLeftMargin(0.15);
    pad1->cd(2)->SetBottomMargin(0.15);

    std::vector<TGraph*> ref_vs_en_graph;

    for (int n = 0; n < angles.size(); n++) {
        string grname = "agr" + IntegerToString(n);
        TGraph* gr = new TGraph();
        gr->SetName(grname.c_str());
        for (double e = eRange[0]; e <= eRange[1]; e += (eRange[1] - eRange[0]) / 100.) {
            gr->AddPoint(e, GetReflectivity(angles[n], e));
        }
        gr->SetLineColor(49 - n * 3);
        gr->SetLineWidth(3);
        ref_vs_en_graph.push_back(gr);
    }

    ref_vs_en_graph[0]->GetXaxis()->SetLimits(eRange[0], eRange[1]);
    ref_vs_en_graph[0]->GetHistogram()->SetMaximum(1);
    ref_vs_en_graph[0]->GetHistogram()->SetMinimum(lowRange);

    ref_vs_en_graph[0]->GetXaxis()->SetTitle("Energy [keV]");
    ref_vs_en_graph[0]->GetXaxis()->SetTitleSize(0.05);
    ref_vs_en_graph[0]->GetXaxis()->SetLabelSize(0.05);
    ref_vs_en_graph[0]->GetYaxis()->SetTitle("Reflectivity");
    ref_vs_en_graph[0]->GetYaxis()->SetTitleOffset(1.5);
    ref_vs_en_graph[0]->GetYaxis()->SetTitleSize(0.05);
    ref_vs_en_graph[0]->GetYaxis()->SetLabelSize(0.05);
    pad1->cd(2)->SetLogy();
    ref_vs_en_graph[0]->Draw("AL");
    for (int n = 1; n < angles.size(); n++) ref_vs_en_graph[n]->Draw("L");

    TLegend* legendA = new TLegend(0.6, 0.75, 0.9, 0.95);
    legendA->SetTextSize(0.03);
    legendA->SetHeader("Angles", "C");  // option "C" allows to center the header
    for (int n = 0; n < angles.size(); n++) {
        std::string lname = "agr" + IntegerToString(n);
        std::string ltitle = DoubleToString(angles[n]) + " degrees";

        legendA->AddEntry(lname.c_str(), ltitle.c_str(), "l");
    }
    legendA->Draw();

    return fCanvas;
}
