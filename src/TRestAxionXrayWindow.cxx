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
/// TRestAxionXrayWindow implements parameters that define the window
/// properties, such as material and thickness. This class will load the
/// transmission data to calculate the transmission for a photon of a
/// given energy and position. The window might be defined as a uniform
/// foil or using a particular structure.
///
/// For the moment, the window geometry is fixed to be a circular window.
/// Therefore, the following are the most basic parameters of any window
/// construction:
///
/// * **material**: The material used by the window. The name must match
/// the name of one of the files found at `data/transmission/` removing
/// the extension.
/// * **thickness**: The depth of the window that defines finally the total
/// x-ray photon absorption.
/// * **radius**: The radius of the window where photons will be accepted.
///
/// On top of that, we might define different types of windows that will
/// allow us to create different structures, such as a strong back where
/// photon opacity is only defined at a regular grid or stripped pattern.
/// We may define the window type using:
///
/// - **type**: We may define different options that define the regions
/// of space where the window material will be effective.
///   -# **foil**: It defines a homogeneous region of material, all the
///    region inside the window radius will present the same efficiency.
///   -# **stripped**: It defines a stripped pattern where the absorption
///   will be effective, otherwise, if the photon is inside the window
///   but it doesn't hit the stripped pattern, the transmission will be
///   equal to 1.
///   -# **grid**: It defines a grid pattern where the absorption
///   will be effective, otherwise, if the photon is inside the window
///   but it doesn't hit the grid pattern, the transmission will be
///   equal to 1.
///
/// In the case that the window type is defined to be *stripped* or
/// *grid*. Then, few additional parameters are necessary to define the
/// pattern structure, made of masking strips on one direction (stripped),
/// or masking strips on both directions (grid).
///
/// * **patternGap**: The distance between 2 masking structures, or
/// periodicity.
/// * **patternWidth**: The width of the masking structure.
/// * **patternOffset**: If this parameter is 0, the first strip will
/// centered at the origin. If not, the pattern will be shifted. In the
/// case of the grid, it will be used to displace the grid on both
/// directions, X and Y.
///
/// Once all the parameters have been defined inside an instance of this
/// class, we will be able to recover the transmission at any given point
/// inside the window using the method TRestAxionXrayWindow::GetEfficiecy.
///
/// The corresponding RML section for initialization through a configuration
/// file would be as follows.
///
/// \code
///	<TRestAxionXrayWindow name="windowTest" verboseLevel="warning" >
///		<parameter name="center" value="(0,0,0)mm" />
///
///		<parameter name="type" value="foil" />
///		<parameter name="thickness" value="0.3um" />
///		<parameter name="material" value="Si3N4" />
///		<parameter name="radius" value="8mm" />
///	</TRestAxionXrayWindow>
/// \endcode
///
/// The pipeline example found at `pipeline/transmission/windowPlot.py` will
/// use a definition with 3 layers to generate the following plot.
///
/// \htmlonly <style>div.image img[src="windowsTransmission.png"]{width:800px;}</style> \endhtmlonly
///
/// ![Windows trasmission at energy ranges (0-5)keV, (5-10)keV and (10-15)keV.](windowsTransmission.png)
///
///--------------------------------------------------------------------------
///
/// RESTsoft - Software for Rare Event Searches with TPCs
///
/// History of developments:
///
/// 2019-March: First concept and implementation of TRestAxionXrayWindow class.
///             Javier Galan
///
/// \class      TRestAxionXrayWindow
/// \author     Javier Galan
///
/// <hr>
///

#include "TRestAxionXrayWindow.h"
using namespace std;

#include "TRestSystemOfUnits.h"
using namespace REST_Units;

ClassImp(TRestAxionXrayWindow);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionXrayWindow::TRestAxionXrayWindow() : TRestMetadata() { Initialize(); }

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
/// corresponding TRestGeant4Metadata section inside the RML.
///
TRestAxionXrayWindow::TRestAxionXrayWindow(const char* cfgFileName, string name)
    : TRestMetadata(cfgFileName) {
    RESTDebug << "Entering TRestAxionXrayWindow constructor( cfgFileName, name )" << RESTendl;

    Initialize();

    LoadConfigFromFile(fConfigFileName, name);
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionXrayWindow::~TRestAxionXrayWindow() {}

///////////////////////////////////////////////
/// \brief Initialization of TRestAxionXrayWindow members. It removes all gases.
///
void TRestAxionXrayWindow::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    fEnergy.clear();
    fTransmission.clear();

    if (fWindowType != "foil" && fPatternGap == 0) {
        RESTError << "TRestAxionXrayWindow::Initialize. fPatternGap cannot be zero!" << RESTendl;
        fPatternGap = 1;
    }

    if (fPatternGap < 0) fPatternGap = -fPatternGap;
}

///////////////////////////////////////////////
/// \brief It reads the data files from the corresponding material that
/// needs to be found in the axiolib database with .sol extension.
/// Usually placed under `data/axion/transmission/`
///
void TRestAxionXrayWindow::ReadMaterial() {
    std::string materialFileName = SearchFile(fMaterial + ".sol");

    RESTDebug << "TRestAxionXrayWindow::ReadMaterial. Reading file : " << materialFileName << RESTendl;

    if (!TRestTools::fileExists(materialFileName)) {
        RESTError << "TRestAxionXrayWindow::ReadMaterial( )" << RESTendl;
        RESTError << "Material file not found : " << materialFileName << RESTendl;
        exit(1);
    }

    FILE* fin = fopen(materialFileName.c_str(), "rt");

    fEnergy.clear();
    fTransmission.clear();
    double en, value;
    while (fscanf(fin, "%lf\t%lf\n", &en, &value) != EOF) {
        RESTDebug << "Energy : " << en << "eV -- Abs : " << value << RESTendl;

        fEnergy.push_back(en / 1000.);
        fTransmission.push_back(TMath::Power(value, fThickness * units("um")));
    }

    RESTDebug << "Items read : " << fEnergy.size() << RESTendl;

    fclose(fin);
}

///////////////////////////////////////////////
/// \brief It returns the window transmission probability for the given energy (in keV)
/// and window position, using energy linear interpolation.
///
/// For the case of patterned window (stripped or grid), it will return 1 if the strip is
/// not hitted.
///
Double_t TRestAxionXrayWindow::GetTransmission(Double_t energy, Double_t x, Double_t y) {
    if (fEnergy.size() == 0) ReadMaterial();

    if ((x - fCenter.X()) * (x - fCenter.X()) + (y - fCenter.Y()) * (y - fCenter.Y()) > fRadius * fRadius)
        return 0;

    if (!HitsPattern(x, y)) return 1.;

    Double_t energyIndex = GetEnergyIndex(energy);

    if (energyIndex < 0) {
        RESTWarning << "Energy : " << energy << " keV is out of range!" << RESTendl;
        return 0;
    }

    // Transmission
    double y2 = fTransmission[energyIndex + 1];
    double y1 = fTransmission[energyIndex];

    // Normalized field
    double x2 = fEnergy[energyIndex + 1];
    double x1 = fEnergy[energyIndex];

    double m = (y2 - y1) / (x2 - x1);
    double n = y1 - m * x1;

    if (m * energy + n < 0) {
        RESTError << "TRestAxionXrayWindow::GetAbsorptionCoefficient. Negative coefficient!" << RESTendl;
        cout << "y2 : " << y2 << " y1 : " << y1 << endl;
        cout << "x2 : " << x2 << " x1 : " << x1 << endl;
        cout << "m : " << m << " n : " << n << endl;
        cout << "E : " << energy << " bin : " << energyIndex << endl;
        return 0;
    }

    return (m * energy + n);
}

///////////////////////////////////////////////
/// \brief It returns true if the window pattern is hitted. False otherwise.
///
Bool_t TRestAxionXrayWindow::HitsPattern(Double_t x, Double_t y) {
    if (fWindowType == "stripped") {
        Double_t xEval = fPatternWidth / 2. + x - fPatternOffset;

        if (xEval > 0) {
            while (xEval > fPatternGap) xEval -= fPatternGap;
        } else {
            while (xEval < 0) xEval += fPatternGap;
        }

        if (xEval > fPatternWidth) {
            return false;
        }
    } else if (fWindowType == "grid") {
        Double_t xEval = fPatternWidth / 2. + x - fPatternOffset;

        if (xEval > 0) {
            while (xEval > fPatternGap) xEval -= fPatternGap;
        } else {
            while (xEval < 0) xEval += fPatternGap;
        }

        if (xEval < fPatternWidth) return true;

        Double_t yEval = fPatternWidth / 2. + y - fPatternOffset;

        if (yEval > 0) {
            while (yEval > fPatternGap) yEval -= fPatternGap;
        } else {
            while (yEval < 0) yEval += fPatternGap;
        }

        if (yEval < fPatternWidth) return true;
        return false;
    }

    return true;
}

///////////////////////////////////////////////
/// \brief It returns the vector element index, from `fEnergy`, that is just below the given input energy.
///
Int_t TRestAxionXrayWindow::GetEnergyIndex(Double_t energy) {
    for (int n = 0; n < fEnergy.size(); n++)
        if (energy < fEnergy[n]) return n - 1;
    return -1;
}

///////////////////////////////////////////////
/// \brief Prints out the transmission probability curve loaded in memory.
/// for debugging pourposes.
///
void TRestAxionXrayWindow::PrintTransmissionData() {
    if (fEnergy.size() == 0) ReadMaterial();
    for (unsigned int n = 0; n < fEnergy.size(); n++)
        cout << "Energy : " << fEnergy[n] << " Transmission : " << fTransmission[n] << endl;
}

///////////////////////////////////////////////
/// \brief Prints on screen the information about the metadata members of TRestAxionXrayWindow
///
void TRestAxionXrayWindow::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    RESTMetadata << "X-ray window type: " << fWindowType << RESTendl;
    RESTMetadata << "Window center: ( " << fCenter.X() << ", " << fCenter.Y() << ", " << fCenter.Z() << ") mm"
             << RESTendl;
    RESTMetadata << "Thickness: " << fThickness * units("um") << " um" << RESTendl;
    RESTMetadata << "Material: " << fMaterial << RESTendl;
    RESTMetadata << "Window radius: " << fRadius << " mm" << RESTendl;
    if (fWindowType != "foil") {
        RESTMetadata << "------" << RESTendl;
        RESTMetadata << "Pattern periodicity: " << fPatternGap << " mm" << RESTendl;
        RESTMetadata << "Pattern width: " << fPatternWidth << " mm" << RESTendl;
        RESTMetadata << "Pattern offset: " << fPatternOffset << " mm" << RESTendl;
    }

    RESTMetadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << RESTendl;
}
