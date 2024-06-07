/*************************************************************************
 * This file is part of the REST software framework.                     *
 *                                                                       *
 * Copyright (C) 2016 GIFNA/TREX (University of Zaragoza)                *
 * For more information see https://gifna.unizar.es/trex                 *
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
 * If not, see https://www.gnu.org/licenses/.                            *
 * For the list of contributors see $REST_PATH/CREDITS.                  *
 *************************************************************************/

/////////////////////////////////////////////////////////////////////////
/// This class describes an axion helioscope signal using different
/// aproximations, calculating an expected rate as a function of energy
/// (independent of position).
///
///    \code
///			<TRestAxionHelioscopeSignal name="BabyIAXO" nature="signal"
///							conversionType="IAXO" bores="2"
///							magnetRadius="35cm" magnetLength="10m"
/// magnetStrength="2T" 							opticsEfficiency="0.3"
/// windowEfficiency="0.8">
///
///				<!-- TRestComponent common fields -->
///				<parameter name="parameterizationNodes" value="{0.0001,0.001,0.01,0.1,1}" />
///				<cVariable name="energy" range="(0,10)keV" bins="20" />
///
///				<!-- Solar flux -->
///				<TRestAxionSolarQCDFlux name="LennertHoofPrimakoff" verboseLevel="warning" >
///					<parameter name="couplingType" value="g_ag"/>
///					<parameter name="couplingStrength" value="1.e-10"/>
///					<parameter name="fluxDataFile"
/// value="Primakoff_LennertHoof_202203.dat"/>
///
///					<parameter name="seed" value="137" />
///				</TRestAxionSolarQCDFlux>
///
///				<!-- Buffer gas -->
///    			<TRestAxionBufferGas name="helium" verboseLevel="warning">
///        			<gas name="He" density="0.0025e-6g/cm^3"/>
///    			</TRestAxionBufferGas>
///
///				<!-- Detector response -->
///				<TRestResponse name="XenonNeon" variable="energy">
///					<parameter name="filename" value="XenonNeon_50Pct_1.4bar.N150f" />
///				</TRestResponse>
///
///			</TRestAxionHelioscopeSignal>
///    </TRestAxionHelioscopeSignal>
///    \endcode
///
///----------------------------------------------------------------------
///
/// REST-for-Physics - Software for Rare Event Searches Toolkit
///
/// History of developments:
///
/// 2023-December: First implementation of TRestAxionHelioscopeSignal
/// Javier Galan
///
/// \class TRestAxionHelioscopeSignal
/// \author: Javier Galan (javier.galan.lacarra@cern.ch)
///
/// <hr>
///
#include "TRestAxionHelioscopeSignal.h"

#include <numeric>

#include "TKey.h"

ClassImp(TRestAxionHelioscopeSignal);

///////////////////////////////////////////////
/// \brief Default constructor
///
TRestAxionHelioscopeSignal::TRestAxionHelioscopeSignal() { Initialize(); }

/////////////////////////////////////////////
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
TRestAxionHelioscopeSignal::TRestAxionHelioscopeSignal(const char* cfgFileName, const std::string& name)
    : TRestComponent(cfgFileName) {
    Initialize();

    LoadConfigFromFile(fConfigFileName, name);

    if (GetVerboseLevel() >= TRestStringOutput::REST_Verbose_Level::REST_Info) PrintMetadata();
}

///////////////////////////////////////////////
/// \brief Default destructor
///
TRestAxionHelioscopeSignal::~TRestAxionHelioscopeSignal() {
    if (fField) delete fField;
}

///////////////////////////////////////////////
/// \brief It will initialize the data frame with the filelist and column names
/// (or observables) that have been defined by the user.
///
void TRestAxionHelioscopeSignal::Initialize() {
    TRestComponent::Initialize();

    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);

    if (fField == nullptr) fField = new TRestAxionField();
}

///////////////////////////////////////////////
/// \brief It returns the intensity/rate (in seconds) corresponding to the
/// class defined helioscope configuration evaluated at the position of the
/// parameter space given by point, which by now is simply the energy.
///
/// The size of the point vector must have the same dimension as the dimensions
/// of the variables of the distribution. Right now is 1-dimension.
///
Double_t TRestAxionHelioscopeSignal::GetSignalRate(std::vector<Double_t> point, Double_t mass) {
    if (GetDimensions() != point.size()) {
        RESTError << "Point should have same dimensions as number of variables!" << RESTendl;
        return 0;
    }

    if (GetDimensions() != 1) {
        RESTError << "Point should have only 1-dimension! Energy" << RESTendl;
        return 0;
    }

    Double_t flux = fFlux->GetFluxAtEnergy(point[0], mass);  // cm-2 s-1 keV-1

    Double_t probability = 0;
    if (ToLower(fConversionType) == "iaxo") {
        probability =
            fOpticsEfficiency * fWindowEfficiency * fField->GammaTransmissionProbability(point[0], mass);

        // We assume all flux ends up inside the spot. No XY dependency of signal.
        Double_t apertureArea = TMath::Pi() * fMagnetRadius * units("cm") * fMagnetRadius * units("cm");

        flux *= fBores * apertureArea;
    }

    if (ToLower(fConversionType) == "amelie") {
        probability = fField->AxionAbsorptionProbability(point[0], mass);

        Double_t apertureArea = TMath::Pi() * fMagnetRadius * units("cm") * fMagnetRadius * units("cm");

        flux *= fBores * apertureArea;
    }

    Double_t signal = flux * probability;

    Double_t normFactor = 1;
    /// There should be only 1-dimension. Rate is integrated to the particular bin size.
    for (size_t n = 0; n < GetDimensions(); n++) {
        normFactor *= (fRanges[n].Y() - fRanges[n].X()) / fNbins[n];
    }

    return normFactor * signal;  // s-1
}

///////////////////////////////////////////////
/// \brief It returns the intensity/rate (in seconds) corresponding to the
/// class defined helioscope configuration integrated in the energy range
/// given by argument.
///
Double_t TRestAxionHelioscopeSignal::GetSignalRate(Double_t mass, Double_t Eo, Double_t Ef) {
    Double_t dE = 0.5;

    Double_t signal = 0;
    for (Double_t en = Eo; en < Ef; en += dE) {
        Double_t flux = fFlux->GetFluxAtEnergy(en, mass);  // cm-2 s-1 keV-1

        /// Below is copy/paste from previous method. Sorry for doing this. (Just replaced point[0] by en)
        Double_t probability = 0;
        if (ToLower(fConversionType) == "iaxo") {
            probability =
                fOpticsEfficiency * fWindowEfficiency * fField->GammaTransmissionProbability(en, mass);

            // We assume all flux ends up inside the spot. No XY dependency of signal.
            Double_t apertureArea = TMath::Pi() * fMagnetRadius * units("cm") * fMagnetRadius * units("cm");
            flux *= fBores * apertureArea;
        }

        if (ToLower(fConversionType) == "amelie") {
            probability = fField->AxionAbsorptionProbability(en, mass);

            Double_t apertureArea = TMath::Pi() * fMagnetRadius * units("cm") * fMagnetRadius * units("cm");

            flux *= fBores * apertureArea;
        }
        /// Above is a copy/paste from previous method. Sorry for doing this.
        signal += flux * probability;
    }

    return signal * dE;  // s-1
}

/////////////////////////////////////////////
/// \brief It will produce a histogram with the distribution using the formula
/// contributions.
///
/// TODO: The histogram is filled just by evaluating the formula, but it would
/// be more realistic that we fill the histograms with a number N of entries
/// that mimic a MC generation scheme similar to TRestComponentDataSet.
///
void TRestAxionHelioscopeSignal::FillHistograms() {
    if (!HasNodes() || !fFlux) return;
    fNodeDensity.clear();

    RESTInfo << "Generating N-dim histogram for " << GetName() << RESTendl;

    int nIndex = 0;
    for (const auto& node : fParameterizationNodes) {
        TString hName = fParameter + "_" + DoubleToString(node);

        Int_t* bins = new Int_t[fNbins.size()];
        Double_t* xlow = new Double_t[fNbins.size()];
        Double_t* xhigh = new Double_t[fNbins.size()];

        for (size_t n = 0; n < fNbins.size(); n++) {
            bins[n] = fNbins[n];
            xlow[n] = fRanges[n].X();
            xhigh[n] = fRanges[n].Y();
        }

        THnD* hNd = new THnD(hName, hName, fNbins.size(), bins, xlow, xhigh);

        // Calculate the bin width in each dimension
        std::vector<double> binWidths;
        for (size_t i = 0; i < fNbins.size(); ++i) {
            double width = static_cast<double>(xhigh[i] - xlow[i]) / bins[i];
            binWidths.push_back(width);
        }

        // Nested loop to iterate over each bin and print its center
        std::vector<int> binIndices(fNbins.size(), 0);  // Initialize bin indices to 0 in each dimension

        bool carry = false;
        while (!carry) {
            // Calculate the center of the current bin in each dimension
            std::vector<double> binCenter;
            for (size_t i = 0; i < fNbins.size(); ++i)
                binCenter.push_back(xlow[i] + (binIndices[i] + 0.5) * binWidths[i]);

            hNd->Fill(binCenter.data(), GetSignalRate(binCenter, node));

            // Update bin indices for the next iteration
            carry = true;
            for (size_t i = 0; i < fNbins.size(); ++i) {
                binIndices[i]++;
                if (binIndices[i] < bins[i]) {
                    carry = false;
                    break;
                }
                binIndices[i] = 0;
            }
        }

        fNodeDensity.push_back(hNd);
        fActiveNode = nIndex;
        nIndex++;
    }
}

/////////////////////////////////////////////
/// \brief Prints on screen the information about the metadata members of TRestAxionSolarFlux
///
void TRestAxionHelioscopeSignal::PrintMetadata() {
    TRestComponent::PrintMetadata();

    RESTMetadata << "Conversion type: " << fConversionType << RESTendl;
    RESTMetadata << " " << RESTendl;

    RESTMetadata << "Magnet bores : " << fBores << RESTendl;
    RESTMetadata << "Magnet radius : " << fMagnetRadius * units("cm") << " cm" << RESTendl;
    RESTMetadata << "Magnet length : " << fMagnetLength * units("m") << " m" << RESTendl;
    RESTMetadata << "Magnet field : " << fMagnetStrength * units("T") << " T" << RESTendl;
    RESTMetadata << " " << RESTendl;

    RESTMetadata << "Optics efficiency : " << fOpticsEfficiency << RESTendl;
    RESTMetadata << "Window efficiency : " << fWindowEfficiency << RESTendl;

    RESTMetadata << "----" << RESTendl;
}

/////////////////////////////////////////////
/// \brief It customizes the retrieval of XML data values of this class
///
void TRestAxionHelioscopeSignal::InitFromConfigFile() {
    TRestComponent::InitFromConfigFile();

    if (fVariables.size() != 1) {
        RESTError << "TRestAxionHelioscopeSignal::InitFromConfigFile."
                  << " Signal should be build with just 1-variable. Energy." << RESTendl;
    } else if (fVariables[0] != "energy") {
        RESTError << "The first variable should be energy. We recommend to name that variable as \"energy\"."
                  << RESTendl;
        RESTError << "Please, double-check the variables definition inside "
                     "TRestAxionHelioscopeSignal::TRestComponent."
                  << RESTendl;
    }

    if (fFlux) {
        delete fFlux;
        fFlux = nullptr;
    }
    fFlux = (TRestAxionSolarFlux*)this->InstantiateChildMetadata("TRestAxionSolarQCDFlux");
	if( fFlux ) fFlux->Initialize();

    if (fGas) {
        delete fGas;
        fGas = nullptr;
    }
    fGas = (TRestAxionBufferGas*)this->InstantiateChildMetadata("TRestAxionBufferGas");

    if (fField == nullptr) fField = new TRestAxionField();

    fField->SetMagneticField(GetMagnetStrength());
    if (ToLower(fConversionType) == "iaxo")
        fField->SetCoherenceLength(GetMagnetLength());
    else if (ToLower(fConversionType) == "amelie")
        fField->SetCoherenceLength(GetMagnetLength());

    fField->AssignBufferGas(fGas);

    TRestAxionHelioscopeSignal::Initialize();
}
