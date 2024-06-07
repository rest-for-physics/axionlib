#include "TCanvas.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TRestAxionHelioscopeSignal.h"

//*******************************************************************************************************
//*** Description: This macro produces an ultimate accesibility zone for the IAXO sensivity in case we
//*** focus the total exposure time in a specific mass. We assume zero background so the coupling is
//*** calculated according to: (Log(20)/Ngamma)^1/4 [10^-10 GeV-1].
//*** Equation (8.14) in https://arxiv.org/pdf/1102.1406
//***
//*** The present macro will scan the mass range applying the best Helium density for each mass. The
//*** vacuum sensitivity will be also calculated and the limit will register the lower between them.
//***
//*** --------------
//*** Usage:
//***
//*** restRoot
//*** [0] .L TransmissionSensitivity.C
//*** TransmissionSensitivity( "IAXO.rml", "BabyIAXOSignal", 1.5, "/tmp/BabyIAXO.txt" );
//***
//*** --------------
//*** Author: Javier Galan
//*******************************************************************************************************

Double_t HeDensity = 0.178e-3;  // g/cm3 at 1bar

Double_t photonEnergy = 4.2;  // keV
Double_t fromEnergy = 0.5;    // keV
Double_t toEnergy = 10;       // keV

Double_t fromMass = 0.001;
Double_t toMass = 100;
Double_t step = 1.01;

int TransmissionSensitivity(std::string rmlFile, std::string section, Double_t years,
                            std::string outputFile = "/tmp/IAXOlimit.txt") {
    Double_t tExposure = years * 365 * 12 * 3600;

    // We recover helioscope settings from RML, but we will redefine the gas internally
    TRestAxionBufferGas gas(rmlFile.c_str());
    TRestAxionSolarQCDFlux flux(rmlFile.c_str());
    flux.Initialize();

    TRestAxionHelioscopeSignal helioscopeGas(rmlFile.c_str(), section.c_str());
    helioscopeGas.AssignBufferGas(&gas);
    helioscopeGas.AssignFlux(&flux);

    TRestAxionHelioscopeSignal helioscopeVacuum(rmlFile.c_str(), section.c_str());
    helioscopeVacuum.AssignBufferGas(nullptr);
    helioscopeVacuum.AssignFlux(&flux);

    FILE* f = fopen(outputFile.c_str(), "wt");
    for (double mass = fromMass; mass <= toMass; mass *= step) {
        helioscopeGas.GetGas()->SetGasDensity("He", gas.GetDensityForMass(mass));

        Double_t gag_vacuum =
            TMath::Power(
                TMath::Log(20) / helioscopeVacuum.GetSignalRate(mass, fromEnergy, toEnergy) / tExposure,
                0.25) *
            1.e-10;
        Double_t gag_gas =
            TMath::Power(TMath::Log(20) / helioscopeGas.GetSignalRate(mass, fromEnergy, toEnergy) / tExposure,
                         0.25) *
            1.e-10;
        if (gag_gas < gag_vacuum)
            fprintf(f, "%lf\t%e\n", mass, gag_gas);
        else
            fprintf(f, "%lf\t%e\n", mass, gag_vacuum);
    }
    fclose(f);

    return 0;
}
