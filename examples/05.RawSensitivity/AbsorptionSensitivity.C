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
//*** The present macro will scan the mass range applying the best Hydrogen-Neon mixture for each mass.
//***
//*** --------------
//*** Usage: 
//***
//*** restRoot
//*** [0] .L AbsorptionSensitivity.C
//*** AbsorptionSensitivity( "CMS.rml", "cmsSignal", 5, "/tmp/CMS.txt" );
//***
//*** --------------
//*** Author: Javier Galan
//*******************************************************************************************************

/// Liquid densities
Double_t HDensity = 0.07085; // g/cm^3"
Double_t NeDensity = 1.207; // g/cm^3"
							 
Double_t photonEnergy = 4.2; //keV
Double_t fromEnergy = 0.5; //keV
Double_t toEnergy = 10; //keV
						

Double_t fromMass = 0.1;
Double_t toMass = 1000;
Double_t step = 1.01;

int AbsorptionSensitivity( std::string rmlFile, std::string section, Double_t years, std::string outputFile = "/tmp/CMSlimit.txt")
{
	Double_t tExposure = years * 365 * 24 * 3600;

	TRestAxionHelioscopeSignal helioscope( rmlFile.c_str(), section.c_str() );

	/// Pure hydrogen
	helioscope.GetGas()->SetGasDensity("H", HDensity / units("g/cm^3") );
	helioscope.GetGas()->SetGasDensity("Ne", 0 / units("g/cm^3") );
	
	// starting mass with pure hydrogen
	Double_t phMass = helioscope.GetGas()->GetPhotonMass(photonEnergy);

	// Vacuum region sensitivity is better with 20% Neon
	helioscope.GetGas()->SetGasDensity("H", 0.8 * HDensity / units("g/cm^3") );
	helioscope.GetGas()->SetGasDensity("Ne", 0.2 * NeDensity / units("g/cm^3") );

	FILE *f = fopen( "/tmp/limit.txt", "wt");
	for( double mass = fromMass; mass <= phMass; mass *= step )
	{
		Double_t gag = TMath::Power(TMath::Log(20)/helioscope.GetSignalRate(mass, fromEnergy, toEnergy)/tExposure, 0.25) * 1.e-10;
		fprintf( f, "%lf\t%e\n", mass, gag );
	}

	for( double NeFraction = 0.0; NeFraction <= 1.0; NeFraction += 0.01 )
	{
		helioscope.GetGas()->SetGasDensity("H", (1-NeFraction) * HDensity / units("g/cm^3") );
		helioscope.GetGas()->SetGasDensity("Ne", NeFraction * NeDensity / units("g/cm^3") );
		Double_t mass = helioscope.GetGas()->GetPhotonMass(photonEnergy);

		Double_t gag = TMath::Power(TMath::Log(20)/helioscope.GetSignalRate(mass, fromEnergy, toEnergy)/tExposure, 0.25) * 1.e-10;
		fprintf( f, "%lf\t%e\n", mass, gag );
	}

	phMass = helioscope.GetGas()->GetPhotonMass(photonEnergy);
	for( double mass = phMass; mass <= toMass; mass *= step )
	{
		Double_t gag = TMath::Power(TMath::Log(20)/helioscope.GetSignalRate(mass, fromEnergy, toEnergy)/tExposure, 0.25) * 1.e-10;
		fprintf( f, "%lf\t%e\n", mass, gag );
	}
	fclose(f);

	return 0;
}
