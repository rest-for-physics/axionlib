#include "TCanvas.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TRestAxionBufferGas.h"
#include "TRestAxionField.h"

Double_t fromEnergy = 0.5;
Double_t toEnergy = 10;

Double_t incidenceAngle = 0.2;
Double_t deltaE = 0.1;

//*******************************************************************************************************
//*** Description: This script will launch the integration of the axion-field with given parameters.
//*** It allows to test different magnetic field cell sizes, for a given mass that can be off-resonance
//*** for dm different from zero, and a given maximum tolerance or error for the integration routine.
//***
//*** The macro sets the TRestAxionField under debug mode to print the different results on screen.
//***
//*** --------------
//*** Usage: restManager FieldIntegrationTests [sX=10] [sX=10] [sZ=10] [dm=0.01] [tolerance=0.1] [Ea=4.2]
//*** --------------
//***
//*** Author: Javier Galan
//*******************************************************************************************************
int REST_Axion_AccurateEfficiencies( std::string fluxFile = "fluxes.rml", std::string fluxName = "LennertHoofPrimakoff", std::string opticsFile = "xmmTrueWolter.rml", std::string opticsName = "xmm" ) {
	
	TRestAxionTrueWolterOptics optics( opticsFile.c_str(), opticsName.c_str() );
	TRestAxionOpticsMirror *mirror = optics.GetMirrorProperties();

	TRestAxionSolarQCDFlux flux( fluxFile.c_str(), fluxName.c_str());
	flux.Initialize();

	Double_t R2sum = 0;
	Double_t fluxSum = 0;
	for( Double_t energy = fromEnergy; energy < toEnergy; energy += deltaE )
	{
		Double_t R = mirror->GetReflectivity(incidenceAngle, energy );
		fluxSum += flux.GetFluxAtEnergy(energy, 0);
		R2sum += flux.GetFluxAtEnergy(energy, 0) * R * R; 
	}

	Double_t R2eff = R2sum/fluxSum;
	std::cout << "R2eff: " << R2eff << std::endl;

	Double_t Rout = optics.GetMaxEntranceRadius();
	Double_t Rin = optics.GetMinEntranceRadius();

	TRestSpiderMask *sMask = optics.GetSpiderMask();
	Double_t na = sMask->GetNumberOfArms();
	Double_t wa = sMask->GetArmsWidth();

	std::vector <Double_t> r = optics.GetR1();
	std::vector <Double_t> th = optics.GetThickness();

	Double_t Aeff = (TMath::Pi() - 0.5*na*wa )*(Rout*Rout-Rin*Rin);
	for( size_t n = 0; n < r.size(); n++ )
		Aeff -= (2*TMath::Pi()-na*wa)*r[n]*th[n];

	std::cout << "Aeff: " << Aeff/Rout/Rout/TMath::Pi() << std::endl;

	TCanvas c;
	c.SetCanvasSize(500, 1000);
	c.SetWindowSize(500, 1000);
	optics.DrawParticleTracks();

	c.Print("tracks.png");


    return 0;
}
