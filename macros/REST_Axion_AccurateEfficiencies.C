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
int REST_Axion_AccurateEfficiencies(std::string fluxFile = "fluxes.rml",
                                    std::string fluxName = "LennertHoofPrimakoff",
                                    std::string opticsFile = "xmmTrueWolter.rml",
                                    std::string opticsName = "xmm") {
    TRestAxionTrueWolterOptics optics(opticsFile.c_str(), opticsName.c_str());
    TRestAxionOpticsMirror* mirror = optics.GetMirrorProperties();

    TRestAxionSolarQCDFlux flux(fluxFile.c_str(), fluxName.c_str());
    flux.Initialize();

	Double_t R2sum = 0;
	Double_t fluxSum = 0;
	Double_t maxFlux = 0;
	for( Double_t energy = fromEnergy; energy < toEnergy; energy += deltaE )
	{
		Double_t R = mirror->GetReflectivity(incidenceAngle, energy );
		fluxSum += flux.GetFluxAtEnergy(energy, 0);
		if( flux.GetFluxAtEnergy(energy, 0) > maxFlux ) maxFlux = flux.GetFluxAtEnergy(energy, 0);
		R2sum += flux.GetFluxAtEnergy(energy, 0) * R * R; 
	}

    Double_t R2eff = R2sum / fluxSum;
    std::cout << "R2eff: " << R2eff << std::endl;

    Double_t Rout = optics.GetMaxEntranceRadius();
    Double_t Rin = optics.GetMinEntranceRadius();

    TRestSpiderMask* sMask = optics.GetSpiderMask();
    Double_t na = sMask->GetNumberOfArms();
    Double_t wa = sMask->GetArmsWidth();

    std::vector<Double_t> r = optics.GetR1();
    std::vector<Double_t> th = optics.GetThickness();

    Double_t Aeff = (TMath::Pi() - 0.5 * na * wa) * (Rout * Rout - Rin * Rin);
    for (size_t n = 0; n < r.size(); n++) Aeff -= (2 * TMath::Pi() - na * wa) * r[n] * th[n];

	std::cout << "Aeff (optics): " << Aeff/Rout/Rout/TMath::Pi() << std::endl;

	TCanvas c;
	c.SetCanvasSize(2400, 1800);
	c.SetWindowSize(2400, 1800);
	c.Divide(2,1);

	c.cd(1);
	optics.GetMirrorProperties()->DrawOpticsPropertiesLinear();

	c.cd(2);
	optics.DrawParticleTracks();

	c.Print("optics.pdf");

	/// Extracted from x-ray window MicromegasStrongBack
	na = 8;
	wa = 2.64*TMath::Pi()/180.;
	Rout = 8.5;
	Rin = 4.55;
	Double_t Ro = 4.25;
	
	Aeff = (TMath::Pi() - 0.5*na*wa )*(Rout*Rout-Rin*Rin) + TMath::Pi()*Ro*Ro;

	std::cout << "Aeff (window): " << Aeff/Rout/Rout/TMath::Pi() << std::endl;

	TRestAxionXrayWindow strongBack("windows.rml", "MicromegasStrongBack");
	TRestAxionXrayWindow mylar("windows.rml", "MicromegasMylar");
	TRestAxionXrayWindow aluminum("windows.rml", "MicromegasAluminumFoil");
	
    TGraph* mylarGraph = new TGraph();
	mylarGraph->SetName( "Mylar" );
	mylarGraph->SetLineColor(49);
	mylarGraph->SetLineWidth(2);

    TGraph* aluminumGraph = new TGraph();
	aluminumGraph->SetName("Aluminum");
	aluminumGraph->SetLineColor(46);
	aluminumGraph->SetLineWidth(2);

    TGraph* solarGraph = new TGraph();
	solarGraph->SetName( "SolarFlux" );
	solarGraph->SetLineColor(43);
	solarGraph->SetLineWidth(2);
	
	Double_t WeffSum = 0;
	Double_t AlSum = 0;
	Double_t MySum = 0;
	for( Double_t energy = deltaE; energy < toEnergy; energy += deltaE )
	{
		Double_t tMy = mylar.GetTransmission( energy, 0 , 0 );
		Double_t tAl = aluminum.GetTransmission( energy, 0 , 0 );
		
		mylarGraph->SetPoint(mylarGraph->GetN(), energy, tMy);
		aluminumGraph->SetPoint(aluminumGraph->GetN(), energy, tAl);
		solarGraph->SetPoint(solarGraph->GetN(), energy, flux.GetFluxAtEnergy(energy,0)/maxFlux);

		WeffSum += flux.GetFluxAtEnergy(energy, 0) * tMy * tAl; 
		AlSum += flux.GetFluxAtEnergy(energy, 0) * tAl; 
		MySum += flux.GetFluxAtEnergy(energy, 0) * tMy; 
	}

	WeffSum = WeffSum/fluxSum;
	AlSum = AlSum/fluxSum;
	MySum = MySum/fluxSum;
	std::cout << "AlSum: " << AlSum << std::endl;
	std::cout << "MySum: " << MySum << std::endl;
	std::cout << "WeffSum: " << WeffSum << std::endl;
	
	TCanvas c2;
	c2.SetCanvasSize(1200, 900);
	c2.SetWindowSize(1200, 900);
	//c2.SetLogy();

    TPad* pad2 = new TPad("pad1", "This is pad1", 0.01, 0.02, 0.99, 0.97);
 //   pad1->Divide(2, 2);
	pad2->SetLogy();
    pad2->Draw();

    ////// Drawing reflectivity versus angle
    pad2->SetRightMargin(0.09);
    pad2->SetLeftMargin(0.25);
    pad2->SetBottomMargin(0.15);

    mylarGraph->GetXaxis()->SetLimits(0, 10);
 //   mylarGraph->GetHistogram()->SetMaximum(1);
    mylarGraph->GetHistogram()->SetMinimum(0);

	mylarGraph->GetXaxis()->SetTitle("Energy [keV]");
	mylarGraph->GetXaxis()->SetTitleSize(0.04);
	mylarGraph->GetXaxis()->SetLabelSize(0.04);
	mylarGraph->GetYaxis()->SetTitle("Transmission");
	mylarGraph->GetYaxis()->SetTitleOffset(1.2);
	mylarGraph->GetYaxis()->SetTitleSize(0.04);
	mylarGraph->GetYaxis()->SetLabelSize(0.04);
    mylarGraph->Draw("AL");
    aluminumGraph->Draw("L");
    //solarGraph->Draw("L");

    Double_t lx1 = 0.6, ly1 = 0.55, lx2 = 0.8, ly2 = 0.75;
    TLegend* legend = new TLegend(lx1, ly1, lx2, ly2);
	legend->SetTextSize(0.03);
	//   legend->SetHeader("Widnows", "C");  // option "C" allows to center the header
	legend->AddEntry( "Mylar", "Mylar", "l");
	legend->AddEntry( "Aluminum", "Aluminum", "l");
	legend->Draw();

	c2.Print("windows.pdf");

    c.Print("tracks.png");

    return 0;
}
