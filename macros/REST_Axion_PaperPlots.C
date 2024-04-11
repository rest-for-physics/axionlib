#include <TRestTools.h>
//*******************************************************************************************************
//*** Description: 
//***
//*** Arguments by default are (in order):
//***
//*** The generated plots are the results from `TRestAxionField::GetMassDensityScanning`,
//*** `TRestAxionField::GammaTransmissionFWHM` and `TRestAxionBufferGas::GetMassDensity`.
//***
//*** --------------
//*** Usage: restManager PlotResonances [options] [ma_max=0.1] [ma_min=0] [Ea=4.2] [Bmag=2.5] [Lcoh=10000]
//*** [gasName=He]
//*** [cutoff=5] [n_ma=10000]
//***
//*** Being all of them optional arguments.
//*** --------------
//*** Author: Javier Galan
//*******************************************************************************************************

const Int_t samples = 1000000;
std::string fluxesString = "LennertHoofPrimakoff:LennertHoofABC";

int SolarPlots()
{
    std::vector<std::string> fluxes = TRestTools::GetOptions(fluxesString);

	TCanvas *c1 = new TCanvas("c1", "My canvas", 800, 350);
	c1->GetFrame()->SetBorderSize(6);
	c1->GetFrame()->SetBorderMode(-1);

	TPad *pad1 = new TPad("pad1", "This is pad1", 0.01, 0.02, 0.99, 0.97);
	pad1->Divide(2, 1);
	pad1->Draw();

	TRestAxionSolarQCDFlux *solarFlux1 = new TRestAxionSolarQCDFlux("fluxes.rml", fluxes[0]);
	solarFlux1->Initialize();
	solarFlux1->PrintMetadata();

	TRestAxionSolarQCDFlux *solarFlux2 = new TRestAxionSolarQCDFlux("fluxes.rml", fluxes[1]);
	solarFlux2->Initialize();
	solarFlux2->PrintMetadata();

	Int_t enBins = 2000;
	Double_t enMax = 20;

	//// Producing histograms ///
	TH2D *EvsRHisto1 = new TH2D("FluxEvsR1", "", enBins, 0, enMax, 100, 0, 1);
	for (int n = 0; n < samples; n++ )
	{
		std::pair<Double_t,Double_t> x = solarFlux1->GetRandomEnergyAndRadius();
		EvsRHisto1->Fill(x.first, x.second);
	}

	TH2D *EvsRHisto2 = new TH2D("FluxEvsR2", "", enBins, 0, enMax, 100, 0, 1);
	for (int n = 0; n < samples; n++ )
	{
		std::pair<Double_t,Double_t> x = solarFlux2->GetRandomEnergyAndRadius();
		EvsRHisto2->Fill(x.first, x.second);
	}

	TRandom3 *rnd = new TRandom3(0);
	TH2D *solarDisk = new TH2D("solar_disk", "SolarDisk", 120, -1.2, 1.2, 120, -1.2, 1.2);
	for (int n = 0; n < samples; n++ )
	{
		Double_t angle = rnd->Rndm() * 2 * TMath::Pi();
		std::pair<Double_t,Double_t> x = solarFlux1->GetRandomEnergyAndRadius();

		solarDisk->Fill(x.second * TMath::Cos(angle), x.second * TMath::Sin(angle));
	}

	TH1D *enSpt1 = EvsRHisto1->ProjectionX();
	TH1D *rSpt1 = EvsRHisto1->ProjectionY();

	TH1D *enSpt2 = EvsRHisto2->ProjectionX();
	TH1D *rSpt2 = EvsRHisto2->ProjectionY();

	pad1->cd(1);
	pad1->cd(1)->SetRightMargin(0.12);
	pad1->cd(1)->SetLeftMargin(0.15);
	pad1->cd(1)->SetBottomMargin(0.15);

	EvsRHisto1->SetStats(0);
	EvsRHisto1->GetXaxis()->SetTitle("Energy [keV]");
	EvsRHisto1->GetXaxis()->SetTitleSize(0.05);
	EvsRHisto1->GetXaxis()->SetLabelSize(0.05);
	EvsRHisto1->GetYaxis()->SetTitle("Solar radius");
	EvsRHisto1->GetYaxis()->SetTitleSize(0.05);
	EvsRHisto1->GetYaxis()->SetLabelSize(0.05);
	EvsRHisto1->Draw("colz");

	pad1->cd(2);
	//pad1->cd(2)->SetLogy();
	pad1->cd(2)->SetRightMargin(0.09);
	pad1->cd(2)->SetLeftMargin(0.15);
	pad1->cd(2)->SetBottomMargin(0.15);

	enSpt1->GetYaxis()->SetTitleSize(0.05);
	enSpt1->SetStats(0);
	enSpt1->GetXaxis()->SetTitle("Energy [keV]");
	enSpt1->GetXaxis()->SetTitleSize(0.05);
	enSpt1->GetXaxis()->SetLabelSize(0.05);
	enSpt1->GetYaxis()->SetTitle("Counts [keV^{-1}]");
	enSpt1->GetYaxis()->SetTitleSize(0.05);
	enSpt1->GetYaxis()->SetLabelSize(0.05);
	//enSpt1->SetFillStyle(4050);
	//enSpt1->SetFillColor(kBlue - 9);
	enSpt1->SetFillColorAlpha(kBlue,0.15);
	enSpt1->SetLineColor(kBlack);
	enSpt1->Scale(enBins/enMax);
    enSpt1->GetYaxis()->SetRangeUser(0., 500000);
	enSpt1->Draw("HIST FL");
	//enSpt2->SetFillColor(kRed - 9);
	enSpt2->SetFillColorAlpha(kRed,0.15);
	enSpt2->SetLineColor(kBlack);
	enSpt2->Scale(enBins/enMax);
	enSpt2->Draw("HIST SAME");

	TLegend* legend = new TLegend(0.87, 0.75, 0.55, 0.9);
	legend->AddEntry(enSpt1, "Primakoff", "lf");
	legend->AddEntry(enSpt2, "ABC", "lf");
	legend->Draw("same");

    c1->Print("/tmp/solarFlux.pdf");
	return 0;
}

int MagnetPlots()
{
	TRestAxionMagneticField *field = new TRestAxionMagneticField( "fields.rml", "babyIAXO_2024_cutoff" );
	TCanvas *c = field->DrawTracks( TVector3(0,0,8000), 100 );
	c->Print("/tmp/tracks.pdf");
	return 0;
}

int OpticsPlots()
{
	TCanvas *c1 = new TCanvas("c1", "My canvas", 800, 350);
	TRestAxionTrueWolterOptics *optics = new TRestAxionTrueWolterOptics("xmmTrueWolter.rml", "xmm" );
	TPad *pad = optics->DrawDensityMaps(7000, 4.2, 0.0, samples, 7500, 9, true);
	c1->Print("/tmp/optics.pdf");
	
	TCanvas *c2 = new TCanvas("c2", "My canvas", 600, 800);
	TPad *pad2 = optics->DrawParticleTracks(0.001, 40, true);
	c2->Print("/tmp/optics2.pdf");

	return 0;
}

int MirrorPlots()
{
	TRestAxionOpticsMirror *mirror = new TRestAxionOpticsMirror("mirrors.rml", "default");
	TCanvas *c = mirror->DrawOpticsProperties("[4,7,10](0,1){0.73,0.79,0.9,0.96}:[0.25,0.5,0.75,1](0,10){0.73,0.79,0.9,0.96}", 1.e-5, 1.e-3, true);
	c->Print("/tmp/mirrors.pdf");
	return 0;
}

int REST_Axion_PaperPlots(std::string plotsString = "magnetPlots" )
{
	if( plotsString == "magnetPlots" )
		MagnetPlots();
	if( plotsString == "opticsPlots" )
		OpticsPlots();
	if( plotsString == "mirrorPlots" )
		MirrorPlots();
	return 0;
}
