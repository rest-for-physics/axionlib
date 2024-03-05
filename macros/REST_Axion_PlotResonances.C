#include "TCanvas.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TRestAxionBufferGas.h"
#include "TRestAxionField.h"
//*******************************************************************************************************
//*** Description: This script plots the transmission probability as a function of the axion mass for a given
// gas and an
//*** axion energy until a maximum axion mass.
//***
//*** Arguments by default are (in order):
//*** - *ma_max* = 0.1: Maximum axion mass (in eV) to be plotted.
//*** - *ma_min* = 0: Minimum axion mass (in eV) to be plotted.
//*** - *Ea* = 4.2: Axion energy (in keV).
//*** - *Bmag* = 2.5: Magnetic field (in T).
//*** - *Lcoh* = 10000: Coherence length (in mm).
//*** - *gasName* = "He": Gas name.
//*** - *cutoff* = It defines a cut-off for the number of resonances to be drawn.
//*** - *n_ma* = 10000: Number of points to be plotted.
//***
//*** The generated plots are the results from `TRestAxionField::GetMassDensityScanning`,
//*** `TRestAxionField::GammaTransmissionFWHM` and `TRestAxionBufferGas::GetMassDensity`.
//***
//*** --------------
//*** Usage: restManager PlotResonances [options] [ma_max=0.1] [ma_min=0] [Ea=4.2] [Bmag=2.5] [Lcoh=10000] [gasName=He]
//*** [cutoff=5] [n_ma=10000]
//***
//*** Being all of them optional arguments.
//*** --------------
//*** Author: Fran Candón
//*******************************************************************************************************
//
const Double_t fProbMax = 3e-18;
const Double_t fNGammaMax = 6000;

const TVector2 fEnergyRange = TVector2(0.5,10);
const Double_t fExposureTime = 30 * 3600. * 18; // s
const Double_t fArea = TMath::Pi() * 35 * 35; // cm2
const Double_t fReBinning = 100; // Transforms 0.001 keV into 0.1keV											  

int REST_Axion_PlotResonances(std::string optionString = "", double ma_max = 0.1, double ma_min = 0, double Ea = 4.2, double Bmag = 2., double Lcoh = 10000, std::string gasName = "He", int cutoff = 5, int n_ma = 1000) {

	gStyle->SetPadRightMargin(0.13);
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);

	std::vector<std::string> options = TRestTools::GetOptions(optionString);

	bool title = std::find(options.begin(), options.end(), "title") != options.end();
	bool vacuum = std::find(options.begin(), options.end(), "vacuum") != options.end();
	bool drawLine = std::find(options.begin(), options.end(), "drawLine") != options.end();
	bool sumProb = std::find(options.begin(), options.end(), "sumProb") != options.end();
	bool legend = std::find(options.begin(), options.end(), "legend") != options.end();
	bool nGamma = std::find(options.begin(), options.end(), "nGamma") != options.end();

    TRestAxionField* ax = new TRestAxionField();
    ax->SetMagneticField(Bmag);
    ax->SetCoherenceLength(Lcoh);
    ax->SetAxionEnergy(Ea);

	/// Could be stored as a data member of TRestAxionPlotResonances for further access
	std::vector<std::pair<Double_t, Double_t>> Psettings = ax->GetMassDensityScanning(gasName, ma_max, 20);

	if( Psettings.size() > cutoff ) Psettings.resize(cutoff);

    // Creates the vector of axion masses (This could be a data member of TRestAxionPlotResonances).
    std::vector<double> m_a;
    double ma_step = (ma_max - ma_min) / n_ma;
	m_a.push_back( ma_min ); // The first mass drawn will get zero probability for the curve filling.
							 // So we repeat the value

    std::vector<double> sum_prob;
	sum_prob.push_back(0);
    for (int i = 0; i < n_ma; i++) {
        m_a.push_back(ma_min + i * ma_step);
        sum_prob.push_back(0);
    }

    // Computes the Vacuum probability
	ax->AssignBufferGas(nullptr);
    std ::vector<double> prob_vac;

	// Probability is not zero, but we introduce an artifact (virtual point) to make proper TGraph filling 
	prob_vac.push_back(0);
    for (size_t j = 1; j < m_a.size(); j++) {
		prob_vac.push_back(ax->GammaTransmissionProbability(m_a[j]));
    }

    // Computes the sum of all the probabilities (Adding vacuum)
    if (vacuum == true) {
        for (size_t i = 0; i < m_a.size(); i++) {
            sum_prob[i] += prob_vac[i];
        }
    }

    std::vector<TGraph*> grp;
    TRestAxionBufferGas* gas = new TRestAxionBufferGas();

    for (const auto& p : Psettings) {
        // Creates the gas and the axion field
        gas->SetGasDensity(gasName, p.second);
        ax->AssignBufferGas(gas);

        // Obtain the probability for each axion mass
        std::vector<double> prob;
		// Probability is not zero, but we introduce an artifact to make proper curve filling
		prob.push_back(0);
        for (size_t j = 1; j < m_a.size(); j++) {
			prob.push_back(ax->GammaTransmissionProbability(m_a[j]));

            sum_prob[j] += prob[j];
        }

        TGraph* gr = new TGraph(m_a.size(), &m_a[0], &prob[0]);
        grp.push_back(gr);
    }


    ///// PLOTS /////

    // Plot the density scan
    grp[0]->SetLineColor(kBlue - 3);
	if( title ) grp[0]->SetTitle("Transmission probability as a function of the axion mass");
	else grp[0]->SetTitle("");
    grp[0]->GetXaxis()->SetTitle("m_{a} [eV]");
    grp[0]->GetYaxis()->SetTitle("P_{a#gamma}");
    grp[0]->GetXaxis()->SetLimits(ma_min+0.00001, ma_max);
	grp[0]->GetXaxis()->SetLabelSize(0.04);
    grp[0]->GetXaxis()->SetLabelFont(42);
	grp[0]->GetXaxis()->SetTitleSize(0.04);
    grp[0]->GetXaxis()->SetTitleFont(42);
    grp[0]->GetYaxis()->SetRangeUser(0, fProbMax);
	grp[0]->GetYaxis()->SetLabelSize(0.04);
    grp[0]->GetYaxis()->SetLabelFont(42);
	grp[0]->GetYaxis()->SetTitleSize(0.04);
    grp[0]->GetYaxis()->SetTitleFont(42);

	grp[0]->SetLineColor(kBlack);
	grp[0]->SetFillColorAlpha(kRed,0.5);
	grp[0]->SetLineWidth(2);
    grp[0]->Draw("AFL");

    for (const auto g : grp) {
        g->SetLineColor(kBlack);
        g->SetLineWidth(2);
		g->SetFillColorAlpha(kRed,0.15);
        g->Draw("FL SAME");
    }

    // PLot of the sum of all the probabilities
	TGraph *sum;
	if( sumProb )
	{
		sum = new TGraph(m_a.size(), &m_a[0], &sum_prob[0]);
		sum->SetLineColor(kBlue + 3);
		sum->SetLineWidth(3);
		sum->Draw("SAME");
	}

    // Plot of the vacuum probability
    TGraph* vac = new TGraph(m_a.size(), &m_a[0], &prob_vac[0]);
    vac->SetLineColor(kBlack);
    vac->SetLineWidth(2);
	vac->SetFillColorAlpha(kBlue,0.25);
    vac->Draw("FL SAME");

	TGraph *gama;
	if( nGamma ) // The following code could be a method of a future TRestAxionPlotResonances::DrawNGamma()
	{
		std ::vector<double> NGamma;

		TRestAxionSolarQCDFlux* sFlux = new TRestAxionSolarQCDFlux("fluxes.rml", "LennertHoofPrimakoff");
		sFlux->Initialize();
		TH1F* spec = sFlux->GetEnergySpectrum();

		/// The original gets 1eV binning to describe monocromatic lines.
		/// We do not need such a high resolution 
		TH1F *rebinned = dynamic_cast<TH1F*>(spec->Rebin(fReBinning,"rebinned"));
		/// We need to renormalize since the bin contents are already renormalized to cm-2 s-1 keV-1
		rebinned->Scale(1./fReBinning);
		Double_t deltaE = rebinned->GetBinCenter(2) - rebinned->GetBinCenter(1);

		if(!rebinned) {
			std::cout << "Energy spectrum is nullptr!" << std::endl;
		}
		else
		{
			Double_t maxNGamma = 0;
			// Obtain the Ngamma for each axion mass
			for (size_t j = 0; j < m_a.size(); j++) {
				std::cout << "Producing Ngamma for mass : " << m_a[j] << std::endl;

				Double_t ng = 0;

				// Vacuum
				ax->AssignBufferGas(nullptr);
				for( int n = 1; n <= rebinned->GetNbinsX(); n++ )
				{
					Double_t energy = rebinned->GetBinCenter(n);
					if ( energy < fEnergyRange.X() || energy > fEnergyRange.Y() )
						continue;

					ax->SetAxionEnergy( rebinned->GetBinCenter(n) );
					ng += rebinned->GetBinContent(n) * ax->GammaTransmissionProbability(m_a[j]);
				}

				/// Gas phase
				for (const auto& p : Psettings) {
					// Creates the gas and assignss it to the axion field
					gas->SetGasDensity(gasName, p.second);
					ax->AssignBufferGas(gas);

					for( int n = 1; n <= rebinned->GetNbinsX(); n++ )
					{
						Double_t energy = rebinned->GetBinCenter(n);
						if ( energy < fEnergyRange.X() || energy > fEnergyRange.Y() )
							continue;

						ax->SetAxionEnergy( rebinned->GetBinCenter(n) );
						ng += rebinned->GetBinContent(n) * ax->GammaTransmissionProbability(m_a[j]);
					}
				}
				ng = ng * deltaE * fArea * fExposureTime;
				if( ng > maxNGamma ) maxNGamma = ng;
				NGamma.push_back( ng );
			}

			gama = new TGraph(m_a.size(), &m_a[0], &NGamma[0]);
			gama->SetLineColor(kRed + 3);
			gama->SetLineWidth(3);
			gama->GetYaxis()->SetRangeUser(0, fNGammaMax);
			gama->Scale(fProbMax/fNGammaMax);
			gama->Draw("SAME");

			std::cout << "Max: " << maxNGamma << std::endl;

			/// 0.001 is to remove the first tick, which is zero, and overlaps with 0.1eV
			TGaxis *A1 = new TGaxis(0.1,0,0.1,fProbMax,0.001,fNGammaMax,510,"+L");
			A1->SetTitle("N_{#gamma}");
			A1->SetTitleOffset(1.5);
			A1->SetLabelSize(0.04);
			A1->SetLabelFont(42);
			A1->SetTextFont(42);
			A1->Draw();
		}
	}

    // Plot of the vertical line
	TLine* verticalLine;
	if( drawLine )
	{
		verticalLine = new TLine(Psettings[0].first, c1->GetUymin(), Psettings[0].first, fProbMax);
		verticalLine->SetLineColor(kGreen - 3);
		verticalLine->SetLineWidth(2);
		verticalLine->Draw("same");
	}

    // Plot of the legend
	if( legend )
	{
		TLegend* legend = new TLegend(0.87, 0.7, 0.47, 0.9);
		legend->AddEntry(grp[0], "P_{a#gamma} for each P-step", "lf");
		legend->AddEntry(vac, "P_{a#gamma} for vacuum", "lf");
		if( sumProb ) legend->AddEntry(sum, "#Sigma P_{a#gamma}", "l");
		if( nGamma ) legend->AddEntry(gama, "N_{#gamma}", "l");
		if( drawLine ) legend->AddEntry(verticalLine, "m_{a} where P_{a#gamma}^{vac} = max(P_{a#gamma}^{vac}/2) ", "l");
		legend->Draw("same");
		c1->Draw();
	}

    c1->Print("/tmp/resonances.pdf");

    return 0;
}
