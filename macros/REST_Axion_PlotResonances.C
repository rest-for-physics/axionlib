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
//*** - *vacuum* = true: If true, the vacuum probability is added to the sum of all the probabilities.
//*** - *n_ma* = 10000: Number of points to be plotted.
//***
//*** The generated plots are the results from `TRestAxionField::GetMassDensityScanning`,
//*** `TRestAxionField::GammaTransmissionFWHM` and `TRestAxionBufferGas::GetMassDensity`.
//***
//*** --------------
//*** Usage: restManager PlotResonances [ma_max=0.1] [ma_min=0] [Ea=4.2] [Bmag=2.5] [Lcoh=10000] [gasName=He]
//[vacuum=true] [n_ma=7000]
//***
//*** Being all of them optional arguments.
//*** --------------
//*** Author: Fran Candón
//*******************************************************************************************************

int REST_Axion_PlotResonances(double ma_max = 0.1, double ma_min = 0, double Ea = 4.2, double Bmag = 2.5,
                              double Lcoh = 10000, std::string gasName = "He", Bool_t vacuum = true,
                              int n_ma = 10000) {
    TRestAxionField* ax = new TRestAxionField();
    ax->SetMagneticField(Bmag);
    ax->SetCoherenceLength(Lcoh);
    ax->SetAxionEnergy(Ea);

    vector<std::pair<Double_t, Double_t>> pair = ax->GetMassDensityScanning(gasName, ma_max, Ea, 3);
    std::vector<double> m_a;
    std::vector<double> sum_prob;

    // Creates the vector of axion masses
    double ma_step = (ma_max - ma_min) / n_ma;
    for (int i = 0; i < n_ma; i++) {
        m_a.push_back(ma_min + i * ma_step);
        sum_prob.push_back(0);
    }

    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    std::vector<TGraph*> grp;

    TRestAxionBufferGas* gas = new TRestAxionBufferGas();

    for (const auto& p : pair) {
        //   for (size_t i = 0; i < pair.size(); i++) {
        // Creates the gas and the axion field
        gas->SetGasDensity(gasName, p.second);
        ax->AssignBufferGas(gas);

        // Obtain the probability for each axion mass
        std::vector<double> prob;
        for (int j = 0; j < n_ma; j++) {
            prob.push_back(ax->GammaTransmissionProbability(m_a[j]));
            sum_prob[j] += prob[j];
        }

        TGraph* gr = new TGraph(n_ma, &m_a[0], &prob[0]);
        grp.push_back(gr);
    }

    // Computes the Vacuum probability
    TRestAxionField* ax_vac = new TRestAxionField();
    std ::vector<double> prob_vac;
    for (int j = 0; j < n_ma; j++) {
        prob_vac.push_back(ax_vac->GammaTransmissionProbability(m_a[j]));
    }

    // Computes the sum of all the probabilities
    if (vacuum == true) {
        for (int i = 0; i < n_ma; i++) {
            sum_prob[i] += prob_vac[i];
        }
    }

    ///// PLOTS /////

    // Plot the density scan
    grp[0]->SetLineColor(kBlue - 3);
    grp[0]->SetTitle("Transmission probability as a function of the axion mass");
    grp[0]->GetXaxis()->SetTitle("m_{a} [eV]");
    grp[0]->GetYaxis()->SetTitle("P_{ag}");
    grp[0]->GetXaxis()->SetLimits(ma_min, ma_max);
    double ylim = 3.5e-18;
    grp[0]->GetYaxis()->SetRangeUser(0, ylim);
    grp[0]->Draw("AL");

    for (const auto g : grp) {
        g->SetLineColor(kBlue - 3);
        g->Draw("SAME");
    }

    // PLot of the sum of all the probabilities
    TGraph* sum = new TGraph(n_ma, &m_a[0], &sum_prob[0]);
    sum->SetLineColor(kBlue + 3);
    sum->SetLineWidth(3);
    sum->Draw("SAME");

    // Plot of the vacuum probability
    TGraph* vac = new TGraph(n_ma, &m_a[0], &prob_vac[0]);
    vac->SetLineColor(kCyan + 1);
    vac->SetLineWidth(2);
    vac->Draw("SAME");

    // Plot of the vertical line
    TLine* verticalLine = new TLine(pair[0].second, c1->GetUymin(), pair[0].second, ylim);
    verticalLine->SetLineColor(kGreen - 3);
    verticalLine->SetLineWidth(2);
    verticalLine->Draw("same");

    // Plot of the legend
    TLegend* legend = new TLegend(0.9, 0.7, 0.48, 0.9);
    legend->AddEntry(grp[0], "P_{ag} for each mass", "l");
    legend->AddEntry(vac, "P_{ag} for vacuum", "l");
    legend->AddEntry(sum, "Sum of all the probs.", "l");
    legend->AddEntry(verticalLine, "m_{a} where P_{ag}^{vac} = max(P_{ag}^{vac}/2) ", "l");
    legend->Draw("same");
    c1->Draw();

    c1->Print("test.png");

    return 0;
}
