///////////////////////////////////////////////////////
/// \file      AxionPlotResonances.C
/// \brief     This script plots the transmission probability as a function of the axion mass for a given gas
/// and axion energy until a maximum axion mass.
// Arguments by default are (in order):
// -ma_max = 0.1 eV (axion mass)
// -Ea = 4.2 keV (axion energy)
//  gasName = "He" (Gas name)
//  vacuum = true (Vacuum added to the sum or not)
//  n_ma = 100000 (Number of points to be plotted)
//
// Example of usage:
// restRoot -l
// .L AxionPlotResonances.C
// AxionPlotResonances(0.1, 4.2, "He", true, 100000)
//
///////////////////////////////////////////////////
/// \date       11/12/2023
/// \author     Francisco R. Candon
///////////////////////////////////////////////////

int REST_Axion_PlotResonances(double ma_max = 0.1, double Ea = 4.2, std::string gasName = "He",
                              Bool_t vacuum = true, int n_ma = 100000) {
    TRestAxionField* ax = new TRestAxionField();

    // m_a min to be considered
    double ma_min = 0;  // eV

    std::pair<std::vector<double>, std::vector<double>> pair =
        ax->GetMassDensityScanning(gasName, ma_max, Ea);
    std::vector<double> photonMass = pair.first;
    std::vector<double> density = pair.second;
    std::vector<double> m_a;
    std::vector<double> sum_prob;

    // Creates the vector of axion masses
    double ma_step = (ma_max - ma_min) / n_ma;
    for (int i = 0; i < n_ma; i++) {
        m_a.push_back(ma_min + i * ma_step);
        sum_prob.push_back(0);
    }

    int size = density.size();
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    TGraph* grp[size];

    for (int i = 0; i < size; i++) {
        // Creates the gas and the axion field
        TRestAxionBufferGas* gas = new TRestAxionBufferGas();
        gas->SetGasDensity(gasName, density[i]);
        TRestAxionField* ax = new TRestAxionField();
        ax->AssignBufferGas(gas);
        double photonMass = gas->GetPhotonMass(Ea);

        // Obatain the maximum probability for the axion mass
        double max_prob = ax->GammaTransmissionProbability(2.5, 10000, Ea, photonMass);

        // Obtain the probability for each axion mass
        std::vector<double> prob;
        for (int j = 0; j < n_ma; j++) {
            prob.push_back(ax->GammaTransmissionProbability(2.5, 10000, Ea, m_a[j]));
            sum_prob[j] += prob[j];
        }

        grp[i] = new TGraph(n_ma, &m_a[0], &prob[0]);

        delete ax;
        delete gas;
    }

    // Computes the Vacuum probability
    TRestAxionField* ax_vac = new TRestAxionField();
    std ::vector<double> prob_vac;
    for (int j = 0; j < n_ma; j++) {
        prob_vac.push_back(ax_vac->GammaTransmissionProbability(2.5, 10000, Ea, m_a[j]));
    }

    // Computes the sum of all the probabilities
    if (vacuum == true) {
        for (int i = 0; i < n_ma; i++) {
            sum_prob[i] += prob_vac[i];
        }
    }

    ///// PLOTS /////

    grp[0]->SetLineColor(kBlue - 3);
    grp[0]->SetTitle("Transmission probability as a function of the axion mass");
    grp[0]->GetXaxis()->SetTitle("m_a [eV]");
    grp[0]->GetYaxis()->SetTitle("P_{ag}");
    grp[0]->GetXaxis()->SetLimits(0, ma_max);
    double ylim = 3.5e-18;
    grp[0]->GetYaxis()->SetRangeUser(0, ylim);
    grp[0]->Draw("AL");
    for (int k = 1; k < size; k++) {
        grp[k]->SetLineColor(kBlue - 3);
        grp[k]->Draw("SAME");
    }
    TGraph* sum = new TGraph(n_ma, &m_a[0], &sum_prob[0]);
    sum->SetLineColor(kBlue + 3);
    sum->SetLineWidth(3);
    sum->Draw("SAME");

    TGraph* vac = new TGraph(n_ma, &m_a[0], &prob_vac[0]);
    vac->SetLineColor(kCyan + 1);
    vac->SetLineWidth(2);
    vac->Draw("SAME");

    TLine* verticalLine = new TLine(photonMass[0], c1->GetUymin(), photonMass[0], ylim);
    verticalLine->SetLineColor(kGreen - 3);  // 2 is red
    verticalLine->SetLineWidth(2);
    verticalLine->Draw("same");

    TLegend* legend = new TLegend(0.9, 0.7, 0.48, 0.9);
    legend->AddEntry(grp[0], "P_{ag} for each mass", "l");
    legend->AddEntry(vac, "P_{ag} for vacuum", "l");
    legend->AddEntry(sum, "Sum of all the probs.", "l");
    legend->AddEntry(verticalLine, "m_{a} where P_{ag}^{vac} = max(P_{ag}^{vac}/2) ", "l");
    legend->Draw("same");
    c1->Draw();

    return 0;
}
