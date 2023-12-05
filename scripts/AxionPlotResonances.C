int AxionPlotResonances() {
    cout << "Hello world!" << endl;
    TRestAxionField* ax = new TRestAxionField();
    double ma_max = 0.1;
    double Ea = 4.2;
    int n_ma = 100000;
    int vacio = 1;  // 1 if you want to include the vacuum probability, 0 if not
    std::string gasName = "He";
    std::pair<std::vector<double>, std::vector<double>> pair = ax->GetMassDensityScanning(gasName,ma_max,Ea);
    std::vector<double> photonMass = pair.first;
    std::vector<double> density = pair.second;
    std::vector<double> m_a;
    std::vector<double> sum_prob;
    double ma_step = (ma_max - 0.00025) / n_ma;
    for (int i = 0; i < n_ma; i++) {
        m_a.push_back(0.00025 + i * ma_step);
        sum_prob.push_back(0);
    }

    int size = density.size();
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    TGraph* grp[size];
    for (int i = 0; i < size; i++) {
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
        grp[i] = new TGraph(size, &m_a[0], &prob[0]);
        delete ax;
        delete gas;
    }
    // Computes the Vacuum probability
    TRestAxionField* ax_vac = new TRestAxionField();
    std ::vector<double> prob_vac;
    for (int j = 0; j < n_ma; j++) {
        prob_vac.push_back(ax_vac->GammaTransmissionProbability(2.5, 10000, Ea, m_a[j]));
    }
    if (vacio == 1) {
        for (int i = 0; i < n_ma; i++) {
            sum_prob[i] += prob_vac[i];
        }
    }

    grp[0]->SetLineColor(kRed);
    grp[0]->SetTitle("Transmission probability as a function of the axion mass");
    grp[0]->GetXaxis()->SetTitle("m_a [eV]");
    grp[0]->GetYaxis()->SetTitle("Transmission probability");
    grp[0]->GetXaxis()->SetLimits(0, ma_max);
    grp[0]->GetYaxis()->SetRangeUser(0, 3.5e-18);
    grp[0]->Draw("AL");
    for (int k = 1; k < size; k++) {
        grp[k]->SetLineColor(kRed);
        grp[k]->Draw("SAME");
    }
    TGraph* gr_vac = new TGraph(n_ma, &m_a[0], &sum_prob[0]);
    gr_vac->SetLineColor(kBlack);
    gr_vac->SetLineWidth(4);
    gr_vac->Draw("SAME");

    TGraph* grp_vac = new TGraph(n_ma, &m_a[0], &prob_vac[0]);
    grp_vac->SetLineColor(kBlue);
    grp_vac->SetLineWidth(2);
    grp_vac->Draw("SAME");

    TLine* verticalLine = new TLine(photonMass[0], c1->GetUymin(), photonMass[0], c1->GetUymax());
    verticalLine->SetLineColor(2);  // 2 is red
    verticalLine->SetLineWidth(2);
    verticalLine->Draw("same");
    c1->Draw();
    return 0;
}
