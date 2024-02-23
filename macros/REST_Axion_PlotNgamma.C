#include "TRestAxionField.h"
#include "TRestAxionSolarFlux.h"
#include "TRestAxionSolarModel.h"
#include "TRestAxionSolarQCDFlux.h"
#include "TMath.h"
//*******************************************************************************************************
//*** Description: It computes the number of photons for each gas step for each axion mass

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
//*** - *E_a_min* = 0.0: Minimum energy (in keV) to be plotted.
//*** - *E_a_max* = 20.0: Maximum energy (in keV) to be plotted.
//***
//*** The generated plots are the results from `TRestAxionField::GetMassDensityScanning`,
//*** `TRestAxionField::GammaTransmissionFWHM`, `TRestAxionBufferGas::GetMassDensity`,
//*** `TRestAxionField::GammaTransmissionProbability` and `TRestAxionSolarQCDFlux`.

//***
//*** --------------
//*** Usage: restManager PlotNgamma [ma_max=0.1] [ma_min=0] [Ea=4.2] [Bmag=2.5] [Lcoh=10000] [gasName=He]

//*** [vacuum=true] [n_ma=10000] [E_a_min=0.0] [E_a_max=20.0]
//***
//*** Being all of them optional arguments.
//*** --------------
//*** Author: Fran Candón
//*******************************************************************************************************

int REST_Axion_PlotNgamma(double ma_max = 0.1, double ma_min = 0, double Ea = 4.2, double Bmag = 2.5,

                          double Lcoh = 10000, std::string gasName = "He", Bool_t vacuum = true,
                          int n_ma = 100, double E_a_min = 0.0, double E_a_max = 20.0) {

    // Factors for the INTEGRAL
    double exp_time = 60 * 60; // seconds
    double detect_eff = 0.2;
    double optic_eff = 0.8;
    double area = TMath::Pi()* 35; // cm^2

    // Creates the vector of axion masses
    std::vector<double> m_a;
    std::vector<double> n_photons;
    std::vector<double> n_photons_gas;

    double ma_step = (ma_max - ma_min) / n_ma;
    for (int i = 0; i < n_ma; i++) {
        m_a.push_back(ma_min + i * ma_step);
        n_photons.push_back(0);
        n_photons_gas.push_back(0);
    }

    // Creates the solar axion flux
    TRestAxionSolarQCDFlux* sFlux = new TRestAxionSolarQCDFlux("fluxes.rml", "Gianotti");
    sFlux->Initialize();
    TH1F* spec = sFlux->GetEnergySpectrum();
    spec->Rebin(20000 / n_ma);
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    spec->Draw();
    c1->Draw();

    // ******************* VACUUM PHASE ********************
    // Creates the axion field and the vacuum
    TRestAxionField* ax_vac = new TRestAxionField();
    ax_vac->SetMagneticField(Bmag);
    ax_vac->SetCoherenceLength(Lcoh);

    // Create a vector for saving the integrals that corresponds to the Number of photones

    double energy_step = (spec->GetBinCenter(1) - spec->GetBinCenter(0));

    TH1* convoluted = (TH1F*) spec->Clone();
    convoluted ->SetName("convoluted");

    for (int i = 0; i <= m_a.size(); ++i) {
        for (int j = 0; j <= spec->GetNbinsX(); ++j) {
            double en = convoluted->GetBinCenter(j+1);
            ax_vac->SetAxionEnergy(en);
            double Paxion = ax_vac->GammaTransmissionProbability(m_a[i]);
            convoluted->SetBinContent(j+1, spec->GetBinContent(j+1) * Paxion );
        }
        double integral = convoluted->Integral(convoluted->FindBin(0),convoluted ->FindBin(20));
        n_photons[i] = integral*optic_eff*exp_time*area*detect_eff; // *optEff*detectEff*area*time
    }
    delete ax_vac;
    delete convoluted;

    // ******************* GAS PHASE ********************

    // Creates the axion field and the gas
    TRestAxionField* ax = new TRestAxionField();
    TRestAxionBufferGas* gas = new TRestAxionBufferGas();
    ax->SetMagneticField(Bmag);
    ax->SetCoherenceLength(Lcoh);
    ax->SetAxionEnergy(Ea);

    // Obtain the pressure steps
    vector<std::pair<Double_t, Double_t>> pareja = ax->GetMassDensityScanning(gasName, ma_max, 20);
    std::vector<TGraph*> grp;

    // Loop that computes N_gamma for each gas pressure
    TH1* convolutedgas = (TH1F*) spec->Clone();
    convolutedgas ->SetName("convolutedgas");

    for (const auto& p : pareja) {
        // Creates the gas and the axion field
        cout << "Density: " << p.second << endl;
        gas->SetGasDensity(gasName, p.second);
        ax->AssignBufferGas(gas);

        for (int i = 0; i < n_ma; i++) {
            for (int j = 0; j <= spec->GetNbinsX(); ++j) {
            double en = convoluted->GetBinCenter(j+1);
            ax->SetAxionEnergy(en);
            double Paxion = ax->GammaTransmissionProbability(m_a[i]);
            convoluted->SetBinContent(j+1, spec->GetBinContent(j+1) * Paxion );
        }
            double integral = convolutedgas->Integral(convolutedgas->FindBin(0),convolutedgas ->FindBin(20));
            n_photons_gas[i] = integral*optic_eff*exp_time*area*detect_eff; // *optEff*detectEff*area*time
        }
        TGraph* gr_gas = new TGraph(n_ma, &m_a[0], &n_photons_gas[0]);
        grp.push_back(gr_gas);
    }

    // ******************* PLOTS ********************
    delete ax_vac;
    cout << "SIZE:" << grp.size() << endl;
    TCanvas* c3 = new TCanvas("c3", "c3", 800, 600);
    TGraph* gr = new TGraph(n_ma, &m_a[0], &n_photons[0]);

    gr->SetLineColor(kRed);

    gr->SetTitle("Number of Photons vs Axion Mass");
    gr->GetXaxis()->SetTitle("Axion Mass (eV)");
    gr->GetYaxis()->SetTitle("Number of Photons");

    // Set the y axis to log scale
    c3->SetLogy();
    c3->SetLogx();
    gr->SetLineWidth(3);
    gr->Draw("AC");

    // Plot of the gas steps
    for (const auto g : grp) {
        // Set a line of a different color for each interation
        g->SetLineColor(kBlue);
        g->Draw("SAME");
    }
    c3->Draw();

    return 0;
}
