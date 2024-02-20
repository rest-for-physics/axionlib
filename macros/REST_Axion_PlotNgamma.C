#include "TRestAxionSolarQCDFlux.h"
#include "TRestAxionSolarFlux.h"
#include "TRestAxionSolarModel.h"
#include "TRestAxionField.h"
//*******************************************************************************************************
//*** Description: It computes the number of photones for each gass step for eah axion mass
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
//*** `TRestAxionField::GammaTransmissionFWHM`, `TRestAxionBufferGas::GetMassDensity`, `TRestAxionField::GammaTransmissionProbability` and `TRestAxionSolarQCDFlux`.	
//***
//*** --------------
//*** Usage: restManager PlotResonances [ma_max=0.1] [ma_min=0] [Ea=4.2] [Bmag=2.5] [Lcoh=10000] [gasName=He]
//*** [vacuum=true] [n_ma=10000] [E_a_min=0.0] [E_a_max=20.0]
//***
//*** Being all of them optional arguments.
//*** --------------
//*** Author: Fran Candón
//*******************************************************************************************************

int REST_Axion_PlotNgamma2(double ma_max = 0.1, double ma_min = 0, double Ea = 4.2, double Bmag = 2.5,
                              double Lcoh = 10000, std::string gasName = "He", Bool_t vacuum = true,
                              int n_ma = 100,double E_a_min = 0.0, double E_a_max = 20.0){

    
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
    TRestAxionSolarQCDFlux *sFlux = new TRestAxionSolarQCDFlux("fluxes.rml", "Gianotti");
    sFlux->Initialize();
    TH1F *spec1 = sFlux->GetEnergySpectrum();
    spec1->Rebin(20000/ n_ma);
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    spec1->Draw();
    c1->Draw();
    
    // ******************* VACUUM PHASE ********************
    // Creates the axion field and the vacuum
    TRestAxionField* ax_vac = new TRestAxionField();
    ax_vac->SetMagneticField(Bmag);
    ax_vac->SetCoherenceLength(Lcoh);

    // Create a vector for saving the integrals that corresponds to the Number of photones
    
    double energy_step = (E_a_max - E_a_min) / n_ma;

    for (int i = 0; i <= m_a.size(); ++i) {
        // TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
        TH1F* vac = new TH1F("histogram", "histogram", n_ma, E_a_min, E_a_max);
        cout << m_a[i] << endl;
        double E = 0;
        for (int j = 0; j <= m_a.size(); ++j) {
            E = E_a_min + ((j) * energy_step - energy_step / 2);
            ax_vac->SetAxionEnergy(E);
            vac->Fill(E, ax_vac->GammaTransmissionProbability(m_a[i]));
            //cout << "j: " << j << endl;
        }
        TH1F *ng = new TH1F();
        *ng = (*spec1)*(*vac); 
        double integral = ng->Integral(ng->FindBin(0), ng->FindBin(20));
        n_photons[i] = integral*0.2*0.8*14*60*60;    
        delete vac;
        delete ng;
    }
    delete ax_vac;

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
    for (const auto& p : pareja) {
        // Creates the gas and the axion field
        cout << "Density: " << p.second << endl;
        gas->SetGasDensity(gasName, p.second);
        ax->AssignBufferGas(gas);

        for (int j = 0; j < n_ma; j++) {
            TH1F* gas = new TH1F("histogram", "histogram", n_ma, E_a_min, E_a_max);
            double E = 0;
            for (int k = 10; k <= m_a.size(); ++k) {
                E = E_a_min + ((k) * energy_step);
                ax->SetAxionEnergy(E);
                gas->Fill(E, ax->GammaTransmissionProbability(m_a[j]));
            }
            TH1F *ng = new TH1F();
            *ng = (*spec1)*(*gas);
            double integral = ng->Integral(ng->FindBin(0), ng->FindBin(20));
            n_photons_gas[j] = integral*0.2*0.8*14*60*60;
            delete gas;
            delete ng;
        }
        TGraph* gr_gas = new TGraph(n_ma, &m_a[0], &n_photons_gas[0]);
        grp.push_back(gr_gas);
    }




    // ******************* PLOTS ********************
    delete ax_vac;
    cout << "SIZE:" <<grp.size() << endl;
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
