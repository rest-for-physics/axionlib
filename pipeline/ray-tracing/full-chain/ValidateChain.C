#include <TH1D.h>
#include <TRestRun.h>

Int_t ValidateChain(std::string fname) {
    std::cout << "Filename : " << fname << std::endl;
    TRestRun* run = new TRestRun(fname);

    if (run->GetEntries() != 100) {
        std::cout << "Error. Number of entries is not 100!" << std::endl;

        return 1;
    }

    run->GetAnalysisTree()->Draw("axionPhoton_probability", "axionPhoton_probability");

    run->GetEntry(10);
    run->GetAnalysisTree()->PrintObservables();

    TH1D* h = (TH1D*)run->GetAnalysisTree()->GetHistogram();

    if (h == nullptr) {
        std::cout << "Problems generating axionPhoton histogram" << std::endl;
        return 2;
    }

    Double_t integral = h->Integral() / run->GetEntries();
    std::cout << "Average axion-photon probability : " << integral << std::endl;

    run->GetAnalysisTree()->Draw("optics_efficiency", "optics_efficiency");
    TH1D* l = (TH1D*)run->GetAnalysisTree()->GetHistogram();
    if (l == nullptr) {
        std::cout << "Problems generating optics_efficiency histogram" << std::endl;
        return 6;
    }

    Double_t integral3 = l->Integral() / run->GetEntries();
    std::cout << "Average optics efficiency : " << integral3 << std::endl;

    run->GetAnalysisTree()->Draw("window_transmission", "window_transmission");
    TH1D* g = (TH1D*)run->GetAnalysisTree()->GetHistogram();
    if (g == nullptr) {
        std::cout << "Problems generating window_transmission histogram" << std::endl;
        return 4;
    }

    Double_t integral2 = g->Integral() / run->GetEntries();
    std::cout << "Average window transmission : " << integral2 << std::endl;

    if (integral < 9.71489e-24 || integral > 9.91489e-24) {
        std::cout << "Axion-photon probability is not within the expected range!" << std::endl;
        return 3;
    }

    if (integral2 < 0.125 || integral2 > 0.155) {
        std::cout << "Average window transmission is not within the expected range!" << std::endl;
        return 5;
    }

    if (integral3 < 0.045 || integral3 > 0.055) {
        std::cout << "Optics efficiency is not within the expected range!" << std::endl;
        return 7;
    }
    delete run;

    return 0;
}
