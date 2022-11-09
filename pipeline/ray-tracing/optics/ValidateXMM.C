#include <TH1D.h>
#include <TRestRun.h>

Int_t ValidateXMM(std::string fname) {
    std::cout << "Filename : " << fname << std::endl;
    TRestRun* run = new TRestRun(fname);

    if (run->GetEntries() != 1000) {
        std::cout << "Error. Number of entries is not 1000!" << std::endl;

        return 1;
    }

    run->GetAnalysisTree()->Draw("optics_efficiency", "optics_efficiency");

    run->GetAnalysisTree()->PrintObservables();

    TH1D* h = (TH1D*)run->GetAnalysisTree()->GetHistogram();

    if (h == nullptr) {
        std::cout << "Problems generating histogram" << std::endl;
        return 2;
    }

    Double_t integral = h->Integral();
    std::cout << "Efficiency : " << integral / run->GetEntries() << std::endl;
    if (integral < 150 || integral > 250) {
        std::cout << "Optics efficiency is not within the expected range!" << std::endl;
        return 3;
    }

    delete run;

    return 0;
}
