#include <TH1D.h>
#include <TRestRun.h>

Int_t ValidateTransmission(std::string fname, Double_t lowEff = 0.07, Double_t highEff = 0.09) {
    std::cout << "Filename : " << fname << std::endl;
    TRestRun* run = new TRestRun(fname);

    if (run->GetEntries() != 10000) {
        std::cout << "Error. Number of entries is not 10000!" << std::endl;

        return 1;
    }

    run->GetAnalysisTree()->Draw("bore_transmission", "bore_transmission");

    TH1D* h = (TH1D*)run->GetAnalysisTree()->GetHistogram();

    if (h == nullptr) {
        std::cout << "Problems generating histogram" << std::endl;
        return 2;
    }

    Double_t integral = h->Integral();
    std::cout << "Photon detection efficiency : " << integral / run->GetEntries() << std::endl;
    Double_t efficiency = integral / run->GetEntries();
    if (efficiency < lowEff || efficiency > highEff) {
        std::cout << "Wrong number of photons detected!" << std::endl;
        return 3;
    }

    delete run;

    return 0;
}
