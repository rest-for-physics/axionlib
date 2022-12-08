#include <TH1D.h>
#include <TRestRun.h>

Double_t probTolerance = 1.e-21;
Double_t fieldTolerance = 0.01;

Int_t Validate(Double_t prob = 5.9978e-19, Double_t fieldAverage = 1.631) {
    TRestRun* run = new TRestRun("AxionPhotonProbability.root");

    if (run->GetEntries() != 100) {
        std::cout << "Error. Number of entries is not 10000!" << std::endl;

        return 1;
    }

    run->GetAnalysisTree()->Draw("axionPhoton_probability", "axionPhoton_probability");

    TH1D* h = (TH1D*)run->GetAnalysisTree()->GetHistogram();
    if (h == nullptr) {
        std::cout << "Problems generating probability histogram" << std::endl;
        return 2;
    }

    Double_t probability = h->Integral() / run->GetEntries();
    std::cout << "Average axion-photon probability: " << probability << std::endl;

    if (probability < prob - probTolerance || probability > prob + probTolerance) {
        std::cout << "Wrong axion-photon probability!" << std::endl;
        return 3;
    }

    run->GetAnalysisTree()->Draw("axionPhoton_fieldAverage", "axionPhoton_fieldAverage");

    TH1D* h2 = (TH1D*)run->GetAnalysisTree()->GetHistogram();
    if (h2 == nullptr) {
        std::cout << "Problems generating field average histogram" << std::endl;
        return 4;
    }

    Double_t field = h2->Integral() / run->GetEntries();
    std::cout << "Average magnetic field: " << field << " T" << std::endl;

    if (field < fieldAverage - fieldTolerance || field > fieldAverage + fieldTolerance) {
        std::cout << "Wrong average field!" << std::endl;
        return 5;
    }

    delete run;

    return 0;
}