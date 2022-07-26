#include <TRestRun.h>

Int_t ValidateXMM(std::string fname) {
    std::cout << "Filename : " << fname << std::endl;
    TRestRun* run = new TRestRun(fname);

    if (run->GetEntries() != 1000) {
        std::cout << "Error. Number of entries is not 1000!" << std::endl;

        return 1;
    }

    delete run;

    return 0;
}
