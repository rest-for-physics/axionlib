#include <iostream>

#include "TRestAxionSolarModel.h"

int REST_Axion_Testing(string configFile = "example.rml") {
    // Initialise TRestAxionSolarModel class from xml file.
    TRestAxionSolarModel* test = new TRestAxionSolarModel(configFile.c_str(), "FirstTest");
    // Print metadata; however:
    // This is currently returning default E1, E2 etc. They only get set after calling GetSolarAxionFlux.
    // Redefine constructor?
    test->PrintMetadata();
    Double_t res1 = test->GetSolarAxionFlux(1.0, 8.0, 1.0, 1.0);
    // Manually change model. Probably need to revise this.
    test->SetSolarAxionSolarModel("arXiv_1302.6283_Primakoff");
    test->PrintMetadata();
    Double_t res2 = test->GetSolarAxionFlux(1.0, 8.0, 1.0, 1.0);

    std::cout << "Compare flux from two models: " << res1 << " vs " << res2 << "." << std::endl;

    delete test;
    // Re-initilise the class, using Maurizios table (currently only possible via destructor).
    // Currently, the results will be written automatically from the GetSolarAxionFlux(...) function. Probably
    // not ideal for large files. Adjust?
    test = new TRestAxionSolarModel(configFile.c_str(), "SecondTest");
    test->PrintMetadata();
    // Currrently, we CANNOT use GetSolarAxionFlux(...) with tables.
    // res1 = test.GetSolarAxionFlux(1.0, 8.0, 1.0, 1.0);
    // std::cout << "Flux from Maurizio's table: " << res1 << "." << std::endl;

    return 0;
}
