#include "TRestAxionSolarModel.h"

int REST_Axion_Demo(std::string configFile = "examples/axionlib.xml") {
    // Initialise TRestAxionSolarModel class from xml file.
    TRestAxionSolarModel test (configFile.c_str(), "SolarAxionLibTest");

    // Define some energies
    std::vector<double> ergs = {0.1,1.0,2.5,5.0,7.5,10.0};
    // Calculate the Primakoff (over the whole Sun) and axion-electron-interaction (only for 50%
    // of the Solar disc) spectra. Print out the results.
    std::vector<double> res = test.GetSolarAxionFluxGAGamma(ergs);
    std::cout << "[ "; for (auto x = res.begin(); x != res.end(); ++x) { std::cout << *x << " "; }; std::cout << "]" << std::endl;
    res = test.GetSolarAxionFluxGAE(ergs, 0.5);
    std::cout << "[ "; for (auto x = res.begin(); x != res.end(); ++x) { std::cout << *x << " "; }; std::cout << "]" << std::endl;

    return 0;
}
