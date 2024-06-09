Double_t Eo = 0.5;  // keV
Double_t Ef = 10;   // keV

Int_t GenerateSignalComponents(std::string rmlFile, std::string name,
                               Double_t totalYears = 1.5, size_t firstSteps = 5, Double_t firstYears = 0.5) {

	Double_t totalExposureTime = totalYears * 365. * 12 * 3600;
	Double_t firstStepsExposure = firstYears * 365. * 12 * 3600;

    TRestAxionField field;
    std::vector<std::pair<Double_t, Double_t>> scanSteps =
        field.GetMassDensityScanning("He", 0.25, 20);  // Up to 0.25 eV

    TRestAxionHelioscopeSignal gasPhase(rmlFile.c_str(), name.c_str());
    std::filesystem::path filePath = rmlFile;
    std::filesystem::path newExtension = ".settings";
    std::string settingsFile = filePath.replace_extension(newExtension);
    std::filesystem::path rootExtension = ".root";
    std::string componentsFile = "Signals" + (std::string)filePath.replace_extension(rootExtension);

    std::filesystem::create_directories("output");
    settingsFile = "output/" + settingsFile;
    componentsFile = "output/" + componentsFile;

    TFile* f = TFile::Open((TString)componentsFile, "RECREATE");
    std::vector<Double_t> exposureTimes;
    for (size_t n = 0; n < scanSteps.size(); n++) {
        Double_t mass = scanSteps[n].first;
        Double_t density = scanSteps[n].second;

        gasPhase.SetName("P" + (TString)IntegerToString(n + 1));
        gasPhase.GetGas()->SetGasDensity("He", density);
        gasPhase.RegenerateHistograms();

        Double_t nGamma = gasPhase.GetSignalRate(mass, Eo, Ef);
        Double_t ksvzFactor = 3.75523 * mass;
        Double_t exposureTime = 0;
        if (n < firstSteps)
            exposureTime = firstStepsExposure / firstSteps;
        else
            exposureTime = TMath::Log(20.) / TMath::Power(ksvzFactor, 4) / nGamma;
        exposureTimes.push_back(exposureTime);

        gasPhase.Write("P" + (TString)IntegerToString(n + 1));
    }
    f->Close();

    std::cout << "Generated signals file: " << componentsFile << std::endl;

    Double_t generatedExposureTime = 0;
    for (size_t n = firstSteps; n < scanSteps.size(); n++) generatedExposureTime += exposureTimes[n];

    FILE* g = fopen(settingsFile.c_str(), "wt");
    for (size_t n = 0; n < firstSteps; n++) fprintf(g, "%lf\tP%d\n", exposureTimes[n], (int)(n + 1));

    for (size_t n = firstSteps; n < scanSteps.size(); n++)
        fprintf(g, "%lf\tP%d\n",
                (totalExposureTime - firstStepsExposure) * exposureTimes[n] / generatedExposureTime,
                (int)(n + 1));
    fclose(g);

    std::cout << "Generated settings file: " << settingsFile << std::endl;

    return 0;
}
