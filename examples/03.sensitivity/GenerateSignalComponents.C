Int_t GenerateSignalComponents(std::string rmlFile, std::string name) {
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

    FILE* g = fopen(settingsFile.c_str(), "wt");
    TFile* f = TFile::Open((TString)componentsFile, "RECREATE");
    for (size_t n = 0; n < scanSteps.size(); n++) {
        gasPhase.SetName("P" + (TString)IntegerToString(n + 1));
        gasPhase.GetGas()->SetGasDensity("He", scanSteps[n].second);
        gasPhase.RegenerateHistograms();
        gasPhase.Write("P" + (TString)IntegerToString(n + 1));
        fprintf(g, "P%d\n", (int)(n + 1));
    }
    fclose(g);

    f->Close();

    std::cout << "Generated signals file: " << componentsFile << std::endl;

    return 0;
}
