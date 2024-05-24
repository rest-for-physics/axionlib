Int_t GenerateSignalComponents() {
    TRestAxionField field;
    std::vector<std::pair<Double_t, Double_t>> scanSteps =
        field.GetMassDensityScanning("He", 0.25, 20);  // Up to 0.25 eV

    TRestAxionHelioscopeSignal gasPhase("BabyIAXO.rml", "GasSignal");
    gasPhase.RegenerateParametricNodes(0.001, 10, 1.02, true);

    TFile* f = TFile::Open("SignalComponents.root", "RECREATE");

    vacuumPhase.Write("Vacuum");

    FILE* g = fopen("GasPhase.settings", "wt");
    for (size_t n = 0; n < scanSteps.size(); n++) {
        // std::cout << "Resonant mass : " << scanSteps[n].first << std::endl;
        // std::cout << "Density : " << scanSteps[n].second << std::endl;

        gasPhase.SetName("P" + (TString)IntegerToString(n + 1));
        gasPhase.GetGas()->SetGasDensity("He", scanSteps[n].second);
        gasPhase.RegenerateHistograms();
        gasPhase.Write("P" + (TString)IntegerToString(n + 1));
        fprintf(g, "P%d\n", (int)(n + 1));
    }
    fclose(g);

    f->Close();

    return 0;
}
