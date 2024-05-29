<<<<<<< Updated upstream
Int_t GenerateSignalComponents(std::string rmlFile, std::string name) {
    TRestAxionField field;
    std::vector<std::pair<Double_t, Double_t>> scanSteps =
        field.GetMassDensityScanning("He", 0.25, 20);  // Up to 0.25 eV

    TRestAxionHelioscopeSignal gasPhase(rmlFile.c_str(), name.c_str());
=======
Double_t Eo = 0.5; //keV
Double_t Ef = 10; //keV

Int_t GenerateSignalComponents(std::string rmlFile, std::string name, Double_t totalExposureTime=3*300*12*3600, size_t skipSteps = 1)
{
	TRestAxionField field;
	std::vector<std::pair<Double_t, Double_t>> scanSteps = field.GetMassDensityScanning( "He", 0.25, 20 ); // Up to 0.25 eV

	TRestAxionHelioscopeSignal gasPhase( rmlFile.c_str(), name.c_str());
	//gasPhase.RegenerateParametricNodes( startMass, upperMass, step, true );
>>>>>>> Stashed changes

    std::filesystem::path filePath = rmlFile;
    std::filesystem::path newExtension = ".settings";
    std::string settingsFile = filePath.replace_extension(newExtension);
    std::filesystem::path rootExtension = ".root";
    std::string componentsFile = "Signals" + (std::string)filePath.replace_extension(rootExtension);

    std::filesystem::create_directories("output");
    settingsFile = "output/" + settingsFile;
    componentsFile = "output/" + componentsFile;

<<<<<<< Updated upstream
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
=======
	TFile *f = TFile::Open( (TString) componentsFile, "RECREATE" );
	std::vector <Double_t> exposureTimes;
	for( size_t n = skipSteps; n < scanSteps.size(); n++ )
	{
		Double_t mass = scanSteps[n].first;
		std::cout << "Mass: " << mass << std::endl;
		Double_t density = scanSteps[n].second;
		std::cout << "Setting density: " << density << std::endl;

		gasPhase.SetName( "P" + (TString) IntegerToString(n+1) );
		gasPhase.GetGas()->SetGasDensity("He", density);
		gasPhase.RegenerateHistograms();

		Double_t nGamma = gasPhase.GetSignalRate(mass, Eo, Ef);
		std::cout << "nGamma: " << nGamma << std::endl;
		std::cout << "Mass: " << mass << std::endl;
		Double_t ksvzFactor = 3.75523 * mass;
		std::cout << "Log(20)/ksvzFactor : " << TMath::Log(20)/TMath::Power( ksvzFactor, 4 ) << std::endl;
		Double_t exposureTime = TMath::Log(20.) / TMath::Power(ksvzFactor, 4) / nGamma;
		exposureTimes.push_back( exposureTime );
		std::cout << "Exposure time: " << exposureTime << std::endl;

		gasPhase.Write( "P" + (TString) IntegerToString(n+1) );
	}
	f->Close();
>>>>>>> Stashed changes

    std::cout << "Generated signals file: " << componentsFile << std::endl;

<<<<<<< Updated upstream
    return 0;
=======
	Double_t generatedExposureTime = 0;
	for( size_t n = skipSteps; n < scanSteps.size(); n++ )
		generatedExposureTime += exposureTimes[n-skipSteps];


	FILE *g = fopen(settingsFile.c_str(), "wt");
	for( size_t n = skipSteps; n < scanSteps.size(); n++ )
		fprintf( g, "%lf\tP%d\n", totalExposureTime * exposureTimes[n-skipSteps] / generatedExposureTime, (int) (n+1) );
	fclose(g);

	std::cout << "Generated settings file: " << settingsFile << std::endl;

	return 0;
>>>>>>> Stashed changes
}
