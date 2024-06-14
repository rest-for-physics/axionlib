#include <filesystem>
#include <iostream>
#include <vector>

#include "TGeoManager.h"
#include "TRestTask.h"
#include "TSystem.h"
using namespace std;
namespace fs = std::filesystem;

#ifndef RESTTask_ComputingTimesPerMass
#define RESTTask_ComputingTimesPerMass

//*******************************************************************************************************
//***
//*** Description: 
//***
//*** --------------
//***
//*** Bla bla bla
//*** --------------
//*** Usage: restManager ComputingTimesPerMass "/full/path/file_*pattern*.root"
//*** --------------
//***
//*******************************************************************************************************
Int_t REST_Axion_ComputingTimesPerMass(std::string namePattern) {
	std::map<double,int> statistics;
	std::map<double,double> computingTime;

    std::vector<std::string> b = TRestTools::GetFilesMatchingPattern(namePattern, true);

    Double_t totalTime = 0;
    int cont = 0;
    for (int i = 0; i < b.size(); i++) {
        string filename = b[i];
		std::cout << filename << std::endl;
        cont++;
        TRestRun* run = new TRestRun();
        TRestAxionGeneratorProcess* gen;

        TFile* f = new TFile(filename.c_str());

        /////////////////////////////
        // Reading metadata classes

        TIter nextkey(f->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)nextkey())) {
            string className = key->GetClassName();
            if (className == "TRestRun") {
                run = (TRestRun*)f->Get(key->GetName());
            }
            if (className == "TRestAxionGeneratorProcess") {
                gen = (TRestAxionGeneratorProcess*)f->Get(key->GetName());
            }
        }

		double mass = StringTo2DVector( gen->GetDataMemberValue("fAxionMassRange") ).X();


        if (run->GetEndTimestamp() == 0 || run->GetRunLength() < 0 || run->GetEntries() == 0) {
			continue;
        } else if (run->GetRunLength() > 0)
		{
			mass = mass * 1000.;
			statistics[mass] += run->GetEntries();
			computingTime[mass] += run->GetRunLength() / 3600.;
		}

        delete run;

        f->Close();
    }

	std::cout << "Contents of the stats map:" << std::endl;
	for (const auto& [mass, stats] : statistics) {
		std::cout << "Key: " << mass << ", Value: " << stats << std::endl;
	}
	std::cout << "Contents of the computing map:" << std::endl;
	for (const auto& [mass, cTime] : computingTime) {
		std::cout << "Key: " << mass << ", Value: " << cTime << std::endl;
	}


    return 0;
}
#endif
