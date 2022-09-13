#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TRestAxionWolterOptics.h"
#include "TRestTask.h"
#include "TRestTools.h"
#include "rapidcsv.h"

#ifndef RestTask_Axion_XMMAngleEff
#define RestTask_Axion_XMMAngleEff

//*******************************************************************************************************
//*** Description: It creates an efficiency curve from different runs under different angles. 
//*** It's also possible to compare data from a csv file with this curve.
//*** 
//***
//*** --------------
//*** Usage: restManager XMMAngleEffCSV [pathAndPattern = "/trueWolter/OpticsBench_Yaw_*_Dev_*_BabyIAXO_Run*.root", csvFile = ""]
//***
//*** Where the optional arguments are:
//***
//*** - The path and pattern name of the files that should be read in. The * is in the places where there is a difference in the file names.
//***
//*** - The path and name of the csv file with data to compare to the data from the prvious files. If not given it will be left out.
//*******************************************************************************************************
Int_t REST_Axion_XMMAngleEffCSV(TString pathAndPattern = "./trueWolter/OpticsBench_Yaw_*_Dev_*_BabyIAXO_Run*.root", TString csvFile = "") {
    vector<string> files = TRestTools::GetFilesMatchingPattern((string)pathAndPattern);
    
    if (files.size() == 0) {
        RESTError << "Files not found!" << RESTendl;
        return -1;
    }

    TRestRun* run = new TRestRun();
    std::vector<float> efficiencies(files.size());
    std::vector<float> angles(files.size());
    std::vector<float> deviations(files.size());
    Double_t eff0 = 0;
    std::vector<float> totalCounts(files.size());

    for (unsigned int n = 0; n < files.size(); n++) {
        run->OpenInputFile(files[n]);

        run->PrintMetadata();

        Int_t obsID = run->GetAnalysisTree()->GetObservableID("optics_efficiency");
        if (obsID == -1) {
            RESTError << "No observable \"" << "optics_efficiency" << "\" in file " << files[n] << RESTendl;
            continue;
        }
        Double_t eff = 0;
        Double_t counts = 0;
        Double_t angle = std::stod(run->GetMetadataMember("optics::fYaw")) * 180. * 60. / 3.1415;
        Double_t dev = std::stod(run->GetMetadataMember("deviation::fDevAngle")) * 180. / 3.1415;

        for (int i = 0; i < run->GetEntries(); i++) {
            run->GetAnalysisTree()->GetBranch((TString)"optics_efficiency")->GetEntry(i);
            Double_t value = run->GetAnalysisTree()->GetDblObservableValue(obsID);
            //std::cout << value << " for i " << i << std::endl;
            counts++;
            eff += value;
        }
        std::cout.precision(10);
        std::cout << "eff: " << eff << " angle: " << angle << std::endl;
        efficiencies[n] = eff;
        angles[n] = angle;
        totalCounts[n] = counts;
        if (angle == 0 && dev == 0) eff0 = eff; // and dev 0
    }
    Double_t angleMax = *max_element(angles.begin(), angles.end());


    delete run;
    std::vector<float> effError;
    for (int n = 0; n < files.size(); n++) {
        std::cout << "eff: " << efficiencies[n] << " total " << totalCounts[n] << std::endl;
        effError.push_back(TMath::Sqrt(efficiencies[n] + efficiencies[n] / eff0) / eff0);
        std::cout << "eff: " << efficiencies[n] << " error " << (TMath::Sqrt(efficiencies[n] + efficiencies[n] / eff0) / eff0) << std::endl;
        efficiencies[n] = efficiencies[n] / eff0;   
    }
    

    TCanvas c("", "", 1200, 800);

    TGraphErrors effGraph((Int_t)files.size(), &angles[0], &efficiencies[0], nullptr, &effError[0]);
    effGraph.SetName("effGraph");
    effGraph.SetTitle("Relative Efficiency Comparison");
    effGraph.GetXaxis()->SetTitle("Angle to optical axis [arcmin]");
    effGraph.GetYaxis()->SetTitle("Efficiency");

    effGraph.SetMarkerColor(kRed);
    effGraph.SetMarkerStyle(21);
    effGraph.Draw("AP");
    
    /// Add the csv file if given
    if (csvFile != ""){
        rapidcsv::Document doc(csvFile);

        std::vector<float> anglesExp = doc.GetColumn<float>("angle[arcmin]");
        std::vector<float> efficienciesExp = doc.GetColumn<float>("effectiveArea");
        Double_t effExp0 = 0;

        for (unsigned int n = 0; n < anglesExp.size(); n++) {
            if (anglesExp[n] < 0.01 && anglesExp[n] > 0) effExp0 = efficienciesExp[n];
        }
        std::cout << "effExp: " << effExp0 << std::endl;
        for (unsigned int n = 0; n < anglesExp.size(); n++) {
            efficienciesExp[n] = efficienciesExp[n] / effExp0;
        }

        TGraph expGraph(anglesExp.size(), &anglesExp[0], &efficienciesExp[0]);
        expGraph.SetName("expGraph");
        expGraph.SetMarkerColor(kBlue);
        expGraph.SetMarkerStyle(2);
        expGraph.Draw("L");
    }

    effGraph.GetYaxis()->SetLimits(0.0, angleMax);
    effGraph.GetHistogram()->SetMaximum(1);
    effGraph.GetHistogram()->SetMinimum(0);


    /*auto legend = new TLegend(0.1, 0.7, 0.2, 0.4);
    legend->AddEntry("effGraph", "expGraph", "lep");
    legend->Draw();*/

    c.Print("./trueWolter/XMM_Efficiency.png");

    return 0;
}
#endif
