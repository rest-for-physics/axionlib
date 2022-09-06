#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TRestAxionWolterOptics.h"
#include "TRestTask.h"
#include "TRestTools.h"

#ifndef RestTask_Axion_XMMAngleEff
#define RestTask_Axion_XMMAngleEff

//*******************************************************************************************************
//*** Description: It creates an efficiency curve from different runs under different angles
//***
//*** 
//***
//*** --------------
//*** Usage: restManager XMMAngleEff [pathAndPattern = "/trueWolter/OpticsBench_Yaw_*_Dev_*_BabyIAXO_Run*.root"]
//***
//*** Where the optional arguments are:
//***
//*** - 
//***
//*******************************************************************************************************
Int_t REST_Axion_XMMAngleEff(TString pathAndPattern = "./trueWolter/OpticsBench_Yaw_*_Dev_*_BabyIAXO_Run*.root") {
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

    for (unsigned int n = 0; n < files.size(); n++) {
        run->OpenInputFile(files[n]);

        run->PrintMetadata();

        Int_t obsID = run->GetAnalysisTree()->GetObservableID("optics_efficiency");
        if (obsID == -1) {
            RESTError << "No observable \"" << "optics_efficiency" << "\" in file " << files[n] << RESTendl;
            continue;
        }
        Double_t eff = 0;
        Double_t angle = std::stod(run->GetMetadataMember("optics::fYaw")) * 180. / 3.1415;
        Double_t dev = std::stod(run->GetMetadataMember("deviation::fDevAngle")) * 180. / 3.1415;
        
        for (int i = 0; i < run->GetEntries(); i++) {
            run->GetAnalysisTree()->GetBranch((TString)"optics_efficiency")->GetEntry(i);
            Double_t value = run->GetAnalysisTree()->GetDblObservableValue(obsID);
            //std::cout << value << " for i " << i << std::endl;
            eff += value;
        }
        std::cout.precision(10);
        std::cout << "eff: " << eff << " angle: " << angle << std::endl;
        efficiencies[n] = eff;
        angles[n] = angle;
        if (angle == 0 && dev == 0) eff0 = eff; // and dev 0
    }

    for (unsigned int n = 0; n < files.size(); n++) {
        efficiencies[n] = efficiencies[n] / eff0;
    }


    std::cout << "eff: " << efficiencies[0] << " " << eff0 << std::endl;
    
    delete run;
    /*  std::vector<Double_t> angle, entranceError, exitError, finalError;
      for (int n = 0; n < divs; n++) {
          angle.push_back((((Double_t)n) + 0.5) * angleMax / divs);
          std::cout << "n : " << totalCounts[n] << " - " << entranceEfficiency[n] / totalCounts[n] << " - "
                    << exitEfficiency[n] / totalCounts[n] << " - " << finalEfficiency[n] / totalCounts[n]
                    << std::endl;
          entranceError.push_back(TMath::Sqrt(entranceEfficiency[n]) / totalCounts[n]);
          entranceEfficiency[n] = entranceEfficiency[n] / totalCounts[n];

          exitError.push_back(TMath::Sqrt(exitEfficiency[n]) / totalCounts[n]);
          exitEfficiency[n] = exitEfficiency[n] / totalCounts[n];

          finalError.push_back(TMath::Sqrt(finalEfficiency[n]) / totalCounts[n]);
          finalEfficiency[n] = finalEfficiency[n] / totalCounts[n];
      }*/

      TCanvas c("", "", 1200, 800);

      TGraph effGraph(files.size(), &angles[0], &efficiencies[0]);
      effGraph.SetName("effGraph");
      effGraph.SetTitle("Geometrical Acceptance");
      effGraph.GetXaxis()->SetTitle("Angle to optical axis [degrees]");
      effGraph.GetYaxis()->SetTitle("Efficiency");
      effGraph.SetMarkerColor(kRed);
      effGraph.SetMarkerStyle(21);
      effGraph.Draw("AP");

      /*effGraph.GetXaxis()->SetLimits(0, angleMax);
      effGraph.GetHistogram()->SetMaximum(1);
      effGraph.GetHistogram()->SetMinimum(0);


      auto legend = new TLegend(0.1, 0.7, 0.48, 0.9);
      legend->AddEntry("effGraph", "REST true Wolter", "lep");
      legend->Draw();*/

      c.Print("./trueWolter/XMM_Efficiency.png");

    return 0;
}
#endif
