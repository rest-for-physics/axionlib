#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TRestAxionWolterOptics.h"
#include "TRestTask.h"

#ifndef RestTask_Axion_XMMEfficiency
#define RestTask_Axion_XMMEfficiency

//*******************************************************************************************************
//*** Description: It creates an efficiency curve as a function of angle respect to axis.
//***
//*** It determines the fraction of photons that:
//***
//*** 1. Enter through the entrance plane (that does not hit the mirror thickness or spider structure),
//*** only photons inside `GetRadialLimits` will be considered in the calculation,
//*** 2. Once they enter the optics, they manage to go though the exit plane, and thus reach the
//*** spot,
//*** 3. and the fraction of photons reaching the spot convoluted with the reflectivity.
//***
//*** --------------
//*** Usage: restManager XMMEfficiency [maxDev=0]
//***
//*** Where the optional arguments are:
//***
//*** - *maxDev*: Maximum deviation (in radians) respect to optical axis to be given to the MC-photons.
//***
//*******************************************************************************************************
Int_t REST_Axion_XMMEfficiency(Double_t maxDev = 0.0005) {
    Int_t nPhotons = 100000;
    Double_t energy = 4.5;

    TRestAxionWolterOptics opt("xmm.rml");

    Double_t angleMax = maxDev * units("deg");
    Int_t divs = 50;
    Double_t angleStep = angleMax / divs;  // degrees

    std::vector<Int_t> totalCounts;            // Counts that enter inside the outer optics radial limit
    std::vector<Double_t> entranceEfficiency;  // Counts that enter inside the optics
    std::vector<Double_t> exitEfficiency;      // Counts that get outside the optics
    std::vector<Double_t> finalEfficiency;     // Counts weighted with the efficiency

    for (int n = 0; n < divs; n++) {
        totalCounts.push_back(0);
        entranceEfficiency.push_back(0);
        exitEfficiency.push_back(0);
        finalEfficiency.push_back(0);
    }

    Int_t photons = 0;
    while (photons < nPhotons) {
        photons++;

        Double_t eff = opt.PropagateMonteCarloPhoton(energy, maxDev);

        /// Check if it is in the area defined by maximum optical radius
        TVector3 inPos = opt.GetEntrancePosition();
        Double_t Rmax = opt.GetRadialLimits().second;
        if (inPos.X() * inPos.X() + inPos.Y() * inPos.Y() > Rmax * Rmax) continue;

        Double_t entryAngle = opt.GetEntranceAngle() * units("deg");
        Int_t angleIndex = (Int_t)(divs * entryAngle / angleMax);
        if (angleIndex >= divs) angleIndex = divs - 1;

        totalCounts[angleIndex]++;

        TVector3 middle = opt.GetMiddlePosition();
        Double_t rM = middle.X() * middle.X() + middle.Y() * middle.Y();

        /// If it hits the entrance mask, the middle position will not be updated.
        /// But if it goes through the Middle position will be updated and will be != 0
        if (rM > 0) {
            entranceEfficiency[angleIndex] += 1.0;

            if (eff > 0) {
                exitEfficiency[angleIndex] += 1.0;
                finalEfficiency[angleIndex] += eff;
            }
        }
    }

    std::vector<Double_t> angle, entranceError, exitError, finalError;
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
    }

    TCanvas c("", "", 1200, 800);

    TGraphErrors entranceGraph(divs, &angle[0], &entranceEfficiency[0], nullptr, &entranceError[0]);
    entranceGraph.SetName("entranceGraph");
    entranceGraph.SetTitle("Geometrical acceptance and reflectivity");
    entranceGraph.GetXaxis()->SetTitle("Angle to axis [degrees]");
    entranceGraph.GetYaxis()->SetTitle("Efficiency");
    entranceGraph.SetMarkerColor(kRed);
    entranceGraph.SetMarkerStyle(21);
    entranceGraph.Draw("ALP");

    entranceGraph.GetXaxis()->SetLimits(0, angleMax);
    entranceGraph.GetHistogram()->SetMaximum(1);
    entranceGraph.GetHistogram()->SetMinimum(0);

    TGraphErrors exitGraph(divs, &angle[0], &exitEfficiency[0], nullptr, &exitError[0]);
    exitGraph.SetName("exitGraph");
    exitGraph.SetMarkerColor(kBlue);
    exitGraph.SetMarkerStyle(20);
    exitGraph.Draw("LP");

    TGraphErrors finalGraph(divs, &angle[0], &finalEfficiency[0], nullptr, &finalError[0]);
    finalGraph.SetName("finalGraph");
    finalGraph.SetMarkerColor(kGreen);
    finalGraph.SetMarkerStyle(20);
    finalGraph.Draw("LP");

    auto legend = new TLegend(0.1, 0.7, 0.48, 0.9);
    legend->AddEntry("entranceGraph", "Photon enters optics", "lep");
    legend->AddEntry("exitGraph", "Photon exits optics", "lep");
    legend->AddEntry("finalGraph", "Includes reflectivity", "lep");
    legend->Draw();

    c.Print("/tmp/XMM_Efficiency.png");

    return 0;
}
#endif
