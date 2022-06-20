#include "TCanvas.h"
#include "TLatex.h"
#include "TRestAxionWolterOptics.h"
#include "TRestTask.h"

#ifndef RestTask_Axion_XMMPlots
#define RestTask_Axion_XMMPlots

//*******************************************************************************************************
//***
//*** Your HELP is needed to verify, validate and document this macro
//*** This macro might need update/revision.
//***
//*******************************************************************************************************
Int_t REST_Axion_XMMPlots(Double_t z = 7000, Double_t eMax = 0, Double_t deviation = 0.0005) {
    Int_t nPhotons = 50000;

    TCanvas c("", "", 1200, 1200);

    TRestAxionWolterOptics opt("xmm.rml");
    opt.DrawScatterMaps(z, eMax, deviation, nPhotons);

    /* Not working
c.cd();
TLatex* texxt = new TLatex(0.11, 0.45, title.c_str());
texxt->SetTextColor(1);
texxt->SetTextSize(5);
texxt->Draw("same");
    */

    std::string fname =
        "/tmp/XMM_SMaps_dev_" + DoubleToString(deviation, "%6.4lf") + "_z_" + DoubleToString(z) + ".png";
    c.Print(fname.c_str());

    TCanvas c2("", "", 1200, 1200);
    opt.DrawDensityMaps(z, eMax, deviation, nPhotons);

    fname = "/tmp/XMM_DMaps_dev_" + DoubleToString(deviation, "%6.4lf") + "_z_" + DoubleToString(z) + ".png";
    c2.Print(fname.c_str());

    TCanvas c1("", "", 1000, 1600);
    opt.DrawParticleTracks(deviation, 30);
    fname = "/tmp/XMM_Tracks_dev_" + DoubleToString(deviation, "%6.4lf") + ".png";
    c1.Print(fname.c_str());

    return 0;
}
#endif
