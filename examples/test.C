void test() {
    TFile f("Run_babyIAXO_1Volume_XY.root");
    // TFile f("Run_babyIAXO_1Volume_sphereIn.root");
    // TFile f("Run_babyIAXO_1Volume_sphereOut.root");
    TRestAnalysisTree* anaTree;
    f.GetObject("AnalysisTree", anaTree);
    // anaTree->Draw("axGen2_posY:axGen2_posZ");
    // anaTree->Scan("axGen2_posX");
    // anaTree->Scan("axGen2_probability");
    anaTree->Draw("axGen2_probability:axGen2_posY:axGen2_posX", "", "", 10000,
                  0);  // probability for the plan
    // anaTree->Draw("axGen2_probability>0:axGen2_posY:axGen2_posX","","",10000,0); // To get the shape of the
    // geometry
    // anaTree->Draw("axGen2_probability:TMath::ATan(axGen2_posZ/axGen2_posY):TMath::ACos(axGen2_posX/15000)","","",10000,0);
    // // probability for the sphere
}
