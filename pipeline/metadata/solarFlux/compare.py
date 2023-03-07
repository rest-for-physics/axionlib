#!/usr/bin/python3

totalSamples = 100000

outfname = "flux_combined.png"

import math
import ROOT
from ROOT import (
    TChain,
    TFile,
    TTree,
    TCanvas,
    TPad,
    TRandom3,
    TH1D,
    TH2D,
    TH3D,
    TProfile,
    TProfile2D,
    TProfile3D,
    TGraph,
    TGraph2D,
    TF1,
    TF2,
    TF3,
    TFormula,
    TLorentzVector,
    TVector3,
)

ROOT.gSystem.Load("libRestFramework.so")
ROOT.gSystem.Load("libRestAxion.so")

c1 = TCanvas("c1", "My canvas", 1200, 600)
c1.GetFrame().SetBorderSize(6)
c1.GetFrame().SetBorderMode(-1)

pad1 = TPad("pad1", "This is pad1", 0.01, 0.02, 0.99, 0.97)
pad1.Divide(2, 1)
pad1.Draw()

primakoffLH = ROOT.TRestAxionSolarFlux("fluxes.rml", "LennertHoofPrimakoff")
primakoffG = ROOT.TRestAxionSolarFlux("fluxes.rml", "Gianotti")

if primakoffLH.GetError():
    print(primakoffLH.GetErrorMessage())
    print("\nSolar flux initialization failed! Exit code : 101")
    exit(101)

LH_spt = TH2D("LennertHoof_2D", "Energy versus solar radius", 200, 0, 20, 100, 0, 1)
for x in range(totalSamples):
    x = primakoffLH.GetRandomEnergyAndRadius()
    LH_spt.Fill(x[0], x[1])

G_spt = TH2D("Gianotti_2D", "Energy versus solar radius", 200, 0, 20, 100, 0, 1)
for x in range(totalSamples):
    x = primakoffG.GetRandomEnergyAndRadius()
    G_spt.Fill(x[0], x[1])

pad1.cd(1)
enSpt_G = G_spt.ProjectionX()
enSpt_G.Scale(primakoffG.GetTotalFlux() / 100.0)  # 100eV binsize
enSpt_G.SetTitle("Energy spectrum")
enSpt_G.GetYaxis().SetTitleSize(0.05)
enSpt_G.SetStats(0)
enSpt_G.GetXaxis().SetTitle("Energy [keV]")
enSpt_G.GetXaxis().SetTitleSize(0.05)
enSpt_G.GetXaxis().SetLabelSize(0.05)
enSpt_G.GetYaxis().SetTitle("Flux [cm-2 s-1 keV-1]")
enSpt_G.GetYaxis().SetTitleSize(0.05)
enSpt_G.GetYaxis().SetLabelSize(0.05)
enSpt_G.SetFillStyle(4050)
enSpt_G.SetFillColorAlpha(ROOT.kBlue, 0.1)
enSpt_G.SetLineColor(ROOT.kBlue + 2)
enSpt_G.Draw("hist")

enSpt_LH = LH_spt.ProjectionX()
enSpt_LH.Scale(primakoffLH.GetTotalFlux() / 100.0)  # 100eV binsize
enSpt_LH.SetFillStyle(4050)
enSpt_LH.SetFillColorAlpha(ROOT.kOrange, 0.1)
enSpt_LH.SetLineColor(ROOT.kOrange + 3)
enSpt_LH.Draw("hist SAME")

legend = ROOT.TLegend(0.45, 0.85, 0.88, 0.76)
legend.AddEntry(enSpt_G, "Primakoff Gianotti", "f")
legend.AddEntry(enSpt_LH, "Primakoff Lennert&Hoof", "f")
legend.SetTextSize(0.03)
legend.Draw()

pad1.cd(2)
rSpt_G = G_spt.ProjectionY()
rSpt_G.SetTitle("Radial distribution")
rSpt_G.GetXaxis().SetTitleSize(0.05)
rSpt_G.GetXaxis().SetTitle("Solar radius")
rSpt_G.GetYaxis().SetTitleSize(0.05)
rSpt_G.GetYaxis().SetLabelSize(0.05)
rSpt_G.SetFillStyle(4050)
rSpt_G.SetFillColorAlpha(ROOT.kBlue, 0.1)
rSpt_G.SetLineColor(ROOT.kBlue + 2)
rSpt_G.SetStats(0)
rSpt_G.Draw()

rSpt_LH = LH_spt.ProjectionY()
rSpt_LH.SetFillStyle(4050)
rSpt_LH.SetFillColorAlpha(ROOT.kOrange, 0.1)
rSpt_LH.SetLineColor(ROOT.kOrange + 3)
rSpt_LH.Draw("SAME")

legend.Draw()

c1.Print(outfname)

print("File: " + str(outfname) + " was generated [\033[92m OK \x1b[0m]")

print("")

exit(0)
