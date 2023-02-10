#!/usr/bin/python3

### This script is now obsolete. It should be updated to use a specific optics inherited class.

import ROOT, math

outfname = "opticsBasic.png"

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

### Creating a canvas and pad for drawing
c1 = TCanvas("c1", "My canvas", 800, 700)
c1.GetFrame().SetBorderSize(6)
c1.GetFrame().SetBorderMode(-1)

pad1 = TPad("pad1", "This is pad1", 0.01, 0.02, 0.99, 0.97)
pad1.Divide(2, 2)
pad1.Draw()

totalSamples = 10000
genSize = 80

basicOptics = ROOT.TRestAxionOptics("optics.rml", "basic")
spiderOptics = ROOT.TRestAxionOptics("optics.rml", "basic_spider")

rings = basicOptics.GetNumberOfRings()
print("Number of rings (no-spider): " + str(rings))

maxRingRadius = basicOptics.GetMaxRingRadius()
print("Max ring radius (no-spider): " + str(maxRingRadius))

rings_sp = spiderOptics.GetNumberOfRings()
print("Number of rings (spider): " + str(rings_sp))

maxRingRadius_sp = spiderOptics.GetMaxRingRadius()
print("Max ring radius (spider): " + str(maxRingRadius_sp))

# OP: Optics plane GP: Generator plane OPS: Optics plane spider
graphsOP = []
for n in range(rings):
    gr = ROOT.TGraph()
    gr.SetMarkerStyle(20)
    gr.SetMarkerSize(0.5)
    gr.SetMarkerColor(38 + n)
    gr.SetTitle("Optics plane (no-Spider)")
    gr.GetXaxis().SetTitle("X [mm]")
    gr.GetXaxis().SetTitleSize(0.05)
    gr.GetXaxis().SetLabelSize(0.05)
    gr.GetYaxis().SetTitle("Y [mm]")
    gr.GetYaxis().SetTitleOffset(1)
    gr.GetYaxis().SetTitleSize(0.05)
    gr.GetYaxis().SetLabelSize(0.05)
    graphsOP.append(gr)

graphsGP = []
for n in range(rings):
    gr = ROOT.TGraph()
    gr.SetMarkerStyle(20)
    gr.SetMarkerSize(0.5)
    gr.SetMarkerColor(38 + n)
    gr.SetTitle("Generator plane (no-Spider)")
    gr.GetXaxis().SetTitle("X [mm]")
    gr.GetXaxis().SetTitleSize(0.05)
    gr.GetXaxis().SetLabelSize(0.05)
    gr.GetYaxis().SetTitle("Y [mm]")
    gr.GetYaxis().SetTitleOffset(1)
    gr.GetYaxis().SetTitleSize(0.05)
    gr.GetYaxis().SetLabelSize(0.05)
    graphsGP.append(gr)

graphsOPS = []
for n in range(rings):
    gr = ROOT.TGraph()
    gr.SetMarkerStyle(20)
    gr.SetMarkerSize(0.5)
    gr.SetMarkerColor(38 + n)
    gr.SetTitle("Optics plane (Spider)")
    gr.GetXaxis().SetTitle("X [mm]")
    gr.GetXaxis().SetTitleSize(0.05)
    gr.GetXaxis().SetLabelSize(0.05)
    gr.GetYaxis().SetTitle("Y [mm]")
    gr.GetYaxis().SetTitleOffset(1)
    gr.GetYaxis().SetTitleSize(0.05)
    gr.GetYaxis().SetLabelSize(0.05)
    graphsOPS.append(gr)

graphsGPS = []
for n in range(rings):
    gr = ROOT.TGraph()
    gr.SetMarkerStyle(20)
    gr.SetMarkerSize(0.5)
    gr.SetMarkerColor(38 + n)
    gr.SetTitle("Generator plane (Spider)")
    gr.GetXaxis().SetTitle("X [mm]")
    gr.GetXaxis().SetTitleSize(0.05)
    gr.GetXaxis().SetLabelSize(0.05)
    gr.GetYaxis().SetTitle("Y [mm]")
    gr.GetYaxis().SetTitleOffset(1)
    gr.GetYaxis().SetTitleSize(0.05)
    gr.GetYaxis().SetLabelSize(0.05)
    graphsGPS.append(gr)

rnd = TRandom3(0)
for n in range(totalSamples):
    x = (10 + maxRingRadius) * (rnd.Rndm() - 0.5) * 2
    y = (10 + maxRingRadius) * (rnd.Rndm() - 0.5) * 2

    pos = TVector3(x, y, 0)
    direction = TVector3((rnd.Rndm() - 0.5) * 0.01, (rnd.Rndm() - 0.5) * 0.01, 1)

    posEntrance = basicOptics.GetPositionAtEntrance(pos, direction)

    ring = basicOptics.GetEntranceRing(pos, direction)
    if ring >= 0:
        graphsGP[ring].SetPoint(graphsGP[ring].GetN(), x, y)
        graphsOP[ring].SetPoint(graphsOP[ring].GetN(), posEntrance.X(), posEntrance.Y())

pad1.cd(2)
graphsOP[0].GetXaxis().SetLimits(-maxRingRadius - 20, maxRingRadius + 20)
graphsOP[0].GetHistogram().SetMaximum(maxRingRadius + 20)
graphsOP[0].GetHistogram().SetMinimum(-maxRingRadius - 20)
graphsOP[0].Draw("AP")
for n in range(1, rings):
    graphsOP[n].Draw("P")

pad1.cd(1)
graphsGP[0].GetXaxis().SetLimits(-maxRingRadius - 20, maxRingRadius + 20)
graphsGP[0].GetHistogram().SetMaximum(maxRingRadius + 20)
graphsGP[0].GetHistogram().SetMinimum(-maxRingRadius - 20)
graphsGP[0].Draw("AP")
for n in range(1, rings):
    graphsGP[n].Draw("P")

for n in range(totalSamples):
    x = (10 + maxRingRadius) * (rnd.Rndm() - 0.5) * 2
    y = (10 + maxRingRadius) * (rnd.Rndm() - 0.5) * 2

    pos = TVector3(x, y, 0)
    direction = TVector3((rnd.Rndm() - 0.5) * 0.01, (rnd.Rndm() - 0.5) * 0.01, 1)

    posEntrance = spiderOptics.GetPositionAtEntrance(pos, direction)

    ring = spiderOptics.GetEntranceRing(pos, direction)
    if ring >= 0:
        graphsGPS[ring].SetPoint(graphsGPS[ring].GetN(), x, y)
        graphsOPS[ring].SetPoint(
            graphsOPS[ring].GetN(), posEntrance.X(), posEntrance.Y()
        )

pad1.cd(4)
graphsOPS[0].GetXaxis().SetLimits(-maxRingRadius_sp - 20, maxRingRadius_sp + 20)
graphsOPS[0].GetHistogram().SetMaximum(maxRingRadius_sp + 20)
graphsOPS[0].GetHistogram().SetMinimum(-maxRingRadius_sp - 20)
graphsOPS[0].Draw("AP")
for n in range(1, rings):
    graphsOPS[n].Draw("P")

pad1.cd(3)
graphsGPS[0].GetXaxis().SetLimits(-maxRingRadius_sp - 20, maxRingRadius_sp + 20)
graphsGPS[0].GetHistogram().SetMaximum(maxRingRadius_sp + 20)
graphsGPS[0].GetHistogram().SetMinimum(-maxRingRadius_sp - 20)
graphsGPS[0].Draw("AP")
for n in range(1, rings):
    graphsGPS[n].Draw("P")

c1.Print(outfname)

print("All tests passed!  [\033[92m OK \x1b[0m]")

print("")

exit(0)
