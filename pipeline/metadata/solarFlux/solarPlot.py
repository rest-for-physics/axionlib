#!/usr/bin/python3

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

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--rml", dest="rmlfile", type=str, help="The input RML file .rml")
parser.add_argument(
    "--out",
    dest="outfname",
    type=str,
    help="The output filename in png,pdf,C or any other ROOT accepted format",
)
parser.add_argument(
    "--fluxname",
    dest="fluxname",
    type=str,
    help="The name of the flux definition to be chosen from the RML",
)
parser.add_argument(
    "--N", dest="samples", type=int, help="The number of generated particles"
)
args = parser.parse_args()

rmlfile = "fluxes.rml"
if args.rmlfile != None:
    rmlfile = args.rmlfile

outfname = "flux.png"
if args.outfname != None:
    outfname = args.outfname

fluxname = "combined"
if args.fluxname != None:
    fluxname = args.fluxname


samples = 20000
if args.samples != None:
    samples = args.samples

validation = True
if (
    rmlfile == "fluxes.rml"
    and fluxname == "combined"
    and outfname == "flux.png"
    and samples == 20000
):
    validation = True

ROOT.gSystem.Load("libRestFramework.so")
ROOT.gSystem.Load("libRestAxion.so")

c1 = TCanvas("c1", "My canvas", 800, 700)
c1.GetFrame().SetBorderSize(6)
c1.GetFrame().SetBorderMode(-1)

pad1 = TPad("pad1", "This is pad1", 0.01, 0.02, 0.99, 0.97)
pad1.Divide(2, 2)
pad1.Draw()

combinedFlux = ROOT.TRestAxionSolarQCDFlux(rmlfile, fluxname)
combinedFlux.Initialize()
combinedFlux.PrintMetadata()

if combinedFlux.GetError():
    print(combinedFlux.GetErrorMessage())
    print("\nSolar flux initialization failed! Exit code : 101")
    exit(101)

comb_spt = TH2D("comb_spt", "Energy versus solar radius", 20000, 0, 20, 100, 0, 1)
for x in range(samples):
    x = combinedFlux.GetRandomEnergyAndRadius()
    comb_spt.Fill(x[0], x[1])

rnd = TRandom3(0)
solarDisk = TH2D("solar_disk", "SolarDisk", 120, -1.2, 1.2, 120, -1.2, 1.2)
for x in range(samples):
    angle = rnd.Rndm() * 2 * math.pi
    x = combinedFlux.GetRandomEnergyAndRadius()

    solarDisk.Fill(x[1] * math.cos(angle), x[1] * math.sin(angle))

pad1.cd(1)
pad1.cd(1).SetRightMargin(0.09)
pad1.cd(1).SetLeftMargin(0.15)
pad1.cd(1).SetBottomMargin(0.15)

comb_spt.SetStats(0)
comb_spt.GetXaxis().SetTitle("Energy [keV]")
comb_spt.GetXaxis().SetTitleSize(0.05)
comb_spt.GetXaxis().SetLabelSize(0.05)
comb_spt.GetYaxis().SetTitle("Solar radius")
comb_spt.GetYaxis().SetTitleSize(0.05)
comb_spt.GetYaxis().SetLabelSize(0.05)
comb_spt.Draw("colz")

pad1.cd(2)
pad1.cd(2).SetLogy()
pad1.cd(2).SetRightMargin(0.09)
pad1.cd(2).SetLeftMargin(0.15)
pad1.cd(2).SetBottomMargin(0.15)
enSpt = comb_spt.ProjectionX()
enSpt.SetTitle("Energy spectrum")
enSpt.GetYaxis().SetTitleSize(0.05)
enSpt.SetStats(0)
enSpt.SetFillStyle(4050)
enSpt.SetFillColor(ROOT.kBlue - 9)
enSpt.SetLineColor(ROOT.kBlack)
enSpt.Draw()

pad1.cd(3)
pad1.cd(3).SetRightMargin(0.09)
pad1.cd(3).SetLeftMargin(0.15)
pad1.cd(3).SetBottomMargin(0.15)
rSpt = comb_spt.ProjectionY()
rSpt.SetTitle("Radial distribution")
rSpt.GetYaxis().SetTitleSize(0.05)
rSpt.SetStats(0)
rSpt.SetFillStyle(4050)
rSpt.SetFillColor(ROOT.kBlue - 9)
rSpt.SetLineColor(ROOT.kBlack)
rSpt.Draw()

pad1.cd(4)
pad1.cd(4).SetRightMargin(0.09)
pad1.cd(4).SetLeftMargin(0.15)
pad1.cd(4).SetBottomMargin(0.15)
solarDisk.SetStats(0)
solarDisk.GetXaxis().SetTitle("X")
solarDisk.GetXaxis().SetTitleSize(0.05)
solarDisk.GetXaxis().SetLabelSize(0.05)
solarDisk.GetYaxis().SetTitle("Y")
solarDisk.GetYaxis().SetTitleOffset(1)
solarDisk.GetYaxis().SetTitleSize(0.05)
solarDisk.GetYaxis().SetLabelSize(0.05)
solarDisk.Draw("colz")

c1.Print(outfname)
print("Generated file : " + outfname)

print("\nMaximum energy bin is " + str(enSpt.GetMaximumBin()))
if validation:
    if enSpt.GetMaximumBin() != 8001:
        print("\nMaximum Bin is not the expected one (8001)! Exit code : 1")
        exit(1)

print("\nMaximum radius bin is " + str(rSpt.GetMaximumBin()))

if validation:
    if rSpt.GetMaximumBin() != 25:
        print("\nMaximum Bin is not the expected one (25)! Exit code : 2")
        exit(2)

exit(0)
