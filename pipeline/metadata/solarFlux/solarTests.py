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


ROOT.gSystem.Load("libRestFramework.so")
ROOT.gSystem.Load("libRestAxion.so")

monoFlux = ROOT.TRestAxionSolarFlux("fluxes.rml", "mono")
monoFlux.LoadTables()
monoFlux.PrintMetadata()

if monoFlux.GetError():
    print(monoFlux.GetErrorMessage())
    print("\nSolar flux initialization failed! Exit code : 101")
    exit(101)

x = monoFlux.GetRandomEnergyAndRadius()
print(x[0])
print(x[1])
if int(x[0] * 100) != 800 or int(x[1] * 100) != 83:
    print("\nMonochromatic flux values seem to be wrong! Exit code : 201")
    exit(201)

print("[\033[92m OK \x1b[0m]")

continuumFlux = ROOT.TRestAxionSolarFlux("fluxes.rml", "Gianotti")
continuumFlux.LoadTables()
continuumFlux.PrintMetadata()

if continuumFlux.GetError():
    print(continuumFlux.GetErrorMessage())
    print("\nSolar flux initialization failed! Exit code : 101")
    exit(101)

x = monoFlux.GetRandomEnergyAndRadius()
if int(x[0] * 100) != 400 or int(x[1] * 100) != 24:
    print("\nContinuum flux values seem to be wrong! Exit code : 202")
    exit(202)

print("All tests passed!  [\033[92m OK \x1b[0m]")

print("")

exit(0)
