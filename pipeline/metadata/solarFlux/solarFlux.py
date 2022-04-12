#!/usr/bin/python3

totalSamples = 20000

outfname = "flux.png"

import math
import ROOT
from ROOT import (
    TChain, TFile, TTree, TCanvas, TPad, TRandom3,
    TH1D, TH2D, TH3D,
    TProfile, TProfile2D, TProfile3D,
    TGraph, TGraph2D,
    TF1, TF2, TF3, TFormula,
    TLorentzVector, TVector3)

ROOT.gSystem.Load("libRestFramework.so")
ROOT.gSystem.Load("libRestAxion.so")

c1 = TCanvas( 'c1', 'My canvas', 800,700 )
c1.GetFrame().SetBorderSize( 6 )
c1.GetFrame().SetBorderMode( -1 )

pad1 = TPad("pad1","This is pad1",0.01,0.02,0.99,0.97);
pad1.Divide(2,2)
pad1.Draw()

combinedFlux = ROOT.TRestAxionSolarFlux("fluxes.rml", "combined")
combinedFlux.LoadTables()
combinedFlux.PrintMetadata()

if combinedFlux.GetError():
    print ( combinedFlux.GetErrorMessage() )
    print ( "\nSolar flux initialization failed! Exit code : 101" )
    exit(101)

comb_spt = TH2D("comb_spt", "Energy versus solar radius", 200, 0, 20, 100, 0, 1 )
for x in range(totalSamples):
    x = combinedFlux.GetRandomEnergyAndRadius()
    comb_spt.Fill( x[0], x[1] )

rnd = TRandom3(0)
solarDisk = TH2D("solar_disk", "SolarDisk", 120, -1.2, 1.2, 120, -1.2, 1.2 )
for x in range(totalSamples):
    angle = rnd.Rndm() * 2 * math.pi
    x = combinedFlux.GetRandomEnergyAndRadius()

    solarDisk.Fill( x[1]*math.cos(angle), x[1]*math.sin(angle) )

pad1.cd(1)
comb_spt.SetStats(0)
comb_spt.GetXaxis().SetTitle("Energy [keV]")
comb_spt.GetXaxis().SetTitleSize(0.05);
comb_spt.GetXaxis().SetLabelSize(0.05);
comb_spt.GetYaxis().SetTitle("Solar radius")
comb_spt.GetYaxis().SetTitleSize(0.05);
comb_spt.GetYaxis().SetLabelSize(0.05);
comb_spt.Draw("box")

pad1.cd(2)
enSpt = comb_spt.ProjectionX()
enSpt.SetTitle("Energy spectrum")
enSpt.GetYaxis().SetTitleSize(0.05);
enSpt.SetStats(0)
enSpt.Draw()
if ( enSpt.GetMaximumBin() != 41 ):
    exit(1)

pad1.cd(3)
rSpt = comb_spt.ProjectionY()
rSpt.SetTitle("Radial distribution")
rSpt.GetYaxis().SetTitleSize(0.05);
rSpt.SetStats(0)
rSpt.Draw()
if ( rSpt.GetMaximumBin() != 25 ):
    exit(2)

pad1.cd(4)
solarDisk.SetStats(0)
solarDisk.GetXaxis().SetTitle("X")
solarDisk.GetXaxis().SetTitleSize(0.05);
solarDisk.GetXaxis().SetLabelSize(0.05);
solarDisk.GetYaxis().SetTitle("Y")
solarDisk.GetYaxis().SetTitleOffset(1)
solarDisk.GetYaxis().SetTitleSize(0.05);
solarDisk.GetYaxis().SetLabelSize(0.05);
solarDisk.Draw()

c1.Print(outfname)

monoFlux = ROOT.TRestAxionSolarFlux("fluxes.rml", "mono")
monoFlux.LoadTables()
monoFlux.PrintMetadata()

if monoFlux.GetError():
    print ( monoFlux.GetErrorMessage() )
    print ( "\nSolar flux initialization failed! Exit code : 101" )
    exit(101)

x = monoFlux.GetRandomEnergyAndRadius()
if ( int( x[0]*100 ) != 800 or int( x[1]*100 ) != 58 ):
    print ( "\nMonochromatic flux values seem to be wrong! Exit code : 201" )
    exit(201)
    
print ("[\033[92m OK \x1b[0m]")

continuumFlux = ROOT.TRestAxionSolarFlux("fluxes.rml", "Gianotti")
continuumFlux.LoadTables()
continuumFlux.PrintMetadata()

if continuumFlux.GetError():
    print ( continuumFlux.GetErrorMessage() )
    print ( "\nSolar flux initialization failed! Exit code : 101" )
    exit(101)

x = monoFlux.GetRandomEnergyAndRadius()
if ( int( x[0]*100 ) != 400 or int( x[1]*100 ) != 24 ):
    print ( "\nMonochromatic flux values seem to be wrong! Exit code : 202" )
    exit(202)

print ("All tests passed!  [\033[92m OK \x1b[0m]")

print ("")

exit(0)
