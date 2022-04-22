#!/usr/bin/python3

import ROOT, math

outfname = "windowsTransmission.png"

from ROOT import (
     TChain, TFile, TTree, TCanvas, TPad, TRandom3,
     TH1D, TH2D, TH3D,
     TProfile, TProfile2D, TProfile3D,
     TGraph, TGraph2D,
     TF1, TF2, TF3, TFormula,
     TLorentzVector, TVector3)

ROOT.gSystem.Load("libRestFramework.so")
ROOT.gSystem.Load("libRestAxion.so")

### Creating a canvas and pad for drawing
c1 = TCanvas( 'c1', 'My canvas', 1200, 400 )
c1.GetFrame().SetBorderSize( 6 )
c1.GetFrame().SetBorderMode( -1 )

pad1 = TPad( "pad1", "This is pad1", 0.01, 0.02, 0.99, 0.97 );
pad1.Divide( 3, 1 )
pad1.Draw()

totalSamples = 100000

cathode = ROOT.TRestAxionXrayWindow("windows.rml", "cathode")
strongBack = ROOT.TRestAxionXrayWindow("windows.rml", "strongBack")
siFoil = ROOT.TRestAxionXrayWindow("windows.rml", "siliconFoil")

radius = cathode.GetWindowRadius()
print ( "\nGetting window radius" )
if radius != 8:
    print ("\nThe window radius is not as expected! Exit code : 102")
    exit(102)
print ( "[\033[92m OK \x1b[0m]" )

# L: Low energy (0-5)keV M: Medium energy H: High energy
histL = ROOT.TH2D("Low", "Low energy", 80, -radius-2, radius+2, 80, -radius-2, radius+2 )
histM = ROOT.TH2D("Medium", "Medium energy", 80, -radius-2, radius+2, 80, -radius-2, radius+2 )
histH = ROOT.TH2D("High", "High energy", 80, -radius-2, radius+2, 80, -radius-2, radius+2 )

#for n in range(rings):
#    gr = ROOT.TGraph()
#    gr.SetMarkerStyle(20)
#    gr.SetMarkerSize(0.5)
#    gr.SetMarkerColor(38+n)
#    gr.SetTitle("Optics plane (no-Spider)")
#    gr.GetXaxis().SetTitle("X [mm]")
#    gr.GetXaxis().SetTitleSize(0.05);
#    gr.GetXaxis().SetLabelSize(0.05);
#    gr.GetYaxis().SetTitle("Y [mm]")
#    gr.GetYaxis().SetTitleOffset(1)
#    gr.GetYaxis().SetTitleSize(0.05);
#    gr.GetYaxis().SetLabelSize(0.05);
#    graphsOP.append(gr)


rnd = TRandom3(0)
for n in range(totalSamples):
     x = (radius+2) * (rnd.Rndm() - 0.5) * 2
     y = (radius+2) * (rnd.Rndm() - 0.5) * 2
     en = 15. * rnd.Rndm()
 
     cth = cathode.GetTransmission( en, x, y)
     stb = strongBack.GetTransmission( en, x, y)
     sil = siFoil.GetTransmission( en, x, y)

     if( cth == 0 ):
         continue
     #print ("x: " + str(x) + " y: " + str(y) + " en: " + str(en) )
     #print ("cathode: " + str(cth) + " sBack: " + str(stb) + " siFoil: " + str(sil) )
     tr = cth * stb * sil

     if( en < 5 ):
         histL.Fill( x, y, tr )
     elif( en < 10):
         histM.Fill( x, y, tr )
     else:
         histH.Fill( x, y, tr )

## This is to work out a graph with the MC obtained efficiencies (TODO)
grCathode = ROOT.TGraph()
grCathode.SetMarkerStyle(20)
grCathode.SetMarkerSize(1.5)
grCathode.SetMarkerColor(38)
grCathode.SetTitle("Aluminum cathode 20nm")
grCathode.GetXaxis().SetTitle("Energy [keV]")
grCathode.GetXaxis().SetTitleSize(0.05);
grCathode.GetXaxis().SetLabelSize(0.05);
grCathode.GetYaxis().SetTitle(" ")
grCathode.GetYaxis().SetTitleOffset(1)
grCathode.GetYaxis().SetTitleSize(0.05);
grCathode.GetYaxis().SetLabelSize(0.05);

grSBack = ROOT.TGraph()
grSBack.SetMarkerStyle(20)
grSBack.SetMarkerSize(1.5)
grSBack.SetMarkerColor(39)
grSBack.SetTitle("Silicon strong-back 200um")
grSBack.GetXaxis().SetTitle("Energy [keV]")
grSBack.GetXaxis().SetTitleSize(0.05);
grSBack.GetXaxis().SetLabelSize(0.05);
grSBack.GetYaxis().SetTitle(" ")
grSBack.GetYaxis().SetTitleOffset(1)
grSBack.GetYaxis().SetTitleSize(0.05);
grSBack.GetYaxis().SetLabelSize(0.05);

grSiFoil = ROOT.TGraph()
grSiFoil.SetMarkerStyle(20)
grSiFoil.SetMarkerSize(1.5)
grSiFoil.SetMarkerColor(40)
grSiFoil.SetTitle("Silicon foil 500nm")
grSiFoil.GetXaxis().SetTitle("Energy [keV]")
grSiFoil.GetXaxis().SetTitleSize(0.05);
grSiFoil.GetXaxis().SetLabelSize(0.05);
grSiFoil.GetYaxis().SetTitle(" ")
grSiFoil.GetYaxis().SetTitleOffset(1)
grSiFoil.GetYaxis().SetTitleSize(0.05);
grSiFoil.GetYaxis().SetLabelSize(0.05);

pad1.cd(1)
histL.SetStats(0)
histL.Draw("colz")

pad1.cd(2)
histM.SetStats(0)
histM.Draw("colz")

pad1.cd(3)
histH.SetStats(0)
histH.Draw("colz")
#graphsOP[0].GetXaxis().SetLimits(-maxRingRadius-20,maxRingRadius+20);
#graphsOP[0].GetHistogram().SetMaximum(maxRingRadius+20);
#graphsOP[0].GetHistogram().SetMinimum(-maxRingRadius-20);
#graphsOP[0].Draw("AP")
#for n in range(1,rings):
#    graphsOP[n].Draw("P")

c1.Print(outfname)

print ("All tests passed!  [\033[92m OK \x1b[0m]")

print ("")

exit(0)
