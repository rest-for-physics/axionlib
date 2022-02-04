#!/usr/bin/python3

import ROOT

ROOT.gSystem.Load("libRestFramework.so")
ROOT.gSystem.Load("libRestAxion.so")

mcplOptics_1 = ROOT.TRestAxionMagneticField("setups.rml", "mcpl")
mcplOptics_1.PrintMetadata()

mcplOptics_2 = ROOT.TRestAxionMagneticField("setups.rml", "mcpl2")
mcplOptics_2.PrintMetadata()

mcplOptics_3 = ROOT.TRestAxionMagneticField("setups.rml", "mcpl3")
mcplOptics_3.PrintMetadata()


print ( "Setups where loaded properly" )
print ("[\033[92m OK \x1b[0m]")

print ("All tests passed!")

exit(0)
