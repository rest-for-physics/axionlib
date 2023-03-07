#!/usr/bin/python3

import ROOT

ROOT.gSystem.Load("libRestFramework.so")
ROOT.gSystem.Load("libRestAxion.so")

mcplOptics_1 = ROOT.TRestAxionMCPLOptics("optics.rml", "mcpl")
mcplOptics_1.PrintMetadata()

## MCPL class is not implemented. We just check that it is properly prototyped

# rings  = mcplOptics_1.GetNumberOfRings()
# print( "Number of rings: " + str(rings) )
# if( rings != 5 ):
#    print( "\nError! Number of rings is not 5!" )
#    exit(1)
# print (" [\033[92m OK \x1b[0m]")
#

print("All tests passed!")

exit(0)
