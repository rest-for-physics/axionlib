#!/usr/bin/python3

import ROOT

ROOT.gSystem.Load("libRestFramework.so")
ROOT.gSystem.Load("libRestAxion.so")

mcplOptics_1 = ROOT.TRestAxionMCPLOptics("setups.rml", "mcpl")
mcplOptics_1.PrintMetadata()

rings  = mcplOptics_1.GetNumberOfRings()
print( "Number of rings: " + str(rings) )
if( rings != 5 ):
    print( "\nError! Number of rings is not 5!" )
    exit(1)
print (" [\033[92m OK \x1b[0m]")

#mcplOptics_2 = ROOT.TRestAxionMagneticField("setups.rml", "mcpl2")
#mcplOptics_2.PrintMetadata()

#mcplOptics_3 = ROOT.TRestAxionMagneticField("setups.rml", "mcpl3")
#mcplOptics_3.PrintMetadata()



print ("All tests passed!")

exit(0)
