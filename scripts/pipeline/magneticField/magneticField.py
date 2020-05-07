#!/usr/bin/python

import ROOT

ROOT.gSystem.Load("libGdml.so")
ROOT.gSystem.Load("libRestCore.so")
ROOT.gSystem.Load("libRestEvents.so")
ROOT.gSystem.Load("libRestMetadata.so")
ROOT.gSystem.Load("libRestAxion.so")

myField = ROOT.TRestAxionMagneticField("fields.rml", "bFieldBabyIAXO")

if myField.GetError():
    print myField.GetErrorMessage()
    print "\nMagnetic field initialization failed! Exit code : 101"
    exit(101)

p1 = ROOT.TVector3(0,0,0)
p2 = ROOT.TVector3(0,0,-2500)
p3 = ROOT.TVector3(0,0,2500)
p4 = ROOT.TVector3(0,0,4500)
d1 = ROOT.TVector3(0,0,1)

b1 = int( 1000 * myField.GetTransversalComponent( p1, d1))
b2 = int( 1000 * myField.GetTransversalComponent( p2, d1))
b3 = int( 1000 * myField.GetTransversalComponent( p3, d1))
b4 = int( 1000 * myField.GetTransversalComponent( p4, d1))

print "\nEvaluating field volume Bykovskiy_201906.dat centered at (0,0,0)",
if( b1 != 2007 or b2 != 1998 or b3 != 1998 or b4 != 1119 ):
    print "\nEvaluation of field failed! Exit code : 102"
    exit(102)
print "[\033[92m OK \x1b[0m]"
print ""
print "All tests passed!"

exit(0)
