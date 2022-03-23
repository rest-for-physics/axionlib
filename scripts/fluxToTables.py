#!/usr/bin/python

# This basic script will serve to create the .dat and .spt files required by TRestAxionSolarFlux.
# For the moment it will just integrate the overall flux into 100eV steps and 0.01 radii to generate
# the .dat files.
#
# TODO upgrade the script so that it extracts peaks from a high definition flux file.
# TODO identify the bin size and adjust accordingly.

import sys, os
import ROOT

from ROOT import TH1D, TH2D
import pandas as pd 


import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--file', dest='fname', type=str, help='Input filename .flux')
parser.add_argument('--out', dest='outfname', type=str, help='Output filename .dat')
parser.add_argument('--binsize', dest='binsize', type=float, help='The energy bin size used at the .flux file [in eV].')
parser.add_argument('--skiprows', dest='skiprows', type=int, help='The number of header rows to be skipped')
args = parser.parse_args()

# This is the default for SolarAxionFlux library
skiprows = 3
if args.skiprows != None:
    skiprows = args.skiprows

# This is the default for Primakoff_LennertHoof.flux
binsize = 10
if args.skiprows != None:
    binsize = args.binsize

if args.fname == None:
    parser.print_usage()
    sys.exit(1)

print (args.fname)
print (args.fname.find("."))

outfname = os.path.basename(args.fname)
outfname = outfname[0:outfname.find(".")] + ".dat"
if args.outfname != None:
    outfname = args.outfname + ".dat"
print (outfname)


data = pd.read_csv(args.fname, sep=" ", skiprows=skiprows, header=None)

print (data)

fluxTable = TH2D( "flux", "Flux table", 100, 0., 1., 200, 0., 20. )

#fluxData = data.as_matrix(columns=data.columns[0:])

for n in range(len(data)):
    r = 0.005 + data[0][n]
    en = data[1][n] - 0.005
    flux = data[2][n] * binsize/100.

    fluxTable.Fill( r, en, flux )

    #print( "R: " + str(r) + " En: " + str(en) + " Flux: " + str(flux) + "\n" )

print (outfname)
outf = open( outfname, "w")
print( str(fluxTable.GetNbinsX()) + " " + str(fluxTable.GetNbinsY()) )
for n in range(0,fluxTable.GetNbinsX()):
    for m in range(0,fluxTable.GetNbinsY()):
        outf.write( str(fluxTable.GetBinContent( n+1, m+1 )) )
        if m == fluxTable.GetNbinsY() - 1:
            outf.write("\n")
        else:
            outf.write("\t")
outf.close()
