#!/usr/bin/python

# This basic script serves to create a plain binary file from the XLS files provided by Bykovsky 
# for the Baby-IAXO magnet. We will only upload the resulting bin file, if the original XLS file
# would be required at some point, do not hesitate to contact me at javier.galan@unizar.es.
#
# I upload this script in case it is necessary to convert future field maps. Some revision will be
# needed in order to adapt to other input/output schemes.
#

import sys
import struct
from pandas import DataFrame, read_csv
import matplotlib.pyplot as plt
import pandas as pd 

# Field map is provided only for the top magnetic bore. We will recenter the field map.
xCenter = 0
yCenter = 0
zCenter = 0.7975

print "Starting to read"
# loading the data into a matrix (xyzBdata)
file = r'../../../data/magneticField/Bykovskiy_201906.xls'
df = pd.read_excel(file)

print  df[1:5]

print "Translating to matrix"
xyzBdata = df.as_matrix(columns=df.columns[0:])

print xyzBdata[0][0:6]

print len(xyzBdata)

# Printing out some info about data
xMax = -100000
xMin = 100000
yMax = -100000
yMin = 100000
zMax = -100000
zMin = 100000

xBMax = -100000
xBMin = 100000
yBMax = -100000
yBMin = 100000
zBMax = -100000
zBMin = 100000

f = open("output.txt", "wt")

for x in xyzBdata:
    f.write( str(x[0]) + "\t" +  str(x[1]) + "\t" +  str(x[2]) + "\t" +  str(x[3]) + "\t" +  str(x[4]) + "\t" +  str(x[5]) + "\n" )
    if xMax < x[0]:
        xMax = x[0]
    if yMax < x[1]:
        yMax = x[1]
    if zMax < x[2]:
        zMax = x[2]

    if xMin > x[0]:
        xMin = x[0]
    if yMin > x[1]:
        yMin = x[1]
    if zMin > x[2]:
        zMin = x[2]

    if xBMax < x[3]:
        xBMax = x[3]
    if yBMax < x[4]:
        yBMax = x[4]
    if zBMax < x[5]:
        zBMax = x[5]

    if xBMin > x[3]:
        xBMin = x[3]
    if yBMin > x[4]:
        yBMin = x[4]
    if zBMin > x[5]:
        zBMin = x[5]

f.close()

print "X max : " + str(xMax)
print "Y max : " + str(yMax)
print "Z max : " + str(zMax) + "\n"

print "X min : " + str(xMin)
print "Y min : " + str(yMin)
print "Z min : " + str(zMin) + "\n"

print "BX max : " + str(xBMax)
print "BY max : " + str(yBMax)
print "BZ max : " + str(zBMax) + "\n"

print "BX min : " + str(xBMin)
print "BY min : " + str(yBMin)
print "BZ min : " + str(zBMin) + "\n"

# Creating the binary file
fbin = open('output.bin', 'wb')
count = 0
for x in xyzBdata:
    # We recenter the volume and redefine axis (x becomes z, y becomes x and z becomes y)
    y = [x[1]-yCenter,x[2]-zCenter,x[0]-xCenter,x[4],x[5],x[3]]

    # Baby-IAXO is symmetric respect to z (direction along axis) and x (direction parallel to the ground).
    # I.e. x and y in the previous axis definition.
    # The y-axis symmetry (in the new axis definition) affects the second bore that is not described in this map.
    # Therefore, we register everything in the binary file for y-axis, while for x and z we register 
    # only positive values to apply later on the corresponding symmetry.
    #
    if y[0] >= 0 and y[2] >= 0:
        count = count + 1
        fbin.write(struct.pack('<%df' % len(y), *y))
        if( count < 6 ):
            print len(y)
            print  x[0:6]
            print  y[0:6]
fbin.close()
print "Lines writen : " + str(count)
