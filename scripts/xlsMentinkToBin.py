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
# Since this is a condition for TRestAxionMagneticField.
xCenter = 0
yCenter = 0.8
zCenter = 0

print("Starting to read")
# loading the data into a matrix (xyzBdata)
file = r"../data/magneticField/Mentink_202401.txt"
df = pd.read_excel(file)

print(df[1:5])

print("Translating to matrix")
# xyzBdata = df.values(columns=df.columns[0:])
xyzBdata = df[df.columns[0:]].values

print(xyzBdata[0][0:6])

print(len(xyzBdata))

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
    f.write(
        str(x[0])
        + "\t"
        + str(x[1])
        + "\t"
        + str(x[2])
        + "\t"
        + str(x[3])
        + "\t"
        + str(x[4])
        + "\t"
        + str(x[5])
        + "\n"
    )
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

print("X max : " + str(xMax))
print("Y max : " + str(yMax))
print("Z max : " + str(zMax) + "\n")

print("X min : " + str(xMin))
print("Y min : " + str(yMin))
print("Z min : " + str(zMin) + "\n")

print("BX max : " + str(xBMax))
print("BY max : " + str(yBMax))
print("BZ max : " + str(zBMax) + "\n")

print("BX min : " + str(xBMin))
print("BY min : " + str(yBMin))
print("BZ min : " + str(zBMin) + "\n")

# Creating the binary file
fbin = open("output.bin", "wb")
count = 0
for x in xyzBdata:
    # We recenter the volume and redefine axis (x becomes z, y becomes x and z becomes y)
    # XLS file distances are expressed in m. We translate to mm.
    y = [
        1000 * (x[0] - xCenter),
        1000 * (x[1] - yCenter),
        1000 * (x[2] - zCenter),
        x[3],
        x[4],
        x[5],
    ]

    # Baby-IAXO is symmetric respect to z (direction along axis) and x (direction parallel to the ground).
    # I.e. x and y in the previous axis definition.
    # The y-axis symmetry (in the new axis definition) affects the second bore that is not described in this map.
    # We apply the corresponding symmetries to define the full map of one magnetic bore in the output binary file.

    count = count + 1
    fbin.write(struct.pack("<%df" % len(y), *y))
    if count < 6:
        print(len(y))
        print(x[0:6])
        print(y[0:6])
    # The original file was only missing the x-axis, when we change the sign of x-axis we must change the sign of Bx.
    if y[0] < 0:
        y[0] = -y[0]
        y[3] = -y[3]
        count = count + 1
        fbin.write(struct.pack("<%df" % len(y), *y))
        if count < 6:
            print(len(y))
            print(x[0:6])
            print(y[0:6])
fbin.close()
print("Lines writen : " + str(count))
