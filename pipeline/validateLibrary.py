#!/usr/bin/python
# -*- coding: iso-8859-15 -*-

# This script generates the version header to be used in REST installation
# J.  Galan - Javier.Galan.Lacarra@cern.ch
# 23 - Dec - 2019

#from __future__ import print_function
import os
import sys
import re
import filecmp
import subprocess


def validateClass(className):
    print ""
    print "++++ Validating class : " + className
    with open(className, 'r') as file:
        data = file.read()

        data = data[data.find("::Initialize"):]
        data = getMethodDefinition(data)
        data = removeCppComment(data)

        #print (data)
        #print data.find("SETLIBRARYVERSION(LIBRARY_VERSION);")
        if data.find("SETLIBRARYVERSION(LIBRARY_VERSION);") >= 0:
            print "OK"
            return
        else:
            print( "Problem found at class : " + className )
            print( "SetLibraryVersion was NOT found at Initialization!" )
            sys.exit(1)
    return

def getObservablePositions(data):
    obsposes = {}
    pos = 0
    str = "SETOBSERVABLEVALUE(\""
    while(pos < len(data)):
        pos1 = data.find(str,pos)
        if(pos1 == -1):
           break
        pos1 += len(str)
        pos2 = data.find("\"",pos1)
        if(pos2 == -1):
           break

        name = data[pos1:pos2]
        if(not obsposes.has_key(name)):
            obsposes[name] = pos1

        pos = pos2 + 1
    return obsposes


def getMethodDefinition(text):
    initPos = text.find("{")

    counter = 1
    start = initPos + 1
    while(counter > 0):
        pos1 = text.find("{",start)
        pos2 = text.find("}",start)

        endPosition = pos2 + 1

        if(pos1 != -1 and pos2 != -1):
            if(pos1 < pos2):
                counter = counter + 1
                start = pos1 + 1
            if(pos2 < pos1):
                counter = counter - 1
                start = pos2 + 1
        elif pos1 != -1:
                print "Big error!!"
        else:
                counter = counter - 1
                start = pos2 + 1

    return text[initPos:endPosition].upper()


def removeCppComment(strInput) :
        state = 0;
        strOutput = ''
        strRemoved = ''

        for c in strInput :
            if state == 0 and c == '/' :        # ex. [/]
                state = 1
            elif state == 1 and c == '*' :      # ex. [/*]
                state = 2
            elif state == 1 and c == '/' :      # ex. [#]
                state = 4
            elif state == 1 :                   # ex. [<secure/_stdio.h> or 5/3]
                state = 0

            elif state == 3 and c == '*':       # ex. [/*he**]
                state = 3
            elif state == 2 and c == '*':       # ex. [/*he*]
                state = 3
            elif state == 2:                    # ex. [/*heh]
                state = 2

            elif state == 3 and c == '/':       # ex. [/*heh*/]
                state = 0
            elif state == 3:                    # ex. [/*heh*e]
                state = 2

            elif state == 4 and c == '\\':      # ex. [//hehe\]
                state = 9
            elif state == 9 and c == '\\':      # ex. [//hehe\\\\\]
                state = 9
            elif state == 9:                    # ex. [//hehe\<enter> or //hehe\a]
                state = 4
            elif state == 4 and c == '\n':      # ex. [//hehe<enter>]
                state = 0

            elif state == 0 and c == '\'':      # ex. [']
                state = 5
            elif state == 5 and c == '\\':      # ex. ['\]
                state = 6
            elif state == 6:                    # ex. ['\n or '\' or '\t etc.]
                state = 5
            elif state == 5 and c == '\'':      # ex. ['\n' or '\'' or '\t' ect.]
                state = 0

            elif state == 0 and c == '\"':      # ex. ["]
                state = 7
            elif state == 7 and c == '\\':      # ex. ["\]
                state = 8
            elif state == 8:                    # ex. ["\n or "\" or "\t ect.]
                state = 7
            elif state == 7 and c == '\"':      # ex. ["\n" or "\"" or "\t" ect.]
                state = 0

            if (state == 0 and c != '/') or state == 5 or\
                state == 6 or state == 7 or state == 8 :
                strOutput += c
            else:
                # removed chareters
                strRemoved += c

        return strOutput


files = []

# r=root, d=directories, f = files
for r, d, f in os.walk(sys.argv[1]):
        for file in f:
            validate = 0
            if '.cxx' in file:
#                print ( file )
                with open(os.path.join(r, file)) as fin:
                    if '::InitFromConfigFile' in fin.read():
                        validate = 1
                with open(os.path.join(r, file)) as fin:
                    if '::LoadDefaultConfig' in fin.read():
                        validate = 1
                with open(os.path.join(r, file)) as fin:
                    if '::Initialize' in fin.read():
                        validate = validate + 1
            if validate == 2:
                files.append(os.path.join(r, file))
                validateClass(os.path.join(r, file))

#for f in files:
#    print(f)

#validateProcess(sys.argv[1]);
sys.exit(0)
