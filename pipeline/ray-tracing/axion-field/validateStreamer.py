import os, sys

os.system("restRoot -b -q StreamerOutput.C'(\"AxionPhotonProbability.root\")' | grep Process | grep TRestAxion | grep version | wc -l > output.log 2>&1")

with open('output.log') as f:
    lines = f.readlines()

for line in lines:
    if( line.find("6") == 0):
        print ("The number of processes inside the event data chain is 6. Succeed!")
        sys.exit(0)
    else:
        print ("The number of processes inside the event data chain is NOT 6! Fail!")
        sys.exit(1)

sys.exit(0)
