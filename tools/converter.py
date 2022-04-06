import glob
import paraview.simple as pvs
import sys
import time
import re
import os
from amrex2tec import amrex2tec

# Q: How to fix "no module named paraview" error?
# A: Add paraview python package location to $PYTHONPATH. For example:
#    export PYTHONPATH=$PYTHONPATH:/usr/lib/python3.10/site-packages


def clean_dat(datPath):
    # In the *.dat file, if 'VARIABLES'contain '[' or ']',
    # Paraview will fail. So, remove these characters from *.dat file
    f = open(datPath, "r+")
    # Assume 'VARIABLES' is in the first 10 lines
    for i in range(10):
        p0 = f.tell()
        line = f.readline()
        if(line.find("VARIABLES") != -1):
            p1 = f.tell()
            break
    # Remove the unit.
    lineNew = re.sub(r"\s*?\[(.*?)\]", r"", line)
    # Remove '\n'
    lineNew = re.sub(r"\n", r"", lineNew)
    # Padding space so that lineNew's length is the same as line's
    lineNew += (len(line)-len(lineNew)-1)*" " + "\n"

    f.seek(p0)
    f.write(lineNew)
    f.close()

    return line


def tec2vtk(datPath):
    # The performance of converting a 6G .dat file
    # None: fastest write.        speed: fail    file size: fail
    # LZ4: fastest compressed.    speed: 60s        file size: 1.4G
    # ZLib: balanced              speed: 96s        file size: 1.0G
    # LZMA: smallest file size    speed: 455s       file size: 0.84G
    compressor = "ZLib"

    tStart = time.time()

    line = clean_dat(datPath)

    datFile = datPath.split('/')[-1]
    dat = pvs.TecplotReader(registrationName=datFile, FileNames=datPath)

    vtkPath = datPath.replace('.dat', '_vtk')+"/"

    pvs.SaveData(vtkPath+datFile.split(".")[0]+".vtm", proxy=dat,
                 PointDataArrays=dat.PointData.keys(),
                 CompressorType=compressor, CompressionLevel='5')

    with open(vtkPath+"unit.txt", 'w') as f:
        f.write(line)

    tEnd = time.time()
    print("tec2vtk conversion takes {0} s".format(tEnd-tStart))
    return vtkPath

######################################################################

def amrex2vtk(datPath):
    tStart = time.time()

    tecPath = amrex2tec(datPath)
    vtkPath = tec2vtk(tecPath)
    os.system("rm -rf "+tecPath)

    tEnd = time.time()
    print("amrex2vtk conversion takes {0} s".format(tEnd-tStart))


if __name__ == "__main__":
    files = sys.argv[1:]
    for f in files:
        if f.find(".dat") >= 0:
            print("\n Converting ", f)
            tec2vtk(f)
        else:
            print("\nConverting ", f)
            amrex2vtk(f)
