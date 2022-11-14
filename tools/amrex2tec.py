import os, time, sys

def amrex2tec(datPath):
    tStart = time.time()

    fleksDir = __file__.replace("/tools/amrex2tec.py", "")
    binPath = fleksDir+"/bin/converter.exe"
    if not os.path.exists(binPath):
        os.system("cd "+fleksDir+"; make CONVERTER > /dev/null")

    tecPath = datPath+".dat"
    pltPath = datPath+".plt"

    if (not os.path.exists(tecPath) and not os.path.exists(pltPath) and
        os.system(binPath + " " + datPath + " > /dev/null") == -1):
        print("amrex2tec conversion failed for file ", datPath)

    if not os.path.exists(pltPath) and os.system("preplot " + tecPath) == -1:
        print("tecplot preplot failed for file ", tecPath)

    os.system("rm " + tecPath)

    tEnd = time.time()
    print("amrex2tec conversion takes {0} s".format(tEnd-tStart))

    return tecPath

if __name__ == "__main__":
    files = sys.argv[1:]
    for f in files:
        print("\nConverting ", f)
        amrex2tec(f)
