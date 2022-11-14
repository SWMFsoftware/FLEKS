import os, time, sys

def amrex2tec(datPath):
    tStart = time.time()

    fleksDir = __file__.replace("/tools/amrex2tec.py", "")
    binPath = fleksDir+"/bin/converter.exe"
    if not os.path.exists(binPath):
        os.system("cd "+fleksDir+"; make CONVERTER > /dev/null")

    tecPath = datPath+".dat"

    if not os.path.exists(tecPath) and os.system(binPath + " " + datPath + " > /dev/null") == -1:
        print("amrex2tec conversion failed for file ", datPath)

    tEnd = time.time()
    print("amrex2tec conversion takes {0} s".format(tEnd-tStart))

    return tecPath

if __name__ == "__main__":
    files = sys.argv[1:]
    for f in files:
        print("\nConverting ", f)
        amrex2tec(f)
