import os
import shutil
import subprocess
import sys
import time


def amrex2tec(datPath, savePlt=False):
    tStart = time.time()

    fleksDir = __file__.replace("/tools/amrex2tec.py", "")
    binPath = fleksDir + "/bin/converter.exe"
    if not os.path.exists(binPath):
        subprocess.run(["make", "-C", fleksDir, "CONVERTER"], stdout=subprocess.DEVNULL)

    tecPath = datPath + ".dat"
    pltPath = datPath + ".plt"

    if not os.path.exists(tecPath) and not os.path.exists(pltPath):
        res = subprocess.run([binPath, datPath], stdout=subprocess.DEVNULL)
        if res.returncode != 0:
            print("amrex2tec conversion failed for file ", datPath, flush=True)

    if savePlt:
        if not os.path.exists(pltPath):
            res = subprocess.run(["preplot", tecPath], stdout=subprocess.DEVNULL)
            if res.returncode != 0:
                print("tecplot preplot failed for file ", tecPath, flush=True)

        if os.path.isfile(tecPath) or os.path.islink(tecPath):
            os.remove(tecPath)
        elif os.path.isdir(tecPath):
            shutil.rmtree(tecPath)

    tEnd = time.time()
    print("amrex2tec conversion takes {0} s".format(tEnd - tStart), flush=True)

    return tecPath


if __name__ == "__main__":
    flag = int(sys.argv[1])
    files = sys.argv[2:]
    for f in files:
        print("\nConverting ", f, flush=True)
        amrex2tec(f, flag != 0)
