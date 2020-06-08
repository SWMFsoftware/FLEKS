import struct
import os
import glob
import numpy as np
import matplotlib.pyplot as plt


class FLEKSTP:
    it_ = 0
    ix_ = 1
    iy_ = 2
    iz_ = 3
    iu_ = 4
    iv_ = 5
    iw_ = 6

    def __init__(self, outputDirs, iDomain=0, iSpecies=0):
        if type(outputDirs) == str:
            outputDirs = [outputDirs]

        self.plistfiles = list()
        self.pfiles = list()
        print(outputDirs)
        for outputDir in outputDirs:
            self.plistfiles = self.plistfiles+glob.glob(outputDir+"/FLEKS" +
                                             str(0)+"_particle_list_species_"+str(0)+"_*")

            self.pfiles = self.pfiles + glob.glob(outputDir+"/FLEKS" +
                                         str(0)+"_particle_species_"+str(0)+"_*")

        self.plistfiles.sort()
        self.pfiles.sort()

        print(self.plistfiles)

        self.plists = []
        for fileName in self.plistfiles:
            self.plists.append(self.read_particle_list(fileName))

        self.pset = set()
        for plist in self.plists:
            self.pset.update(plist.keys())

    def read_particle_list(self, fileName):
        # 2 integers + 1 unsigned long long
        listUnitSize = 2*4+8
        nByte = os.path.getsize(fileName)
        nPart = int(nByte/listUnitSize)
        plist = {}
        with open(fileName, 'rb') as f:
            for _ in range(nPart):
                binaryData = f.read(listUnitSize)
                (cpu, id, loc) = struct.unpack('iiQ', binaryData)
                plist.update({(cpu, id): loc})
        return plist

    def read_particle_trajectory(self, partID):
        dataList = list()
        unitSize = 7
        for fileName, plist in zip(self.pfiles, self.plists):
            if partID in plist:
                ploc = plist[partID]
                with open(fileName, 'rb') as f:
                    f.seek(ploc)
                    binaryData = f.read(4*4)
                    (cpu, idtmp, nRecord, weight) = struct.unpack(
                        'iiif', binaryData)
                    binaryData = f.read(4*unitSize*nRecord)
                    dataList = dataList + \
                        list(struct.unpack('f'*nRecord*unitSize, binaryData))

        nRecord = int(len(dataList)/unitSize)
        return np.sort(np.array(dataList).reshape(nRecord, unitSize), axis=FLEKSTP.it_)

    def plot(self, data):
        plt.ion()
        t = data[:, FLEKSTP.it_]

        tNorm = (t-t[0])/(t[-1]-t[0])

        f, ax = plt.subplots(2, 2, figsize=(12, 8))
        ax[0, 0].plot(data[:, FLEKSTP.ix_], data[:, FLEKSTP.iy_], 'k')
        ax[0, 0].scatter(data[:, FLEKSTP.ix_], data[:, FLEKSTP.iy_], c=plt.cm.winter(
            tNorm), edgecolor='none', marker='o', s=20)
        ax[0, 0].set_xlabel('x')
        ax[0, 0].set_ylabel('y')

        ax[0, 1].plot(data[:, FLEKSTP.ix_], data[:, FLEKSTP.iz_], 'k')
        ax[0, 1].scatter(data[:, FLEKSTP.ix_], data[:, FLEKSTP.iz_], c=plt.cm.winter(
            tNorm), edgecolor='none', marker='o', s=20)
        ax[0, 1].set_xlabel('x')
        ax[0, 1].set_ylabel('z')
        ax[1, 0].plot(data[:, FLEKSTP.iy_], data[:, FLEKSTP.iz_], 'k')
        ax[1, 0].scatter(data[:, FLEKSTP.iy_], data[:, FLEKSTP.iz_], c=plt.cm.winter(
            tNorm), edgecolor='none', marker='o', s=20)
        ax[1, 0].set_xlabel('y')
        ax[1, 0].set_ylabel('z')

        v = np.sqrt(data[:, FLEKSTP.iu_]**2 +
                    data[:, FLEKSTP.iv_]**2 + data[:, FLEKSTP.iw_]**2)
        ax[1, 1].plot(data[:, FLEKSTP.it_], v, label='velocity')
        ax[1, 1].scatter(t, v, c=plt.cm.winter(tNorm),
                         edgecolor='none', marker='o', s=20)
        ax[1, 1].set_xlabel('time')
        ax[1, 1].set_ylabel('velocity')


# outputDir = ["/home/yuxichen/dev/SWMF/run_test/RESULTS/2step/PC/test_particles",
#              "/home/yuxichen/dev/SWMF/run_test/RESULTS/2restart/PC/test_particles"]
# tp = FLEKSTP(outputDir, 0, 1)
# d1 = tp.read_particle_trajectory((0, 852))
# tp.plot(d1)
