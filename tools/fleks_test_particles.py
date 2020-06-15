import struct
import os
import glob
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
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
                                             str(iDomain)+"_particle_list_species_"+str(iSpecies)+"_*")

            self.pfiles = self.pfiles + glob.glob(outputDir+"/FLEKS" +
                                         str(iDomain)+"_particle_species_"+str(iSpecies)+"_*")

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
        return np.array(dataList).reshape(nRecord, unitSize)

    def read_initial_loc_with_ID(self, partID):
        unitSize = 7
        for fileName, plist in zip(self.pfiles, self.plists):
            if partID in plist:
                ploc = plist[partID]
                with open(fileName, 'rb') as f:
                    f.seek(ploc)
                    binaryData = f.read(4*4)
                    (cpu, idtmp, nRecord, weight) = struct.unpack(
                        'iiif', binaryData)
                    nRead = 1
                    binaryData = f.read(4*unitSize*nRead)
                    dataList = list(struct.unpack('f'*nRead*unitSize, binaryData))
        return dataList

    def select_particles(self, fSelect):
        selected={}
        icount = 0        
        for pid in self.pset:
            pdata = self.read_initial_loc_with_ID(pid)
            if(fSelect(pid, pdata)):
                selected.update({pid:pdata})
                icount = icount + 1
        return selected


    def plot(self, data):
        plt.ion()
        t = data[:, FLEKSTP.it_]

        tNorm = (t-t[0])/(t[-1]-t[0])
        
        f = plt.figure(figsize=(12, 6))

        nrow = 2
        ncol = 4
        isub = 1
        ax = f.add_subplot(nrow, ncol, isub)
        ax.plot(data[:, FLEKSTP.ix_], data[:, FLEKSTP.iy_], 'k')
        ax.scatter(data[:, FLEKSTP.ix_], data[:, FLEKSTP.iy_], c=plt.cm.winter(
            tNorm), edgecolor='none', marker='o', s=20)
        ax.set_xlabel('x')
        ax.set_ylabel('y')

        isub = isub + 1
        ax = f.add_subplot(nrow, ncol, isub)
        ax.plot(data[:, FLEKSTP.ix_], data[:, FLEKSTP.iz_], 'k')
        ax.scatter(data[:, FLEKSTP.ix_], data[:, FLEKSTP.iz_], c=plt.cm.winter(
            tNorm), edgecolor='none', marker='o', s=20)
        ax.set_xlabel('x')
        ax.set_ylabel('z')

        isub = isub + 1
        ax = f.add_subplot(nrow, ncol, isub)
        ax.plot(data[:, FLEKSTP.iy_], data[:, FLEKSTP.iz_], 'k')
        ax.scatter(data[:, FLEKSTP.iy_], data[:, FLEKSTP.iz_], c=plt.cm.winter(
            tNorm), edgecolor='none', marker='o', s=20)
        ax.set_xlabel('y')
        ax.set_ylabel('z')

        isub = isub + 1
        ax = f.add_subplot(nrow, ncol, isub, projection='3d')
        ax.plot3D(data[:, FLEKSTP.ix_], data[:, FLEKSTP.iy_], data[:, FLEKSTP.iz_])
        ax.scatter(data[:, FLEKSTP.ix_], data[:, FLEKSTP.iy_], data[:, FLEKSTP.iz_], c=plt.cm.winter(tNorm), marker='o',s=3)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')


        isub = isub + 1
        ax = f.add_subplot(nrow, ncol, isub)
        ax.plot(t, data[:, FLEKSTP.iu_], label='Vx')
        ax.scatter(t,data[:, FLEKSTP.iu_] , c=plt.cm.winter(tNorm),
                         edgecolor='none', marker='o', s=20)
        ax.set_xlabel('time')
        ax.set_ylabel('Vx')

        isub = isub + 1
        ax = f.add_subplot(nrow, ncol, isub)
        ax.plot(t, data[:, FLEKSTP.iv_], label='Vy')
        ax.scatter(t, data[:, FLEKSTP.iv_], c=plt.cm.winter(tNorm),
                         edgecolor='none', marker='o', s=20)
        ax.set_xlabel('time')
        ax.set_ylabel('Vy')

        isub = isub + 1
        ax = f.add_subplot(nrow, ncol, isub)
        ax.plot(t, data[:, FLEKSTP.iw_], label='Vz')
        ax.scatter(t, data[:, FLEKSTP.iw_], c=plt.cm.winter(tNorm),
                         edgecolor='none', marker='o', s=20)
        ax.set_xlabel('time')
        ax.set_ylabel('Vz')


        v = np.sqrt(data[:, FLEKSTP.iu_]**2 +
                    data[:, FLEKSTP.iv_]**2 + data[:, FLEKSTP.iw_]**2)
        isub = isub + 1
        ax = f.add_subplot(nrow, ncol, isub)
        ax.plot(data[:, FLEKSTP.it_], v, label='velocity')
        ax.scatter(t, v, c=plt.cm.winter(tNorm),
                         edgecolor='none', marker='o', s=20)
        ax.set_xlabel('time')
        ax.set_ylabel('|V|')

        return f

    def plot_particle(self, partID):
        pData = self.read_particle_trajectory(partID)
        return self.plot(pData)

    


# outputDir = ["/home/yuxichen/test_particles"]
# tp = FLEKSTP(outputDir)
# # d1 = tp.read_particle_trajectory((0, 852))
# f = tp.plot_particle((2141,45))
