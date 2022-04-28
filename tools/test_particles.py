import matplotlib.pyplot as plt
import os
import numpy as np
import glob
import struct


class FLEKSTP(object):
    r"""
    A class that is used to read and plot test particles.    

    Parameters
    -----------
    outputDirs: String
        The path to the test particle dataset. 

    Examples
    ----------
    >>> tp = FLEKSTP("res/run1/PC/test_particles", iSpecies=1)
    >>> pIDs = list(tp.IDs())
    >>> tp.plot_trajectory(pIDs[3])
    >>> tp.save_trajectory_to_csv(pIDs[5])
    >>> pData = tp.read_particles_at_time(6500.8, doSave=True)    
    >>> f = tp.plot_loc(pData)
    """

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

        self.iSpecies = iSpecies
        self.plistfiles = list()
        self.pfiles = list()
        for outputDir in outputDirs:
            self.plistfiles = self.plistfiles+glob.glob(outputDir+"/FLEKS" +
                                                        str(iDomain)+"_particle_list_species_"+str(iSpecies)+"_*")

            self.pfiles = self.pfiles + glob.glob(outputDir+"/FLEKS" +
                                                  str(iDomain)+"_particle_species_"+str(iSpecies)+"_*")

        self.plistfiles.sort()
        self.pfiles.sort()

        self.plists = []
        for fileName in self.plistfiles:
            self.plists.append(self.read_particle_list(fileName))

        self.pset = set()
        for plist in self.plists:
            self.pset.update(plist.keys())

        self.file_time = []
        for filename in self.pfiles:
            record = self._read_the_first_record(filename)
            self.file_time.append(record[FLEKSTP.it_])

        print('Particles of species ', self.iSpecies,
              ' are read from ', outputDirs)
        print('Number of particles: ', len(self.pset))

    def read_particle_list(self, fileName):
        r"""
        Read and return a list of the particle IDs.     
        """
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

    def _read_the_first_record(self, fileName):
        r"""
        Get the first record stored in one file. 
        """
        dataList = list()
        unitSize = 7
        with open(fileName, 'rb') as f:
            while True:
                binaryData = f.read(4*4)

                if not binaryData:
                    break  # EOF

                (cpu, idtmp, nRecord, weight) = struct.unpack(
                    'iiif', binaryData)
                if nRecord > 0:
                    binaryData = f.read(4*unitSize)
                    dataList = dataList + \
                        list(struct.unpack('f'*unitSize, binaryData))
                    return dataList

    def read_particles_at_time(self, time, doSave=False):
        r"""
        Get the information of all the particles at a given time, 
        and save to a csv file with the name "particles_t***.csv" 
        to the current directory if doSave is True. The location of 
        these particles (the csv file) can be visualized with Paraview.

        Parameters
        -----------
        ids: a numpy array of tuples contains the particle IDs.
        pData: a numpy array real number. Contains the particle weight, 
            location and velocity


        Examples
        ----------
        >>> ids, pData = pt.read_particles_at_time(3700, doSave=True)
        """

        nFile = len(self.pfiles)
        for iFile in range(nFile):
            if iFile == 0 and time < self.file_time[iFile]:
                raise Exception(
                    "Error: There is no particle at the given time")

            if iFile == nFile - 1:
                break
            if time > self.file_time[iFile] and time < self.file_time[iFile+1]:
                break

        fileName = self.pfiles[iFile]

        unitSize = 7
        dataList = []
        idList = []
        with open(fileName, 'rb') as f:
            while True:
                binaryData = f.read(4*4)
                if not binaryData:
                    break  # EOF

                (cpu, idtmp, nRecord, weight) = struct.unpack(
                    'iiif', binaryData)
                binaryData = f.read(4*unitSize*nRecord)
                allRecords = list(struct.unpack(
                    'f'*nRecord*unitSize, binaryData))
                for i in range(nRecord):
                    if(allRecords[unitSize*i + FLEKSTP.it_] > time or i == nRecord-1):
                        dataList.append(allRecords[unitSize*i:unitSize*(i+1)])
                        idList.append((cpu, idtmp))
                        break

        npData = np.array(dataList)
        idData = np.array(idList, dtype="i,i")
        if doSave:
            fileName = "particles_t"+str(time)+".csv"
            np.savetxt(fileName, npData, delimiter=",",
                       header="time, x, y, z, ux, uy, uz", comments="")

        return idData, npData

    def IDs(self):
        return self.pset

    def save_trajectory_to_csv(self, partID, fileName=None):
        r""" 
        Save the trajectory of a particle to a csv file.  

        Example
        -----------------
        >>> tp.save_trajectory_to_csv((3,15))
        """
        pData = self.read_particle_trajectory(partID)
        if fileName == None:
            fileName = "trajectory_"+str(partID[0])+"_"+str(partID[1])+".csv"
        np.savetxt(fileName, pData, delimiter=",",
                   header="time, x, y, z, ux, uy, uz", comments="")

    def read_particle_trajectory(self, partID):
        r"""
        Read and return the trajector of a particle. 

        Parameters
        ----------
        partID: particle ID

        Examples
        ----------
        >>> trajectory = tp.read_particle_trajectory((66,888))
        """
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
        r"""
        Read and return the initial location of a test particle
        """

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
                    dataList = list(struct.unpack(
                        'f'*nRead*unitSize, binaryData))
                return dataList

    def select_particles(self, fSelect=None):
        r"""
        Select and return the particles that satisfy the requirement set by the 
        user defined function fSelect. The first argument of fSelect is the particle
        ID, and the second argument is the initial record (one record contains time,
        location, velocity and weight of a particle) of a particle.

        Examples
        ----------
        >>> def fselect(pid, pdata):
        >>>     return pdata[FLEKSTP.it_] < 3601 and pdata[FLEKSTP.ix_] > 12 and abs(pdata[FLEKSTP.iy_])<1 and abs(pdata[FLEKSTP.iz_])<1

        >>> pselected=tp.select_particles(fselect)  
        >>> tp.plot_trajectory(list(pselected.keys())[1])
        """

        selected = {}
        icount = 0

        if fSelect == None:
            def fSelect(id, data): return True

        for pid in self.pset:
            pdata = self.read_initial_loc_with_ID(pid)
            if(fSelect(pid, pdata)):
                selected.update({pid: pdata})
                icount = icount + 1
        return selected

    def plot_trajectory(self, partID):
        r""" 
        Plots the trajectory and velocities of the particle partID. 

        Example
        -----------------
        >>> tp.plot_trajectory((3,15))
        """

        data = self.read_particle_trajectory(partID)

        plt.ion()
        t = data[:, FLEKSTP.it_]

        tNorm = (t-t[0])/(t[-1]-t[0])

        f = plt.figure(figsize=(12, 6))

        nrow = 3
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
        ax.plot3D(data[:, FLEKSTP.ix_],
                  data[:, FLEKSTP.iy_], data[:, FLEKSTP.iz_])
        ax.scatter(data[:, FLEKSTP.ix_], data[:, FLEKSTP.iy_],
                   data[:, FLEKSTP.iz_], c=plt.cm.winter(tNorm), marker='o', s=3)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        isub = isub + 1
        ax = f.add_subplot(nrow, ncol, isub)
        ax.plot(t, data[:, FLEKSTP.ix_], label='x')
        ax.scatter(t, data[:, FLEKSTP.ix_], c=plt.cm.winter(tNorm),
                   edgecolor='none', marker='o', s=20)
        ax.set_xlabel('time')
        ax.set_ylabel('x')

        isub = isub + 1
        ax = f.add_subplot(nrow, ncol, isub)
        ax.plot(t, data[:, FLEKSTP.iy_], label='y')
        ax.scatter(t, data[:, FLEKSTP.iy_], c=plt.cm.winter(tNorm),
                   edgecolor='none', marker='o', s=20)
        ax.set_xlabel('time')
        ax.set_ylabel('y')

        isub = isub + 1
        ax = f.add_subplot(nrow, ncol, isub)
        ax.plot(t, data[:, FLEKSTP.iz_], label='z')
        ax.scatter(t, data[:, FLEKSTP.iz_], c=plt.cm.winter(tNorm),
                   edgecolor='none', marker='o', s=20)
        ax.set_xlabel('time')
        ax.set_ylabel('z')

        isub = isub + 1

        isub = isub + 1
        ax = f.add_subplot(nrow, ncol, isub)
        ax.plot(t, data[:, FLEKSTP.iu_], label='Vx')
        ax.scatter(t, data[:, FLEKSTP.iu_], c=plt.cm.winter(tNorm),
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

        plt.tight_layout()

        return f

    def plot_loc(self, pData):
        r"""
        Plot the location of particles pData

        Examples
        ----------
        >>> pData = tp.read_particles_at_time(3700, doSave=True)
        >>> f = tp.plot_loc(pData)
        """

        px = pData[:, FLEKSTP.ix_]
        py = pData[:, FLEKSTP.iy_]
        pz = pData[:, FLEKSTP.iz_]

        f = plt.figure(figsize=(12, 12))

        nrow = 2
        ncol = 2
        isub = 1
        ax = f.add_subplot(nrow, ncol, isub)
        ax.scatter(px, py, s=1)
        ax.set_xlabel('x')
        ax.set_ylabel('y')

        isub = 2
        ax = f.add_subplot(nrow, ncol, isub)
        ax.scatter(px, pz, s=1)
        ax.set_xlabel('x')
        ax.set_ylabel('z')

        isub = 3
        ax = f.add_subplot(nrow, ncol, isub)
        ax.scatter(py, pz, s=1)
        ax.set_xlabel('y')
        ax.set_ylabel('z')

        isub = 4
        ax = f.add_subplot(nrow, ncol, isub, projection='3d')
        ax.scatter(px, py, pz, s=1)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_ylabel('z')

        return f
