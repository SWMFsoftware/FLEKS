from typing import List, Tuple, Dict, Union, Callable

import matplotlib.pyplot as plt
import os
import numpy as np
import glob
import struct


class FLEKSTP(object):
    r"""
    A class that is used to read and plot test particles. Each particle ID consists of
    a CPU index, a particle index on each CPU, and a location index.
    By default, 7 real numbers saved for each step: time + position + velocity.

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
    >>> ids, pData = tp.read_particles_at_time(6500.8, doSave=True)
    >>> f = tp.plot_loc(pData)
    """

    it_ = 0
    ix_ = 1
    iy_ = 2
    iz_ = 3
    iu_ = 4
    iv_ = 5
    iw_ = 6
    iBx_ = 7
    iBy_ = 8
    iBz_ = 9
    iEx_ = 10
    iEy_ = 11
    iEz_ = 12

    def __init__(
        self,
        outputDirs: Union[str, List[str]],
        iDomain: int = 0,
        iSpecies: int = 0,
        iListStart: int = 0,
        iListEnd: int = -1,
        readAllFiles: bool = False,
    ):
        if type(outputDirs) == str:
            outputDirs = [outputDirs]

        header = outputDirs[0] + "/Header"
        if os.path.exists(header):
            with open(header, "r") as f:
                self.nReal = int(f.readline())
        else:
            raise FileNotFoundError(f"Header file not found in {outputDirs[0]}")

        self.iSpecies = iSpecies
        self.plistfiles = list()
        self.pfiles = list()

        for outputDir in outputDirs:
            self.plistfiles = self.plistfiles + glob.glob(
                f"{outputDir}/FLEKS{iDomain}_particle_list_species_{iSpecies}_*"
            )

            self.pfiles = self.pfiles + glob.glob(
                f"{outputDir}/FLEKS{iDomain}_particle_species_{iSpecies}_*"
            )

        self.plistfiles.sort()
        self.pfiles.sort()

        self.list_index_to_time = []
        if readAllFiles:
            for filename in self.pfiles:
                record = self._read_the_first_record(filename)
                if record == None:
                    continue
                self.list_index_to_time.append(record[FLEKSTP.it_])

        if iListEnd == -1:
            iListEnd = len(self.plistfiles)
        self.plistfiles = self.plistfiles[iListStart:iListEnd]
        self.pfiles = self.pfiles[iListStart:iListEnd]

        self.plists: List[Dict[Tuple[int, int], int]] = []
        for fileName in self.plistfiles:
            self.plists.append(self.read_particle_list(fileName))

        self.IDs = set()
        for plist in self.plists:
            self.IDs.update(plist.keys())

        self.file_time = []
        for filename in self.pfiles:
            record = self._read_the_first_record(filename)
            if record == None:
                continue
            self.file_time.append(record[FLEKSTP.it_])

        print(f"Particles of species {self.iSpecies} are read from {outputDirs}")
        print(f"Number of particles: {len(self.IDs)}")

    def getIDs(self):
        return list(sorted(self.IDs))

    def get_index_to_time(self) -> List:
        r"""
        Getter method for accessing list_index_to_time.
        """
        if len(self.list_index_to_time) == 0:
            print("Index to time mapping was not initialized")
        return self.list_index_to_time

    def read_particle_list(self, fileName: str) -> Dict[Tuple[int, int], int]:
        r"""
        Read and return a list of the particle IDs.
        """
        # 2 integers + 1 unsigned long long
        listUnitSize = 2 * 4 + 8
        nByte = os.path.getsize(fileName)
        nPart = int(nByte / listUnitSize)
        plist = {}
        with open(fileName, "rb") as f:
            for _ in range(nPart):
                binaryData = f.read(listUnitSize)
                (cpu, id, loc) = struct.unpack("iiQ", binaryData)
                plist.update({(cpu, id): loc})
        return plist

    def _read_the_first_record(self, fileName: str) -> Union[List[float], None]:
        r"""
        Get the first record stored in one file.
        """
        dataList = list()
        with open(fileName, "rb") as f:
            while True:
                binaryData = f.read(4 * 4)

                if not binaryData:
                    break  # EOF

                (cpu, idtmp, nRecord, weight) = struct.unpack("iiif", binaryData)
                if nRecord > 0:
                    binaryData = f.read(4 * self.nReal)
                    dataList = dataList + list(
                        struct.unpack("f" * self.nReal, binaryData)
                    )
                    return dataList

    def read_particles_at_time(
        self, time: float, doSave: bool = False
    ) -> Tuple[np.ndarray, np.ndarray]:
        r"""
        Get the information of all the particles at a given time, and save to a csv file
        with the name "particles_t***.csv" in the current directory if doSave is True.

        Returns
        -----------
        ids: a numpy array of tuples contains the particle IDs.
        pData: a numpy real array with the particle weight, location and velocity.

        Examples
        ----------
        >>> ids, pData = pt.read_particles_at_time(3700, doSave=True)
        """

        nFile = len(self.pfiles)
        for iFile in range(nFile):
            if iFile == 0 and time < self.file_time[iFile]:
                raise Exception("Error: There is no particle at the given time")

            if iFile == nFile - 1:
                break
            if time >= self.file_time[iFile] and time < self.file_time[iFile + 1]:
                break

        fileName = self.pfiles[iFile]

        dataList = []
        idList = []
        with open(fileName, "rb") as f:
            while True:
                binaryData = f.read(4 * 4)
                if not binaryData:
                    break  # EOF

                (cpu, idtmp, nRecord, weight) = struct.unpack("iiif", binaryData)
                binaryData = f.read(4 * self.nReal * nRecord)
                allRecords = list(struct.unpack("f" * nRecord * self.nReal, binaryData))
                for i in range(nRecord):
                    if (
                        allRecords[self.nReal * i + FLEKSTP.it_] >= time
                        or i == nRecord - 1
                    ):
                        dataList.append(
                            allRecords[self.nReal * i : self.nReal * (i + 1)]
                        )
                        idList.append((cpu, idtmp))
                        break

        npData = np.array(dataList)
        idData = np.array(idList, dtype="i,i")
        if doSave:
            fileName = f"particles_t{time}.csv"
            header = "cpu,iid,time,x,y,z,ux,uy,uz"
            if self.nReal == 10:
                header += ",bx,by,bz"
            elif self.nReal == 13:
                header += ",bx,by,bz,ex,ey,ez"

            with open(fileName, "w") as f:
                f.write(header + "\n")
                for id_row, data_row in zip(idData, npData):
                    f.write(
                        f"{id_row[0]},{id_row[1]},{','.join(str(x) for x in data_row)}\n"
                    )

        return idData, npData

    def save_trajectory_to_csv(
        self,
        pID: Tuple[int, int],
        fileName: str = None,
        shiftTime: bool = False,
        scaleTime: bool = False,
    ) -> None:
        r"""
        Save the trajectory of a particle to a csv file.

        Parameters
        ----------
        pID: particle ID.
        shiftTime: If set to True, set the initial time to be 0.
        scaleTime: If set to True, scale the time into [0,1] range, only scale time if
                    shiftTime = True.

        Example
        -----------------
        >>> tp.save_trajectory_to_csv((3,15))
        """
        pData = self.read_particle_trajectory(pID)
        if fileName == None:
            fileName = "trajectory_" + str(pID[0]) + "_" + str(pID[1]) + ".csv"
        header = "time [s], X [R], Y [R], Z [R], U_x [km/s], U_y [km/s], U_z [km/s]"
        if self.nReal == 10:
            header += ", B_x [nT], B_y [nT], B_z [nT]"
        if self.nReal == 13:
            header += (
                ", B_x [nT], B_y [nT], B_z [nT], E_x [uV/m], E_y [uV/m], E_z [uV/m]"
            )
        if shiftTime:
            pData[:, 0] -= pData[0, 0]
            if scaleTime:
                pData[:, 0] /= pData[-1, 0]
        np.savetxt(fileName, pData, delimiter=",", header=header, comments="")

    def read_particle_trajectory(self, pID: Tuple[int, int]):
        r"""
        Read and return the trajectory of a particle.

        Parameters
        ----------
        pID: particle ID

        Examples
        ----------
        >>> trajectory = tp.read_particle_trajectory((66,888))
        """
        dataList = list()
        for fileName, plist in zip(self.pfiles, self.plists):
            if pID in plist:
                ploc = plist[pID]
                with open(fileName, "rb") as f:
                    f.seek(ploc)
                    binaryData = f.read(4 * 4)
                    (cpu, idtmp, nRecord, weight) = struct.unpack("iiif", binaryData)
                    binaryData = f.read(4 * self.nReal * nRecord)
                    dataList = dataList + list(
                        struct.unpack("f" * nRecord * self.nReal, binaryData)
                    )

        nRecord = int(len(dataList) / self.nReal)
        return np.array(dataList).reshape(nRecord, self.nReal)

    def read_initial_location(self, pID):
        r"""
        Read and return the initial location of a test particle
        """

        for fileName, plist in zip(self.pfiles, self.plists):
            if pID in plist:
                ploc = plist[pID]
                with open(fileName, "rb") as f:
                    f.seek(ploc)
                    binaryData = f.read(4 * 4)
                    (cpu, idtmp, nRecord, weight) = struct.unpack("iiif", binaryData)
                    nRead = 1
                    binaryData = f.read(4 * self.nReal * nRead)
                    dataList = list(struct.unpack("f" * nRead * self.nReal, binaryData))
                return dataList

    def select_particles(
        self, fSelect: Callable[[Tuple[int, int], List[float]], bool] = None
    ) -> List[Tuple[int, int]]:
        r"""
        Select and return the particles whose initial condition satisfy the requirement
        set by the user defined function fSelect. The first argument of fSelect is the
        particle ID, and the second argument is the initial record (time, location,
        velocity and weight of a particle) of a particle.

        Examples
        ----------
        >>> def fselect(tp, pid):
        >>>     pdata = tp.read_initial_loc_with_ID(pid)
        >>>     intime = pdata[FLEKSTP.it_] < 3601
        >>>     inregion = pdata[FLEKSTP.ix_] > 20
        >>>     return intime and inregion
        >>>
        >>> pselected = tp.select_particles(fselect)
        >>> tp.plot_trajectory(list(pselected.keys())[1])
        """

        if fSelect == None:
            def fSelect(id, data):
                return True

        pselected = list(filter(fSelect, self.IDs))

        return pselected

    def plot_trajectory(self, pID: Tuple[int, int]):
        r"""
        Plots the trajectory and velocities of the particle pID.

        Example
        -----------------
        >>> tp.plot_trajectory((3,15))
        """

        data = self.read_particle_trajectory(pID)
        t = data[:, FLEKSTP.it_]
        tNorm = (t - t[0]) / (t[-1] - t[0])

        ncol = 3
        nrow = 3  # Default for X, V
        if self.nReal == 10:  # additional B field
            nrow = 4
        elif self.nReal == 13:  # additional B and E field
            nrow = 5

        f, axs = plt.subplots(nrow, ncol, figsize=(12, 6), constrained_layout=True)

        # Plot trajectories
        for i, ax in enumerate(axs[0, :]):
            x_id = FLEKSTP.ix_ if i < 2 else FLEKSTP.iy_
            y_id = FLEKSTP.iy_ if i == 0 else FLEKSTP.iz_
            ax.plot(data[:, x_id], data[:, y_id], "k")
            ax.scatter(
                data[:, x_id],
                data[:, y_id],
                c=plt.cm.winter(tNorm),
                edgecolor="none",
                marker="o",
                s=10,
            )
            ax.set_xlabel("x" if i < 2 else "y")
            ax.set_ylabel("y" if i == 0 else "z")

        def plot_data(dd, label, irow, icol):
            axs[irow, icol].plot(t, dd, label=label)
            axs[irow, icol].scatter(
                t, dd, c=plt.cm.winter(tNorm), edgecolor="none", marker="o", s=10
            )
            axs[irow, icol].set_xlabel("time")
            axs[irow, icol].set_ylabel(label)

        def plot_vector(idx, labels, irow):
            for i, (id, label) in enumerate(zip(idx, labels)):
                plot_data(data[:, id], label, irow, i)

        plot_vector([FLEKSTP.ix_, FLEKSTP.iy_, FLEKSTP.iz_], ["x", "y", "z"], 1)
        plot_vector([FLEKSTP.iu_, FLEKSTP.iv_, FLEKSTP.iw_], ["Vx", "Vy", "Vz"], 2)

        if self.nReal > FLEKSTP.iBx_:
            plot_vector(
                [FLEKSTP.iBx_, FLEKSTP.iBy_, FLEKSTP.iBz_], ["Bx", "By", "Bz"], 3
            )

        if self.nReal > FLEKSTP.iEx_:
            plot_vector(
                [FLEKSTP.iEx_, FLEKSTP.iEy_, FLEKSTP.iEz_], ["Ex", "Ey", "Ez"], 4
            )

        return

    def plot_loc(self, pData: np.ndarray):
        r"""
        Plot the location of particles pData.

        Examples
        ----------
        >>> ids, pData = tp.read_particles_at_time(3700, doSave=True)
        >>> f = tp.plot_loc(pData)
        """

        px = pData[:, FLEKSTP.ix_]
        py = pData[:, FLEKSTP.iy_]
        pz = pData[:, FLEKSTP.iz_]

        # create subplot mosaic with different keyword arguments
        skeys = ["A", "B", "C", "D"]
        f, axs = plt.subplot_mosaic(
            "AB;CD",
            per_subplot_kw={("D"): {"projection": "3d"}},
            gridspec_kw={"width_ratios": [1, 1], "wspace": 0.1, "hspace": 0.1},
            figsize=(10, 10),
            constrained_layout=True,
        )

        # Create 2D scatter plots
        for i, (x, y, labels) in enumerate(
            zip([px, px, py], [py, pz, pz], [("x", "y"), ("x", "z"), ("y", "z")])
        ):
            axs[skeys[i]].scatter(x, y, s=1)
            axs[skeys[i]].set_xlabel(labels[0])
            axs[skeys[i]].set_ylabel(labels[1])

        # Create 3D scatter plot
        axs[skeys[3]].scatter(px, py, pz, s=1)
        axs[skeys[3]].set_xlabel("x")
        axs[skeys[3]].set_ylabel("y")
        axs[skeys[3]].set_zlabel("z")

        return
