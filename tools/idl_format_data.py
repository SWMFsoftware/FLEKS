import glob
import numpy as np
import matplotlib.pyplot as plt
import struct
import math
import utilities
import yt
import data_container

x_ = 0
y_ = 1
z_ = 2

class selector:
    def __getitem__(self, keys):
        self.indices = list(keys)
        if len(self.indices) < 3:
            self.indices += [0]*(3-len(self.indices))
        return self.indices

class dataframe:
    def __init__(self):
        self.array = None
        self.name = None
        # Selector for the spatial indices
        self.cut = selector()
        self.cut[:, :, :]

    def setData(self, dataIn, nameIn):
        assert (nameIn.size == dataIn.shape[0])
        assert (dataIn.ndim <= 4)
        shape = list(dataIn.shape) + [1]*(4-dataIn.ndim)
        self.array = np.reshape(dataIn, shape)
        self.name = tuple(nameIn)

    def fixDataSize(self):
        self.name = tuple(self.name)
        assert (len(self.name) == self.array.shape[0])
        assert (self.array.ndim <= 4)
        shape = list(self.array.shape) + [1]*(4-self.array.ndim)
        self.array = np.reshape(self.array, shape)

    def __getitem__(self, keys):
        '''Example: 
            d["varname",3:,1:4:2,3]
        '''
        if type(keys) is str:
            # If the spatial indices are not specified, use self.cut
            keys = [keys] + self.cut.indices
        else:
            keys = list(keys) + [slice(None, None, None)]*(4-len(keys))
        ivar = self.name.index(keys[0])
        return np.squeeze(self.array[ivar, keys[1], keys[2], keys[3]])


class IDLDataSet(object):
    r"""
    A class used to handle the *.out format data. 

    Example:
    >>> ds = IDLDataSet('3d.out')
    >>> dc2d = ds.get_slice('y',1)
    >>> dc3d = ds.get_domain()
    """

    def __init__(self, filename="none"):
        fileList = glob.glob(filename)
        nfiles = len(fileList)
        assert nfiles > 0, 'Error: can not file file!'
        if nfiles > 1:
            fileList.sort()
            print("nfiles = ", nfiles)
            ifile = input("ifile= ? start from 1             ")
            while ifile > nfiles or ifile < 1:
                print("ifile = ", ifile, " is a bad value. Input again")
                ifile = input("ifile= ? start from 1             ")
        else:
            ifile = 1

        self.filename = fileList[ifile - 1]
        self.isOuts = self.filename[-4:] == "outs"
        self.data = dataframe()        
        self.nInstance = None if self.isOuts else 1
        self.npict = 1            
        self.fileformat = None        
        self.variables = None
        self.unit = None        
        self.iter = None
        self.runtime = None
        self.ndim = None
        self.gencoord = None
        self.grid = None

        self.read_data()

    def __post_process_param__(self):

        planet_radius = 1.0 
        
        # Not always correct. 
        for var, val in zip(self.param_name, self.para):
            if var == "xSI":
                planet_radius = float(100*val)

        self.registry = yt.units.unit_registry.UnitRegistry()
        self.registry.add("Planet_Radius", planet_radius,
                          yt.units.dimensions.length)

    def info(self):
        print("\n-----------------------------")
        print("filename    : ", self.filename)
        print("variables   : ", self.variables)
        print("unit        : ", self.unit)
        print("nInstance   : ", self.nInstance)
        print("npict       : ", self.npict)
        print("time        : ", self.runtime)
        print("nIter       : ", self.iter)
        print("ndim        : ", self.ndim)
        print("gencoord    : ", self.gencoord)
        print("grid        : ", self.grid)
        print("-----------------------------")

    def __repr__(self):
        self.info()
        return "\n"

    def get_domain(self):
        r""" 
        This methods returns a dataContainer3D object that contains all the 3D data. 
        """
        dataSets = {}
        for varname in self.data.name:
            idx = self.data.name.index(varname)
            unit = utilities.get_unit(varname, self.unit)
            dataSets[varname] = yt.YTArray(
                np.squeeze(self.data.array[idx, :, :, :]), unit, registry=self.registry)

        labels = ['', '', '']
        axes = [None, None, None]
        for idim in range(self.ndim):
            name = self.variables[idim]
            idx = self.data.name.index(name)
            labels[idim] = name.upper()
            unit = utilities.get_unit('X', self.unit)
            name = name.upper()
            if name in ["X", "Y", "Z"]:
                axes[idim] = yt.YTArray(
                    self.data.array[idx, :, :, :], unit, registry=self.registry)

        if self.ndim == 1:
            dc = data_container.dataContainer1D(
                dataSets, np.squeeze(axes[0]), labels[0],
                step=self.iter, time=self.runtime, filename=self.filename)
        elif self.ndim == 2:
            dc = data_container.dataContainer2D(
                dataSets, np.squeeze(axes[0])[:, 0], np.squeeze(axes[1])[0, :],
                labels[0], labels[1], step=self.iter, time=self.runtime, filename=self.filename)
        else:
            dc = data_container.dataContainer3D(
                dataSets, axes[0][:, 0, 0], axes[1][0, :, 0], axes[2][0, 0, :],
                step=self.iter, time=self.runtime, filename=self.filename)

        return dc

    def save_data(self, saveName, saveFormat="ascii"):
        # Only support ascii so far.
        self.save_ascii_instance(saveName)

    def save_ascii_instance(self, saveName):
        with open(saveName, 'w') as f:
            f.write(self.unit+"\n")
            f.write("{:d}\t{:e}\t{:d}\t{:d}\t{:d}\n".format(
                self.iter, self.runtime, self.ndim, self.nparam, self.nvar))
            [f.write("{:d}\t".format(i)) for i in self.grid]
            f.write('\n')
            if self.nparam > 0:
                [f.write("{:e}\t".format(i)) for i in self.para]
                f.write('\n')
            [f.write(i+" ") for i in self.variables]
            f.write('\n')

            nk = self.grid[2] if self.ndim > 2 else 1
            nj = self.grid[1] if self.ndim > 1 else 1
            ni = self.grid[0]
            for kk in range(nk):
                for jj in range(nj):
                    for ii in range(ni):
                        [f.write("{:e}\t".format(i))
                         for i in self.data.array[:, ii, jj, kk]]
                        f.write('\n')

    def read_data(self):
        if self.fileformat is None:
            with open(self.filename, 'rb') as f:
                EndChar = '<'  # Endian marker (default: little.)
                RecLenRaw = f.read(4)
                RecLen = (struct.unpack(EndChar+'l', RecLenRaw))[0]
                if RecLen != 79 and RecLen != 500:
                    self.fileformat = "ascii"
                else:
                    self.fileformat = "binary"

        if self.fileformat == "ascii":
            self.read_ascii()
        elif self.fileformat == "binary":
            self.read_binary()
        else:
            print("Unknown format = ", self.fileformat)

        ndim = self.ndim
        nvar = self.nvar
        self.data.name = tuple(self.variables)[0:ndim+nvar]
        self.data.fixDataSize()
        self.param_name = self.variables[ndim+nvar:]
        self.__post_process_param__()

    def read_ascii(self):
        if self.nInstance is None:
            # Count how many instances are there.
            with open(self.filename, 'r') as f:
                for i, l in enumerate(f):
                    pass
                nLineFile = i+1

            with open(self.filename, 'r') as f:
                self.nInstanceLength = self.read_ascii_instance(f)

            self.nInstance = round(nLineFile/self.nInstanceLength)

        nLineSkip = (self.npict)*self.nInstanceLength if self.isOuts else 0
        with open(self.filename, 'r') as f:
            if nLineSkip > 0:
                for i, line in enumerate(f):
                    if i == nLineSkip-1:
                        break
            self.read_ascii_instance(f)

    def read_ascii_instance(self, infile):
        nline = 0
        # Read the top header line:
        headline = infile.readline().strip()
        nline += 1
        self.unit = headline.split()[0]

        # Read & convert iters, runtime, etc. from next line:
        parts = infile.readline().split()
        nline += 1
        self.iter = int(parts[0])
        self.runtime = float(parts[1])
        self.ndim = int(parts[2])
        self.gencoord = self.ndim < 0
        self.ndim = abs(self.ndim)
        self.nparam = int(parts[3])
        self.nvar = int(parts[4])

        # Read & convert grid dimensions.
        grid = [int(x) for x in infile.readline().split()]
        nline += 1
        self.grid = np.array(grid)
        self.npoints = abs(self.grid.prod())

        # Quick ref vars:
        time = self.runtime
        npts = self.npoints
        ndim = self.grid.size
        nvar = self.nvar
        npar = self.nparam

        # Read parameters stored in file.
        self.para = np.zeros(npar)
        if npar > 0:
            self.para[:] = infile.readline().split()
            nline += 1

        # Read variable names.
        names = infile.readline().split()
        nline += 1

        # Save grid names (e.g. 'x' or 'r') and save associated params.
        self.dims = names[0:ndim]
        self.variables = np.array(names)

        # Create string representation of time.
        self.strtime = '%4.4ih%2.2im%06.3fs' %\
            (np.floor(time/3600.), np.floor(time % 3600. / 60.0),
             time % 60.0)

        nRow = ndim + nvar
        nCol = self.npoints
        self.data.array = np.zeros((nRow, nCol))

        for i, line in enumerate(infile.readlines()):
            parts = line.split()

            if i >= self.npoints:
                break

            for j, p in enumerate(parts):
                self.data.array[j][i] = float(p)

        nline += self.npoints
        shapeNew = np.append([nRow], self.grid)
        self.data.array = np.reshape(self.data.array, shapeNew, order='F')

        return nline

    def read_binary(self):
        if self.nInstance is None:
            with open(self.filename, 'rb') as f:
                self.read_binary_instance(f)
                self.nInstanceLength = f.tell()
                f.seek(0, 2)
                endPos = f.tell()
            self.nInstance = round(endPos/self.nInstanceLength)

        with open(self.filename, 'rb') as f:
            if self.isOuts:
                f.seek((self.npict)*self.nInstanceLength, 0)
            self.read_binary_instance(f)

    def read_binary_instance(self, infile):
        # On the first try, we may fail because of wrong-endianess.
        # If that is the case, swap that endian and try again.
        EndChar = '<'  # Endian marker (default: little.)
        self.endian = 'little'
        RecLenRaw = infile.read(4)

        RecLen = (struct.unpack(EndChar+'l', RecLenRaw))[0]
        if (RecLen > 10000) or (RecLen < 0):
            EndChar = '>'
            self.endian = 'big'
            RecLen = (struct.unpack(EndChar+'l', RecLenRaw))[0]

        headline = (struct.unpack('{0}{1}s'.format(EndChar, RecLen),
                                  infile.read(RecLen)))[0].strip()
        if str is not bytes:
            headline = headline.decode()
        self.unit = headline.split()[0]

        (OldLen, RecLen) = struct.unpack(EndChar+'2l', infile.read(8))
        pformat = 'f'
        # parse rest of header; detect double-precision file.
        if RecLen > 20:
            pformat = 'd'
        (self.iter, self.runtime,
         self.ndim, self.nparam, self.nvar) = \
            struct.unpack('{0}l{1}3l'.format(
                EndChar, pformat), infile.read(RecLen))
        self.gencoord = self.ndim < 0                 
        self.ndim = abs(self.ndim)
        # Get gridsize
        (OldLen, RecLen) = struct.unpack(EndChar+'2l', infile.read(8))

        self.grid = np.array(struct.unpack('{0}{1}l'.format(EndChar,
                                                            abs(self.ndim)), infile.read(RecLen)))
        self.npoints = abs(self.grid.prod())

        # Quick ref vars:
        time = self.runtime
        npts = self.npoints
        ndim = self.grid.size
        nvar = self.nvar
        npar = self.nparam

        # Read parameters stored in file.
        self.para = np.zeros(npar)
        if npar > 0:
            (OldLen, RecLen) = struct.unpack(EndChar+'2l', infile.read(8))
            self.para[:] = struct.unpack('{0}{1}{2}'.format(EndChar, npar, pformat),
                                         infile.read(RecLen))

        (OldLen, RecLen) = struct.unpack(EndChar+'2l', infile.read(8))
        names = (struct.unpack('{0}{1}s'.format(EndChar, RecLen),
                               infile.read(RecLen)))[0]
        if str is not bytes:
            names = names.decode()

        names.strip()
        names = names.split()

        # Save grid names (e.g. 'x' or 'r') and save associated params.
        self.dims = names[0:ndim]
        self.variables = np.array(names)

        self.strtime = '{0:04d}h{1:02d}m{2:06.3f}s'.format(
            int(time//3600), int(time % 3600//60), time % 60)

        nRow = ndim + nvar
        nCol = self.npoints
        self.data.array = np.zeros((nRow, nCol), dtype=np.float32)

        # Get the grid points...
        (OldLen, RecLen) = struct.unpack(EndChar+'2l', infile.read(8))
        # prod = [1] + pbdat['grid.cumprod().tolist()
        for i in range(0, ndim):
            # Read the data into a temporary grid.
            tempgrid = np.array(struct.unpack('{0}{1}{2}'.format(EndChar, npts, pformat),
                                              infile.read(int(RecLen//ndim))))
            self.data.array[i, :] = tempgrid

            # Get the actual data and sort.
        for i in range(ndim, nvar+ndim):
            (OldLen, RecLen) = struct.unpack(EndChar+'2l', infile.read(8))
            tmp = np.array(struct.unpack('{0}{1}{2}'.format(EndChar,
                                                            npts, pformat), infile.read(RecLen)))
            self.data.array[i, :] = tmp
        # Consume the last record length
        infile.read(4)

        shapeNew = np.append([nRow], self.grid)
        self.data.array = np.reshape(self.data.array, shapeNew, order='F')

    def plot(self, iv1name, iv2name, *dvname):
        ''' 
        dvname: dependent variables name 
        iv1name: the first independent variable name
        iv2name: the second independent variable name
        '''
        x = self.data[iv1name]
        y = self.data[iv2name]

        nvar = len(dvname)
        nRow = int(round(np.sqrt(nvar)))
        nCol = math.ceil(nvar/nRow)
        print('nvar = ', nvar, ' nRow = ', nRow, ' nCol = ', nCol)

        f, axes = plt.subplots(nRow, nCol)
        axes = np.array(axes)  # in case nRow = nCol = 1
        aspect = (y.max() - y.min())/(x.max() - x.min())
        axes = axes.reshape(-1)

        for isub, ax in zip(range(nvar), axes):
            w = self.data[dvname[isub]]
            cs = ax.contourf(x, y, w, levels=100, cmap="rainbow")
            cb = f.colorbar(cs, ax=ax, shrink=aspect*0.6)
            cb.ax.set_yticks([w.min(), w.max()])
            #print('type f = ', type(f))
            #print('type cs = ', type(cs))
            #f.set_ticks([w.min(), w.max()])
            ax.set_aspect(aspect)
            ax.set_xlabel(iv1name)
            ax.set_ylabel(iv2name)
            ax.set_title(dvname[isub])

        # Delete axes for empty subplots
        for ax in axes[nvar:nRow*nCol]:
            f.delaxes(ax)

    def extract_data(self, sat):
        satData = None
        if type(sat) == np.ndarray and sat.ndim == 2 and sat.shape[1] >= self.ndim:
            nVar = self.nvar + self.ndim
            nPoint = sat.shape[0]
            satData = np.zeros((nPoint, nVar))
            for i in range(nPoint):
                satData[i, :] = self.get_data(sat[i, :])
        return satData

    def get_data(self, loc):
        i1, j1, k1 = 0, 0, 0
        while self.data['x'][i1, 0, 0] < loc[x_]:
            i1 = i1 + 1
        while self.data['y'][0, j1, 0] < loc[y_]:
            j1 = j1 + 1
        while self.data['z'][0, 0, k1] < loc[z_]:
            k1 = k1 + 1

        i0 = i1 - 1
        j0 = j1 - 1
        k0 = k1 - 1

        wx0 = (self.data['x'][i1, 0, 0] - loc[x_]) / \
            (self.data['x'][i1, 0, 0] - self.data['x'][i0, 0, 0])
        wy0 = (self.data['y'][0, j1, 0] - loc[y_]) / \
            (self.data['y'][0, j1, 0] - self.data['y'][0, j0, 0])
        wz0 = (self.data['z'][0, 0, k1] - loc[z_]) / \
            (self.data['z'][0, 0, k1] - self.data['z'][0, 0, k0])

        wx1 = 1.0 - wx0
        wy1 = 1.0 - wy0
        wz1 = 1.0 - wz0

        w = np.zeros((2, 2, 2))
        w[0, 0, 0] = wx0*wy0*wz0
        w[0, 0, 1] = wx0*wy0*wz1
        w[0, 1, 0] = wx0*wy1*wz0
        w[0, 1, 1] = wx0*wy1*wz1
        w[1, 0, 0] = wx1*wy0*wz0
        w[1, 0, 1] = wx1*wy0*wz1
        w[1, 1, 0] = wx1*wy1*wz0
        w[1, 1, 1] = wx1*wy1*wz1

        res = np.zeros((self.nvar + self.ndim))

        for i in range(2):
            for j in range(2):
                for k in range(2):
                    res = res + w[i, j, k]*self.data.array[:, i0+i, j0+j, k0+k]

        return res
