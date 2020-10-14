import yt
import os
import numpy as np
import glob
import struct
import data_container
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from random import random

from yt.funcs import setdefaultattr
from yt.frontends.boxlib.api import BoxlibHierarchy, BoxlibDataset
from yt.fields.field_info_container import FieldInfoContainer

from utilities import plot_unit_planet, plot_unit_si, get_unit

class FLEKSFieldInfo(FieldInfoContainer):
    l_units = "code_length"
    v_units = "code_velocity"
    p_units = "code_pressure"
    b_units = "code_magnetic"
    e_units = "code_magnetic * code_velocity"
    rho_units = "code_density"
    mass_units = "code_mass"

    # TODO: find a way to avoid repeating s0, s1...
    known_other_fields = (
        ("Bx", (b_units,   ["magnetic_field_x"], r"B_x")),
        ("By", (b_units,   ["magnetic_field_y"], r"B_y")),
        ("Bz", (b_units,   ["magnetic_field_z"], r"B_z")),
        ("Ex", (e_units, [], r"E_x")),
        ("Ey", (e_units, [], r"E_y")),
        ("Ez", (e_units, [], r"E_z")),
        ("X", (l_units, [], r"X")),
        ("Y", (l_units, [], r"Y")),
        ("Z", (l_units, [], r"Z")),
        ("rhos0", (rho_units, [], r"\rho")),
        ("uxs0", (v_units, [], r"u_x")),
        ("uys0", (v_units, [], r"u_y")),
        ("uzs0", (v_units, [], r"u_z")),
        ("pxxs0", (p_units, [], r"P_{xx}")),
        ("pyys0", (p_units, [], r"P_{yy}")),
        ("pzzs0", (p_units, [], r"P_{zz}")),
        ("pxys0", (p_units, [], r"P_{xy}")),
        ("pxzs0", (p_units, [], r"P_{xz}")),
        ("pyzs0", (p_units, [], r"P_{yz}")),
        ("rhos1", (rho_units, [], r"\rho")),
        ("uxs1", (v_units, [], r"u_x")),
        ("uys1", (v_units, [], r"u_y")),
        ("uzs1", (v_units, [], r"u_z")),
        ("pxxs1", (p_units, [], r"P_{xx}")),
        ("pyys1", (p_units, [], r"P_{yy}")),
        ("pzzs1", (p_units, [], r"P_{zz}")),
        ("pxys1", (p_units, [], r"P_{xy}")),
        ("pxzs1", (p_units, [], r"P_{xz}")),
        ("pyzs1", (p_units, [], r"P_{yz}")),
    )

    known_particle_fields = (
        ("particle_weight", (mass_units, ["p_w"], r"weight")),
        ("particle_position_x", (l_units, ["p_x"], "x")),
        ("particle_position_y", (l_units, ["p_y"], "y")),
        ("particle_position_z", (l_units, ["p_z"], "z")),
        ("particle_velocity_x", (v_units, ["p_ux"], r"u_x")),
        ("particle_velocity_y", (v_units, ["p_uy"], r"u_y")),
        ("particle_velocity_z", (v_units, ["p_uz"], r"u_z")),
    )

    extra_union_fields = (
        (mass_units, "particle_mass"),
    )

    def __init__(self, ds, field_list):
        super(FLEKSFieldInfo, self).__init__(ds, field_list)

        # setup nodal flag information
        for field in ds.index.raw_fields:
            finfo = self.__getitem__(('raw', field))
            finfo.nodal_flag = ds.nodal_flags[field]

    def setup_fluid_fields(self):
        from yt.fields.magnetic_field import \
            setup_magnetic_field_aliases

        for field in self.known_other_fields:
            fname = field[0]
            self.alias(("mesh", fname), ('boxlib', fname))

        # TODO: I do not know the purpose of the following function call. --Yuxi
        setup_magnetic_field_aliases(
            self, "FLEKS", ["B%s" % ax for ax in "xyz"])

    def setup_fluid_aliases(self):
        super(FLEKSFieldInfo, self).setup_fluid_aliases("mesh")

    def setup_particle_fields(self, ptype):
        super(FLEKSFieldInfo, self).setup_particle_fields(ptype)


class FLEKSHierarchy(BoxlibHierarchy):

    def __init__(self, ds, dataset_type="boxlib_native"):
        super(FLEKSHierarchy, self).__init__(ds, dataset_type)

        is_checkpoint = True
        for ptype in self.ds.particle_types:
            self._read_particles(ptype, is_checkpoint)

    def _detect_output_fields(self):
        super(FLEKSHierarchy, self)._detect_output_fields()

        # now detect the optional, non-cell-centered fields
        self.raw_file = self.ds.output_dir + "/raw_fields/"
        self.raw_fields = []
        self.field_list += [('raw', f) for f in self.raw_fields]
        self.raw_field_map = {}
        self.ds.nodal_flags = {}
        self.raw_field_nghost = {}


class FLEKSDataset(BoxlibDataset):
    r"""
    Read the plot the AMReX format output from FLEKS.

    Parameters
    ----------
    output_dir : String
        The path to the data.

    Examples
    --------

    >>> import fleks
    >>> ds = fleks.FLEKSDataset("res/run1/PC/3d_particle*n00004750_amrex")
    """

    _index_class = FLEKSHierarchy
    _field_info_class = FLEKSFieldInfo

    def __init__(self, output_dir,
                 cparam_filename=None,
                 fparam_filename=None,
                 dataset_type='boxlib_native',
                 storage_filename=None,
                 units_override=None,
                 unit_system="mks"):
        self.default_fluid_type = "mesh"
        self.default_field = ("mesh", "density")
        self.fluid_types = ("mesh", "index", "raw")

        super(FLEKSDataset, self).__init__(output_dir,
                                           cparam_filename,
                                           fparam_filename,
                                           dataset_type,
                                           storage_filename,
                                           units_override,
                                           unit_system)

    def _parse_parameter_file(self):
        super(FLEKSDataset, self)._parse_parameter_file()

        fleks_header = os.path.join(self.output_dir, "FLEKSHeader")
        with open(fleks_header, "r") as f:
            plot_string = f.readline().lower()
            self.radius = float(f.readline())  # should be in unit [m]

        # It seems the second argument should be in the unit of cm.
        self.unit_registry.add("Planet_Radius", 100*self.radius,
                               yt.units.dimensions.length)

        if plot_string.find("si") != -1:
            self.parameters["fleks_unit"] = "si"
        elif plot_string.find("planet") != -1:
            self.parameters["fleks_unit"] = "planet"
        elif plot_string.find("pic") != -1:
            self.parameters["fleks_unit"] = "pic"
        else:
            self.parameters["fleks_unit"] = "unknown"

        particle_types = glob.glob(self.output_dir + "/*/Header")
        particle_types = [cpt.split(os.sep)[-2] for cpt in particle_types]
        if len(particle_types) > 0:
            self.parameters["particles"] = 1
            self.particle_types = tuple(particle_types)
            self.particle_types_raw = self.particle_types
        else:
            self.particle_types = ()
            self.particle_types_raw = ()

    def _set_code_unit_attributes(self):
        unit = self.parameters["fleks_unit"]

        setdefaultattr(self, 'time_unit', self.quan(1, get_unit("time", unit)))
        setdefaultattr(self, 'length_unit', self.quan(1, get_unit("X", unit)))
        setdefaultattr(self, 'mass_unit', self.quan(1, get_unit("mass", unit)))
        setdefaultattr(self, 'velocity_unit',
                       self.quan(1, get_unit("u", unit)))
        setdefaultattr(self, 'magnetic_unit',
                       self.quan(1, get_unit("B", unit)))
        setdefaultattr(self, 'density_unit',
                       self.quan(1, get_unit("rho", unit)))
        setdefaultattr(self, 'pressure_unit',
                       self.quan(1, get_unit("p", unit)))

    def get_slice(self, norm, cut_loc):
        r""" 
        This method returns a dataContainer2D object that contains a slice along
        the 'norm' direction at 'cut_loc' 

        Parameters
        ---------------------
        norm: String 
        'x', 'y' or 'z' 

        cut_log: Float         
        """ 

        axDir = {'X': 0, 'Y': 1, 'Z': 2}
        idir = axDir[norm.upper()]

        if type(cut_loc) != yt.units.yt_array.YTArray:
            cut_loc = self.arr(cut_loc, 'code_length')

        # Define the slice range -------------------
        slice_dimension = self.domain_dimensions
        slice_dimension[idir] = 1

        left_edge = self.domain_left_edge
        right_edge = self.domain_right_edge

        dd = (right_edge[idir] - left_edge[idir])*1e-6
        left_edge[idir] = cut_loc - dd
        right_edge[idir] = cut_loc + dd
        # ----------------------------------------------

        abArr = self.arbitrary_grid(left_edge, right_edge, slice_dimension)

        dataSets = {}
        for _, var in self.field_list:
            dataSets[var] = np.squeeze(abArr[var])

        axLabes = {0: ('Y', 'Z'), 1: ('X', 'Z'), 2: ('X', 'Y')}

        axes = []
        for axis_label in axLabes[idir]:
            ax_dir = axDir[axis_label]
            axes.append(np.linspace(self.domain_left_edge[ax_dir],
                                    self.domain_right_edge[ax_dir], self.domain_dimensions[ax_dir]))

        return data_container.dataContainer2D(
            dataSets, axes[0], axes[1], axLabes[idir][0], axLabes[idir][1], norm, cut_loc)

    def get_domain(self):
        r"""
        Reading in all the simulation data into a 3D box. It returns a dataContainer3D object. 
        """ 
        domain = self.covering_grid(
            level=0, left_edge=self.domain_left_edge, dims=self.domain_dimensions)

        dataSets = {}
        for _, var in self.field_list:
            dataSets[var] = domain[var]

        axes = []
        for idim in range(self.dimensionality):
            axes.append(np.linspace(self.domain_left_edge[idim],
                                    self.domain_right_edge[idim], self.domain_dimensions[idim]))

        return data_container.dataContainer3D(dataSets, axes[0], axes[1], axes[2])

    def plot_slice(self, norm, cut_loc, vars, unit_type="planet",
                   *args, **kwargs):
        r"""Plot 2D slice

        Parameters
        ----------
        norm : string, one of 'x', 'y' or 'z'. Norm direction of the slice.

        cut_loc : float. The location of the slice.

        vars : a list or string of plotting variables. 
            Example: "Bx rhos0" or ["Bx", "rhos0"]

        unit_type : The unit system of the plots. "planet" or "si".

        Examples
        --------

        >>> import fleks
        >>> vars = ["rhos0", "uzs0", "Bz", "pxxs1", "Ex"]
        >>> splt = ds.plot_slice("y", 0.0, vars)
        >>> splt.display()
        """

        if type(vars) == str:
            vars = vars.split()

        center = self.domain_center
        idir = "xyz".find(norm.lower())
        center[idir] = cut_loc
        splt = yt.SlicePlot(self, norm, fields=vars, center=center,
                            origin='native', *args, **kwargs)
        for var in vars:
            splt.set_log(var, False)
            splt.set_unit(var, get_unit(var, unit_type))

        splt.set_axes_unit(get_unit("X", unit_type))

        return splt

    def plot_phase(self, left_edge, right_edge, x_field, y_field, z_field,
                   unit_type="planet", x_bins=128, y_bins=128):
        r"""Plot phase space distribution of particle

        Parameters
        ----------
        left_edge & right_edge : a list of float numbers.
            Define the box region to plot.
            Example: left_edge=(9, -0.1, 0.2), right_edge=(10, 1, 2)

        x_field & y_field: string
            The x- y- axes. Potential input: "p_ux", "p_uy", "p_uz".
            Eventhough this is assumed to be a 'phase' plot, the x- and y- axex
            can also be set as 'p_x', 'p_y' or 'p_z' and it will show the 
            spatial distribution of the particles.             

        z_field: string
            It is usually the particle wegith: "p_w".

        unit_type : string
            The unit system of the plots. "planet" or "si".

        Examples
        --------
        >>> phase = ds.plot_phase([8.75, -1, -1], [9.25, 0, 0], 
                                "p_ux", "p_uy", "p_w")
        >>> phase.show()
        """

        dd = self.box(left_edge, right_edge)
        var_type = 'particle'
        plot = yt.PhasePlot(dd, (var_type, x_field), (var_type, y_field),
                            (var_type, z_field), weight_field=None,
                            x_bins=x_bins, y_bins=y_bins)
        plot.set_unit((var_type, x_field), get_unit(x_field, unit_type))
        plot.set_unit((var_type, y_field), get_unit(y_field, unit_type))
        plot.set_unit((var_type, z_field), get_unit(z_field, unit_type))

        return plot

    def plot_particles(self, left_edge, right_edge, x_field, y_field, z_field,
                       unit_type="planet", x_bins=128, y_bins=128):
        r"""Plot the particle position

        Parameters
        ----------
        left_edge & right_edge : a list of float numbers.
            Define the box region to plot.
            Example: left_edge=(9, -0.1, 0.2), right_edge=(10, 1, 2)

        x_field & y_field: string
            The x- y- axes. Potential input: "p_x", "p_y", "p_z".

        z_field: string
            It is usually the particle wegith: "p_w".

        unit_type : string
            The unit system of the plots. "planet" or "si".

        Examples
        --------
        >>> phase = ds.plot_particles([8, -1, -1], [10, 0, 0], "p_x",
                     "p_y", "p_w", unit_type="planet")        
        >>> phase.show()
        """
        dd = self.box(left_edge, right_edge)
        var_type = 'particle'

        nmap = {"p_x": "particle_position_x",
                "p_y": "particle_position_y", "p_z": "particle_position_z"}
        plot = yt.ParticlePlot(self, (var_type, nmap[x_field]), (var_type, nmap[y_field]),
                               (var_type, z_field), data_source=dd)
        plot.set_axes_unit((get_unit(x_field, unit_type),
                            get_unit(y_field, unit_type)))
        plot.set_unit((var_type, z_field), get_unit(z_field, unit_type))

        return plot


class FLEKSTP(object):
    r"""
    A class that is used to read and plot test particles.

    Parameters
    -----------
    outputDirs: String
        The path to the test particle dataset. 

    Examples
    ----------
    >>> tp = FLEKSTP("res/run1/PC/test_particles")
    >>> pIDs = tp.IDs()
    >>> tp.plot((2141,45))
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

    def IDs(self):
        return self.pset

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

    def plot_data(self, data):
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

    def plot(self, partID):
        r""" 
        Plots the trajectory and velocities of a particle. 

        Example
        -----------------
        >>> tp.plot((3,15))
        """        
        pData = self.read_particle_trajectory(partID)
        return self.plot_data(pData)
