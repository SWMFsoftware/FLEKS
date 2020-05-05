import yt
import os
import numpy as np
import glob

from yt.funcs import setdefaultattr
from yt.frontends.boxlib.api import BoxlibHierarchy, BoxlibDataset
from yt.fields.field_info_container import FieldInfoContainer

plot_unit_planet = {
    "time": "s",
    "t": "s",
    "mass": "amu",
    "rho": "amu/cm**3",
    "u": "km/s",
    "ux": "km/s",
    "uy": "km/s",
    "uz": "km/s",
    "p": "nPa",
    "pxx": "nPa",
    "pxy": "nPa",
    "pxz": "nPa",
    "pyy": "nPa",
    "pyz": "nPa",
    "pzz": "nPa",
    "B": "nT",
    "Bx": "nT",
    "By": "nT",
    "Bz": "nT",
    "E": "nT*km/s",
    "Ex": "nT*km/s",
    "Ey": "nT*km/s",
    "Ez": "nT*km/s",
    "X": "Planet_Radius",
    "Y": "Planet_Radius",
    "Z": "Planet_Radius",
    "p_x": "Planet_Radius",
    "p_y": "Planet_Radius",
    "p_z": "Planet_Radius",
    "p_ux": "km/s",
    "p_uy": "km/s",
    "p_uz": "km/s",
    "p_w": "amu"
}

plot_unit_si = {
    "time": "s",
    "t": "s",
    "mass": "kg",
    "rho": "kg/m**3",
    "u": "m/s",
    "ux": "m/s",
    "uy": "m/s",
    "uz": "m/s",
    "p": "Pa",
    "pxx": "Pa",
    "pxy": "Pa",
    "pxz": "Pa",
    "pyy": "Pa",
    "pyz": "Pa",
    "pzz": "Pa",
    "B": "T",
    "Bx": "T",
    "By": "T",
    "Bz": "T",
    "E": "T*m/s",
    "Ex": "T*m/s",
    "Ey": "T*m/s",
    "Ez": "T*m/s",
    "X": "m",
    "Y": "m",
    "Z": "m",
    "p_x": "m",
    "p_y": "m",
    "p_z": "m",
    "p_ux": "m/s",
    "p_uy": "m/s",
    "p_uz": "m/s",
    "p_w": "kg"
}


def get_unit(var, unit_type="planet"):
    if var[-1].isdigit():
        # Example: pxxs0 -> pxx
        var = var[0:-2]

    if unit_type == "planet":
        return plot_unit_planet[var]
    elif unit_type == "si":
        return plot_unit_si[var]
    else:
        return "unitary"


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
        ("Bx", (b_units,   [], r"B_x")),
        ("By", (b_units,   [], r"B_y")),
        ("Bz", (b_units,   [], r"B_z")),
        ("Ex", (e_units, [], r"E_x")),
        ("Ey", (e_units, [], r"E_y")),
        ("Ez", (e_units, [], r"E_z")),
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
        output_dir = glob.glob(output_dir)[0]
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

    def plot_slice(self, norm, cut_loc, vars, unit_type="planet",
                   *args, **kwargs):
        r"""Plot 2D slice

        Parameters
        ----------
        norm : string, one of 'x', 'y' or 'z'. Norm direction of the slice.

        cut_loc : float. The location of the slice.

        vars : a list of plotting variables. Example: ["Bx", "rhos0"]

        unit_type : The unit system of the plots. "planet" or "si".

        Examples
        --------

        >>> import fleks
        >>> vars = ["rhos0", "uzs0", "Bz", "pxxs1", "Ex"]
        >>> splt = ds.plot_slice("y", 0.0, vars)
        >>> splt.display()
        """

        center = self.domain_center
        idir = "xyz".find(norm.lower())
        center[idir] = cut_loc
        splt = yt.SlicePlot(self, norm, fields=vars, center=center,
                            origin='native', *args, **kwargs)
        for var in vars:
            splt.set_log(var, False)
            splt.set_unit(var, get_unit(var, unit_type))

        splt.set_axes_unit(get_unit("X", unit_type))
        # splt.display()
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
        plot = yt.ParticlePlot(self, (var_type, x_field), (var_type, y_field),
                               (var_type, z_field), data_source=dd)
        plot.set_axes_unit((get_unit(x_field, unit_type),
                            get_unit(y_field, unit_type)))
        plot.set_unit((var_type, z_field), get_unit(z_field, unit_type))

        return plot
