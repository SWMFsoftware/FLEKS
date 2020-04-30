import yt
import os
import numpy as np
import glob

from yt.funcs import setdefaultattr
from yt.frontends.boxlib.api import BoxlibHierarchy, BoxlibDataset
from yt.fields.field_info_container import FieldInfoContainer

l_units = "code_length"
v_units = "code_velocity"
p_units = "code_pressure"
b_units = "code_magnetic"
e_units = "code_magnetic * code_velocity"
rho_units = "code_density"
mass_units = "code_mass"

class FLEKSFieldInfo(FieldInfoContainer):

    known_other_fields = (
        ("Bx", (b_units,   [], None)),
        ("By", (b_units,   [], None)),
        ("Bz", (b_units,   [], None)),
        ("Ex", (e_units, [], None)),
        ("Ey", (e_units, [], None)),
        ("Ez", (e_units, [], None)),
        ("rhos0", (rho_units, [], r"\rho")),
        ("uxs0", (v_units, [], r"U_x")),
        ("uys0", (v_units, [], r"U_y")),
        ("uzs0", (v_units, [], r"U_z")),
        ("pxxs0", (p_units, [], r"P_{xx}")),
        ("pyys0", (p_units, [], r"P_{yy}")),
        ("pzzs0", (p_units, [], r"P_{zz}")),
        ("pxys0", (p_units, [], r"P_{xy}")),
        ("pxzs0", (p_units, [], r"P_{xz}")),
        ("pyzs0", (p_units, [], r"P_{yz}")),
        ("rhos1", (rho_units, [], r"\rho")),
        ("uxs1", (v_units, [], r"U_x")),
        ("uys1", (v_units, [], r"U_y")),
        ("uzs1", (v_units, [], r"U_z")),
        ("pxxs1", (p_units, [], r"P_{xx}")),
        ("pyys1", (p_units, [], r"P_{yy}")),
        ("pzzs1", (p_units, [], r"P_{zz}")),
        ("pxys1", (p_units, [], r"P_{xy}")),
        ("pxzs1", (p_units, [], r"P_{xz}")),
        ("pyzs1", (p_units, [], r"P_{yz}")),    
    )

    known_particle_fields = (
        ("particle_weight", (mass_units, [], None)),
        ("particle_position_x", (l_units, [], None)),
        ("particle_position_y", (l_units, [], None)),
        ("particle_position_z", (l_units, [], None)),
        ("particle_velocity_x", (v_units, [], None)),
        ("particle_velocity_y", (v_units, [], None)),
        ("particle_velocity_z", (v_units, [], None)),
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

        #TODO: I do not know the purpose of the following function call. --Yuxi 
        setup_magnetic_field_aliases(self, "FLEKS", ["B%s" % ax for ax in "xyz"])

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


def _skip_line(line):
    if len(line) == 0:
        return True
    if line[0] == '\n':
        return True
    if line[0] == "=":
        return True
    if line[0] == ' ':
        return True


class FLEKSDataset(BoxlibDataset):

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
            self.radius = float(f.readline())

        self.code_length = 1 # m
        self.code_velocity = 1 # m/s
        self.code_mass = 1 # kg 
        self.code_density = 1 # kg/m^3 
        self.code_magnetic = 1 # T
        #self.code_electric = 1 # v/m 
        self.code_pressure = 1 # Pa

        if plot_string.find("si") != -1:
            self.parameters["fleks_unit"] = "si"
        elif plot_string.find("planet") != -1:
            self.parameters["fleks_unit"] = "planet"
            self.code_length = self.radius
            self.code_velocity = 1e3 # km/s
            self.code_mass =  1.66054e-27 # amu
            self.code_density = 1.66054e-21 # amu/cc 
            self.code_magnetic = 1e-9 # nT
            #self.code_electric = self.code_magnetic * self.code_velocity # v/m  
            self.code_pressure = 1e-9 # nPa           
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
        setdefaultattr(self, 'length_unit', self.quan(self.code_length, "m"))
        setdefaultattr(self, 'mass_unit', self.quan(self.code_mass, "kg"))
        setdefaultattr(self, 'time_unit', self.quan(1, "s"))
        setdefaultattr(self, 'velocity_unit', self.quan(self.code_velocity, "m/s"))
        setdefaultattr(self, 'magnetic_unit', self.quan(self.code_magnetic, "T"))
        setdefaultattr(self, 'density_unit', self.quan(self.code_density, "kg/m**3"))
        setdefaultattr(self, 'pressure_unit', self.quan(self.code_pressure, "Pa"))

ds = FLEKSDataset("res/run4/PC/3d_fluid*0_amrex")
print("bx_min = ", ds.r['Bx'].v.min())
print("bx_max = ", ds.r['Bx'].v.max())
splt=yt.SlicePlot(ds,"z","pxxs0")   
splt.set_log("pxxs0", False)
splt.set_unit('pxxs0', 'nPa')
splt.set_axes_unit("code_length")
splt.set_xlabel("X")
# splt.set_zlim("pxxs0", 7e-5, 9e-5)
splt.save()

if False:
    filename = glob.glob('PC/plots/3d_particle_*0000_amrex')[0]
    ds = FLEKSDataset(filename)
    ds.field_list
    dd = ds.box([-0.001, -0.001, -0.005], [0, 0, 0])
    #dd = ds.all_data()
    position = yt.ParticlePlot(ds, ('particle', 'particle_position_x'),
                               ('particle', 'particle_position_y'), ('particle', 'particle_weight'), data_source=dd)
    # position.set_axes_unit("km")
    # save result
    position.save("position.png")

    position = yt.PhasePlot(dd, ('particle', 'particle_velocity_x'),
                            ('particle', 'particle_velocity_y'), ('particle', 'particle_weight'), weight_field=None, x_bins=128, y_bins=128)
    # save result
    position.save("phase.png")
