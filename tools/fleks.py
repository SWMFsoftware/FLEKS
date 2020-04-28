import yt
import os
import numpy as np
import glob

from yt.funcs import setdefaultattr
from yt.frontends.boxlib.api import WarpXFieldInfo, BoxlibHierarchy, BoxlibDataset

FLEKSFieldInfo = WarpXFieldInfo

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
        setdefaultattr(self, 'length_unit', self.quan(1.0, "m"))
        setdefaultattr(self, 'mass_unit', self.quan(1.0, "kg"))
        setdefaultattr(self, 'time_unit', self.quan(1.0, "s"))
        setdefaultattr(self, 'velocity_unit', self.quan(1.0, "m/s"))


filename = glob.glob('res/run2/PC/3d_particle_*0001_amrex')[0]
ds = FLEKSDataset(filename)
ds.field_list
dd = ds.box([-0.003,-0.002, -0.001],[0,0,0])
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
