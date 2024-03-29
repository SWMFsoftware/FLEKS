import glob
import os
from idl_format_data import IDLDataSet
from yt_interface import FLEKSDataset
from test_particles import FLEKSTP
from data_container import compare


def load(filename, iDomain=0, iSpecies=0, readFieldData=False):
    files = glob.glob(filename)
    if len(files) == 0:
        raise Exception('Error: can not find the file/directory!')
    filename = files[0]

    basename = os.path.basename(os.path.normpath(filename))

    if basename == 'test_particles':
        return FLEKSTP(filename, iDomain=iDomain, iSpecies=iSpecies)
    elif basename.find('.') != -1 and basename.split('.')[-1] in ['out', 'outs']:
        return IDLDataSet(filename)
    elif basename[-6:] == '_amrex':
        return FLEKSDataset(filename, readFieldData)
    else:
        raise Exception('Error: unknown file format!')
