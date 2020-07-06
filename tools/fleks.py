import glob
from idl_format_data import IDLDataSet
from yt_interface import FLEKSDataset, FLEKSTP

def load(filename):
    files = glob.glob(filename)
    if len(files) == 0:
        raise Exception('Error: can not find the file/director!')
    filename = files[0]
    
    tmp = filename.split('/')[-1]
    if tmp == 'test_particles':
        return FLEKSTP(filename)
    elif tmp.find('.') and tmp.split('.')[-1] == 'out':
        return IDLDataSet(filename)
    elif tmp[-6:] == '_amrex':
        return FLEKSDataset(filename)
    else:
        raise Exception('Error: unknown file format!')
