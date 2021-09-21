import sys, os
import ase

def add_ipi_paths(base=os.path.expanduser("~")+"/i-pi/"):
    """ Adds to system paths so that one can run from the jupyter 
    notebooks without any settings after having just cloned the i-PI repo
    and compiled the FORTRAN driver. """
    
    sys.path.append(base)
    os.environ['PATH'] += (":"+base+"/bin/")
    
def read_ipi_xyz(filename):
    """ Reads a file in xyz i-PI format and returns it in ASE format. """
    
    from ipi.utils.io import read_file
    file_handle = open(filename, "r")
    frames = []
    while True:
        try:
            ret = read_file("xyz", file_handle)
            frames.append(ase.Atoms(ret["atoms"].names, 
                                    positions=ret["atoms"].q.reshape((-1,3))*0.529177, 
                                    cell=ret["cell"].h.T*0.529177))
        except EOFError:
            break
        except:
            raise
    return frames