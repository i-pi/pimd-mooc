import sys, os
import ase
import numpy as np
import re

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
                                    cell=ret["cell"].h.T*0.529177, pbc=True))
        except EOFError:
            break
        except:
            raise
    return frames

def read_ipi_output(filename):
    """ Reads an i-PI output file and returns a dictionary with the properties in a tidy order. """
    
    f = open(filename, "r")
    
    regex = re.compile(".*column *([0-9]*) *--> ([^ {]*)")
    
    fields = []; cols = []
    for line in f:
        if line[0] == "#":
            match = regex.match(line)
            if match is None:
                print("Malformed comment line: ", line)
                raise ValueError()
            fields.append(match.group(2))
            cols.append(slice(int(match.group(1))-1,int(match.group(1))))
        else:
            break # done with header
    f.close()
    
    columns = {}
    raw = np.loadtxt(filename)
    for i, c in enumerate(fields):
        columns[c] = raw[:,cols[i]].T
        if columns[c].shape[0] == 1:
            columns[c].shape = columns[c].shape[1]
    return columns

def correlate(x, y, xbar=None, ybar=None, normalize=True):
    """ Computes the correlation function of two quantities. 
       It can be given the exact averages as parameters."""
    if xbar is None:
        xbar = x.mean()
    if ybar is None:
        ybar = y.mean()
        
    cf = np.correlate(x - xbar, y - ybar, mode='same')
    return cf[len(x)//2:]/(((x-xbar)*(y-ybar)).sum() if normalize else 1)

def autocorrelate(x, xbar=None, normalize=True):
    """ Computes the autocorrelation function of a trajectory. 
    It can be given the exact average as a parameter"""
    
    if xbar is None:
        xbar = x.mean()
    acf = np.correlate(x - xbar, x-xbar, mode='same')
    return acf[len(x)//2:]/(((x-xbar)*(x-xbar)).sum() if normalize else 1)