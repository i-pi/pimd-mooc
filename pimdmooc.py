import sys, os

def add_ipi_paths(base="~/i-pi/"):
    """ Adds to system paths so that one can run from the jupyter 
    notebooks without any settings after having just cloned the i-PI repo
    and compiled the FORTRAN driver. """
    
    sys.path.append(base)
    os.environ['PATH'] += (":"+base+"/bin/")