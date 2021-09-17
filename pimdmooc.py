import sys, os





def add_ipi_paths(base="~/i-pi/"):
    sys.path.append(base)
    os.environ['PATH'] += (base+"/bin/")