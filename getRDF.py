#!/usr/local/bin/python

__author__ = 'Igor Poltavsky'
__version__ = '1.0'

""" get_rdf.py
The script reads simulation time, potential energy, positions and forces from
standard i-PI output files and computes a conventional and PPI radial distribution function (RDF) estimators.
The output is saved to two files which are created in the folder which contains the input files.
The results are printed out in the format: "distance", "RDF".

The script assumes that the input files are in 'xyz' format, with prefix.pos_*.xyz (positions) and
prefix.for_*.xyz (forces) naming scheme.
This would require the following lines in input.xml file:
<trajectory filename='pos' stride='n' format='xyz' cell_units='angstrom'> positions </trajectory>
where n is the same integer number.

Syntax:
   python rdf_ppi.py "prefix" "simulation temperature (in Kelvin)" "element A" "element B" "number of bins for RDF"
   "minimum distance (in Angstroms)" "maximum distance (in Angstroms) "number of time frames to skip in the beginning
   of each file (default 0)"

WARNING:
   Since the python code for computing RDF is inefficient the fortran function in f90 folder must be compiled to
   compute RDFs.
"""

import numpy as np
import sys
import glob
import os
from ipi.utils.units import unit_to_internal, unit_to_user, Constants, Elements
from ipi.utils.io import read_file


def RDF(prefix, filetype, A, B, nbins, r_min, r_max, width, ss=0, unit='angstrom'):

    # Adding fortran functions (when exist)
    try:
        import fortran
    except:
        print('WARNING: No compiled fortran module for fast calculations have been found.\n'
              'Proceeding the calculations is not possible.')
        sys.exit(0)

    skipSteps = int(ss)                                                  # steps to skip
    nbins = int(nbins)

    fns_pos = sorted(glob.glob(prefix + ".pos*"))
    #print(fns_pos)
    fn_out_rdf = prefix + '.' + A + B + ".rdf.dat"
    fn_out_rdf2 = prefix + '.' + A + B + ".rdf2.dat" 
    fn_out_rdf3 = prefix + '.' + A + B + ".rdf3.dat" 
    fn_out_rdf4 = prefix + '.' + A + B + ".rdf4.dat" 
    fn_out_rdf5 = prefix + '.' + A + B + ".rdf5.dat" 

    # check that we found the same number of positions files
    nbeads = len(fns_pos)

    # open input and output files
    ipos = [open(fn, "r") for fn in fns_pos]

    # Species for RDF
    species = (A, B)
    speciesMass = np.array([Elements.mass(species[0]), Elements.mass(species[1])], order='F')

    r_min = unit_to_internal('length', unit, float(r_min))  # Minimal distance for RDF
    r_max = unit_to_internal('length', unit, float(r_max))  # Maximal distance for RDF
    width = unit_to_internal('length', unit, float(width))  # Maximal distance for RDF
    dr = (r_max - r_min) / nbins  # RDF step
    rdf = np.array([[r_min + (0.5 + i) * dr, 0] for i in range(nbins)], order='F')   # conventional RDF

    # Kernel Density Estimation (use Gaussian instead of spikes for the datapoints)
    forma = np.zeros((nbins,nbins))
    for bin in range(nbins):
        forma[bin] = 1.0/(np.sqrt(2.0)*np.pi)/width*np.exp(-0.5*(rdf[:, 0]-rdf[bin, 0])**2/width**2) 
    rdfG = np.copy(rdf) 

    # RDF auxiliary variables
    cell = None         # simulation cell matrix
    inverseCell = None  # inverse simulation sell matrix
    cellVolume = None   # simulation cell volume
    natomsA = 0  # the total number of A type particles in the system
    natomsB = 0  # the total number of B type particles in the system
    # Here A and B are the types of elements used for RDF calculations

    natoms = 0  # total number of atoms
    ifr = 0     # time frame number
    pos, mass = None, None  # positions, and mass arrays
    noteof = True  # end of file test variable

    while noteof:  # Reading input files and calculating PPI correction

        if ifr % 100 == 0:
            print('\r number of beads = {:d},  Processing frame {:d}'.format(nbeads, ifr), end=" ")
            sys.stdout.flush()

        try:
            for i in range(nbeads):
                ret = read_file(filetype, ipos[i], dimension='length')
                if natoms == 0:
                    mass, natoms = ret["atoms"].m, ret["atoms"].natoms
                    pos = np.zeros((nbeads, 3 * natoms), order='F')
                cell = ret["cell"].h
                inverseCell = ret["cell"].get_ih()
                cellVolume = ret["cell"].get_volume()
                pos[i, :] = ret["atoms"].q
        except EOFError:  # finished reading files
            noteof = False

        if noteof:
            if ifr >= skipSteps:  # RDF calculations

                species_A = [3 * i + j for i in np.where(mass == speciesMass[0])[0] for j in range(3)]
                species_B = [3 * i + j for i in np.where(mass == speciesMass[1])[0] for j in range(3)]
                natomsA = len(species_A)
                natomsB = len(species_B)
                posA = np.zeros((nbeads, natomsA), order='F')
                posB = np.zeros((nbeads, natomsB), order='F')
                for bead in range(nbeads):
                    posA[bead, :] = pos[bead, species_A]
                    posB[bead, :] = pos[bead, species_B]

                fortran.updaterdf(rdf, posA, posB, natomsA / 3, natomsB / 3, nbins, r_min, r_max, cell,
                                   inverseCell, nbeads, speciesMass[0], speciesMass[1])#, forma.T)

                ifr += 1

            else:
                ifr += 1
            

#            if ifr > skipSteps and ifr % 1000 == 0:
#
#                rdfG[:, 1] = np.sum(forma * rdf.T[1].astype(int), axis=1)
#
#                # Some constants
#                if(A==B): 
#                    const = 2./(natomsA*(natomsB-1))
#                else:
#                    const = 1./(natomsA*natomsB)
#                const /= float(ifr - skipSteps)
#
#            #    # Normalization
#                _rdf = np.copy(rdfG)
#                _rdf[:, 1] *= const / nbeads
#                
#                # Creating RDF from N(r)
#                const, dr = cellVolume / (4 * np.pi / 3.0), _rdf[1, 0] - _rdf[0, 0]
#                _rdf[:, 1] = const * _rdf[:, 1] * dr / ((_rdf[:, 0] + 0.5 * dr)**3 - (_rdf[:, 0] - 0.5 * dr)**3)
#                _rdf[:, 0] = unit_to_user('length', unit, _rdf[:, 0])
#                
#                # Writing the results into files
#                np.savetxt(fn_out_rdf, _rdf)

        else:
            # Get the kernel density smoothed RDF by summing the contrbutions.
            rdfG[:, 1] = np.sum(forma * rdf.T[1].astype(int), axis=1)

            # Some constants
            if(A==B): 
                const = 2./(natomsA*(natomsB-1))
            else:
                const = 1./(natomsA*natomsB)
            const /= float(ifr - skipSteps)

         #   # Normalization
            _rdf = np.copy(rdfG)
            _rdf[:, 1] *= const / nbeads
            
            # Creating RDF from N(r)
            const, dr = cellVolume / (4 * np.pi / 3.0), _rdf[1, 0] - _rdf[0, 0]
            _rdf[:, 1] = const * _rdf[:, 1] * dr / ((_rdf[:, 0] + 0.5 * dr)**3 - (_rdf[:, 0] - 0.5 * dr)**3)
            _rdf[:, 0] = unit_to_user('length', unit, _rdf[:, 0])
            
            # Writing the results into files
            np.savetxt(fn_out_rdf, _rdf)

def main(*arg):

    RDF(*arg)

if __name__ == '__main__':

    main(*sys.argv[1:])
