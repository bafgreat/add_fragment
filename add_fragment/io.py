#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke and Dr. Shabhajit Das"
__status__ = "production"
from ase.io import read, write


def load_system(filename):
    '''
    A function to load the system from a file. The function reads the system from a
    file, formats it into a suitable ASE atoms object, and returns it.

    parameters:
    -----------
    filename : str : name of the file from which the system is to be loaded
    returns:
    ase.Atoms : ASE atoms object representing the loaded system
    '''
    # Load the system from the file
    ase_atoms = read(filename)

    return ase_atoms