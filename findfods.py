#!/usr/bin/python3
#Atom Class
# Author: Angel-Emilio Villegas Sanchez

from  globaldata import GlobalData
from Molecule import *
from Analysis import *
from graphing import *
import sys

dat = GlobalData()
if len(sys.argv) == 1:
    print("No arguments were given, please provide XYZ file name")
    exit(1)
elif len(sys.argv) == 2:
    print("One argument passed. Creating FOD Prediction.")
    mol = Molecule(sys.argv[1])
    mol.CreateCLUSTER()
    mol.CreateFRMORB()
    mol.CreateXYZ()
elif len(sys.argv) == 3:
    if sys.argv[1] == "list":
        assert len(sys.argv) == 3, "You did not provide a list of files to analyze."
        print("You provided the 'list' flag. Your input file is expected to have several filenames for comparison.")
        mols = CreateMolecules(sys.argv[2])
        for m in mols:
            m.GeFBEdges()
        EdgeDist_FFODs(mols)
        Angles_Hist(mols)

    elif sys.argv[1] == "check":
        mol = Molecule(sys.argv[2])
        mol._CheckChemValency()

    elif sys.argv[1] == "createsdf":
        mol = Molecule(sys.argv[2])
        mol.CreateSDF()

    else:
        print("Two arguments passed. Reverse Determination of Relaxed FODs.")
        mol = Molecule(sys.argv[1], sys.argv[2])
        mol.CreateCompXYZ()
        mol.GeFBEdges()
        EdgeDist_FFODs([mol])
        #mol._debug_printBFODsXYZ()
        #mol._debug_printAtoms()
