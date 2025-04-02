#!/usr/bin/python3
#Atom Class
# Author: Angel-Emilio Villegas Sanchez

from .Molecule import Molecule
from FODLego.graphing import *
import sys
import logging
logging.basicConfig(format="%(levelname)s:%(filename)s:%(funcName)s(): %(message)s", level=logging.DEBUG)

def main():
    if len(sys.argv) == 1:
        logging.warning("No arguments were given, please provide XYZ file name")
        exit(1)
    elif len(sys.argv) == 2:
        logging.info("One argument passed. Creating FOD Prediction.")
        mol = Molecule(sys.argv[1])
        mol.CreateCLUSTER()
        mol.CreateFRMORB()
        mol.CreateXYZ()
    elif len(sys.argv) == 3:
        if sys.argv[1] == "list":
            assert len(sys.argv) == 3, "You did not provide a list of files to analyze."
            print("You provided the 'list' flag. Your input file is expected to have several filenames for comparison.")
            mols = CreateMolecules(sys.argv[2])
            MolecularSet(mols).saveset('duyen.pickle')
            #EdgeDist_FFODs(mols)
            #Angles_Hist(mols)
            #tables(mols)
            # Histogram_Radii(mols)
            #Histogram_Deviation(mols)
            #ffod_radii()
        elif sys.argv[1] == "open":
            logging.info("Creating an open-shell calculation with alpha/beta FODs")
            mol = Molecule(sys.argv[2], openshell=True)
            mol.CreateCLUSTER()
            mol.CreateFRMORB()
            mol.CreateXYZ()
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
            #mol.GeFBEdges()
main()
