# RDKit for BondPrediction
from rdkit import Chem
from rdkit.Chem import Mol
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import rdmolops
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import GetPeriodicTable
from rdkit import Geometry
# Numpy
from scipy.spatial import distance
# Others
from os import remove
import logging
logger = logging.getLogger(__name__)

# Custom Made library
from FODLego.globaldata import GlobalData
from FODLego.ElementaryClasses import *
from FODLego.Bond import *
from FODLego.FOD import FOD
from FODLego.BFOD import *


class Molecule:
    def __init__(self, source, RelaxedFODs = None, OgXYZ = None, openshell = False) -> None:
        # Molecular Parameters
        self.mAtoms: List[Atom] = []
        self.mComment = ''
        self.mQ = 0
        self.mBonds: List[Bond] = []
        self.mBFODs = set()
        self.mFFODs = set()
        self.mCFODs = set()
        self.mFODs = []
        self.mValidStruct = True
        self.mOpen = openshell

        # Associated files
        self.mSrc: str = source
        self.mOgXYZ = OgXYZ
        self.mTargetFile = RelaxedFODs

        # Load File
        self._LoadSource()
        self.__RD_Bonds()
        self.CheckStericity()
        self.CalculateFODs()
        self._InterFOD_Dist()
        self._CheckChemValency()

        # Reverse Determine Parameters with a target file
        if RelaxedFODs != None:
            self.mRelaxPos = []
            self.__LoadTargetFODs()
            self.ReverseDetermination()
            if OgXYZ != None:
                self.CreateCompXYZ()

    #Public Methods
    def CalculateFODs(self):
        """
        Loop through the atoms and call respective methods to 
        calculate the FOD shells
        """
        # TODO: Find a way to just add them deterministically, since sets mess order
        for atom in self.mAtoms:
            atom.mFODStruct.PrepareShells(self.mAtoms)
            # Add the calculated FODs to the molecule
            for bfod in atom.mFODStruct.mBFODs:
                self.mBFODs.add(bfod)
            for ffod in atom.mFODStruct.mFFODs:
                self.mFFODs.add(ffod)
            for cfod in atom.mFODStruct.mCore:
                self.mCFODs.add(cfod)
        self.mFODs = set.union(self.mBFODs, self.mCFODs, self.mFFODs)

    def _InterFOD_Dist(self):
        # Create a list with all FODs. Make into List for index use.
        allFODs = set.union(self.mBFODs, self.mFFODs)
        allFODs = list(allFODs)
        n = len(allFODs)
        # Double loop without double counting
        for i in range(n):
            fod1 = allFODs[i]
            for j in range(i+1,n):
                fod2 = allFODs[j]
                if dist(fod1.mPos, fod2.mPos) < 0.4:
                    print(f'Valence FOD at {fod1.mPos} is very close to FOD at {fod2.mPos}')

    def CreateSDF(self) -> None:
        """
        Create a SDF file!
        """ 
        writer = Chem.SDWriter(self.mSrc + '.sdf')
        writer.write(self.rdmol)

    def CreateXYZ(self) -> None:
        """
        Create an XYZ file with
        """
        with open("lego.xyz",'w') as output:
            #First 2 lines
            output.write(str(len(self.mAtoms) + len(self.mFODs)) + '\n')
            output.write(self.mComment)
             
            # Write all atoms
            for atom in self.mAtoms:
                atom_coords = ' '.join([f"{x:7.4f}" for x in atom.mPos])
                output.write(f"{atom.mName} {atom_coords}\n")
             
            # Initialize separate strings for each channel
            up_fods = []
            down_fods = []

            for fod in self.mFODs:
                fod_coords = ' '.join([f"{x:7.4f}" for x in fod.mPos])
                if fod.mChannel is True:
                    up_fods.append(f"X {fod_coords}\n")
                else:
                    down_fods.append(f"He {fod_coords}\n")

            # Combine all lines with True channels first, then False channels
            output.write(''.join(up_fods + down_fods))
    
    def CreateCLUSTER(self) -> None:
        """
        Creates a CLUSTER file that will serve as an input file for FLOSIC to begin
        """
        cluster = open("CLUSTER", "w")
        # CLUSTER Preamble
        cluster.write("LDA-PW91*LDA-PW91\n")
        cluster.write("NONE\n")
        cluster.write(f"{len(self.mAtoms)}\n")

        # Loop thorugh each atom for coordinates
        for atom in self.mAtoms:
            coordinate = " ".join( f"{x * GlobalData.ANG2AU:10.5f}" for x in atom.mPos) + " " + str(atom.mZ) + " ALL" + '\n'
            cluster.write(coordinate)
        cluster.write(f"{self.mQ} {0.0}")  # TODO: Make a variable that contains sum of all spins
        cluster.close()
    
    def CreateFRMORB(self) -> None:
        cluster = open("FRMORB", "w")
        # CLUSTER Preamble
        cluster.write(f"{len(GlobalData.mFODs)} 0\n")
        # Loop thorugh each atom for coordinates
        for fod in GlobalData.mFODs:
            coordinate = " ".join( f"{x * GlobalData.ANG2AU:7.4f}" for x in fod.mPos) + '\n'
            cluster.write(coordinate)
        cluster.close()

    def ClosedMol(self) -> bool:
        """
        Checks that all atoms in the molecule have a closed shell
        """
        for atom in self.mAtoms:
            if atom._CheckFullShell() == False:
                return False
        return True

    def CreateCompXYZ(self) -> None:
        """
        A file, with the postfix "legocomp.xyz", with two sets of FODs is written.
        FODLego predictions are written as 'X' atoms.
        The relaxed FODs are written as 'He' atoms.
        """
        file = self.mComment
        prefix = file[:-4]
        with open(f"{prefix}_legocomp.xyz",'w') as output:
            #First 2 lines
            output.write(str(len(self.mAtoms) + 2*len(self.mFODs)) + '\n')
            output.write(self.mComment)
            #Write all atoms
            for atom in self.mAtoms:
                output.write(' '.join([atom.mName,*[str(x) for x in atom.mPos]]) + '\n')

            #Write all FODs
            fods = set.union(self.mBFODs, self.mCFODs, self.mFFODs)
            for bfod in fods:
                xyz = " ".join([str(x) for x in bfod.mPos])   
                output.write(f"X {xyz}\n")

            # Write the Relaxed FODs read from comparison
            for relaxed in self.mRelaxPos:
                pos = " ".join([str(x) for x in relaxed])   
                output.write(f"He {pos}\n")

    def ReverseDetermination(self) -> None:
        """
        This function executes the reversedetermination of paramters for all Target FODs that have
        been associated with the predicted FODs.
        """
        from FODLego.Shells import FODShell
        self.__AssociateTargets()
        
        #Loop through atoms to deterpmine the average and the variance
        for atom in self.mAtoms:
            for shell in atom.mFODStruct.mCoreShells:
                # For a tetra, get average radii
                if shell.mShape == 'point':
                    core_dist = norm(atom.mPos - shell.mfods[0].mPos)
                elif shell.mShape == 'tetra':
                    targets = np.vstack([x.mAssocFOD.mPos for x in shell.mfods])
                    target_R = distance.cdist([atom.mPos], targets, 'euclidean')
                    # Pred. Stats. of Target Shell
                    shell.mTarget_u_R = np.mean(target_R)
                    shell.mTarget_s2_R = np.var(target_R)
                    # Pred. Stats. of Predicted Shell
                    pred_radii = [x.mR for x in shell.mfods]
                    shell.mPred_u_R = np.mean(pred_radii)
                    shell.mTarget_s2_R = np.var(pred_radii)

    # Getters
    def GeFBEdges(self):
        for at in self.mAtoms:
            edges = at.GetAssocEdges_B_F_FOD()
            if len(edges) > 0:
                logger.info(edges)
            else:
                logger.info("No edges found")

    def get_ffod(self, type):
        return [f for f in self.mFFODs if isinstance(f,type)]

    def get_bfod(self, type):
        return [f for f in self.mBFODs if isinstance(f,type)]

    # Logical Expressions
    def has_fodtype(self, type) -> bool:
        """
        Checks if the molecule contains a specific type of FOD.
        """
        has = False
        for fod in self.mFFODs:
            if isinstance(fod,type):
                print(self.mSrc)
                has = True
        return has

    #Private Methods
    def __AssociateTargets(self):
        """
        This function associates the target FOD that is nearest to the predicted FOD (i.e. the output of this program). 
        First it checks for CFODs and removes them from the list (assuming we have the most accuracy on them), and then we start finding the FFODs and BFODs.
        Note: This is a bit messy.
        """
        #TODO: Create a function that loops over the things, instead of doing 3 for loops....
        from FODLego.FFOD import FFOD

        # Create a vector of the distances, vertical
        rlx = np.vstack([x for x in self.mRelaxPos])
        for fod in self.mCFODs:
            # Get the minimum distance to Target FOD
            distances = distance.cdist([fod.mPos],rlx, 'sqeuclidean')
            index = np.argmin(distances[0])

            # Print general information
            # Create appropriate associate fod
            if isinstance(fod, CFOD):
                fod.mAssocFOD = CFOD(fod.mAtom, rlx[index])
                # Exclude the FOD that has been associated from the relaxed list.
                rlx = np.delete(rlx,index,axis=0)
            else:
                print("Invalid classification for associated FOD")

        # BFODs
        for fod in self.mBFODs:
            # Get the minimum distance to Target FOD
            distances = distance.cdist([fod.mPos],rlx, 'sqeuclidean')
            index = np.argmin(distances[0])

            # Print general information
            # Create appropriate associate fod
            if isinstance(fod, BFOD):
                fod.mAssocFOD = BFOD(
                    fod.mBold,
                    fod.mMeek,
                    rlx[index])
                rlx = np.delete(rlx,index,axis=0)

        # FFODs
        for fod in self.mFFODs:
            # Get the minimum distance to Target FOD
            distances = distance.cdist([fod.mPos],rlx,'sqeuclidean')
            index = np.argmin(distances[0])

            if isinstance(fod, FFOD):
                fod.mAssocFOD = FFOD(
                    fod.mAtom,
                    target=rlx[index])
                rlx = np.delete(rlx,index,axis=0)

    def __CreateCoordMap(self, file):
        """
        Create a map for the Embedding to work....
        """
        XYZ = open(file, "r")
        count = int(XYZ.readline()) #Read Size
        self.mComment = XYZ.readline() #Read Size
        map = {}
        # Fill atom information 
        for i in range(count):
            coor = XYZ.readline().split()
            atom_xyz = [float(x) for x in coor[1:4]]
            point = Geometry.Point3D(*atom_xyz)
            map[i] = point
        XYZ.close()
        return map

    def __LoadSMILES(self, tmp=None):
        """
        Creates a 3D conformer of the SMILES structure according to the "working with 3D Molecules" section of the RDKit documentation.
        """
        # Seed a random number
        from random import randint
        seed = randint(0,2000)

        # Prepare SMILES Molecule with rdkit
        if tmp != None:
            core = Chem.MolFromXYZFile(tmp)
            try:
                self.rdmol = AllChem.ConstrainedEmbed(self.rdmol, core, randomseed=seed, maxAttempts=8000)
                AllChem.MMFFOptimizeMolecule(self.rdmol)

                # Load onto FODLego scheme
                for i,atom in enumerate(self.rdmol.GetAtoms()):
                    coor = self.rdmol.GetConformer().GetAtomPosition(i)
                    name = atom.GetSymbol()
                    self.mAtoms.append(Atom(i, name,coor, self)) #Name and Position
            except Exception as e:
                self.mValidStruct = False
                print(f'File with {tmp} not cannot be embeded')
                print(e.args[0])
        else:
            AllChem.EmbedMolecule(self.rdmol, maxAttempts=8000, randomSeed=seed)
            AllChem.MMFFOptimizeMolecule(self.rdmol)

            for i, at in enumerate(self.rdmol.GetAtoms()):
                name = at.GetSymbol()
                pos = self.rdmol.GetConformer().GetAtomPosition(i)
                self.mAtoms.append(Atom(i, name,[pos.x, pos.y, pos.z], self))

    def __LoadTargetFODs(self):
        """
        Load a file of FOD position to reversedetermine their parameters.
        """
        # Read size of UP/DOWN FODs
        if isinstance(self.mTargetFile, str):
            TargetF = open(self.mTargetFile, "r")
            relx = TargetF.readline().split()
            upcount = int(relx[0])
            downcount = int(relx[1])

            # Load positions as ndarrays
            for i in range(upcount):
                coor = TargetF.readline()
                coor = coor.replace('D','E')
                coor = coor.split()
                if self.mTargetFile[-6:] == 'FRMORB':
                    atom_xyz = np.array([float(x)*(GlobalData.AU2ANG) for x in coor[0:3]])
                elif self.mTargetFile[-6:] == 'Target':
                    atom_xyz = np.array([float(x) for x in coor[0:3]])
                else:
                    atom_xyz = np.array([float(x) for x in coor[0:3]])
                self.mRelaxPos.append(atom_xyz)  # Name and Position
            TargetF.close()
        elif isinstance(self.mTargetFile, list):
            print("A list was passed")
            pass

    def __LoadPDB(self):
        """
        Credit to Betelgeuse in stackoverflow describing how to get this:
        https://stackoverflow.com/questions/71915443/rdkit-coordinates-for-atoms-in-a-molecule
        """
        for i, atom in enumerate(self.rdmol.GetAtoms()):
            position = self.rdmol.GetConformer().GetAtomPosition(i)
            pos = [position.x, position.y, position.z]
            self.mAtoms.append(Atom(i, atom.GetSymbol(), pos, self))

    def __RD_Bonds(self):
        """
        Calculate Bonds using the RDKit library.
        This will be used for prototyping  
        """ 
        rdmolops.Kekulize(self.rdmol)
        #Loop through the
        for rdbond in self.rdmol.GetBonds():
            # Get atoms
            i1 = rdbond.GetBeginAtomIdx()   
            i2 = rdbond.GetEndAtomIdx() 
            Atom1 = self.mAtoms[i1]
            Atom2 = self.mAtoms[i2]

            # Set order and create bonds
            order = rdbond.GetBondTypeAsDouble()
            newbond = Bond(Atom1, Atom2, order)
            self.mBonds.append(newbond)
            Atom1.AddBond(Atom2, order)
            Atom2.AddBond(Atom1, order)
        
        # Find out valence of atoms after connectivity
        for at in self.mAtoms:
            at.mValCount = at.FindValence()

    def _CheckChemValency(self) -> None:
        """
        Check that the valency of the atoms is no more than four.
        This assumes max sp3 hybridization of orbitals.
        """
        for atom in self.rdmol.GetAtoms():
            if atom.GetTotalValence() > 4 and self.mValidStruct == True:
                logger.warning(f"Non simple bonding at {atom.GetSymbol()} in {self.mSrc}. Valency above 4.")
                self.mValidStruct = False

    def _CheckRadicals(self):
        if Chem.Descriptors.NumRadicalElectrons(self.rdmol) > 0:
            self.mValidStruct = False
 
    def CheckStericity(self):
        """
        Determine Steric number of all atoms. Currently assumes that the atoms are closed-shell
        TODO: Implement some way to easily add an open-shell calculation, in which there might be 
        open-shells
        """
        for atom in self.mAtoms:
            atom.CalcSteric()
        
    def CountFODs(self):
        """
        Returns the amout of FODs in the molecule
        """
        count = 0
        for atom in self.mAtoms:
            count += len(atom.mFODStruct.mfods)
        return count
    
    def _LoadSource(self):
        """
        Loads the molecule using RDKit. Can take different forms of input to create the molecule.
        Args:
            src: The
            ogXYZ:
        """
        # Parameters
        rdDistGeom.EmbedParameters.enableSequentialRandomSeeds=False
        rdDistGeom.EmbedParameters.maxIterations=1
        AllChem.EmbedParameters.clearConfs = True

        def CLUST2XYZ(input_file, output_file):
            """
            Load a CLUSTER File into an rdmol
            """
            # Function to process the file and convert to .xyz format
            with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
                lines = infile.readlines()

                # Skip the first two lines and get the atom count from the third line
                atom_count = int(lines[2].split()[0])
                outfile.write(f"{atom_count}\n\n")  # XYZ format first line (atom count) and an empty comment line

                # Process the remaining lines
                for line in lines[3:3+atom_count]:
                    parts = line.split()
                    x, y, z = parts[:3]  # First three are the coordinates
                    x = float(x)*GlobalData.AU2ANG
                    y = float(y)*GlobalData.AU2ANG
                    z = float(z)*GlobalData.AU2ANG

                    atomZ = int(parts[3])  # Atomic number at index 3
                    element = GetPeriodicTable().GetElementSymbol(atomZ)

                    # Write the element and coordinates in .xyz format
                    outfile.write(f"{element} {x:10.5f} {y:10.5f} {z:10.6f}\n")

        def read_xyz():
            for i, atom in enumerate(self.rdmol.GetAtoms()):
                position = self.rdmol.GetConformer().GetAtomPosition(i)
                pos = [position.x, position.y, position.z]
                self.mAtoms.append(Atom(i, atom.GetSymbol(), pos, self))

        def create_xyz(file: str) -> None:
            self.mComment = self.mSrc + '\n'
            self.rdmol = Chem.MolFromXYZFile(file)
            read_xyz()
            rdDetermineBonds.DetermineConnectivity(self.rdmol)
            if len(self.rdmol.GetBonds()) > 0:
                rdDetermineBonds.DetermineBondOrders(self.rdmol, charge=self.mQ)

        if self.mOgXYZ == None:
            if self.mSrc[-3:] == "pdb":  # WiP. This is a test
                self.rdmol = Chem.MolFromPDBFile(self.mSrc)
                self.__LoadPDB()

            elif self.mSrc[-3:] == "xyz":
                create_xyz(self.mSrc)

            elif self.mSrc[-7:] == "CLUSTER":
                # Turn CLUSTER into an .xyz file  
                CLUST2XYZ(self.mSrc, "CLUST.xyz")
                create_xyz("CLUST.xyz")
                remove("CLUST.xyz")

            else:  # SMILES
                print('THis is SMILES')
                m = Chem.MolFromSmiles(self.mSrc)
                self.rdmol = Chem.AddHs(m)
                self.mComment = self.mSrc + '\n'
                self.__LoadSMILES()

        else:
                self.rdmol = Chem.MolFromSmiles(self.mSrc)
                self.rdmol = Chem.AddHs(self.rdmol)
                self.__LoadSMILES(self.mOgXYZ)

        # Read RDKit documentation. They said this was in experimentation. It might help later down the line.
        Chem.rdMolDescriptors.CalcOxidationNumbers(self.rdmol)

    #Debugging Methods
    def debug_printTargPred():
        c = np.vstack([x for x in self.mRelaxPos])
        for pfod in GlobalData.mFODs:
            # Get the minimum distance to Target FOD
            distances = distance.cdist([pfod.mPos],c, 'sqeuclidean')
            index = np.argmin(distances[0])
            if __debug__:
                print('-'*50)
                print(f"Index: {index}")
                print(f'Predicted: {pfod.mPos}')
                print(f'{type(pfod)}')
                print(f'Target: {c[index]}')
                print(f'Distance: {sqrt(distances[0][index])}')

    def _debug_printAtoms(self):
        """Print atom names and positions"""
        for atom in self.mAtoms:  
            print("---------------------------------")
            print(atom.mName, "at", atom.mPos)
            at = self.rdmol.GetAtomWithIdx(atom.mI)
            print(f'RDValency: {at.GetTotalValence()}')
            print(f'RDImplicitVal: {at.GetImplicitValence()}')
            print(f'RDExplicitVal: {at.GetExplicitValence()}')
            print(f'RDFormalQ: {at.GetFormalCharge()}')
            print(f'RDNeighbors: {[x.GetSymbol() for x in at.GetNeighbors()]}')
            print(f'RDNeighbors: {at.GetHybridization()}')
            print(f'Valency: {atom.mValCount}')
            print(f'Steric Number: {atom.mSteric}')
            print(f'Free Pairs: {atom.mFreePairs}')
            print("Shell (Core) Structure:",*atom.mFODStruct.mCore)
            print(f'RD Bonds: {at.GetBonds()}')
            print(f'RD Total Degree of atom: {at.GetTotalDegree()}')
            print(f'RD Oxidation: {at.GetProp("OxidationNumber")}')
            
            print('BondedAtoms: ')
            for b in atom.mBonds:
                bonded = b.mAtoms[1].mName
                print(f'-- Bonded to {bonded}({b.mAtoms[1].mI}). Order {b.mOrder}')
            
            closedshell = atom._CheckFullShell()
            print (f'Shell Full: {closedshell}')
            if (closedshell == False): print ("###NONCLOSED SHELL SYSTEM###")

    def _debug_printBondMatrix(self):
        print("##BOND MATRIX##")
        str4 = [[0] * len(self.mAtoms) for _ in range(len(self.mAtoms))]
        for b in self.mBonds:
            str4[b.mAtoms[0]][b.mAtoms[1]] = b.mOrder
        for atom in str4:
            print(atom)    

    def _debug_printBFODs(self):
        from FODLego.ElementaryClasses import Atom
        for atom in GlobalData.mAtoms:
            print(f'In atom {atom.mI}:')
            for bfod in atom.mFODStruct.mBFODs:
                print(bfod)

    def _debug_printBFODsXYZ(self):
        with open(f"lego.xyz",'w') as output:
            #First 2 lines
            output.write(str(len(self.mAtoms) + len(GlobalData.mFODs)) + '\n')
            output.write(self.mComment)
            #Write all atoms
            for atom in self.mAtoms:
                output.write(' '.join([atom.mName,*[str(x) for x in atom.mPos]]) + '\n')
            #Write all FODs
            for bfod in GlobalData.mFODs:
                xyz = " ".join([str(x) for x in bfod.mPos])   
                output.write(f"X {xyz}\n")

    def _debug_CompareTargetFODs(self):
        """
        This function enumerates the number of predicted FODs (generated by FODLego)
        and the number of FODs from your target file.
        """
        print("-"*30)
        print(f'You have {len(GlobalData.mFODs)} Predicted FODs')
        print(f'You have {len(self.mRelaxPos)} Target FODs')

    #String Output
    def __str__(self) -> str:
        print(self.mTargetFile)
        str2 = f"Atoms: {len(self.mAtoms)}\n"
        str3 = f"Bonds: {len(self.mBonds)}\n"
        str4 = f"FODs: {len(self.mFODs)}\n"
        return str2+str3+str4
