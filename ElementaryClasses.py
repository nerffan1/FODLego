#Description: This file contains as set of classes that will implement the FOD Heuristic Solution 
#  following the paradigm of Object-Oriented Programming (OOP). Encapsulation for FOD_Structure, FODs,
#  FOD_Orbital, among other things will be included in this file.
# Roadmap: Use polymorphism for Closed Hybrid Shells (e.g. sp3).
#  -  Add Tetrahedra class and several attributes/methods to manipulate them
#  - In far future, somehow implement the triaugmented triangular prism that corresponds to sp3d5 ( 9 FODs, 18 electrons) 
#Author: Angel-Emilio Villegas S.
from ast import Global
from  globaldata import GlobalData
import math3d
import numpy as np
from numpy import sqrt 
from typing import List
from scipy.spatial.transform import Rotation as R
import csv

################# FOD STRUCTURE #################

class Bond:
    def __init__(self,start,end,order):
        self.mAtoms = (start,end)
        self.mOrder = order
    def __str__(self) -> str:
        return f"From {self.mAtoms[0]} to {self.mAtoms[1]}. Order: {self.mOrder}"
    
class Atom:
    def __init__(self, index: int, Name: str, Pos):
        #Known Attributes
        self.mName = Name
        self.mPos = np.array(Pos) 
        self.mI = index
        self.mZ = GlobalData.GetElementAtt(self.mName, "AtomicNumber")
        self.mPeriod = int(GlobalData.GetZAtt(self.mZ, "Period" ))
        self.mGroup = int(GlobalData.GetZAtt(self.mZ, "Group" ))
        self.mValCount = self._FindValence()
        #Undetermined Attributes
        self.mSteric = 0
        self.mFreePairs = 0
        self.mCharge = 0 #Temporarily zero for all atoms. Ionic atoms can be dealt with in the future
        self.mBonds = []
        self.mFODStruct = FODStructure(self)
        self.mCompleteVal = False
        
    def AddBond(self, atom2: int, order: int):
        self.mBonds.append(Bond(self.mI, atom2,order))

    def CalcSteric_test(self) -> None:
        """
        TODO: Need to account for systems where the valence electrons + bonding FODs
        """
        self.mFreePairs = int((self.mSteric - np.sum([2*x.mOrder for x in self.mBonds ]))/2)
        self.mSteric = self.mFreePairs + len(self.mBonds)
    
    def CalcSteric(self) -> None:
        """
        This function assumes that the system at hand only contains Closed Shell calculations.
        TODO: Need to account for systems where the valence electrons + bonding FODs
        """
        #Electrons involved in the Bond
        bondelec = np.sum([2*bond.mOrder for bond in self.mBonds])
        # The difference between the total electrons and the number of electrons that fill the shell
        self.mCharge = (self.mZ + bondelec) - GlobalData.GetFullElecCount(self.mGroup, self.mPeriod)
        self.mFreePairs = int(GlobalData.mShellCount[self.mPeriod] - bondelec)/2
        self.mSteric = self.mFreePairs + len(self.mBonds)
    
    def _FindValence(self):
        """
        This method finds the number of electrons in the valence shell by finding the 
        difference between the current Group and the last ClosedShell Group. Only works up
        to 5th period.
        TODO: This can be saved in GlobalData
        """
        if self.mGroup < 4:
            return self.mGroup
        else:
            if self.mPeriod < 4:
                return (2 + (self.mGroup - 12))
            elif self.mPeriod < 6:
                return (self.mGroup)
                    
    def _CheckFullShell(self):
        """
        Check that the atom has a closed shell.
        Future: Add a variable that saves the info so that looping every time
        this is called (if called more than once) is unnecesary
        TODO: Alternatively, just check the amount of electrons
        """
        #Determine How many electrons are needed to fill shell
        for ClGrp in GlobalData.mClosedGroups:
            if self.mGroup < ClGrp:
                checkshell = ClGrp - self.mGroup + self.mCharge
                break 
        for bond in self.mBonds:
            checkshell -= bond.mOrder
        if checkshell == 0:
            return True
        else:
            return False
    
    # Additional Functions
    def __str__(self):
        pass    

######################## FOD Structure Class  ########################
class FODStructure:
    def __init__(self, parent: Atom):
        self.mAtom = parent
        self.mCore = [] #A list of FODShells
        self.mValence = [] #A list of FOfDs
        self.mfods = [] #All Finalized FODs
        self.mLastBond = 0
    
    # Setter Functions
    def AddValFOD(self, fods: List[np.ndarray], at2=False,finalize=False):
        """
        This function adds the FODs to the valence of the current atom's FOD Structure.
        It also accepts a secondary atom, at2, in order to add the FOD to its valence. The
        Finalize parameter is True only for the first atom that is writing the FODs; by finalizing
        we mean that the FOD is added to the list of overall FODs.
        """
        # Add to current valence
        for fod in fods:
            #Assertions
            print(fod)
            assert isinstance(fod, np.ndarray), "The variable is not a NumPy ndarray."

            self.mValence.append(fod)
            
            #Add to bonded atom valence, if passed
            if at2 != False:
                at2.mFODStruct.AddValFOD([fod])

            #Add to the final list of atoms, without duplicating
            if finalize:
                if len(self.mfods) == 0:
                    self.mfods = fod  
                else:
                    self.mfods = np.vstack((self.mfods,fod))

    def PrepareShells(self, atoms: List[Atom]):
        """
        This function will determine the creation of Core FODs, those that are not 
        related to bonding. 
        Comment: The scheme is easy in the first 3 periods of the Periodic
        Table, but it will become trickier ahead if Hybridization heuristics don't work.
        Currently it only works for closed shell calculations (V 0.1).
        TODO: This will assume that we are doing up to the n=3 shell, with sp3 hybridization
        TODO: Take into account free pairs of electrons
        TODO: For mOrder=2, there are many schemes to determine direction
        TODO: Find an elegant solution to do exceptions for H bonds
        TODO: Account previous bonds formed, by talling previous FODs, or looking back at mBonds 
        TODO: GlobalData.GetFullElecCount() Can be precalculated ahead of time and placed as a member vatrable
        """
        at1 = self.mAtom
        
        def SingleBond(at2: Atom):
            dom,weak,dir = DominantAtom(at1,at2,True)
            bfod = AxialPoint_Simple(dom, weak, dir)
            self.AddValFOD([dom.mPos + bfod], at2, True)
        
        def DoubleBond(at2: Atom, bond: Bond):
            """ 
            Create the FODs representing the double bond. Currently the FOD filling is unidirectional (sequential)
            and  does not account for the next atom in the iteration to see if we can further accomodate the bonding FODs 
            """
            def D_Bond_Height(dom: Atom, axisproj: float):
                """
                Return radial distance from interatomic (bonding) axis. 
                For heterogenous atoms....
                """
                # Radius Logic
                elecs = GlobalData.GetFullElecCount(dom.mGroup, dom.mPeriod)
                rad = GlobalData.mRadii[elecs][dom.mZ]
                return sqrt(rad**2 - ((axisproj)**2))
                            
            def CrossPerpDir():
                vector_for_cross = []
                for otherb in self.mAtom.mBonds:
                    if otherb != bond:
                        vector_for_cross.append(atoms[otherb.mAtoms[1]].mPos)
                vector_for_cross -= self.mAtom.mPos
                vec = np.cross(*vector_for_cross)
                return vec/np.linalg.norm(vec)
            
            def D_FFOD_Direction():
                if self.mAtom.mFreePairs == 0:
                    #Determine DFFOD Direction
                    if len(self.mAtom.mBonds) == 3: # Planar
                        return  CrossPerpDir()
                    elif len(self.mAtom.mBonds) == 2:
                        if self.mAtom.mBonds.index(bond) == 0:
                            return  RandomPerpDir(dir)
                        else:
                            # Cross product between atom and already-placed FODs
                            vector_for_cross = []
                            for fods in self.mValence:
                                vector_for_cross.append(fods)
                            vector_for_cross -= self.mAtom.mPos
                            return np.cross(*vector_for_cross)   
                else:
                    return RandomPerpDir(dir)

            #Information
            axis2fod = np.ndarray(3)
            dom,sub,fugal = DominantAtom(at1,at2, True)

            if GlobalData.GetFullElecCount(self.mAtom.mGroup,self.mAtom.mPeriod) <= 18:
                #Find perpendicular unit vector
                if self.mAtom.mFreePairs == 0:
                    axis2fod = D_FFOD_Direction()
                  
                #Add both FODs of the Double Bond
                midpoint = AxialPoint_Simple(dom,sub,fugal)
                # Determine Vertical Projection of DFFOD
                axis2fod *= D_Bond_Height(dom, np.linalg.norm(midpoint))

                #Add FODs
                self.AddValFOD([dom.mPos + midpoint + axis2fod],at2,True)
                self.AddValFOD([dom.mPos + midpoint - axis2fod],at2,True)
    
        
        def AddFreeElectron(free: int):
            """
            TODO: Change conditionals as to create more concise code 
            TODO: Create a series of variables for chosen constants, NO magic numbers
            """
            #Variables
            vector_for_cross = []
            fulle = GlobalData.GetFullElecCount(self.mAtom.mGroup, self.mAtom.mPeriod)
            vector_for_cross = self.mAtom.mPos - [fod for fod in self.mValence]
                
            if len(self.mAtom.mBonds) == 1:
                #Useful Information
                at2 = atoms[self.mAtom.mBonds[0].mAtoms[1]]
                free_dir = self.mAtom.mPos - at2.mPos
                free_dir /= np.linalg.norm(free_dir)

            if free == 1:
                # For a single atom
                if len(self.mAtom.mBonds) == 3 :
                    dr = vector_for_cross.sum(0)
                elif len(self.mAtom.mBonds) == 2 :
                    free_dir = AddNormals(vector_for_cross)
                    dr = free_dir*GlobalData.mRadii[fulle][self.mAtom.mZ]
                elif len(self.mAtom.mBonds) == 1:
                    l = GlobalData.mVert[fulle][self.mAtom.mZ]
                    if (self.mAtom.mPeriod < 3):
                        dr = free_dir*l/sqrt(8) #Place at midsphere distance
                    else:
                        dr = free_dir*GlobalData.mRadii[fulle][self.mAtom.mZ]
                self.AddValFOD([self.mAtom.mPos+dr],False,True)

            elif free == 2:
                #Direction away from atom, on plane of other bonds, or neighboring atoms
                axis2fod = np.cross(*vector_for_cross)
                axis2fod /= np.linalg.norm(axis2fod)
                #Get angle to atoms of FFODs (ffod-atom-ffod)
                ab = np.dot(*vector_for_cross)
                ab_abs =  np.linalg.norm(vector_for_cross,axis=1)
                theta = np.deg2rad(220) - np.arccos(ab/(ab_abs[0]*ab_abs[1]))
                theta /= 2

                if len(self.mValence) == 2: #Might be different for FODs or for number of bonds
                    #Determine the free direction
                    dxy  = AddNormals(vector_for_cross)
                    if len(at1.mBonds) == 2:
                        dxy = at1.mPos - [atoms[bonds.mAtoms[1]].mPos for bonds in at1.mBonds]
                        dxy = dxy.sum(0)
                        dxy /= np.linalg.norm(dxy)

                    #TODO: Must determine logic to choose distance, based of bonding as well 
                    Rad = np.linalg.norm(vector_for_cross[0]) # Remember these are the atom-BFOD distances
                    dxy *= Rad*np.cos(theta)
                    axis2fod *= Rad*np.sin(theta)
                elif len(self.mValence) == 1:
                    #Add both FODs of the Double Bond
                    axis2fod *= GlobalData.mVert[fulle][self.mAtom.mZ]/2
                    dxy = self.mAtom.mPos + free_dir*np.linalg.norm(axis2fod)/np.tan(np.deg2rad(54.735))
                
                self.AddValFOD([self.mAtom.mPos + dxy + axis2fod],False,True)
                self.AddValFOD([self.mAtom.mPos + dxy - axis2fod],False,True)

            elif free == 3:
                if self.mAtom.mZ < 10:
                    #Create the starting FOD. Begin with the vertical component
                    R_f = sqrt(3/8)*np.linalg.norm(at1.mPos - at1.mFODStruct.mValence[0])
                elif at1.mPeriod > 2:                                                            
                    R_f = np.linalg.norm(at1.mPos - at1.mFODStruct.mValence[0])
                # Begin First FOD 
                axis2fod = RandomPerpDir(free_dir)
                axis2fod *= np.sin(np.deg2rad(70.5288))*R_f 
                #Create the horizontal component, parallel to the bonding axis
                horizontal = np.cos(np.deg2rad(70.5288))*R_f*free_dir
                #Add vectors, and rotate to create equilateral triangle
                dr = horizontal + axis2fod
                #Create rotations and rotated FODs
                rot_fods = RotatePoints(3, dr, free_dir)

                #Translate to the atom of interest, and add to Valence 
                self.AddValFOD(at1.mPos + rot_fods,False,True)

        def PlaceFODs_Triple(fugal: np.ndarray, a: float, rad: float, at2: Atom, c: float):
            """ This function create an FOD at a certain distance, based of a and rad,
            which are quantities assumed to dominate the interaction.
            """
            equilat_r =a/sqrt(3)
            if at1.mZ == at2.mZ:
                if at1.mPeriod ==2 & at2.mPeriod == 2:
                    theta = np.arccos((c/2)/rad)
                elif at1.mPeriod > 2 & at2.mPeriod > 2:
                    theta = np.arctan(rad/(c/2)) 
            else:
                theta = np.arcsin((a/2)/rad) 
            dr = fugal*rad*np.cos(theta) + RandomPerpDir(fugal)*rad*np.sin(theta)
            return RotatePoints(3, dr, fugal)

        def TripleBond(at2: Atom):
            """
            #TODO: Create a helper funtion for conditional statements
            """
            #Determine dominant atom and direction
            atom, fugal = DominantAtom(at1,at2,False)
            # Variables
            c = np.linalg.norm(fugal) # Distance
            fugal /= np.linalg.norm(fugal)  
            elecs = GlobalData.GetFullElecCount(atom.mGroup, atom.mPeriod)
            rad = GlobalData.mRadii[elecs][atom.mZ] #Radius
            a = GlobalData.mVert[elecs][atom.mZ] #Edge
            #Place FODs
            rot_fods = PlaceFODs_Triple(fugal, a, rad, at2, c)
            self.AddValFOD(atom.mPos + rot_fods,at2,True)

        def AddBFODs():
            for bond in self.mAtom.mBonds:
                bonded_at = atoms[bond.mAtoms[1]] 
                if bonded_at.mI  > self.mAtom.mI:
                    if bond.mOrder == 1:
                        SingleBond(bonded_at)
                    elif bond.mOrder == 2:
                        DoubleBond(bonded_at, bond)
                    elif bond.mOrder == 3:
                        TripleBond(bonded_at)
        
        def AddFFODs():
            if self.mAtom.mFreePairs == 2:
                if self.mAtom.mSteric >= 3:
                    AddFreeElectron(2)
            elif self.mAtom.mFreePairs == 1:
                if self.mAtom.mSteric >= 2:
                    AddFreeElectron(1)
            elif self.mAtom.mFreePairs == 3:
                AddFreeElectron(3)

        def AddCoreElectrons():
            #Count core electrons and
            core_elec = self.mAtom.mZ - self.mAtom.mValCount
            if core_elec != 0:
                for shell in GlobalData.mGeo_Ladder[core_elec]:
                    if shell == 'point':
                        self.mCore.append(self.Point(self.mAtom))
                    elif shell == 'tetra':
                        self.mCore.append(self.Tetrahedron(10,self.mAtom))
        
        #Prepare the valence shell first, since it will help determine the
        # orientation of the inner shells
        AddBFODs()
        AddFFODs()
        AddCoreElectrons()
             
    def FinalizeFODs(self):
        """
        Add all FODs in the FODStructure, to the atom, so that they can be
        printed when creating the XYZ output file
        """
        #Add the core FODs
        for shell in self.mCore:
            shell.mfods += self.mAtom.mPos   
            if len(self.mfods) == 0:
                self.mfods = [shell.mfods]
            else:
                self.mfods = np.vstack((self.mfods,shell.mfods))  ###HOW TO concatenate FODs, easily
    
    class FODShell:
        def __init__(self, shape, fods, owner: Atom):
            self.mOwner = owner
            self.mShape = shape 
            self.mfods = np.array(fods)
        def __str__(self):
            return self.mShape

    class Point(FODShell):
        def __init__(self, owner: Atom):
            super().__init__('point', [0.0,0.0,0.0], owner)

    class Tetrahedron(FODShell):
        """
        Tetrahedron Class: FOD Geometry corresponding to sp3 "hybridization' geometry.
        Roadmap: There will  be different functions to create compound transformations of FODs (e.g. the base, or peak
        of the tetrahedron), and to rotate them in the proper direction as well.
        """
        def __init__(self, core_amount, owner: Atom):
            super().__init__('tetra', GlobalData.mTetraGeo, owner)
            self.mfods *= GlobalData.mRadii[core_amount][self.mOwner.mZ]
        #Class Methods
        def CreateTetra(self):
            """
            This method arranges 4 FOD points in a tetrahedral form. Depending on the bonding, and free electrons,
            the direction vector will be different. 
            """
            pass

        def RotateTetra(self):
            pass
################# ADDITIONAL FUNCTIONS #################

def AddNormals(vectors: list[np.array]) -> np.ndarray:
    """This function normalizes all the given vectors, adds their normal,
    and normalizes the resultant vector.
    
    This function is motivated by the observation that FFODs are closer to the
    bisecting angle in double FFODs, rather than a weighted """
    free_dir_norm = [1/np.linalg.norm(x) for x in vectors]
    free_dir = vectors * np.reshape(free_dir_norm,(len(vectors),1))
    free_dir = free_dir.sum(0)
    free_dir /= np.linalg.norm(free_dir)
    return free_dir

def DominantAtom(at1: Atom, at2: Atom, all=True):
    """
    This Function determines the dominant atom in the bonding and its distance to the weaker atom.
    If the 
    at1: An atom
    at2: An atom bonding to at2
    all: Boolean to return dominant and weak atom. Default only return dominant atom.
    """
    if at1.mPeriod < at2.mPeriod:
        fugal = at2.mPos - at1.mPos
        dom = at1
        sub = at2
    elif at1.mPeriod > at2.mPeriod:
        fugal = at1.mPos - at2.mPos
        dom = at2
        sub = at1
    elif at1.mZ > at2.mZ:
            fugal = at2.mPos - at1.mPos
            dom = at1
            sub = at2
    elif at1.mZ <= at2.mZ:
            fugal = at1.mPos - at2.mPos
            dom = at2
            sub = at1
    # Either return dom and sub, or just the dominant atom. 
    if all:
        return dom, sub, fugal
    else:
        return dom, fugal

def AxialPoint_Simple(dom:Atom, sub:Atom, dir:np.ndarray) -> np.ndarray:
            """
            Return FOD location for a Single FOD representing a Single Bond.
            If the bonding atoms are the same type, then the midpoint is chosen.
            Atom1 is assumed to be bigger than Atom2.  
            TODO: Should just return a length since we can calculate dominant one
            """
            Z1 = dom.mZ
            Z2 = sub.mZ
            if dom.mZ == sub.mZ:
                return dir*0.5
            else:
                Z1 = 0.4 if Z1 == 1 else Z1
                Z2 = 0.4 if Z2 == 1 else Z2
                r = sqrt(Z1/Z2)
                g = r/(1+r)
            if g < 0.5:
                return dir*g
            else:
                return dir*(1-g)
            
def RotatePoints(n:int,fod0:np.ndarray,axis:np.ndarray) -> List[np.ndarray]:
    """
    This function creates n points in a circle, starting from your fod0 (i.e. the first fod in the circle).

    n: Number of points to equally distribute on a circle
    ogp: The Original Point, the first element in your circle
    axis: The axis of rotation
    """
    assert len(axis) == 3, "The array must have 3 dimensions"
    fodsRotated = [fod0]
    step = (2*np.pi)/n
    for i in range(1,n):
        rot = R.from_rotvec(((step*i*2*np.pi)/3)*axis)
        fod = np.matmul(rot.as_matrix(),fod0)
        fodsRotated.append(fod)
    return fodsRotated

def RandomPerpDir(ref: np.ndarray) -> np.ndarray:
    """
    This function returns a random perpendicular direction with respect to your reference direction.

    Ref: Reference direction
    """
    if ref[0] == 0: return np.array([1.0,0.0,0.0])
    elif ref[1] == 0: return np.array([0.0,1.0,0.0])
    elif ref[2] == 0: return np.array([0.0,0.0,1.0])
    else:
        b_z = -(10*ref[0] + 2*ref[1])/dir[2]
        randperp = np.array([10,2,b_z])
        randperp /= np.linalg.norm(randperp)
        return randperp     