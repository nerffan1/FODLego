#Description: This file contstains as set of classes that will implement the FOD Heuristic Solution
#  following the paradigm of Object-Oriented Programming (OOP). Encapsulation for FOD_Structure, FODs,
#  FOD_Orbital, among other things will be included in this file.
# Roadmap: Use polymorphism for Closed Hybrid Shells (e.g. sp3).
#  -  Add Tetrahedra class and several attributes/methods to manipulate them
#  - In far future, somehow implement the triaugmented triangular prism that corresponds to sp3d5 ( 9 FODs, 18 electrons) 
#Author: Angel-Emilio Villegas S.
import logging
logger = logging.getLogger(__name__)
import numpy as np
from numpy import sqrt
from numpy.linalg import norm
from typing import List
from scipy.spatial.transform import Rotation as R
from scipy.spatial.distance import cdist
from rdkit.Chem import GetPeriodicTable as PT
#
# FODLego Dependencies
from FODLego.Bond import *
from FODLego.FOD import *
import FODLego.Shells as Shells
from  .globaldata import GlobalData

################# FOD STRUCTURE #################

class Atom:
    def __init__(self, index: int, Name: str, Pos, owner):
        #Undetermined Attributes
        self.mSteric = 0
        self.mFreePairs = 0
        self.mCharge = 0  # In the future can be changed 
        self.mBonds = []
        self.mGlobalBonds = []
        self.mFODStruct = FODStructure(self)
        self.mCompleteVal = False
        #Known Attributes
        self.mName = Name
        self.mPos = np.array(Pos)
        self.mI = index
        self.mZ = PT().GetAtomicNumber(Name)
        self.mPeriod = GlobalData.GetPeriod(self.mZ)
        self.mGroup = GlobalData.GetRow(self.mZ)
        self.mValCount = 0
        self.mOwner = owner
        

    def GetMonoCovalRad(self): 
        elecs = GlobalData.GetFullElecCount(self.mGroup, self.mPeriod)
        return GlobalData.mRadii[elecs][self.mZ]

    def GetMonoCovalEdge(self):
        """
        Get the FOD edge distance of a monoatomic calculation
        """
        elecs = GlobalData.GetFullElecCount(self.mGroup, self.mPeriod)
        return GlobalData.mVert[elecs][self.mZ]

    def GetLastAtomRadius(self):
        dist = np.linalg.norm(self.mFODStruct.mCore[-1].mPos - self.mPos)
        return dist

    def GetValenceFODs(self):
        return self.mFODStruct.mValence

    def GetBonds(self):
        return self.mBonds

    def GetBFODs(self):
        return self.mFODStruct.mBFODs
    
    def GetFFODs(self):
        return self.mFODStruct.mFFODs

    def AddBond(self, atom2, order: int):
        self.mBonds.append(Bond(self, atom2,order))

    def AddBFOD(self, fod):
        self.mFODStruct.mBFODs.append(fod)
    
    def GetVec2BFODs(self):
        return [ x.mPos - self.mPos for x in self.mFODStruct.mBFODs]

    def GetVectoNeighbors(self):
        """
        This function returns the sum of the vectors that start on
        neighboring (bonded) atoms and end on current atom.
        """
        freedir = []
        for bond in self.mBonds:
            at2 = bond.mAtoms[1]
            freedir.append(at2.mPos - self.mPos)
        return freedir

    def AverageBFODDir(self):
        """
        Returns the average vector of all the of all BFODs in atom. It is of particular use for finding the direction of FFODs.
        If the displacement is too small a planar structure is suggested; thus, a cross-product is taken in order
        TODO: Might require reimplementation for cases where we just want the ATOM-ATOM vector, instead of the FOD vectors?
        TODO: Use magic numbers in some list elsewhere
        TODO: Resolve whether to use Atoms or FODs as reference
        """
        resultant = np.zeros(3)
        bfods = self.mFODStruct.mBFODs
        for bfod in bfods:
            # If the ffod atom is Meek, then the vector points in its direction! So, no issue here!
            if self == bfod.mMeek:
                resultant += bfod.mBondDir
            # When the atom is Bold, the BondDir points away, so you must get the negative
            else:
                resultant -= bfod.mBondDir

        # If the average displacement is too small, then the 3 points are planar
        resultant /= len(self.mFODStruct.mBFODs)

        if np.linalg.norm(resultant) < .2:
            vecs = [bfod.mPos-self.mPos for bfod in bfods]
            vec1 = tofrom(vecs[0],vecs[1])
            vec2 = tofrom(vecs[1],vecs[2])
            cross = np.cross(vec1,vec2)

            # Get average distance from other FODs
            d = np.mean([norm(b) for b in vecs])
            # Set norm and direction
            resultant = normalize(cross)*d

        # Return the average
        return resultant

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
    
    def FindValence(self):
        """
        This method finds the number of electrons in the valence shell by finding the 
        difference between the current Group and the last ClosedShell Group. Only works up
        to 5th period.
        """
        # TODO: Assume in the meantime that bonded atoms have are closed shell
        if len(self.mBonds) > 0:
            if self.mGroup < 4:
                return self.mGroup
            else:
                if self.mPeriod < 4:
                    return (2 + (self.mGroup - 12))
                elif self.mPeriod < 6:
                    return (self.mGroup)
        else:
            return 0

                    
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

    #Additional Functions
    def GetAssocEdges_B_F_FOD(self,typ=FOD):
        """
        Returns the inter FFOD-BFOD distances.
        """
        dists = []
        if len(self.mFODStruct.mFFODs) > 0 and exists(typ,self.mFODStruct.mFFODs):
                ffods = [x.mAssocFOD.mPos for x in self.mFODStruct.mFFODs]
                bfods = [x.mAssocFOD.mPos for x in self.mFODStruct.mBFODs]
                pairdD = cdist(ffods,bfods)
                dd = np.tril(pairdD)
                dd[dd==0] = np.nan
                dists = dd[~np.isnan(dd)]
        return dists

    def get_assoc_radii_b_f_fod(self,typ=FOD):
        """
        Returns the inter FFOD-BFOD distances.
        """
        bfod_radii = []
        ffod_radii = []
        if len(self.mFODStruct.mFFODs) > 0 and exists(typ,self.mFODStruct.mFFODs):
            # Get the radii from the atom using the dominant/schema
            # If the atom owning the FOD is the Bold atom, then the distance
            # to that atom is mBoldR

            # BFOD radii
            for bfod in self.GetBFODs():
                    d = dist(bfod.mAssocFOD.mPos, self.mPos)
                    bfod_radii.append(d)

            # FFOD Radii
            from FFOD import SFFOD
            for ffod in self.GetFFODs():
                ffod_radii.append(ffod.mAssocFOD.mR)
                if ffod.mAssocFOD.mR < .75 and self.mZ == 7 and typ == SFFOD:
                    print(self.mOwner)

        # Return the rad
        return bfod_radii, ffod_radii

    # Additional Functions
    def __str__(self):
        pass    

######################## FOD Structure Class  ########################
#HOW TO concatenate FODs, easily
class FODStructure:
    def __init__(self, parent: Atom):
        self.mAtom = parent
        self.mCore = [] #A list of FODShells
        self.mCoreShells = []
        self.mValence = [] #A list of FODs
        self.mBFODs: List[BFOD] = []
        self.mFFODs: List[FOD] = []
        self.mfods = [] #All Finalized FODs
        self.mLastBond = 0
    
    # Setter Functions
    def AddValFOD(self, fods: List[np.ndarray], at2=False,finalize=False):
        """
        This function adds the FODs to the valence of the current atom's FOD Structure.
        It also accepts a secondary atom, at2, in order to add the FOD to its valence. The
        Finalize parameter is True only for the first atom that is writing the FODs; by finalizing
        we mean that the FOD is added to the list of overall FODs.
        TODO: Currently only takes positions. Do we want to also add it here?
        """
        # Add to current valence
        for fod in fods:
            #Assertions
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

    def _AddCoreShell(self, shell):
        """
        This function adds the Shell object to the FODStructure object
        and it also adds the FODs of that Shell object to the mCore of 
        the FODStructure.

       Args:
            shell: The Shell object that will be added to the FOD 
            Structure.
            ch: Boolean representing the alpha/beta (up/down) channels 
            of the electrons. True for Alpha, False for Beta.
        """
        self.mCoreShells.append(shell)
        # Add individual FODs to the electronic structure
        for fod in shell.mfods:
            GlobalData.mFODs.append(fod)
            self.mCore.append(fod)

    def PrepareShells(self, atoms: List[Atom]):
        """
        This function will determine the creation of Core FODs, those 
        that are not related to bonding. The scheme is easy in the first
        3 periods of the Periodic Table, but it will become trickier 
        ahead if Hybridization heuristics don't work. Currently it only 
        works for closed shell calculations (V 0.1).
        """
        # TODO: This will assume that we are doing up to the n=3 shell,
        #  with sp3 hybridization
        # TODO: Take into account free pairs of electrons
        # TODO: For mOrder=2, there are many schemes to determine
        #  direction
        # TODO: Find an elegant solution to do exceptions for H bonds
        # TODO: Account previous bonds formed, by talling previous FODs,
        #  or looking back at mBonds 
        # TODO: GlobalData.GetFullElecCount() Can be precalculated ahead
        #  of time and placed as a member vatrable

        #Lazy loading in order to 
        from FODLego.BFOD import SBFOD, DBFOD, TBFOD
       
        at1 = self.mAtom

        def SingleBond(at2: Atom, curr_bond: Bond):
           """
           Creates a SBFOD along the axis. The heuristics is found in SBFOD Class
           """
           #Add the BFODs, new way
           boldMeek = BoldMeek(at1,at2)
           newFOD = SBFOD(*boldMeek)
           _AddBFOD(curr_bond, at1, at2, newFOD)

        def DoubleBond(at2: Atom, curr_bond: Bond):
            """
            Create the FODs representing the double bond. Currently the FOD filling is unidirectional (sequential)
            and  does not account for the next atom in the iteration to see if we can further accomodate the bonding FODs 
            """
                            
            def HeightDir_fromNeighborBFODs():
                """
                This returns the unit-vector of the direction that a DBFOD is displaced away from 
                the bonding axis (called 'Height' throughout this code). It is obtained by a series of steps:
                1) The cross-product of the FOD-Atom-FOD is obtained
                2) Measure the angle between the Bonding Axis (BA) and the vector found in (1) 
                3) If beyond a certain threshold, then rotate the direction 
                
                """
                # (1) Obtain the FOD-Atom-FOD Cross Product. The height in a planar structure
                vector_for_cross = []
                for otherb in self.mAtom.mBonds:
                    if otherb != curr_bond:
                        vector_for_cross.append(otherb.mAtoms[1].mPos)
                vector_for_cross -= self.mAtom.mPos
                D = np.cross(*vector_for_cross)
                D = normalize(D) 

                # (2) Measure Angle between D and the Bonding Axis (BA) 
                BA = at2.mPos - at1.mPos
                angle = AngleBetween(BA,D)
                
                #(3) Check for orthogonality
                thres_max = 1.01*(np.pi/2) # +1% Deviation from 100%
                thres_min = .99*(np.pi/2) # -1% Deviation from 100%
                if angle > thres_max or angle < thres_min:
                    #Get axis of rotation
                    axis = np.cross(BA,D)
                    axis = normalize(axis)
                    #Get rotation angle
                    if angle > thres_max: 
                        angle_diff =  np.pi/2 - angle
                    elif angle < thres_min:
                        angle_diff = (np.pi/2) - angle
                    #Rotate
                    D = RotateVec(D, axis*angle_diff)
 
                return D

            def D_BFOD_Direction():
                """
                - Atoms with 3 bonds: Then the DBFODs are based of the other
                  two FODs. It accounts for bent molecules. Prime example: C60
                - Atoms with 2 bonds: The 
                """
                if self.mAtom.mFreePairs == 0:
                    #Determine DFFOD Direction
                    if len(self.mAtom.mBonds) == 3:
                        return  HeightDir_fromNeighborBFODs()
                    elif len(self.mAtom.mBonds) == 2:
                        if self.mAtom.mBonds.index(curr_bond) == 0:
                            return RandomPerpDir(dir)
                        else:
                            # Cross product between atom and already-placed FODs
                            vector_for_cross = []
                            for fods in self.mBFODs:
                                vector_for_cross.append(fods.mPos)
                            vector_for_cross -= self.mAtom.mPos
                            return np.cross(*vector_for_cross)
                elif self.mAtom.mFreePairs == 1:
                    vector_for_cross = []
                    for otherb in self.mAtom.mBonds:
                            vector_for_cross.append(otherb.mAtoms[1].mPos)
                    vector_for_cross -= self.mAtom.mPos
                    D = np.cross(*vector_for_cross)
                    return normalize(D)
                else:
                    return RandomPerpDir(dir)

            #Information
            axis2fod = np.ndarray(3)
            dom,sub= BoldMeek(at1,at2)

            if GlobalData.GetFullElecCount(self.mAtom.mGroup,self.mAtom.mPeriod) <= 18:
                #Find perpendicular unit vector
                    dir = sub.mPos - dom.mPos
                    axis2fod = D_BFOD_Direction()
                    # Create FODs and link
                    f1 = DBFOD(dom,sub,axis2fod)
                    f2 = DBFOD(dom,sub,-axis2fod)
                    _AddBFOD(curr_bond, dom, sub, f1, f2) # Does dom/sub matter, or are at1/at2 fine?

        def _AddBFOD(curr_bond: Bond, at1: Atom, at2: Atom, *fods):
            """
            This function adds a new FOD to the individual atoms, to the list in GlobalData, and to the FODStructure
            """ 
            # TODO: Make Bonds easier to deal with by having one instance instead of one per bond per atom (i.e. there are 2 instances of Bond that are slightly for both atoms in a bond)
            # The main purpose of this function is not to not duplicate the FOD in GlobalData
            # by adding the FODs in each individual atom. This makes this class a type of 
            # FOD manager in addition to constructing the structure.
            # The reason why we don't load FODs directly
            # is because we don't know whether at1 or at2
            # is the self.mAtom

            # Create siblings
            curr_bond.SetFODs(fods)
            if len(fods) == 2:
                fods[0].AddSibling(fods[1])
                fods[1].AddSibling(fods[0])
            elif len(fods) == 3:
                fods[0].AddSibling(fods[1],fods[2])
                fods[1].AddSibling(fods[0],fods[2])
                fods[2].AddSibling(fods[0],fods[1])
            # Add to atoms and globaldata
            for fod in fods:
                at1.AddBFOD(fod)  # Maybe remove this, and instead do a getter function
                at2.AddBFOD(fod)
                GlobalData.mFODs.append(fod)
                GlobalData.mBFODs.append(fod)
        
        def _AddFFOD(*ffods):
            """
            TODO: Put this on a bigger scope
            """
            # Create siblings. Manually seemed the fastest way to implement.
            if len(ffods) == 2:
               ffods[0].AddSibling(ffods[1])
               ffods[1].AddSibling(ffods[0])
            elif len(ffods) == 3:
               ffods[0].AddSibling(ffods[1],ffods[2])
               ffods[1].AddSibling(ffods[0],ffods[2])
               ffods[2].AddSibling(ffods[0],ffods[1])
            for ffod in ffods:
                self.mFFODs.append(ffod)
                GlobalData.mFODs.append(ffod)          

        def AddFreeElectron(free: int):
            """
            TODO: Change conditionals as to create more concise code 
            TODO: Create a series of variables for chosen constants, NO magic numbers
            """

            def ChoosePerpDir():
                """
                Return the direction of the double bond FODs
                """
                heightdir = np.array([0.0,0.0,0.0])
                if len(self.mAtom.mBonds) == 2:
                    heightdir = np.cross(*[fod.mPos - at1.mPos for fod in self.mBFODs],axis=0)
                elif len(self.mAtom.mBonds) == 1:
                    heightdir = np.cross(*[fod.mPos - at1.mPos for fod in self.mBFODs],axis=0)
                elif len(self.mAtom.mBonds) == 3:
                    heightdir = np.cross(*[bond.mAtoms_p[1] - bond.mAtoms_p[0] for bond in self.mAtom.mBonds], axis=0)

                # Add a check in case the cross product gives 0
                if (heightdir == np.array([0.0,0.0,0.0])).all():
                    heightdir = RandomPerpDir(self.mBFODs[0].mPos)
                return normalize(heightdir)

            if len(self.mAtom.mBonds) == 1:
                #Useful Information
                at2 = self.mAtom.mBonds[0].mAtoms[1]
                free_dir = self.mAtom.mPos - at2.mPos

            if free == 2:
                from FODLego.FFOD import DFFOD
                heightdir = ChoosePerpDir()
                a = DFFOD(at1,heightdir)
                b = DFFOD(at1,-heightdir)
                _AddFFOD(a,b)

            elif free == 3:
                from FODLego.FFOD import TFFOD
                dir0 = RandomPerpDir(free_dir)
                norms = RotateNormals(3, dir0, normalize(free_dir)) 
                f1 = TFFOD(at1, norms[0])
                f2 = TFFOD(at1, norms[1])
                f3 = TFFOD(at1, norms[2])
                _AddFFOD(f1,f2,f3)

        def TripleBond(at2: Atom, curr_bond: Bond):
            """
            #TODO: Create a helper funtion for conditional statements
            """
            # Place BFODs, for new implementation
            bonddir = tofrom(at2.mPos,at1.mPos)
            boldmeek = BoldMeek(at1,at2)
            dir0 = RandomPerpDir(bonddir)
            norms = RotateNormals(3, dir0, normalize(bonddir)) 
            # Create BFODs and link
            f1 = TBFOD(*boldmeek, norms[0])
            f2 = TBFOD(*boldmeek, norms[1])
            f3 = TBFOD(*boldmeek, norms[2])
            _AddBFOD(curr_bond, at1, at2, f1, f2, f3)

        def AddBFODs():
            """
            Finish determining BFODs. This is done after initializing all initial BFODs.
            """
            for bond in self.mAtom.mBonds:
                bonded_at = bond.mAtoms[1]
                if bonded_at.mI  > self.mAtom.mI:
                    if bond.mOrder == 1:
                        SingleBond(bonded_at, bond)
                    elif bond.mOrder == 2:
                        DoubleBond(bonded_at, bond)
                    elif bond.mOrder == 3:
                        TripleBond(bonded_at, bond)

        def AddFFODs():
            from FODLego.FFOD import SFFOD, DFFOD, TFFOD
            if self.mAtom.mFreePairs == 2:
                if self.mAtom.mSteric >= 3:
                    AddFreeElectron(2)
            elif self.mAtom.mFreePairs == 1:
                if self.mAtom.mSteric >= 2:
                    _AddFFOD(SFFOD(self.mAtom))
            elif self.mAtom.mFreePairs == 3:
                AddFreeElectron(3)

        def AddCoreElectrons():
            #Count core electrons and
            core_elec = self.mAtom.mZ - self.mAtom.mValCount
            logger.debug(f"\n Atom: {self.mAtom.mName}\n Core electrons: {core_elec}")
            if core_elec != 0:
                for shell in GlobalData.mGeo_Ladder[core_elec]:
                    if shell == 'pt':
                        self._AddCoreShell(Shells.Point(self.mAtom))
                        if self.mAtom.mOwner.mOpen == True:
                            self._AddCoreShell(Shells.Point(self.mAtom, False))
                    elif shell == 'tetra':
                        self._AddCoreShell(Shells.Tetra(self.mAtom, 10))
                        if self.mAtom.mOwner.mOpen == True:
                            self._AddCoreShell(Shells.Tetra(self.mAtom, 10, False))
                    elif shell == 'a_triaug':
                        pass # For future development: Beyond scope


        # Prepare the valence shell first, since it will help determine the
        # orientation of the inner shells
        AddBFODs()
        AddCoreElectrons()
        AddFFODs()

        # Define Valence
        self.mValence = self.mBFODs + self.mFFODs
            
################# ADDITIONAL FUNCTIONS #################
def BoldMeekDir(at1: Atom, at2: Atom, all=True):
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

def BoldMeek(at1: Atom, at2: Atom):
    """
    This Function determines the dominant and meek atom in the bonding.
    at1: An atom
    at2: An atom bonding to at2
    """
    if at1.mPeriod < at2.mPeriod:
        dom = at1
        sub = at2
    elif at1.mPeriod > at2.mPeriod:
        dom = at2
        sub = at1
    elif at1.mZ > at2.mZ:
            dom = at1
            sub = at2
    elif at1.mZ <= at2.mZ:
            dom = at2
            sub = at1
    # Either return dom and sub, or just the dominant atom. 
    return dom, sub 

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

def InverseSqRatio(dom, sub):
    Z1 = dom.mZ
    Z2 = sub.mZ
    if dom.mZ == sub.mZ:
        return 0.5
    else:
        Z1 = 0.4 if Z1 == 1 else Z1
        Z2 = 0.4 if Z2 == 1 else Z2
        r = sqrt(Z1/Z2)
        g = r/(1+r)
        if g < 0.5:
            return g
        else:
            return (1-g)
