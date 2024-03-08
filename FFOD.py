from Funcs import *
from FOD import *
from ElementaryClasses import *
from copy import copy
class FFOD(FOD):
    def __init__(self, atom: Atom, HeightDir = np.zeros(3), target=None):
        super().__init__(copy(atom.mPos))
        self.mAtom = atom
        self.mAngle = 0.0
        self.mR = 0.0
        # Vectors
        self.mHeight = HeightDir
        self.mFreeDir = atom.AverageBFODDir()
        # Reverse Determination parameters
        if isinstance(target, np.ndarray):
            self.mPos = target

class SFFOD(FFOD):
    def __init__(self, atom: Atom):
            super().__init__(atom)
            self.DetermineParameters()

    def DetermineParameters(self) -> None:
        """
        This function determines the distance away from parent atom of the SFFOD. 
        """
        # Leave average of bonding atoms.
        if len(self.mAtom.mBonds) == 3 :
            dr = self.mFreeDir

        #Perhaps add a restriction here later
        elif len(self.mAtom.mBonds) == 2 :
            dr = normalize(self.mFreeDir)*self.mAtom.GetMonoCovalRad()
        
        # The following section has 2 options. In atoms with just a 1s core shell 
        # tend to be closer to the nucleus when their valence is complete.
        elif len(self.mAtom.mBonds) == 1:
            if (self.mAtom.mPeriod < 3):
                dr = normalize(self.mFreeDir)*self.mAtom.GetMonoCovalEdge()/sqrt(8)  #Place at midsphere distance?
            else:
                dr = self.mFreeDir
        # Set position
        self.mPos += dr

    ### Measurement ###
    # Description: This is after associating an FOD point from an optimized output to a "predicted" 
    # FOD (generated by this program). 
    def MeasureR(self) -> float:
        return np.linalg.norm(self.mPos - self.mAtom.mPos)

    def ReverseDetermination(self, targetFOD: np.ndarray):
        """
        Determines the parameters
        """
        atom2ffod = targetFOD - self.mAtom
        self.mR = norm(atom2ffod)
        self.mAngle = 0.0
        self.mHeight = 0.0
        self.mFreeDir = normalize(atom2ffod)

class DFFOD(FFOD):
    def __init__(self, atom: Atom, heightdir: np.ndarray):
        super().__init__(atom, heightdir)
        self.DetermineParameters()

    def DetermineParameters(self):
        #Direction away from atom, on plane of other bonds, or neighboring atoms
        self.mR = self.DetermineR()
        self.mPos += normalize(self.mFreeDir)*self.mR*np.cos(self.mAngle)
        self.mPos += self.mHeight*self.mR*np.sin(self.mAngle)

    def DetermineR(self) -> float:
        # Determine the angle at which the DFFODs open! 220 rule
        neighbors = self.mAtom.GetVectoNeighbors()
        bfod = self.mAtom.GetBFODs()[-1]
        if self.mAtom == bfod.mMeek:
            F = bfod.mMeekR
        else:
            F = bfod.mBoldR

        if len(self.mAtom.mBonds) == 2:
            theta = (np.deg2rad(220) - AngleBetween(*neighbors))/2
        else: 
            theta = (np.deg2rad(220) - AngleBetween(*self.mAtom.GetVec2BFODs()))/2
        self.mAngle = theta
    
        if self.mAtom.mPeriod < 3:
            # FreeDirection and Atom-BFOD Angle, plust lengths needed for determination
            vec2bfod = bfod.mPos - self.mAtom.mPos
            phi = AngleBetween(self.mFreeDir,vec2bfod) 
            E = self.mAtom.GetMonoCovalEdge()
            
            # Finalize 
            return  np.sqrt(E**2 - F**2 + 2*F*E*np.cos(phi)*np.cos(theta))
        else:
            return F

    def ReverseDetermination(self, targetFOD: np.ndarray):
        """
        Determines the parameters from an FOD of choice
        """
        atom2ffod = targetFOD - self.mAtom
        self.mR = norm(atom2ffod)
        self.mAngle = 0.0
        self.mHeight = 0.0
        self.mFreeDir = normalize(atom2ffod)

class TFFOD(FFOD):
    def __init__(self, atom: Atom, heightdir: np.ndarray):
        super().__init__(atom, heightdir)
        self.DetermineParameters()

    def DetermineParameters(self):
        # Determine FreeDir projection
        if self.mAtom.mPeriod < 3:
            h = self.DetermineHeight_Tight()
        else:
            # Get Free_direction shift
            h = self.DetermineHeight()
        # Finalize Position by giving height and bonding direction offsets
        self.mPos += self.mHeight*h
        self.mPos += self.mFreeDir  # This was set in DetermineHeight()
        # Measure angle
        self.mAngle = AngleBetween(self.mFreeDir,self.mPos-self.mAtom.mPos)

    def DetermineHeight_Tight(self) -> float:
        # This length is the radius of a circumscribed circle on a triangle
        # of sides equal to the Atom-SBFOD distance.
        h = self.GetSBFODDist()/sqrt(3)
        self.mFreeDir = normalize(self.mFreeDir)*h/np.tan(np.deg2rad(70.52))
        return h

    def DetermineHeight(self) -> float:
        """
        Determines the "height" of the FFOD: The distance from the bonding axis
        to the FFOD
        """
        # Use the distance of last shell to place this shell
        proj_freedir = self.mAtom.GetLastAtomRadius()
        # Changing the mfreedir to end at projection point
        self.mFreeDir = normalize(self.mFreeDir)*proj_freedir
        # For the FOD-Atom distance use the distance to SBFOD in atom valence
        l = self.GetSBFODDist()
        # Return the height
        return sqrt(l**2 - proj_freedir**2)
    
    def GetSBFODDist(self) -> float:
        if len(self.mAtom.mBonds) == 1:
            bfod = self.mAtom.GetBFODs()[0]
            if bfod.mMeek == self.mAtom:
                return bfod.mMeekR
            else:
                return bfod.mBoldR
