from FODLego.Funcs import *
from FODLego.ElementaryClasses import *
from FODLego.globaldata import *

class BFOD(FOD):
    """
    The BFOD class is in charge of encapsulating the reparametrization of Bonding FODs. This class describes the FODs using the angles, distances, heights, and ratios
    along the bonding axis in order to characterize the FODs. Some of these attributes of this class will tend toward zero in the SBFOD, but we know from some examples
    that they are not always zero (e.g. they might lie slightly off the bonding axis, so there might be an angle).
    """
    def __init__(self, boldAt: Atom, meekAt: Atom, target=None, pol=True):
        super().__init__(ch=pol)
        # Atoms
        self.mBold = boldAt
        self.mMeek = meekAt
        # Angles
        self.mBoldAngle = 0.0
        self.mMeekAngle = 0.0
        # Vectors
        self.mBondDir = meekAt.mPos - boldAt.mPos  # Always in direction away fromBold atom
        self.mHeight = np.zeros(3)
        # Distances
        self.mBondDist = np.linalg.norm(self.mBondDir)
        self.mMeekR = 0.0
        self.mBoldR = 0.0
        # Misc 
        self.mBoldPortion = 0.0
        # Reverse Determination
        if isinstance(target, np.ndarray):
            self.mPos = target
            self.RevDet()

    def Calc_AxisBoldPortion(self, Zbold:int, Zmeek:int) -> float:
            """
            Finds the portion (from 0 to 1) of the bonding distance that the Bold atom covers.
            For example, the FODs might be placed at the 0.33 mark of the line starting
            from the bold atom, ending on the meek atom.   
 
            Parameters
            ----------
            - Zbold : The atomic number of the dominant (bold) atom.
            - Zmeek : The atomic number of the weaker (meek) atom.      
            
            Returns
            -------
            Float
            The portion of the bonding axis covered by the dominant atom. 
            """
            if Zbold == Zmeek:
                return 0.5 #Equal atoms, then return the midpoint
            else:
                Zbold = 0.4 if Zbold == 1 else Zbold
                r = sqrt(Zmeek/Zbold)
                g = r/(1+r)
            # The portion should always be under 0.5 since it is the closer to the
            # dominant atom.
            if g < 0.5:
                return g
            else:
                return (1-g)

    def IsMonoatomic(self) -> bool:
        """
        This function returns whether or not we are dealing with a monoatomic bond
        """
        if self.mMeek.mZ == self.mBold.mZ:
            return True
        else:
            return False

    def RevDet(self):
        """"
        This function uses the target FOD (e.g. an optimized FOD from FLOSIC) and returns the parameters
        """
        # Helper data
        bold2fod = self.mPos -  self.mBold.mPos
        meek2fod = self.mPos - self.mMeek.mPos
        # Angles
        self.mBoldAngle = AngleBetween(bold2fod, self.mBondDir)
        self.mMeekAngle = AngleBetween(meek2fod, -self.mBondDir)
        # Distances
        self.mMeekR = norm(meek2fod)
        self.mBoldR = norm(bold2fod)
        # Misc 
        self.mBoldPortion = (np.cos(self.mBoldAngle)*self.mBoldR)/self.mBondDist
        if self.mBoldPortion > 0.5:
            self.mBoldPortion
        if (self.mBoldPortion > 1):
           print('Reverse Determination issue. Axial projection proportion over 1.')

    def PrintParams(self):
        print(f"{self.mBold.mName}-{self.mMeek.mName}")
        print(f"BoldR: {self.mBoldR}")
        print(f"MeekR: {self.mMeekR}")
        print(f"Bold Theta: {np.rad2deg(self.mBoldAngle)}")
        print(f"Meek Theta: {np.rad2deg(self.mMeekAngle)}")
        print(f"BoldPortion: {self.mBoldPortion}")

class SBFOD(BFOD):
    """
    This is the Single Bonding FOD (SBFOD) class
    """
    def __init__(self, bold: Atom, meek: Atom, ch=True):
        super().__init__(bold,meek, pol=ch)
        self.mBoldPortion = self.Calc_AxisBoldPortion(bold.mZ, meek.mZ)
        self.DetermineParamenters()

    def DetermineParamenters(self):
        """
        Add the single BFOD along the axis. 
        """
        bfod = self.mBondDir*self.mBoldPortion
        self.mPos = self.mBold.mPos + bfod
        self.mBoldR = self.mBondDist*self.mBoldPortion
        self.mMeekR = self.mBondDist*(1-self.mBoldPortion)

class DBFOD(BFOD):
    def __init__(self, bold: Atom, meek: Atom, heightdir: np.ndarray, ch=True):
        super().__init__(bold,meek,pol=ch)
        self.mHeight = heightdir
        self.DetermineParameters()

    def GetHeight(self, r) -> float:
        """
        Return radial distance from interatomic (bonding) axis. 
        For heterogenous atoms....
        """
        rad = self.mBold.GetMonoCovalRad()
        val = np.sqrt(rad**2 - r**2)
        if val < 0:
            print("The DBFOD prediction cannot be simply determined. ")
        return val

    def GetBondAxProj(self) -> float:
        """
        This function determines the ratio between the projection of the FOD-ATOM
        and the Bonding Axis.
        """
        if self.mBold.mZ == self.mMeek.mZ:
            return 0.5
        else:
            if self.mBold.mPeriod == self.mMeek.mPeriod:
                return InverseSqRatio(self.mBold, self.mMeek)
            else:
                self.mBoldAngle = np.deg2rad(54)
                proj = np.cos(self.mBoldAngle)*self.mBold.GetMonoCovalRad()
                return proj/self.mBondDist

    def DetermineParameters(self):
        """
        Determine the Atom-Atom-FOD angles.
        """
        # Set BondProjection
        self.mBoldPortion = self.GetBondAxProj()

        # Set FOD Position
        delta_bond = self.mBondDir*self.mBoldPortion  # TODO: Just use AxialPoint_Simple?
        delta_height = self.mHeight*self.GetHeight(np.linalg.norm(delta_bond))

        # Finalize parameters
        self.mPos = self.mBold.mPos + delta_bond + delta_height

        # Measure Meek Angle-Distance
        toFOD = self.mMeek.mPos - self.mPos
        self.mMeekAngle = AngleBetween(self.mBondDir, toFOD)
        self.mMeekR = np.linalg.norm(toFOD)

        # Measure Bold Angle-Distance
        toFOD = self.mPos - self.mBold.mPos
        self.mBoldAngle = AngleBetween(self.mBondDir, toFOD)
        self.mBoldR = np.linalg.norm(self.mBold.mPos - self.mPos)

class TBFOD(BFOD):
    def __init__(self, bold: Atom, meek: Atom, heightdir: np.ndarray):
        super().__init__(bold,meek)
        self.mHeight = heightdir 
        self.mBoldAngle = np.deg2rad(54.735) 
        self.DetermineParameters()

    def DetermineParameters(self):
        """
        This function determines the distance of the TBFOD away from the bonding
        axis. 
        """
        # Determine the location of the FOD
        rad = self.mBold.GetMonoCovalRad()
        a = self.mBold.GetMonoCovalEdge()
        c = self.mBondDist 
        # This section modifies the height depending on the nature of the monoatomic bond
        # 2nd period elements tend to tighten BFODs.
        if self.IsMonoatomic():
            if self.mBold.mPeriod == 2 and self.mMeek.mPeriod == 2:
                theta = np.arccos((c/2)/rad)
            elif self.mBold.mPeriod > 2 and self.mMeek.mPeriod > 2:
                theta = np.arctan(rad/(c/2)) 
        else:
            #theta = np.arcsin((a/2)/rad) 
            theta = np.deg2rad(54.735)
        self.mBoldAngle = theta
        dr = normalize(self.mBondDir)*rad*np.cos(theta) 
        dl = self.mHeight*rad*np.sin(theta)
        self.mPos = self.mBold.mPos + dr + dl

        # Determine the Meek Atom angle
        self.DetermineMeek()

        # Determine Distances
        self.mBoldR = dist(self.mPos, self.mBold.mPos)
        self.mMeekR = dist(self.mPos, self.mMeek.mPos)
        self.mBoldPortion = (np.cos(self.mBoldAngle)*self.mBoldR)/c

    def DetermineMeek(self):
        """
        Determines the Bond-Meek-FOD angle. I.e., the angle between the bonding axis
        and the Meek-FOD vector.
        """
        self.mMeekAngle = AngleBetween(-self.mBondDir, self.mPos - self.mMeek.mPos)

    def DetermineHeight(self):
        pass
